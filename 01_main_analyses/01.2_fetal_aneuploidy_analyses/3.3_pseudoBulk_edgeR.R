## pb with EdgeR - transcriptional impact of abnormal karyotype on fLiver cells

outDir = '~/lustre_mt22/Aneuploidy/Results/3.3_pseudoBulk_edgeR_withCyclCells/nov24/without_published2n/genoAssay_withX/'
formula = '~ %s'

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)




##-------------------##
##   Libraries     ####
##-------------------##

library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(readxl)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")



fit_model <- function(pb,colDat,formula,geneMap,groupID='group',MDS_groups = c('Genotype','donorID','sex'),coef=NULL){
  
  if(!grepl('%s',formula,fixed=TRUE) || grepl('%s',sub('%s','',formula,fixed=TRUE),fixed=TRUE))
    stop("Invalid formula, must contain one instance of %s")
  #Convert it now in case it is still invalid in more subtle ways
  formula = as.formula(sprintf(formula,groupID))
  
  
  # create an edgeR object with counts and grouping factor
  y <- DGEList(pb, group = colDat$Genotype)
  y$samples = cbind(y$samples,colDat[match(rownames(y$samples),colDat$donorID),!colnames(colDat) %in% colnames(y$samples)])
  
  ## Plot library sizes
  df = y$samples
  df$donorID = rownames(df)
  p = ggplot(df,aes(reorder(donorID,`lib.size`),fill=group,y=`lib.size`))+
    geom_col()+
    scale_fill_manual(values = col25)+
    theme_classic(base_size = 12)+xlab('')+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          axis.line = element_blank(),axis.text = element_text(color='black'),
          panel.border = element_rect(fill=F,colour = 'black'))
  
  
  
  print(p)
  
  
  
  # filter out genes with low counts
  print("Dimensions before subsetting:")
  print(dim(y))
  print("")
  all(y$samples$lib.size>1e5)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  print("Dimensions after subsetting:")
  print(dim(y))
  print("")
  
  
  
  ## Normalization 
  y <- calcNormFactors(y)
  # MDS plot
  cluster<-as.factor(colDat$Genotype[match(rownames(y$samples),colDat$donorID)])
  df = plotMDS(y, pch=16,col=col25[cluster], main="MDS") 
  df = data.frame(x = df$x,y=df$y,donorID = names(df$x))
  df = cbind(df,colDat[match(df$donorID,colDat$donorID),!colnames(colDat) %in% colnames(df)])
  for(f in MDS_groups){
    df$group = df[[f]]
    p = ggplot(df,aes(x,y,col=group))  +
      geom_point()+
      theme_classic(base_size = 10)+
      scale_color_manual(values = col25)+
      theme(axis.line = element_blank(),axis.text = element_text(color='black'),
            panel.border = element_rect(fill=F,colour = 'black'))
    print(p)
  }
  
  
  # create a design matrix: here we have multiple donors so also consider that in the design matrix
  design <- model.matrix(formula,data=y$samples)
  
  # estimate dispersion
  y <- estimateDisp(y, design = design,robust = T)
  # fit the model
  fit <- glmQLFit(y, design)
  plotQLDisp(fit)
  
  # Extract coefficients
  if(is.null(coef)){
    contrast = makeContrasts(groupT21-groupdiploid, levels=colnames(design))  
    qlf<-glmQLFTest(fit, contrast=contrast)
  }else{
    qlf<-glmQLFTest(fit, coef = coef)
  }
  
  
  # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
  tt <- topTags(qlf, n = Inf,p.value = 0.05,)
  tt <- tt$table
  tt = annotateGenes(tt,geneMap = geneMap)
  
  ## Calculate logCPM
  y$logCPM <- cpm(y, log=TRUE, prior.count = 1)
  
  return(list("fit"=fit, "design"=design, "y"=y, 'tt'=tt))
}


##-----------------------##
##        Params          #
##-----------------------##
#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

geneMap = read.delim('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data/Donor15680/liver/fLiver_MY_200531_10043298/filtered_feature_bc_matrix/features.tsv.gz',header = F)
colnames(geneMap) = c('ensID','geneSym','GEX')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))



kMerge=11
keepCylcingCells=T
remove_cyclingCells=F
ageMatch = T


##------------------------------##
##        DGE analysis        ####
##------------------------------##
tissue = 'liver'
for(tissue in c('liver')){
  # Create output directory
  outDir_fp = file.path(outDir,tissue)
  plotDir = outDir_fp
  if(!dir.exists(plotDir)){
    message('Making plot directory')
    dir.create(plotDir,recursive = T)
  }
  
  big.srat = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS')
  big.srat$Sex[big.srat$Sex == 'female'] = 'F'
  big.srat$Sex[big.srat$Sex == 'male'] = 'M'
  big.srat$annot[big.srat$annot == 'Mesenchyme'] = 'Fibroblast'
  big.srat$annot[big.srat$annot == 'Placenta'] = 'Trophoblast'
  big.srat$annot[big.srat$annot == 'NPC'] = 'Neuronal'
  big.srat$finalAnn = as.character(big.srat$annot)
  big.srat$finalAnn[big.srat$finalAnn %in% c('earlyMK')] = 'MK'
  big.srat$finalAnn[big.srat$finalAnn %in% c('promyelocyte','myelocyte')] = 'Myelocyte'
  big.srat$finalAnn[big.srat$finalAnn %in% c('proMono')] = 'Monocyte'
  big.srat$finalAnn[big.srat$finalAnn %in% c('LMPP_ELP','pro.B.cell','pre.B.cell')] = 'B.cell.prog'
  big.srat$finalAnn[big.srat$finalAnn %in% c('ILC.precursor','T.cell','NK_T')] = 'NK.T'
  big.srat$ageGroup = big.srat$gestationalAge
  
  
  # big.srat = filter_sratObj(tissue = tissue,ageMatch=ageMatch,remove_cyclingCells=remove_cyclingCells,
  #                           srat_in_fp='~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_annotated_0424.RDS')
  # big.srat = big.srat[[1]]
  
  
  ###############
  # Exploration #
  ###############
  #Define sample metadata.  The HDBR samples ages are just the initial estimates which are rather inprecise
  metadata = big.srat@meta.data[,c("Genotype",'finalAnn','donorID','gestationalAge','ageGroup','Sex','assay','termination','cellID')]
  sampDat = metadata %>% group_by(donorID,Sex,gestationalAge,ageGroup,assay,termination,Genotype) %>% summarise(nCells = length(unique(cellID))) %>% as.data.frame()
  
  colnames(sampDat) = c('donor','sex','age','ageGroup','tech','type','group','nCells')
  
  # Check that there's no duplicated donorID
  if(n_distinct(sampDat$donor) != nrow(sampDat)){
    stop('Non-unique donorID detected!')
  }
  rownames(sampDat) = sampDat$donor
  
  if(sum(is.na(sampDat))>0){
    sampDat[is.na(sampDat)] = '??'  
  }
  
  
  #Loop first over genotype
  n_deg_perChr_allCT = tibble()
  nCell_perCT_perGeno = data.frame()
  #--------------------------------------
  unique(big.srat$Genotype)
  
  for(geno in c('T21')){
    if(geno %in% c('diploid'))
      next
    print(geno)
    
    ##-----------------------------------------------------------------##
    ##    Define control genes - excluding the trisomic chromosome   ####
    ##-----------------------------------------------------------------##
    
    # Define control genes for Triploid case to be genes on chrX, excluding PAR region and genes which escape X-inactivation
    if(geno == 'complete_trisomy'){
      
      # specify ChrX and chrY pairing regions
      # chrXY_par = GRanges(c('X','X','Y','Y'),
      #                     IRanges(c(10001,155701383,10001,56887903),c(2781479,156030895,2781479,57217415)),
      #                     region = c('PAR1','PAR2','PAR1','PAR1'))
      # parRegion_genes = subsetByOverlaps(gns,chrXY_par)
      # # Remove genes from PAR region from list of control genes
      # controlGenes = controlGenes[!controlGenes %in% parRegion_genes$gene_id]
      # Remove genes which escape X-inactivation, which lies on strata 2-3 mainly
      #Definition from https://academic.oup.com/gbe/article/5/10/1863/522116, which is hg19
      #strata = c(0,2.78,5.04,8.43,30.62,55.78,75.53,99.98,130.82,145.73,155.72)*1e6
      #These are the boundaries lifted over and the PAR definitions from a hg38 source
      strata = c(0,2.781479,5.12,8.46,30.60,55.75,76.31,100.73,131.69,146.65,155.701383)*1e6
      strata = GRanges('X',IRanges(strata,width=c(diff(strata),1e9)))
      names(strata) = paste0('strata',seq_along(strata))
      names(strata)[1] = 'PAR1'
      names(strata)[length(strata)]='PAR2'
      o = findOverlaps(gns,strata)
      gns$strata = 'NA'
      gns$strata[queryHits(o)] = names(strata)[subjectHits(o)]
      
      controlGenes = gns[seqnames(gns) == 'X']
      controlGenes = controlGenes[!controlGenes$strata %in% c('PAR1','PAR2','strata2','strata3')]$gene_id
      
      tgtChrs=c(1:22,'X')
    }else if(geno %in% c('T21','T18','T22','MX')){
      chr = gsub('T|M','',geno)
      controlGenes = gns[!seqnames(gns) %in% c('X','Y',chr)]$gene_id
      if(var %in% c('geno_withX','genoAssay_withX')){
        tgtChrs = c(1:22,'X')
      }else if(var == 'geno'){
        tgtChrs=as.character(c(1:22))
      }
      
    }
    
    
    
    
    ##-----------------------------------##
    ##    Extract relevant metadata    ####
    ##-----------------------------------##
    
    genoDat = sampDat[sampDat$group %in% c('diploid',geno),]
    #Define default comparison
    genoDat$group = factor(genoDat$group,levels = c('diploid',geno))
    
    #Subset seurat object to keep only the relevant genotypes
    #srat.geno = subset(big.srat,subset = cellID %in% big.srat$cellID[big.srat$Genotype %in% c('diploid',geno)])
    
    out = list()
    ## Loop over cell types
    for(tgtIdx in c(1:n_distinct(big.srat@meta.data$finalAnn[big.srat$Genotype == geno]))){
      
      tgtCell = unique(big.srat@meta.data$finalAnn[big.srat$Genotype==geno])[tgtIdx]
      
      if(tgtCell %in% c('?','unknown','others','VCAM1+.EI.macrophage','doublets','lowQual','Trophoblast','Neuron')){
        next
      }
      
      #if(tgtCell %in% c('ME','Fibroblast') & geno == "complete_trisomy"){next}
      if(tgtCell %in% c('pre.B.cell') & geno == "T18"){next}
      if(tgtCell %in% c('pre.B.cell','Mesothelial_cells') & geno == "MX"){next}
      if(tgtCell %in% c('Mesenchyme') & geno == "T21"){next}
      
      message(sprintf("\n\n------- Consider cell type %s from geno %s tissue %s",tgtCell,geno,tissue))
      #srat = subset(srat.geno,subset = finalAnn == tgtCell)
      srat = subset(big.srat,subset = cellID %in% big.srat$cellID[big.srat$Genotype %in% c('diploid',geno) & big.srat$finalAnn == tgtCell])
      
      
      nCellsGroup = table(factor(srat@meta.data$Genotype))
      
      if((!all(nCellsGroup>=50)) | length(nCellsGroup) <2){
        message(sprintf('Low number of cells detected'))
        print(nCellsGroup)
        next
      }
      
      #Check how many from individual donors
      nCells = table(srat@meta.data$donorID)
      if(sum(nCells>30)<3){
        message(sprintf("Too few effective replicates.  Skipping..."))
        print(nCells)
        next
      }
      message("Found the following number of high confidence cells")
      print(nCells)
      
      # remove individuals with < 30 cells
      donorID_toKeep = names(nCells[nCells >= 30])
      donorID_toRemove = unique(big.srat$donorID[big.srat$Genotype == 'diploid' & big.srat$published_ann_2 != 'NA'])
      donorID_toKeep = donorID_toKeep[!donorID_toKeep %in% donorID_toRemove]
      
      # Downsample ME
      if(tgtCell %in% c('ME','LE')){
        message('Downsampling LE/ME')
        max_nCell = max(nCells[names(nCells) %in% donorID_toKeep])
        min_nCell = min(nCells[names(nCells) %in% donorID_toKeep])
        if(min_nCell <= 5000){
          min_nCell = 5000 # Set a lower cap - do not downsample cells to below this min_nCell value
        }
        if(max_nCell > 8000){
          cellID_toKeep = c()
          # Downsample to min_nCell
          set.seed(1237)
          for(d in donorID_toKeep){
            if(nCells[d] > min_nCell){
              cellID = sample(srat$cellID[srat$donorID == d],min_nCell)
            }else{
              cellID = srat$cellID[srat$donorID == d]
            }
            
            cellID_toKeep = c(cellID_toKeep,cellID)
          }
        }
      }
      
      #OK, we're going ahead, create the needed objects
      if(tgtCell %in% c('ME','LE')){
        toc = srat@assays$RNA@counts[,colnames(srat@assays$RNA@counts) %in% cellID_toKeep]
      }else{
        toc = srat@assays$RNA@counts[,colnames(srat@assays$RNA@counts) %in% srat$cellID[srat$donorID %in% donorID_toKeep]]
      }
      
      toc = toc[rownames(toc) %in% geneMap$geneSym,]
      m = match(rownames(toc),geneMap$geneSym)
      sum(is.na(m))
      rownames(toc) = geneMap$ensID[m]
      mDat = data.frame(row.names = colnames(toc),
                        cellID = colnames(toc),
                        donor = srat@meta.data$donorID[match(colnames(toc),srat$cellID)])
      mDat = merge(mDat,genoDat,by='donor')
      mDat = mDat[match(colnames(toc),mDat$cellID),]
      rownames(mDat) = colnames(toc)
      
      # check that rownames(mDat) is in correct order
      if(!all(rownames(mDat) == mDat$cellID)){
        stop(sprintf('Incorrect mDat cellID order for tissue %s cell type %s genotype %s. Please check!',tissue, tgtCell, geno))
      }
      # Remove cellID column from mDat
      mDat = mDat[,colnames(mDat) != 'cellID']
      
      coords = gns[rownames(toc)]
      
      if(n_distinct(mDat$tech) == 2){
        formula_toUse = paste0(formula,' + tech')
      }
      print(formula_toUse)
      
      
      ##------    from here, instead of using the DESeq2 function, we replace it with edgeR pipeline    ---------####
      formula = '~ %s'
      donorID = 'donorID'
      group = 'Genotype'
      
      ## Drop irrelevant genes
      #Uninteresting genes
      w = which(as.character(seqnames(coords)) %in% tgtChrs)
      coords = coords[w]
      seqlevels(coords) = tgtChrs
      toc = toc[coords$gene_id,]
      #Non-expressed genes
      w = which(rowSums(toc>0)>=3)
      coords = coords[w]
      toc = toc[w,]
      
      ##=========================##
      # Get genomic coordinates
      #Order aa and bb by genomic coords
      #If it's positive stranded or unknown use start, else end
      coords$TSS = ifelse(strand(coords)=='-',end(coords),start(coords))
      o = order(seqnames(coords),coords$TSS)
      coords = coords[o]
      toc = toc[o,]
      
      ##=========================##
      # Create pseudobulk 
      pb = do.call(cbind,lapply(split(colnames(toc),mDat[,donorID]),function(e) rowSums(toc[,e,drop=FALSE])))
      pb_byCT_byAgeGroup[[paste0(tgtCell,'_',ageGroup)]] = pb
      
      colDat = mDat[match(colnames(pb),mDat[,donorID]),]
      colDat$Genotype[colDat$Genotype == '2n'] = 'diploid'
      colDat$Genotype = factor(colDat$Genotype,c('diploid','T21'))
      rownames(colDat) = colDat[,donorID]
      
      
      ##=========================##
      # Fit EDGER model
      pdf(file.path(outDir,paste0(tgtCell,'_',ageGroup,'_pbEdgeR_plots.pdf')))
      out[[paste0(tgtCell,'_',ageGroup)]] = fit_model(pb,colDat,formula = '~ 0 + %s',geneMap,groupID='group')
      dev.off()
      
      
      
      
      
      
      
      
      
      
      ## Add number of cells going into each comparison
      tmp = data.frame(nCell_2n = length(srat$cellID[srat$donorID %in% donorID_toKeep & srat$cellID %in% colnames(toc) & srat$Genotype == 'diploid']),
                       nCell_AK = length(srat$cellID[srat$donorID %in% donorID_toKeep & srat$cellID %in% colnames(toc) & srat$Genotype == geno]),
                       ct = tgtCell,
                       geno = geno)
      
      nCell_perCT_perGeno = rbind(nCell_perCT_perGeno,tmp)
    } # Loop over cell types
    
    saveRDS(out,paste0(outDir_fp,'/',geno,'_',tissue,'_pseudoBulkDESeq2out.RDS'))
  } # Loop over genotypes
  
  
  write.csv(nCell_perCT_perGeno,file.path(outDir,'liver_nCell_perCT_perGeno_used_pbDEGs.csv'))
  
}


print('Completed!')
