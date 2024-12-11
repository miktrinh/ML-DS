### DGE analysis by celltype by genotype ###

#outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withoutCyclCells/oct23_ageMatched'
outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jul24/without_published2n/genoAssay'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)



##-------------------##
##   Libraries     ####
##-------------------##

library(tidyverse)
library(Seurat)
library(GenomicFeatures)
library(ComplexHeatmap)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")


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



srat_inDir = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0724.RDS'

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
  
  
  big.srat = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0724.RDS')
  big.srat$Sex[big.srat$Sex == 'female'] = 'F'
  big.srat$Sex[big.srat$Sex == 'male'] = 'M'
  big.srat$finalAnn[big.srat$finalAnn %in% c('earlyMK')] = 'MK'
  big.srat$finalAnn[big.srat$finalAnn %in% c('promyelocyte','myelocyte')] = 'Myelocyte'
  big.srat$finalAnn[big.srat$finalAnn %in% c('proMono')] = 'Monocyte'
  big.srat$finalAnn[big.srat$finalAnn %in% c('LMPP_ELP','pro.B.cell','pre.B.cell')] = 'B.cell.prog'
  big.srat$finalAnn[big.srat$finalAnn %in% c('ILC.precursor','T.cell','NK_T')] = 'NK.T'
  big.srat$ageGroup = big.srat$gestationalAge
  
  # big.srat = filter_sratObj(tissue = tissue,ageMatch=ageMatch,remove_cyclingCells=remove_cyclingCells,
  #                           srat_in_fp='~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0724.RDS')
  # big.srat = big.srat[[1]]
  
  
  ### To get a decent control - we need at least 3 diploid-donors to have at least 3 values for the boxplots
  #sum_CT = sum_CT[(sum_CT$finalAnn %in% sum_CT$finalAnn[sum_CT$Genotype=='diploid' & sum_CT$nDonor_with100cells >=3]) &
  #                  !(sum_CT$finalAnn %in% c('Other','Leukocytes','Erythroblasts')),]
  
  #ct_toRemove = sum_CT[,c('finalAnn','Genotype','nDonor_with100cells')] %>% 
  #  mutate(type = ifelse(Genotype == 'diploid','diploid','triploid')) %>% group_by(finalAnn, type) %>% 
  #  summarise(nDonor = ifelse(type == 'diploid',(unique(nDonor_with100cells)>=3),sum(nDonor_with100cells>=1))) %>% 
  #  distinct() %>% group_by(finalAnn) %>% summarise(toRemove = sum(nDonor ==0)>0)
  
  #sum_CT = sum_CT[sum_CT$finalAnn %in% ct_toRemove$finalAnn[ct_toRemove$toRemove==F],]
  
  #ct_toKeep = unique(sum_CT$finalAnn)
  
  
  ### Subset to keep cells from relevant tissue only (both AK and REF)
  
  #big.srat = subset(big.srat, subset = finalAnn %in% ct_toKeep)
  
  ###############
  # sub-sampling #
  ###############
  #### Perform sub-sampling #####
  # Perform sub-sampling for each donor to reach tgtCells
  # n_tgtCells = 100
  # all_cellIDs_toKeep = c()
  # for(ct in unique(big.srat$finalAnn)){
  #   
  #   
  #   nCells = table(big.srat@meta.data$donorID[big.srat@meta.data$finalAnn == ct])
  #   # Remove donors with fewer cells than tgtCells
  #   donor_toKeep = names(nCells[nCells>as.numeric(n_tgtCells)])
  #   if(length(donor_toKeep) <=1){
  #     next
  #   }
  #   
  #   ## For liver, look at the proportion of cycling cells --> subset cycling cells to match the lowest perc (in AK presumably)
  #   cyclingCells = big.srat@meta.data[big.srat@meta.data$finalAnn == ct & big.srat@meta.data$donorID %in% donor_toKeep,] %>% 
  #     group_by(Genotype,donorID) %>% 
  #     summarise(n_G1 = sum(Phase == 'G1'),
  #               n_G2M = sum(Phase == 'G2M'),
  #               n_S = sum(Phase == 'S'),
  #               n_cell = n(),
  #               frac_G2M = n_G2M/n_cell,
  #               frac_S = n_S/n_cell,
  #               n_G2M_tokeep = ifelse(tissue_toKeep == 'liver',round(min(frac_G2M) * n_tgtCells),0),
  #               n_S_tokeep = ifelse(tissue_toKeep == 'liver',round(min(frac_S) * n_tgtCells),0),
  #               n_G1_tokeep = n_tgtCells - n_G2M_tokeep - n_S_tokeep)
  #   
  #   
  #   set.seed(1230)
  #   mdat = big.srat@meta.data[big.srat@meta.data$finalAnn == ct & big.srat@meta.data$donorID %in% donor_toKeep,]
  #   mdat$n_tgtCells_toKeep = 'NA'
  #   mdat$n_tgtCells_toKeep[mdat$Phase=='G1'] = cyclingCells$n_G1_tokeep[match(mdat$donorID[mdat$Phase=='G1'],cyclingCells$donorID)]
  #   mdat$n_tgtCells_toKeep[mdat$Phase=='G2M'] = cyclingCells$n_G2M_tokeep[match(mdat$donorID[mdat$Phase=='G2M'],cyclingCells$donorID)]
  #   mdat$n_tgtCells_toKeep[mdat$Phase=='S'] = cyclingCells$n_S_tokeep[match(mdat$donorID[mdat$Phase=='S'],cyclingCells$donorID)]
  #   
  #   cells_toKeep = mdat %>% group_by(donorID,Phase) %>% 
  #     mutate(id = seq(1,n()),
  #            selected = ifelse(id %in% sample(1:max(id),unique(n_tgtCells_toKeep)),T,F))  
  #   
  #   
  #   cells_toKeep = cells_toKeep$cellID[cells_toKeep$selected]
  #   all_cellIDs_toKeep = unique(c(all_cellIDs_toKeep,cells_toKeep))
  # }
  # 
  # big.srat.sub = subset(big.srat, subset = cellID %in% all_cellIDs_toKeep)
  
  
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
  unique(big.srat$Genotype[!is.na(big.srat$Genotype)])
  for(geno in c('T21')){
    out=list()
    if(geno %in% c('diploid'))
      next
    print(geno)
    genoDat = sampDat[sampDat$group %in% c('diploid',geno),]
    #Define default comparison
    genoDat$group = factor(genoDat$group,levels = c('diploid',geno))
    
    #Subset seurat object to keep only the relevant genotypes 
    srat.geno = subset(big.srat,subset = Genotype %in% c('diploid',geno))
    
    # Subset for ageMatch
    # if(ageMatch){
    #   ageGroup_toKeep = (rowSums(table(srat.geno$ageGroup,srat.geno$Genotype)>0) == 2)
    #   ageGroup_toKeep = names(ageGroup_toKeep[ageGroup_toKeep==T])
    #   if(geno == 'T22'){
    #     ageGroup_toKeep = c('15pcw','16pcw','14pcw')
    #   }
    #   srat.geno = subset(srat.geno,subset = ageGroup %in% ageGroup_toKeep)
    # }
    
    ## Loop over cell types
    for(tgtIdx in c(1:n_distinct(srat.geno@meta.data$finalAnn[srat.geno$Genotype!='diploid']))){
      
      tgtCell = unique(srat.geno@meta.data$finalAnn[srat.geno$Genotype!='diploid'])[tgtIdx]  
      if(tgtCell %in% c('?','unknown','others','VCAM1+.EI.macrophage')){
        next
      }
      
      if(tgtCell %in% c('ME','Fibroblast') & geno == "complete_trisomy"){next}
      if(tgtCell %in% c('pre.B.cell') & geno == "T18"){next}
      if(tgtCell %in% c('Fibroblast','pre.B.cell') & geno == "MX"){next}
      message(sprintf("\n\n------- Consider cell type %s from geno %s tissue %s",tgtCell,geno,tissue))
      srat = subset(srat.geno,subset = finalAnn == tgtCell)
      #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
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
      
      #OK, we're going ahead, create the needed objects
      toc = srat@assays$RNA@counts[,colnames(srat@assays$RNA@counts) %in% srat$cellID[srat$donorID %in% donorID_toKeep]]
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
      
      #formula = ifelse(geno == 'T22','~ %s + sex')
      formula = '~ %s + sex'
      # if(n_distinct(mDat$sex) == 1){
      #   formula = gsub('+ sex','',formula)
      # }
      # ak_ageGroup = unique(mDat$ageGroup[mDat$group != 'diploid'])
      # dip_ageGroup = unique(mDat$ageGroup[mDat$group == 'diploid'])
      # if(n_distinct(mDat$ageGroup) == 1 | sum(dip_ageGroup %in% ak_ageGroup) == 0){
      #   formula = gsub('+ ageGroup','',formula)
      # }
      # 
      pdf(file.path(plotDir,paste0(geno,'_',tissue,'___unmerge___',tgtCell,'.pdf')))
      out[[paste0(tgtCell,'_unmerged')]] = compareCell(toc = toc,
                                                       mDat = mDat[match(colnames(toc),rownames(mDat)),],
                                                       coords = gns[rownames(toc)],
                                                       cellTypeName=tgtCell,
                                                       formula=formula,
                                                       donorID='donor',groupID='group')
      dev.off()
      
      pdf(file.path(plotDir,paste0(geno,'_',tissue,'___merge',kMerge,'___',tgtCell,'.pdf')))
      out[[paste0(tgtCell,'_merge',kMerge)]] = compareCell(toc,mDat,coords,
                                                          kMerge=kMerge,
                                                          formula = formula,
                                                          cellTypeName=tgtCell)
      dev.off()
    }
    
    
    
    nCell_perCT_perGeno = data.frame()
    #--------------------------------------
    for(geno in unique(big.srat$Genotype[!is.na(big.srat$Genotype)])){
      if(geno %in% c('diploid'))
        next
      print(geno)
      genoDat = sampDat[sampDat$group %in% c('diploid',geno),]
      #Define default comparison
      genoDat$group = factor(genoDat$group,levels = c('diploid',geno))
      
      #Subset seurat object to keep only the relevant genotypes 
      srat.geno = subset(big.srat,subset = Genotype %in% c('diploid',geno))
      
      ## Loop over cell types
      for(tgtIdx in c(1:(length(unique(srat.geno@meta.data$finalAnn[srat.geno$Genotype!='diploid']))))){
        
        tgtCell = unique(srat.geno@meta.data$finalAnn[srat.geno$Genotype!='diploid'])[tgtIdx]  
        if(tgtCell %in% c('?','unknown','others','VCAM1+.EI.macrophage')){
          next
        }
        
        if(tgtCell %in% c('ME','Fibroblast') & geno == "complete_trisomy"){next}
        if(tgtCell %in% c('pre.B.cell') & geno == "T18"){next}
        if(tgtCell %in% c('Fibroblast','pre.B.cell') & geno == "MX"){next}
        message(sprintf("\n\n------- Consider cell type %s from geno %s tissue %s",tgtCell,geno,tissue))
        srat = subset(srat.geno,subset = finalAnn == tgtCell)
        
        #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
        nCellsGroup = table(factor(srat@meta.data$Genotype))
        
        if((!all(nCellsGroup>=50)) | length(nCellsGroup) <2){
          message(sprintf('Low number of cells detected'))
          print(nCellsGroup)
          next
        }
        
        #Check how many from individual donors
        nCells = table(srat@meta.data$donorID)
        if(sum(nCells>20)<3){
          message(sprintf("Too few effective replicates.  Skipping..."))
          print(nCells)
          next
        }
        message("Found the following number of high confidence cells")
        print(nCells)
        
        # remove individuals with < 30 cells
        donorID_toKeep = names(nCells[nCells >= 30])
        
        #OK, we're going ahead, create the needed objects
        tmp = data.frame(nCell_2n = length(srat$cellID[srat$donorID %in% donorID_toKeep & srat$Genotype == 'diploid']),
                         nCell_AK = length(srat$cellID[srat$donorID %in% donorID_toKeep & srat$Genotype == geno]),
                         ct = tgtCell,
                         geno = geno)
        
        nCell_perCT_perGeno = rbind(nCell_perCT_perGeno,tmp)
      }
    }
    
    
    write.csv(nCell_perCT_perGeno,file.path(outDir,'nCell_perCT_perGeno_used_pbDEGs.csv'))
    
    
    
    
    
    
    
    
    saveRDS(out,paste0(outDir_fp,'/',geno,'_',tissue,'_pseudoBulkDESeq2out.RDS'))
    
    
    #############################################
    # Plot median log2FC across celltypes       #
    #############################################
    # Get the median log2FC per chromosome for each cell type --> do a boxplot (x = chromosome) + dots coloured by celltypes
    log2FC = tibble()
    #,paste0('merge',kMerge)
    for(dat in c('unmerged')){
      #for(dat in c('unmerged')){
      results = out[grepl(dat,names(out))]
      celltypes = gsub('_unmerged$','',names(results))
      for(idx in 1:length(results)){
        celltype = celltypes[idx]
        #celltype = sapply(strsplit(names(results)[idx],split='_'),'[',1)
        print(celltype)
        log2FC_perChr = results[[idx]][['log2FC_perChr']]
        
        #boxplot(log2FC_perChr,
        #        outline=FALSE,
        #        xlab='Chromosome',
        #        ylab='logFC'
        #)
        #abline(h=0,col='red')
        
        med = sapply(log2FC_perChr,FUN = function(x){median(x)})
        mu = sapply(log2FC_perChr,FUN = function(x){mean(x)})
        #plot(x=names(med),y=med,pch=19)
        #abline(h=0,col='red')
        
        tmp = data.frame(celltype = celltype, med_log2FC = med,mu_log2FC=mu, chr = names(med),type=dat)
        log2FC = rbind(log2FC,tmp)
        
      }
    }
    
    
    log2FC$chr = factor(log2FC$chr,levels = c(1:22,'X'))
    pdf(file.path(outDir_fp,paste0(geno,'_',tissue,'___',dat,'___medianLog2FCperChr.pdf')),width = 13,height = 10)
    p = ggplot(log2FC,aes(x=chr, y=med_log2FC))+
      geom_boxplot()+
      geom_point(aes(col = celltype))+
      scale_color_manual(values = c(col25,brewer.pal(12,'Paired')))+
      geom_hline(yintercept = 0)+
      theme_bw(base_size = 13)+
      facet_wrap(vars(type),ncol=1,scales = 'free')+
      xlab('Chromosome')+ylab('median log2FC')+
      ggtitle(paste0(geno,'_',tissue))
    
    print(p)
    dev.off()
    
    write.table(log2FC,file.path(outDir_fp,paste0(geno,'_',tissue,'___',dat,'___medianLog2FCperChr.csv')),sep = ',',row.names = F,col.names = T)
    
    
    
    
    
    
    ######################################################
    # Plot Number of DEGs per Chr across celltypes       #
    ######################################################
    # Get theNumber of DEGs per chromosome for each cell type --> do a boxplot (x = chromosome) + dots coloured by celltypes
    #,paste0('merge',kMerge)
    for(dat in c('unmerged')){
      celltypes = gsub('_unmerged$','',names(results))
      results = out[grepl(dat,names(out))]
      
      for(idx in 1:length(results)){
        #celltype = sapply(strsplit(names(results)[idx],split='_'),'[',1)
        celltype = celltypes[idx]
        print(celltype)
        n_deg_perChr = results[[idx]][['deTable']]
        n_deg_perChr$celltype = celltype
        n_deg_perChr$geno = geno
        n_deg_perChr$type = dat
        
        n_deg_perChr_allCT = rbind(n_deg_perChr_allCT,n_deg_perChr)
      }
    }
  }
  
  
  write_csv(n_deg_perChr_allCT,file.path(outDir_fp,'nDEGs_perChr_allCT_allGeno.csv'))
  # df = n_deg_perChr_allCT %>% group_by(geno,celltype) %>% summarise(up = sum(nUp),down=-sum(nDown),nExp = sum(nExpressed))
  # df2 = pivot_longer(df,cols = c(3,4),names_to = 'type',values_to = 'nDE')
  # ggplot(df2,aes(reorder(celltype,nDE/nExp),nDE/nExp,fill=type))+geom_col() + 
  #   theme_bw()+
  #   scale_fill_manual(values = col25)+
  #   facet_wrap(vars(geno))+
  #   theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + xlab('')
  
}

  
  




# ### pseudobulk between 3p and 5p
# 
# #OK, we're going ahead, create the needed objects
# toc = big.srat@assays$RNA@counts
# toc = toc[rownames(toc) %in% geneMap$geneSym,]
# m = match(rownames(toc),geneMap$geneSym)
# sum(is.na(m))
# rownames(toc) = geneMap$ensID[m]
# mDat = data.frame(row.names = colnames(toc),
#                   cellID = colnames(toc),
#                   donor = big.srat@meta.data$donorID,
#                   tech= big.srat$assay,
#                   sex = big.srat$Sex)
# mDat = mDat[match(colnames(toc),mDat$cellID),]
# rownames(mDat) = colnames(toc)
# mDat$tech[mDat$tech == "10X 3'v2"] = 'v2'
# mDat$tech[mDat$tech == "scRNA 3' v3.1"] = 'v1'
# 
# # check that rownames(mDat) is in correct order
# if(!all(rownames(mDat) == mDat$cellID)){
#   stop(sprintf('Incorrect mDat cellID order for tissue %s cell type %s genotype %s. Please check!',tissue, tgtCell, geno))
# }
# # Remove cellID column from mDat
# mDat = mDat[,colnames(mDat) != 'cellID']
# table(mDat$donor,mDat$tech)
# coords = gns[rownames(toc)]
# 
# formula = '~ %s + sex'
# 
# pdf(file.path(plotDir,paste0(geno,'_',tissue,'___unmerge___',tgtCell,'.pdf')))
# batch_comp = compareCell(toc = toc,
#                          mDat = mDat[match(colnames(toc),rownames(mDat)),],
#                          coords = gns[rownames(toc)],
#                          cellTypeName='assay',
#                          formula=formula,
#                          donorID='donor',groupID='tech')
# dev.off()
# 
# batch_comp=list(res=res,log2FC_perChr=log2FC_perChr,de=de,deTable=deTab,mergeMap=mergeMap,dds=dds,rldCor = rldCor)
# 
# batch_deg=batch_comp[['de']]



















##------------------------------##
##     Import the results     ####
##------------------------------##

# Import dataset
tissue = 'liver'
# fLiver = filter_sratObj(tissue,ageMatch = T,remove_cyclingCells = F)
# fLiver = fLiver[[1]]
fLiver = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver_2n/jan24/liver_liverREFmerged_clean_processed_annotated_noUnknowns_noPublished2n_0124.RDS')


# Import DEG output (formula = . ~ geno)
resultDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/apr24/without_published2n/geno/'
allGenes = import_pbDEGresults(outDir = resultDir,tissue = 'liver',allGenes = T)
allGenes$geneSym = geneMap$geneSym[match(allGenes$ensID,geneMap$ensID)]

## Import significant DEGs
allDEGs = import_pbDEGresults(outDir = resultDir,tissue = 'liver')
head(allDEGs)
table(allDEGs$geno,allDEGs$ct)
dim(allDEGs[is.na(allDEGs$geneSym),])

#allDEGs.sub$ct_geno = gsub('complete$','Triploid',allDEGs.sub$ct_geno)
#allDEGs.sub$geno[allDEGs.sub$geno == 'complete'] = 'Triploid'

allDEGs.sub = allDEGs

View(allDEGs.sub[allDEGs.sub$ct == 'MEMP_MEP' & allDEGs.sub$geno == 'T21',])
View(allDEGs.sub[allDEGs.sub$ct == 'HSC_MPP' & allDEGs.sub$geno == 'T21',])


Idents(fLiver) = paste0(fLiver$annot_jan24,':',fLiver$Genotype)
DotPlot(fLiver,idents = unique(Idents(fLiver)[grepl('HSC|MEMP',Idents(fLiver)) ]),
        #group.by = 'annot_jan24',
        scale = T,
        #features = allDEGs.sub[allDEGs.sub$ct == 'Kupffer.cell' & allDEGs.sub$geno == 'T21' & allDEGs.sub$log2FoldChange > 0,]$geneSym
        features = rownames(fLiver)[grepl('^IFN',rownames(fLiver))]
        #features = mlds_tam_shared_up$geneSym[mlds_tam_shared_up$avgExpr >= 0.5][31:60]
        #features = df$geneSym[df$module == 'MLDS.TAM.up'][1:60]
        
)+RotatedAxis()

DotPlot(mlds,idents = unique(Idents(mlds)[grepl('MLDS:|TAM:|HSC_MPP|Mono_CD14|MEP|MK|EE|CMP|B.cell',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT',Idents(mlds))]),
        #group.by = 'donorID',
        scale = T,
        features = c(allDEGs.sub[allDEGs.sub$ct == 'MEMP_MEP' & allDEGs.sub$geno == 'T21' & allDEGs.sub$log2FoldChange > 0,]$geneSym,'IFNAR2','IFNGR1','IRF2','IRF6','STAT1','RUNX1','TNFRSF1A','TNFRSF1B')
        #features = geneList[[2]]
        #features = mlds_tam_shared_up$geneSym[mlds_tam_shared_up$avgExpr >= 0.5][31:60]
        #features = df$geneSym[df$module == 'MLDS.TAM.up'][1:60]
        
)+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11)) + xlab('') + ylab('')




## Perform enrichR on all DEGs
upDEG_enriched <- enrichr(unique(allDEGs.sub$geneSym[allDEGs.sub$geno == 'T21' & allDEGs.sub$direction == 'AK_up']), dbs)
downDEG_enriched <- enrichr(unique(allDEGs.sub$geneSym[allDEGs.sub$geno == 'T18' & allDEGs.sub$direction == 'AK_down']), dbs)

View(upDEG_enriched[[6]][1,])
downDEG_enriched[[5]][6,]
plotEnrich(upDEG_enriched[[9]], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
plotEnrich(downDEG_enriched[[6]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")

figureSupp2_enrichR_DEG_AKfLiver = function(){
  
  if(!file.exists()){
    # Import DEG output (formula = . ~ geno)
    resultDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/'
    allDEGs = import_pbDEGresults(outDir = resultDir,tissue = 'liver')
    
    ## Save as supplementary table
    allDEGs2 = allDEGs[,c('geno','ct','contrast','ensID','geneSym','chr','baseMean','log2FoldChange','lfcSE','pvalue','padj','isTF','isCSM','isCosmic','cosmicTier', 'tumourType', 'isTSG', 'cellFrac_overall', 'cellFrac_g1', 'cellFrac_g2')]
    colnames(allDEGs2)[1:2] = c('genotype','cell_type')
    allDEGs2$contrast = gsub('group_','',allDEGs2$contrast)
    allDEGs2$contrast = gsub('complete_trisomy','triploid',allDEGs2$contrast)
    allDEGs2$genotype = gsub('3n','triploid',allDEGs2$genotype)
    rownames(allDEGs2) = gsub('^complete','triploid',rownames(allDEGs2))
    write_delim(allDEGs2,'~/lustre_mt22/Aneuploidy/manuscriptDraft_0124/TableS2_AK.fLiver_pbDESeq2_sigDEGs.tsv',delim = '\t')
  }else{
    allDEGs = read.delim('~/lustre_mt22/Aneuploidy/manuscriptDraft_0124/TableS2_AK.fLiver_pbDESeq2_sigDEGs.tsv',sep = '\t')
  }
  
  ## Perform enrichR on all DEGs
  upDEG_enriched <- enrichr(unique(allDEGs$geneSym[allDEGs$geno == 'T21' & allDEGs$direction == 'AK_up']), dbs)
  downDEG_enriched <- enrichr(unique(allDEGs$geneSym[allDEGs$geno == 'T21' & allDEGs$direction == 'AK_down']), dbs)
  
  # plotEnrich(upDEG_enriched[[1]], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
  # plotEnrich(downDEG_enriched[[9]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
  
  ## Plot up-MSigDB_Hallmark_2020
  
  plotFun_enrichR_ak.fLiver_DEGs = function(noFrame=FALSE,noPlot=FALSE){
    allTerms = data.frame()
    for(i in c(5,6)){
      df = upDEG_enriched[[i]]
      df = df[df$Adjusted.P.value < 0.01,]
      df = df[order(df$Combined.Score,decreasing = F),]
      df$nGene = as.numeric(gsub('/.*$','',df$Overlap))
      df$overlap2 = as.numeric(gsub('/.*$','',df$Overlap)) / as.numeric(gsub('^.*/','',df$Overlap))
      # only keep terms with at least 5% overlap with gene list
      df=df[df$overlap2 > 0.05,]
      df$db = names(upDEG_enriched)[i]
      if(nrow(df) >= 15){
        df = df[1:15,]
      }
      df$yStart = seq(1:nrow(df))
      allTerms = rbind(allTerms,df)
      
      
      p1 = ggplot(df,aes(Combined.Score,reorder(Term,Combined.Score)))+
        geom_point(aes(size = nGene,col=overlap2))+
        geom_segment(aes(x=0,xend=rev(Combined.Score),y=1:n_distinct(df$Term),yend=1:n_distinct(df$Term)),col=grey(0.7))+
        xlab('Combined Score')+ylab('') + ggtitle(names(upDEG_enriched)[i])+
        scale_color_gradient(low='#eba9a9',high = '#a10505')+
        theme_classic(base_size = 11)+
        theme(panel.border = element_rect(fill=F),axis.line = element_blank())
      #print(p1)
    }
    
    allTerms$Term = factor(allTerms$Term,levels = allTerms$Term)
    p3=ggplot(allTerms,aes(Combined.Score,Term))+
      geom_point(aes(size = nGene,col=overlap2))+
      facet_grid(db~.,scales = 'free_y',space = 'free_y')+
      geom_segment(aes(x=0,xend=(Combined.Score),y=yStart,yend=yStart),col=grey(0.7))+
      xlab('Combined Score')+ylab('') + ggtitle(names(upDEG_enriched)[i])+
      scale_color_gradient(low='#eba9a9',high = '#a10505')+
      theme_classic(base_size = 11)+
      theme(panel.border = element_rect(fill=F),axis.line = element_blank())
    print(p3)
  }
  
  saveFig(file.path(plotDir,'Supp.Fig2E_enrichR_AKfLiver_DESeq2_DEGs'),plotFun_enrichR_ak.fLiver_DEGs,rawData=allTerms,width = 7,height = 7,res = 500,useDingbats = T)
  
}




##---------------------------------------------------##
##     Explore the list of DEGs in T21             ####
##---------------------------------------------------##
View(allDEGs)

fLiver$group = paste0(fLiver$annot_mar24,'.',fLiver$Genotype)
#fLiver$group[grepl('Kupffer',fLiver$group)] = gsub('Kupffer.cell','Macrophage',fLiver$group[grepl('Kupffer',fLiver$group)])
avgExpr = AverageExpression(fLiver,group.by = 'group')
avgExpr = avgExpr[['RNA']]

## Interferon genes
genes1 = c('IFNG','IL1B','IL18')
genes2 = c('IFNAR1','IFNAR2','IFNGR1','IFNGR2')

#mtx = avgExpr[grepl('IFN|^IL|^CCL',rownames(avgExpr)),grepl('HSC|MEMP|Mono|Mac|Kup',colnames(avgExpr))]
#mtx = avgExpr[rownames(avgExpr) %in% c(genes1,genes2),grepl('^Mono|Mac|Kup|NK_T|HSC|MEMP',colnames(avgExpr))]
#mtx = mtx[apply(mtx,1,function(x){max(x) > 0.01}),]
# mtx = mtx[c(genes1,genes2),c(colnames(mtx)[grepl('diploid',colnames(mtx))],
#                              colnames(mtx)[grepl('T21',colnames(mtx))],
#                              colnames(mtx)[grepl('T18',colnames(mtx))],
#                              colnames(mtx)[grepl('T22',colnames(mtx))],
#                              colnames(mtx)[grepl('complete_trisomy',colnames(mtx))],
#                              colnames(mtx)[grepl('MX',colnames(mtx))])]
# dim(mtx)
# mtx1 = mtx[,!grepl('HSC|MEMP',colnames(mtx))]
# col_fun = circlize::colorRamp2(c(-1, 0, 3), c(grey(0.99),grey(0.8),grey(0.1)))
# hm = Heatmap(t(scale(t(mtx1))),show_column_dend = F,show_row_dend = F,name='zScaled expression',
#              cluster_rows = F,cluster_columns = F,
#              col = col_fun,
#              column_split=fLiver$annot_jan24[fLiver$annot_jan24 != 'Kupffer.cell'][match(colnames(mtx1),fLiver$group[fLiver$annot_jan24 != 'Kupffer.cell'])])
# draw(hm)
# 
# mtx2 = mtx[rownames(mtx) %in% genes2,grepl('HSC|MEMP',colnames(mtx))]
# Heatmap(t(scale(t(mtx2))),show_column_dend = F,show_row_dend = F,name='zScaled expression',
#         cluster_rows = F,cluster_columns = F,
#         col = col_fun,
#         column_split=fLiver$annot_jan24[match(colnames(mtx2),fLiver$group)])



## Plot log2FC instead of avgExpr
d = allGenes[allGenes$geneSym %in% c(genes1,genes2) & grepl('^Mono|Mac|Kup|NK.T|HSC|MEMP',allGenes$ct),]
d$group = paste0(d$ct,'.',d$geno)
d$group[grepl('NK.T',d$group)] = gsub('NK.T','NK_T',d$group[grepl('NK.T',d$group)])
mtx = pivot_wider(d,id_cols = c('geneSym'),names_from = 'group',values_from = 'log2FoldChange')
mtx = column_to_rownames(mtx,var = 'geneSym')
mtx = mtx[c(genes1,genes2),c(colnames(mtx)[grepl('diploid',colnames(mtx))],
                             colnames(mtx)[grepl('T21',colnames(mtx))],
                             colnames(mtx)[grepl('T18',colnames(mtx))],
                             colnames(mtx)[grepl('T22',colnames(mtx))],
                             colnames(mtx)[grepl('complete_trisomy',colnames(mtx))],
                             colnames(mtx)[grepl('MX',colnames(mtx))])]
dim(mtx)
mtx1 = mtx[,!grepl('HSC|MEMP',colnames(mtx))]
col_fun = circlize::colorRamp2(c(-0.5, 0, 2.5), c(grey(0.99),grey(0.8),grey(0.1)))
col_fun = circlize::colorRamp2(c(-1, 0, 1), c('#1a4a87','white','#a4282c'))

hm = Heatmap(t(scale(t(mtx1))),show_column_dend = F,show_row_dend = F,name='zScaled expression',
             cluster_rows = F,cluster_columns = F,
             col = col_fun,
             column_split=fLiver$annot_jan24[match(colnames(mtx1),fLiver$group)])
draw(hm)

mtx2 = mtx[rownames(mtx) %in% genes2,grepl('HSC|MEMP',colnames(mtx))]
Heatmap(t(scale(t(mtx2))),show_column_dend = F,show_row_dend = F,name='zScaled expression',
        cluster_rows = F,cluster_columns = F,
        col = col_fun,
        column_split=fLiver$annot_jan24[match(colnames(mtx2),fLiver$group)])




# ##--------------------------------------------------------------##
# ##      Comparing the output of geno vs geno+sex              ####
# ##--------------------------------------------------------------##
# ## Conclusion: Majority of the DEGs are not sex-related. Even then, sex is confounded with gestational age... so it could be an age thing...
# ## We will go forward with just "geno"
# out_wSex = import_pbDEGresults(outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno_sex',tissue = 'liver')
# head(out_wSex)
# dim(out_wSex[is.na(out_wSex$geneSym),])
# 
# allDEGs$in_wSex = NA
# out_wSex$in_woSex = NA
# allDEGs$groupID = paste0(allDEGs$ct_geno,':',allDEGs$direction)
# out_wSex$groupID = paste0(out_wSex$ct_geno,':',out_wSex$direction)
# for(g in unique(allDEGs$groupID)){
#   #allDEGs$in_wSex[allDEGs$groupID == g] = ifelse(allDEGs$geneSym[allDEGs$groupID == g] %in% out_wSex$geneSym[out_wSex$groupID == g],T,F)
#   out_wSex$in_woSex[out_wSex$groupID == g] = ifelse(out_wSex$geneSym[out_wSex$groupID == g] %in% allDEGs$geneSym[allDEGs$groupID == g],T,F)
# }
# 
# table(is.na(allDEGs$in_wSex))
# table(allDEGs$in_wSex,allDEGs$ct_geno)
# table(out_wSex$in_woSex)
# 
# ## Plot expression of these genes in each donor
# Idents(big.srat) = paste0(big.srat$annot_jan24,':',big.srat$Genotype)
# big.srat$tmp = paste0(big.srat$Genotype,':',big.srat$Sex,':',big.srat$donorID)
# DotPlot(big.srat,idents = c('HSC_MPP:T21','HSC_MPP:diploid'),group.by = 'tmp',
#         features = c(allDEGs$geneSym[allDEGs$groupID == 'HSC_MPP:T21:AK_up' & allDEGs$in_wSex == F]
#                      #allDEGs$geneSym[allDEGs$groupID == 'MEMP_MEP:T21:AK_up' & allDEGs$in_wSex == F]
#         )) + RotatedAxis()
# 
# 
# 
# 
# ##-----------------------------------------------------##
# ##      Comparing DEGs pCut 0.05 vs 0.1              ####
# ##-----------------------------------------------------##
# # Not worth investigating, as changing pCut requires changes in the compareCell function, to calculate max.pct.cells
# allDEGs_0.1 = import_pbDEGresults(outDir = resultDir,tissue = 'liver',pCut = 0.1,geneMap = geneMap)
# 
# allDEGs_0.1$in_0.05 = NA
# allDEGs$groupID = paste0(allDEGs$ct_geno,':',allDEGs$direction)
# allDEGs_0.1$groupID = paste0(allDEGs_0.1$ct_geno,':',allDEGs_0.1$direction)
# for(g in unique(allDEGs$groupID)){
#   #allDEGs$in_wSex[allDEGs$groupID == g] = ifelse(allDEGs$geneSym[allDEGs$groupID == g] %in% out_wSex$geneSym[out_wSex$groupID == g],T,F)
#   allDEGs_0.1$in_0.05[allDEGs_0.1$groupID == g] = ifelse(allDEGs_0.1$geneSym[allDEGs_0.1$groupID == g] %in% allDEGs$geneSym[allDEGs$groupID == g],T,F)
# }
# 
# table(is.na(allDEGs_0.1$in_0.05))
# table(allDEGs_0.1$in_0.05,allDEGs_0.1$ct_geno)
# table(out_wSex$in_woSex)
# 
# ## Plot expression of these genes in each donor
# Idents(big.srat) = paste0(big.srat$annot_jan24,':',big.srat$Genotype)
# big.srat$tmp = paste0(big.srat$Genotype,':',big.srat$Sex,':',big.srat$donorID)
# DotPlot(big.srat,idents = c('HSC_MPP:T21','HSC_MPP:diploid'),group.by = 'tmp',
#         features = c(allDEGs$geneSym[allDEGs$groupID == 'HSC_MPP:T21:AK_up' & allDEGs$in_wSex == F]
#                      #allDEGs$geneSym[allDEGs$groupID == 'MEMP_MEP:T21:AK_up' & allDEGs$in_wSex == F]
#         )) + RotatedAxis()

##-----------------------------------------------------------------##
##      Module scoring of top DEGs by geno_celltype              ####
##-----------------------------------------------------------------##

## Since we don't have a lot of DEGs per cell type, let's use all the genes as module and score them in AK, inhouse 2n, and published 2n
library(UCell)

## Module scoring with UCell
for(celltype in unique(out$ct)){
  #if(celltype == 'LE'){next}
  print(celltype)
  geneList = split(out$geneSym[out$ct == celltype],paste0(out$ct_geno[out$ct == celltype],'.',out$direction[out$ct == celltype]))
  big.srat <- UCell::AddModuleScore_UCell(big.srat, features = geneList,ncores = 5)
  
  write.csv(big.srat@meta.data,'~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/liver_UCELLscore_pbDEGs_tmp_2401.csv')
  
}


ucell_result = read.csv('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/liver_UCELLscore_pbDEGs_tmp_2401.csv')
df = ucell_result[,c('cellID','donorID','annot_jan24','Genotype',colnames(ucell_result)[grepl('EE',colnames(ucell_result))])]
df = pivot_longer(df,cols = colnames(df)[grepl('EE',colnames(df))],names_to = 'module',values_to = 'score')
df$group = ifelse(df$donorID %in% c('Hsb35','Hsb32','Hsb31'),'diploid_inhouse',df$Genotype)
ggplot(df[df$annot_jan24 %in% c('B.cell','HSC_MPP','MEMP_MEP','EE','MK','Mast.cell')& 
            grepl('T21.AK_up',df$module),],aes(donorID,score,fill=group))+
  #geom_point(size=0.001)+
  geom_boxplot(outlier.size = 0,outlier.colour = 'white')+
  scale_fill_manual(values = col25)+
  facet_grid(annot_jan24 ~ group,scales = 'free_x') + 
  geom_hline(yintercept = 0.10)+
  theme_bw() + theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0))

# Idents(big.srat) = paste0(big.srat$annot_jan24,':',big.srat$Genotype)
# ct = 'MEMP_MEP'
# DotPlot(big.srat,idents = paste0(ct,':',c('T21','diploid')),group.by = 'donorID',
#         features = out$geneSym[out$ct == ct & out$log2FoldChange > 0]) + RotatedAxis()
# 
# 
# 







##--------------------------------------------------##
# Looks like this step is too memory intesive - haven't been able to get it to work
# I need nCell epxressing each gene for both 2n and all the other geno, not just 2n...
## Or at least is epxressed in 10 cells

if(file.exists('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/oct23_ageMatched/fLiver_2n_nCells_perGene_perCelltype.csv')){
  nCell_perGene_perCTGeno = read.csv('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/oct23_ageMatched/fLiver_2n_nCells_perGene_perCelltype.csv')
}else{
  # Calculate number of cells per Geno per CT epxressing each gene
  mtx = big.srat@assays$RNA@counts#[,big.srat$cellID[big.srat$Phase == 'G1']]
  big.srat$group = paste0(big.srat$finalAnn,':',big.srat$Genotype)
  mtx_list = split(colnames(mtx),big.srat$group[match(colnames(mtx),big.srat$cellID)])
  mtx_list = mtx_list[grepl('diploid$',names(mtx_list))]
  
  nCell_perGene_perCTGeno = data.frame()
  for(i in 1:length(mtx_list)){
    print(i)
    mtx.sub = mtx[,mtx_list[[i]]]
    nCellbyGene = rowSums(as.matrix(mtx.sub >= 0.5))
    tmp = data.frame(geneSym = names(nCellbyGene),
                     nCell = nCellbyGene,
                     ct_geno = names(mtx_list)[i])
    nCell_perGene_perCTGeno = rbind(nCell_perGene_perCTGeno,tmp)
  }
  
  nCell_perGene_perCTGeno$geno = sapply(strsplit(nCell_perGene_perCTGeno$ct_geno,split=':'),'[',2)
  nCell_perGene_perCTGeno$ct = sapply(strsplit(nCell_perGene_perCTGeno$ct_geno,split=':'),'[',1)
  write.csv(nCell_perGene_perCTGeno,'fLiver_2n_nCells_perGene_perCelltype.csv')  
}



colnames(nCell_perGene_perCTGeno)[colnames(nCell_perGene_perCTGeno) == 'nCell'] = 'nCell_expr_2n'
allDEGs.sub = merge(allDEGs.sub,nCell_perGene_perCTGeno[,c('geneSym','nCell_expr_2n','ct')],by=c('geneSym','ct'),all.x=T)

# mtx_byGenoCT = do.call(cbind,lapply(split(colnames(mtx),big.srat$group[match(colnames(mtx),big.srat$cellID)]),function(e) rowSums(as.matrix(mtx[,e] > 0.5))))
# mtx_byGenoCT = as.data.frame(mtx_byGenoCT)
# mtx_byGenoCT$gene = rownames(mtx_byGenoCT)
# mtx_byGenoCT = pivot_longer(mtx_byGenoCT, cols = 1:(ncol(mtx_byGenoCT) - 1),names_to = 'group',values_to = 'nCell_expr')
# mtx_byGenoCT$ct = gsub(':.*$','',mtx_byGenoCT$group)
# mtx_byGenoCT$group = gsub('complete_trisomy','3n',mtx_byGenoCT$group)
# mtx_byGenoCT_2n = mtx_byGenoCT[grepl('2n|diploid',mtx_byGenoCT$group),]
# colnames(mtx_byGenoCT_2n)[3] = 'nCell_expr_2n'

# table(allDEGs.sub$ct_geno %in% nCell_perGene_perCTGeno$ct_geno)
# table(allDEGs.sub$ct %in% mtx_byGenoCT_2n$ct)
#table(mtx_byGenoCT$group[!grepl('2n',mtx_byGenoCT$group) & !mtx_byGenoCT$group %in% allDEGs.sub$ct_geno])
# table(allDEGs.sub$geneSym %in% mtx_byGenoCT$gene)
# Add this info to the previous data frame
# allDEGs.sub = merge(allDEGs.sub,mtx_byGenoCT,by.x = c('geneSym','ct_geno','ct'),by.y = c('gene','group','ct'),all.x = T)
# allDEGs.sub = merge(allDEGs.sub,mtx_byGenoCT_2n[,colnames(mtx_byGenoCT_2n) != 'group'],by.x = c('geneSym','ct'),by.y = c('gene','ct'),all.x = T)

# allDEGs.sub$max_nCell = pmax(allDEGs.sub$nCell_expr,allDEGs.sub$nCell_expr_2n)
# allDEGs.sub = allDEGs.sub[allDEGs.sub$max_nCell > 10,]




allDEGs.sub = allDEGs
##----------------------------------------------------------------------##
##      Investigate DEGs across genotype for B.cell.prog              ####
##----------------------------------------------------------------------##
# It looks like T21 and MX results in more global DEGs than T18 - this is somewhat strange as I'd expect T18 to have much stronger consequences...
# Using Endo as an example
endo = allDEGs.sub[allDEGs.sub$ct == 'Endo',]
table(endo$geno)
endo = endo %>% group_by(geneSym,direction) %>% mutate(n_geno = n_distinct(geno))
Idents(big.srat) = big.srat$finalAnn
big.srat$tmp = paste0(big.srat$Genotype,':',big.srat$donorID)
df = endo[endo$n_geno == 1 & endo$direction == 'AK_up' & endo$geno == 'T21' & 
            endo$log2FoldChange > 0.5 & endo$cellFrac_g2 > endo$cellFrac_g1,]
DotPlot(big.srat,idents = 'Endo',group.by = 'Genotype',
        features = c(#unique(endo$geneSym[endo$n_geno == 4 & endo$direction == 'AK_up'])
                     #unique(endo$geneSym[endo$n_geno == 3 & endo$direction == 'AK_down']),
                     #unique(endo$geneSym[endo$n_geno == 2 & endo$direction == 'AK_up']),
                     unique(df$geneSym)[1:100]
                     )) + RotatedAxis() +
  theme(axis.text.x = element_text(size=8))

##------------------------------------------------------------------------------##
##    Look at the number of shared DEs across genotypes in each cell type     ####
##------------------------------------------------------------------------------##
library(UpSetR)
upsetPlotList = list()
for(celltype in unique(allDEGs.sub$ct)){
  df = allDEGs.sub[allDEGs.sub$ct == celltype,]
  df$group = paste0(df$geno,':',df$direction)
  df_gene_list = split(df$geneSym,df$group)
  
  p = upset(fromList(df_gene_list[grepl('up',names(df_gene_list))]))
  
  
  df = df %>% group_by(geneSym,direction) %>% mutate(nGeno = n_distinct(geno))
  deg_MX = df[df$geno == 'MX',]
  df$altered_inMX = ifelse(df$direction == 'AK_up',(df$geneSym %in% deg_MX$geneSym[deg_MX$direction == 'AK_up']),
                           (df$geneSym %in% deg_MX$geneSym[deg_MX$direction == 'AK_down']))
  
  ggplot(df,aes(cellFrac_g1,cellFrac_g2,col=nGeno))+
    geom_point(aes(shape=geno))+
    facet_wrap(vars(direction))+
    geom_hline(yintercept = 50)+
    geom_vline(xintercept = 50)+
    theme_bw()+theme(panel.grid = element_blank())
    
  # Top genes up-regulated: expressed in >=50% AK but  <= 50% 2n
  topUp = df[df$direction == 'AK_up' & df$cellFrac_g1 < 50 & df$cellFrac_g2 > 50,]
  topUp = topUp[order(topUp$log2FoldChange,decreasing = T),]
  big.srat$finalAnn_tmp = ifelse(big.srat$Phase != 'G1','cycling',big.srat$finalAnn)
  Idents(big.srat) = big.srat$finalAnn_tmp
  big.srat$ann = paste0(big.srat$finalAnn,'/',big.srat$Genotype)
  DotPlot(big.srat,idents = celltype,group.by = 'ann',features = unique(topUp$geneSym[topUp$altered_inMX == F & topUp$])) + RotatedAxis() + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))
  
  DotPlot(big.srat,idents = c('HSC_MPP','MEMP_MEP','CMP_GMP','LMPP_ELP','B.cell.prog'),group.by = 'ann',features = unique(topUp$geneSym)) + RotatedAxis() + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))
  
  DotPlot(big.srat,idents = c('HSC_MPP','MEMP_MEP','CMP_GMP','LMPP_ELP','B.cell.prog'),
          group.by = 'ann',features = unique(topUp$geneSym[topUp$geno == 'T21' & topUp$nGeno == 1])) + RotatedAxis() + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))
  
  
  mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML/AML_mergedMLDS_clean_noMTcells_annotated_0923.RDS')
  Idents(mlds) = mlds$finalAnn
  DotPlot(mlds,idents = 'Tumour',
          group.by = 'donorID',features = unique(topUp$geneSym[topUp$geneSym %in% rownames(mlds)])) + RotatedAxis() + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))
  
  DotPlot(big.srat,idents = celltype,group.by = 'ann',features = unique(topUp$geneSym)) + RotatedAxis() + 
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))
  
  #what are the 11 commonly upregulated across all geno?
  df2 = df %>% group_by(direction,geneSym) %>% summarise(nGeno = n_distinct(geno))
    
}
allDEGs.sub %>% group_by(ct,direction) %>% summarise(nGeno = n_distinct(geno))

## remove genes which are DEGs across all cell types
nCT = allDEGs.sub %>% group_by(geneSym,geno) %>% summarise(n_CT = n_distinct(ct))
nCT = nCT %>% group_by(geneSym) %>% mutate(nGeno = n_distinct(geno))
nCT_toRemove = nCT[nCT$n_CT >=10 & nCT$nGeno == 5,]

allDEGs.sub = allDEGs.sub[!allDEGs.sub$geneSym %in% nCT_toRemove$geneSym,]

DotPlot(big.srat,idents = unique(big.srat$ann[grepl('ME\\/',big.srat$ann)]),
        features = df2$geneSym[df2$nGeno == 3 & df2$direction == 'AK_up'])  + RotatedAxis()






##---------------------------------------##
##    Number of DEGs per cell type     ####
##---------------------------------------##
## For each cell type, determine if the gene is expressed in that cell type
Idents(big.srat) = big.srat$finalAnn
mtx_byGenoCT_2n = nCell_perGene_perCTGeno
all_nDEG_summary_v1 = all_nDEG_summary

allDEGs.sub = allDEGs.sub[allDEGs.sub$padj < 0.05,]
# allDEGs.sub$toKeep = ifelse(allDEGs.sub$direction == 'AK_up' & allDEGs.sub$cellFrac_g2 > allDEGs.sub$cellFrac_g1, T,
#                             ifelse(allDEGs.sub$direction == 'AK_down' & allDEGs.sub$cellFrac_g2 < allDEGs.sub$cellFrac_g1, T,F))
# allDEGs.sub = allDEGs.sub[allDEGs.sub$toKeep == T,]
all_nDEG_summary = data.frame()
allDEGs.sub.filtered = data.frame()
for(celltype in unique(big.srat$finalAnn)){
  if(!celltype %in% unique(allDEGs.sub$ct)){next}
  # Genes which are expressed in 2n cells are
  genesExpressed = mtx_byGenoCT_2n[mtx_byGenoCT_2n$ct == celltype & mtx_byGenoCT_2n$nCell_expr_2n >= 3,]
  message(sprintf('%d genes are expressed in 2n %s',n_distinct(genesExpressed$gene),celltype))
  
  ## Extract list of DEGs to consider
  degs = allDEGs.sub[allDEGs.sub$ct == celltype,]
  tmp = as.data.frame(table(degs$geneSym %in% genesExpressed$gene,degs$direction))
  
  ## If a gene is considered "not expressed" in 2n - it cannot be "down-reg" in AK compared to 2n
  ## However, it can be "up-regulated" in AK compared to 2n
  if(n_distinct(tmp$Var1) == 2){
    to_remove = tmp$Freq[tmp$Var1 == F & tmp$Var2 == 'AK_down']
    if(to_remove > 0){
      message(sprintf('Removing %d DEGs due to "not expressed" in 2n but "down-regulated" in AK',to_remove))  
      # Remove these genes
      genes_toRemove = degs$geneSym[!degs$geneSym %in% genesExpressed$gene & degs$direction == 'AK_down']
      degs = degs[!degs$geneSym %in% genes_toRemove,]
    }  
  }
  
  
  ## Calculate fraction of 'expressed' genes being DE
  degs$nGeneExpressed = n_distinct(genesExpressed$gene)
  allDEGs.sub.filtered = rbind(allDEGs.sub.filtered,degs)
  
  degs = degs %>% group_by(ct,ct_geno,nGeneExpressed,geno,direction) %>% summarise(nDEG = n_distinct(geneSym)) %>% 
    group_by(ct,ct_geno,nGeneExpressed,geno) %>% mutate(total_nDEG = sum(nDEG))
  degs$frac_nDEG = degs$nDEG / degs$nGeneExpressed
  degs$frac_total_nDEG = degs$total_nDEG / degs$nGeneExpressed
  
  all_nDEG_summary = rbind(all_nDEG_summary,degs)
}

write.csv(allDEGs.sub.filtered,file.path(outDir,paste0(tissue,'_allDEGs_filtered.csv')))


# nCell used varies because of removal of different individuals depending on the context of the comparison
#nCell_used_inDEGs = big.srat@meta.data %>% group_by(finalAnn,Genotype) %>% summarise(nCell = n_distinct(cellID))
#all_nDEG_summary = merge(all_nDEG_summary,nCell_used_inDEGs,by.x=c('geno','ct'),by.y = c('Genotype','finalAnn'))

if(file.exists(file.path(outDir,'liver_nCell_perCT_perGeno_used_pbDEGs.csv'))){
  nCell_perCT_perGeno = read.csv(file.path(outDir,'liver_nCell_perCT_perGeno_used_pbDEGs.csv'),row.names = 1)
  nCell_perCT_perGeno = nCell_perCT_perGeno[,colnames(nCell_perCT_perGeno) != 'X']
  nCell_perCT_perGeno$geno[nCell_perCT_perGeno$geno == 'complete_trisomy'] = 'Triploid'
}
all_nDEG_summary$geno[all_nDEG_summary$geno == '3n'] = 'Triploid'
all_nDEG_summary = merge(all_nDEG_summary,nCell_perCT_perGeno,by=c('geno','ct'),all.x = T)
all_nDEG_summary$nCell = all_nDEG_summary$nCell_2n + all_nDEG_summary$nCell_AK
all_nDEG_summary$nCell_AKto2n = all_nDEG_summary$nCell_AK/all_nDEG_summary$nCell_2n

write.csv(all_nDEG_summary,file.path(outDir,paste0(tissue,'_nDEGs_perCTperGeno_summary.csv')))
tissue = 'liver'
all_nDEG_summary = read.csv(paste0(tissue,'_nDEGs_perCTperGeno_summary.csv'))

fig1f_2nAK_Liv_fracDEG = function(){
  all_nDEG_summary = read.csv(file.path(outDir,paste0(tissue,'_nDEGs_perCTperGeno_summary.csv')))  
  
  p1 = ggplot(all_nDEG_summary,aes(log10(nCell),frac_total_nDEG))+
    geom_smooth(method = 'lm',se = F,col = 'grey',lwd=0.4) +
    #geom_point(aes(shape=geno,col=ct),size=3)+
    geom_point(aes(col=ct,size=nCell_AKto2n))+
    scale_size_continuous(range = c(0.5,5))+
    facet_wrap(vars(geno))+
    scale_color_manual(values = sample(col25))+
    theme_bw() + theme(panel.grid = element_blank())
  
  all_nDEG_summary$geno = factor(all_nDEG_summary$geno,c('T21','T18','T22','MX','Triploid'))
  
  cell_cols = c('MEMP_MEP' = '#d10808',
                "B.cell.prog" = '#279be3',"B.cell" = '#12436e','MK'='#8d439c','Hepatocyte'=pal37H[c(26)])
  a = rep(grey(0.8),n_distinct(all_nDEG_summary$ct)-length(cell_cols))
  names(a) = unique(all_nDEG_summary$ct[! all_nDEG_summary$ct %in% names(cell_cols)])
  cell_cols = c(cell_cols,a)
  data = all_nDEG_summary
  
  if(file.exists(file.path(plotDir,'Fig1f_2nAKLiver_fractionExpressedGenesDE_v3_rawData.tsv'))){
    data = read.delim(file.path(plotDir,'Fig1f_2nAKLiver_fractionExpressedGenesDE_v3_rawData.tsv'),sep = '\t')
  }
  plotFun_degFrac_perCT_perGeno = function(noFrame=FALSE,noPlot=FALSE){
    # p1 = ggplot(data,aes(log10(nCell),frac_total_nDEG))+
    #   geom_smooth(method = 'lm',se = F,col = 'grey',size=0.4) +
    #   geom_point(data = data[data$ct %in% names(cell_cols),],aes(col=ct,size=nCell_AKto2n))+
    #   #geom_point(data = data[data$ct %in% names(cell_cols),],aes(col=ct,size=nCell_AKto2n))+
    #   #scale_size_continuous(range = c(1,8),name = 'ratio of AK-to-Diploid cell count')+
    #   facet_wrap(vars(geno),ncol=5)+
    #   scale_color_manual(values = cell_cols,name = 'cell type')+
    #   theme_classic(base_size = 13) + 
    #   #scale_y_log10()+
    #   theme(panel.border = element_rect(fill=F),legend.position = 'bottom',legend.title = element_text(size=13))+
    #   xlab('log10 cell count') + ylab('Fraction of genes being DE')
    # print(p1)
    # 
    # p2 = ggplot(data,aes(log10(nCell),frac_total_nDEG))+
    #   geom_smooth(method = 'lm',se = F,col = 'grey',size=0.4) +
    #   geom_point(data = data[data$ct %in% names(cell_cols),],aes(col=ct,size=nCell_AKto2n))+
    #   #geom_point(data = data[data$ct %in% names(cell_cols),],aes(col=ct,size=nCell_AKto2n))+
    #   #scale_size_continuous(range = c(1,8),name = 'ratio of AK-to-Diploid cell count')+
    #   facet_wrap(vars(geno),ncol=5)+
    #   scale_color_manual(values = cell_cols,name = 'cell type')+
    #   theme_classic(base_size = 13) + 
    #   scale_y_log10()+
    #   theme(panel.border = element_rect(fill=F),legend.position = 'bottom',legend.title = element_text(size=13))+
    #   xlab('log10 cell count') + ylab('log10 Fraction of genes being DE')
    # print(p2)
    # data$ct = as.factor(data$ct)
    
    d_maxFrac = data %>% group_by(ct) %>% summarise(yEnd=max(frac_total_nDEG))
    
    data$lineage = ifelse(data$ct == 'HSC_MPP','HSC_MPP',
                          ifelse(data$ct %in% c('B.cell.prog','B.cell',"NK.T"), 'Lymphoid',
                                 ifelse(data$ct %in% c('MEMP_MEP','MK','Mast.cell','EE','ME','LE'), 'Megk/Ery',
                                        ifelse(data$ct %in% c('CMP_GMP','Monocyte','Macrophage','Kupffer.cell','myelocyte','DC2','pDC'), 'Myeloid',
                                               ifelse(data$ct %in% c('Endo','Fibroblast','Hepatocyte'), 'Stromal','others')))))
    data$lineage = factor(data$lineage,c('HSC_MPP','Megk/Ery','Myeloid','Lymphoid','Stromal'))
    data$ct = factor(as.character(data$ct),c('HSC_MPP',
                                             'MEMP_MEP','MK','Mast.cell','EE','ME','LE',
                                             'CMP_GMP','Monocyte','Macrophage','Kupffer.cell','myelocyte','DC2','pDC',
                                             'Endo','Fibroblast','Hepatocyte',
                                             'B.cell.prog','B.cell',"NK.T"
                                             ))
    d_maxFrac$ct = factor(d_maxFrac$ct,c('HSC_MPP',
                                         'MEMP_MEP','MK','Mast.cell','EE','ME','LE',
                                         'CMP_GMP','Monocyte','Macrophage','Kupffer.cell','myelocyte','DC2','pDC',
                                         'Endo','Fibroblast','Hepatocyte',
                                         'B.cell.prog','B.cell',"NK.T"
    ))
    geno_cols = c('diploid' = grey(0.7),
                  'T21' = '#b18db8', #7a4b82',#dea6af',
                  'T18' = '#3d5dad', # '#89cff0',#pal37H[17],#'#3B87C7',
                  'T22' = '#679551',
                  'T13' = '#526691',
                  'MX' = '#e07d26',
                  'Triploid' = '#93221E')
    d_maxFrac$lineage = data$lineage[match(d_maxFrac$ct,data$ct)]
    data$group = ifelse(data$ct %in% c('HSC_MPP'),'HSC_MPP',
                        ifelse(data$ct %in% c('MEMP_MEP','MK','Mast.cell','EE','ME','LE'),'Meg/Ery/Mast',
                               ifelse(data$ct %in% c('CMP_GMP','Monocyte','Macrophage','Kupffer.cell','myelocyte','DC2','pDC'),'Myeloid',
                                      ifelse(data$ct %in% c('B.cell.prog','B.cell','NK.T'),'Lymphoid',
                                             ifelse(data$ct %in% c('Hepatocyte','Fibroblast','Endo'),'Stromal','others')))))
    data$group = factor(data$group,c('HSC_MPP','Meg/Ery/Mast','Myeloid','Lymphoid','Stromal'))
    
    data$ct_numericalID = NA
    for(g in unique(data$group)){
      d = data[data$group == g,]
      d$ct = factor(d$ct,levels = levels(d$ct)[levels(d$ct) %in% d$ct])
      d$ct_numericalID = as.numeric(d$ct)
      data$ct_numericalID[data$group == g] = d$ct_numericalID[match(data$ct[data$group == g],d$ct)]
    }
    
    data2 = data %>% group_by(group,ct,ct_numericalID) %>% summarise(max_frac = max(frac_total_nDEG))
    
    p3 = ggplot(data,aes(ct,log10(frac_total_nDEG)))+
      facet_grid(.~group,scales = 'free_x',space = 'free_x')+
      geom_point(aes(col=geno,size=log10(nCell)),alpha=0.5)+
      geom_segment(data=data2,aes(x = ct_numericalID,xend=ct_numericalID,y=min(log10(data$frac_total_nDEG)),yend=log10(max_frac)),col=grey(0.7))+
      scale_y_continuous(labels = c(10^(-2),10^(-3),10^(-4)),breaks = c(-2,-3,-4))+
      scale_color_manual(values = geno_cols)+
      theme_classic(base_size = 15)+
      theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
      theme(panel.border = element_rect(fill=F,linewidth = 1.5),axis.line = element_blank())+
      xlab('Cell type') + ylab('Fraction of DEGs') 
    print(p3)
  }
  
  saveFig(file.path(plotDir,'Fig1f_2nAKLiver_fractionExpressedGenesDE_v3_log10yScale'),plotFun_degFrac_perCT_perGeno,rawData=data,width = 11,height = 5.5,res = 500,useDingbats = T)
  #saveFig(file.path(plotDir,'Fig1f_2nAKLiver_fractionExpressedGenesDE_v3'),plotFun_degFrac_perCT_perGeno,rawData=data,width = 11,height = 5.5,res = 500,useDingbats = T)
}




## Plot fraction of DEGs being on affected chromosome

df = allDEGs.sub.filtered[allDEGs.sub.filtered$geno != '3n',]
df$chr_type = ifelse(df$chr == gsub('T|M','',df$geno),'local','global')
df = df %>% group_by(geno,ct,chr_type) %>% summarise(nDEG = n_distinct(geneSym)) %>% 
  group_by(geno,ct) %>% mutate(total_DEG = sum(nDEG),chrType_frac = nDEG/total_DEG)


d = df[df$chr_type == 'local',]
d$ct = factor(d$ct,levels = unique(d$ct))
ggplot(d,aes(ct,y=chrType_frac,col=geno))+
  geom_point()+
  geom_line(aes(x=as.numeric(ct),y = chrType_frac,col=geno),lwd=0.4)+
  theme_classic()
  
d$geno = factor(d$geno,levels = c('T21','T18','T22','MX'))
geno_cols = c('diploid' = grey(0.7),
              'T21' = '#b18db8',#7a4b82',#dea6af',
              'T18' = '#3d5dad',#'#89cff0',#pal37H[17],#'#3B87C7',
              'T22' = '#679551',
              'T13' = '#526691',
              'MX' = '#e07d26',
              'Triploid' = '#93221E')
p1 = ggplot(d,aes(geno,y=chrType_frac,fill=geno))+
  geom_boxplot(col='black',alpha=0.6) +
  geom_jitter(aes(col=ct),size=0.8,width = 0.15)+
  scale_fill_manual(values = geno_cols)+
  scale_color_manual(values = fLiver_celltype_cols)+
  theme_classic(base_size = 11) + theme(legend.key.size = unit(2,'mm'))+
  ylab('Fraction of DEGs on affected chromosome') + xlab('Genotype')
pdf('liver_plots/liver_fracDEGsOnAffectedChr.pdf',width = 4.5,height = 3.5)
print(p1)
dev.off()


# ggplot(all_nDEG_summary,aes(geno,frac_nDEG,fill=direction))+
#   geom_col(position = 'dodge')+
#   geom_point(aes(size=nCell))+
#   facet_wrap(vars(ct))+
#   scale_fill_manual(values = col25)+
#   theme_bw() + theme(panel.grid = element_blank())


















# all_nDEG_summary_withBackground = data.frame()
# for(g in unique(all_nDEG_summary$geno)){
#   main = all_nDEG_summary[all_nDEG_summary$geno == g,]
#   bg = all_nDEG_summary[all_nDEG_summary$geno != g,]
#   bg$ct = 'bg'
#   bg$geno = g
#   
#   all_nDEG_summary_withBackground = rbind(all_nDEG_summary_withBackground,main)
#   all_nDEG_summary_withBackground = rbind(all_nDEG_summary_withBackground,bg)
# }
# 
# all_nDEG_summary_withBackground$geno[all_nDEG_summary_withBackground$geno == '3n'] = 'Triploid'
# all_nDEG_summary_withBackground$geno = factor(all_nDEG_summary_withBackground$geno,c('T21','T18','T22','MX','Triploid'))
# 
# 
# all_nDEG_summary_withBackground = read.delim(file.path(plotDir,'Figxx_2nAKLiver_fractionExpressedGenesDE_rawData.tsv'),sep = '\t')
# 
# cell_cols = c('MEMP_MEP' = '#9e1b32',
#               "B.cell.prog" = '#005579',"B.cell" = '#702963','MK'='#9C6DA5','Hepatocyte'=pal37H[c(26)])
# a = rep(grey(0.8),n_distinct(all_nDEG_summary_withBackground$ct)-3)
# names(a) = unique(all_nDEG_summary_withBackground$ct[! all_nDEG_summary_withBackground$ct %in% names(cell_cols)])
# cell_cols = c(cell_cols,a)
# plotFun_genoContribution_no2n = function(noFrame=FALSE,noPlot=FALSE){
#   p1 = ggplot(all_nDEG_summary_withBackground[all_nDEG_summary_withBackground$ct != 'bg',],aes(log10(nCell),frac_total_nDEG))+
#     geom_smooth(data = all_nDEG_summary_withBackground[all_nDEG_summary_withBackground$ct != 'bg',],method = 'lm',se = F,col = 'grey',size=0.4) +
#     #geom_point(aes(shape=geno,col=ct),size=3)+
#     #geom_point(data = all_nDEG_summary_withBackground[all_nDEG_summary_withBackground$ct != 'bg',],aes(col=ct,size=nCell_AKto2n))+
#     geom_point(data = all_nDEG_summary_withBackground[all_nDEG_summary_withBackground$ct != 'bg' & !all_nDEG_summary_withBackground$ct %in% c('MEMP_MEP','B.cell','B.cell.prog','MK','Hepatocyte'),],aes(col=ct,size=nCell_AKto2n))+
#     geom_point(data = all_nDEG_summary_withBackground[all_nDEG_summary_withBackground$ct != 'bg' & all_nDEG_summary_withBackground$ct %in% c('MEMP_MEP','B.cell','B.cell.prog','MK','Hepatocyte'),],aes(col=ct,size=nCell_AKto2n))+
#     #scale_size_continuous(range = c(1,8),name = 'ratio of AK-to-Diploid cell count')+
#     facet_wrap(vars(geno),ncol=3)+
#     scale_color_manual(values = cell_cols,name = 'cell type')+
#     theme_classic(base_size = 13) + 
#     theme(panel.border = element_rect(fill=F),legend.position = 'bottom',legend.title = element_text(size=13))+
#     xlab('log10 cell count') + ylab('Fraction of genes being DE')
#   print(p1)
# }
# 
# saveFig(file.path(plotDir,'Figxx_2nAKLiver_fractionExpressedGenesDE_v2'),plotFun_genoContribution_no2n,rawData=all_nDEG_summary_withBackground,width = 6,height = 3.8,res = 500,useDingbats = T)
# 
# ggplot(all_nDEG_summary[all_nDEG_summary$ct == 'MEMP_MEP',],aes(log10(nCell_AK),log10(nCell_2n),size=frac_total_nDEG))+
#   geom_point(aes(shape=geno,col=ct))+
#   #geom_smooth(method = 'lm') +
#   #facet_wrap(vars(ct))+
#   scale_color_manual(values = sample(col25))+
#   theme_bw() + theme(panel.grid = element_blank())






## Let's see if these DEGs are associated with Phase
DotPlot(big.srat,group.by = 'Phase',
        features = allDEGs.sub$geneSym[allDEGs.sub$ct == 'MEMP_MEP' & allDEGs.sub$geno == 'T21'][1:50]) + RotatedAxis()






##--------------------------------##
##    T21_MEMP unique DEGs      ####
##--------------------------------##
allDEGs.sub.filtered = read.csv(file.path(outDir,paste0(tissue,'_allDEGs_filtered.csv')))

# Check if T21_MEMP degs are unique to T21
deg_summary = allDEGs.sub.filtered %>% group_by(ct,direction,geneSym) %>% summarise(nGeno = n_distinct(geno))

Idents(srat) = srat$annot_jan24
DotPlot(srat,ident=c('MEMP_MEP'),group.by = 'tmp',
        features = deg_summary$geneSym[deg_summary$ct == 'MEMP_MEP' & 
                                         deg_summary$geneSym %in% allDEGs.sub.filtered$geneSym[allDEGs.sub.filtered$ct == 'MEMP_MEP' & allDEGs.sub.filtered$geno == 'T21' & allDEGs.sub.filtered$direction == 'AK_up'] & 
                                         deg_summary$direction == 'AK_up' & deg_summary$nGeno > 1]
        ) + RotatedAxis()


t21_memp_degs = rbind(allDEGs.sub.filtered[allDEGs.sub.filtered$ct == 'MEMP_MEP' & allDEGs.sub.filtered$geno == 'T21' & allDEGs.sub.filtered$direction == 'AK_up' & 
                                               allDEGs.sub.filtered$geneSym %in% deg_summary$geneSym[deg_summary$ct == 'MEMP_MEP' & 
                                                                                                       deg_summary$direction == 'AK_up' & deg_summary$nGeno == 1],],
                      allDEGs.sub.filtered[allDEGs.sub.filtered$ct == 'MEMP_MEP' & allDEGs.sub.filtered$geno == 'T21' & allDEGs.sub.filtered$direction == 'AK_down' & 
                                             allDEGs.sub.filtered$geneSym %in% deg_summary$geneSym[deg_summary$ct == 'MEMP_MEP' & 
                                                                                                     deg_summary$direction == 'AK_down' & deg_summary$nGeno == 1],])
t21_memp_degs = t21_memp_degs[abs(t21_memp_degs$log2FoldChange) > 0.5,]

write.csv(t21_memp_degs,file.path('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/liver_t21MEMP_uniqueDEGs.csv'))

DotPlot(srat,ident=c('MEMP_MEP'),group.by = 'tmp',
        features = t21_memp_degs$geneSym[t21_memp_degs$direction == 'AK_up']
) + RotatedAxis()



## Module scoring with UCell
library(UCell)
geneList = split(t21_memp_degs$geneSym,paste0(gsub(':','.',t21_memp_degs$ct_geno),'.',t21_memp_degs$direction))
srat <- UCell::AddModuleScore_UCell(srat, features = geneList,ncores = 5)
  
write.csv(srat@meta.data,'~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/liver_UCELLscore_t21MEMP_2401.csv')




##-------------------------------------------------























geneExpr_nCell_summary = data.frame()
geneExpr_nCell_byDonorID_summary = data.frame()

mtx = big.srat@assays$RNA@counts[,big.srat$cellID[big.srat$Genotype == 'diploid']]
geneExpr_byCT_2n = do.call(cbind,lapply(split(colnames(mtx),big.srat$finalAnn[match(colnames(mtx),big.srat$cellID)]),function(e) rowSums(as.matrix(mtx[,e] > 0.5))))


for(celltype in unique(big.srat$finalAnn)){
  print(celltype)
  cellIDs = big.srat$cellID[big.srat$finalAnn == celltype & big.srat$Genotype == 'diploid']
  toc = big.srat@assays$RNA@counts[,cellIDs]
  geneExpr = apply(toc,1,function(x){sum(x>0.55)})
  tmp = data.frame(genesym = names(geneExpr),
                   nCell_Expr = geneExpr,
                   totalCell = length(cellIDs),
                   celltype = celltype)
  geneExpr_nCell_summary = rbind(geneExpr_nCell_summary,tmp)
  
  
  
  
  mDat = big.srat@meta.data[big.srat@meta.data$cellID %in% colnames(toc),c('cellID','donorID')]
  toc_donor = lapply(split(colnames(toc),mDat[,'donorID']),function(e) toc[,e])
  
  geneExpr_byDonor = do.call(rbind,lapply(1:length(toc_donor), function(i){
    mtx = toc_donor[[i]]
    if(is.null(ncol(mtx))){
      geneExpr = as.numeric(mtx>0.55)
      names(geneExpr) = names(mtx)
      tmp = data.frame(genesym = names(geneExpr),
                       nCell_Expr = geneExpr,
                       totalCell = 1,
                       celltype = celltype,
                       donor = names(toc_donor)[i])
    }else{
      geneExpr = apply(mtx,1,function(x){sum(x>0.55)})
      tmp = data.frame(genesym = names(geneExpr),
                       nCell_Expr = geneExpr,
                       totalCell = ncol(toc_donor[[i]]),
                       celltype = celltype,
                       donor = names(toc_donor)[i])
    }
    
    return(tmp)
    
  }))
  
  geneExpr_nCell_byDonorID_summary = rbind(geneExpr_nCell_byDonorID_summary,geneExpr_byDonor)
}

head(geneExpr_nCell_summary)
geneExpr_nCell_summary$frac_Expr = geneExpr_nCell_summary$nCell_Expr / geneExpr_nCell_summary$totalCell
geneExpr_nCell_summary$isExpressed = geneExpr_nCell_summary$frac_Expr >= 0.3

nGeneExpressed_byCT = geneExpr_nCell_summary %>% group_by(celltype) %>% summarise(nExpressed = sum(isExpressed))


degs.filtered$nExpressed = nGeneExpressed_byCT$nExpressed[match(degs.filtered$ct,nGeneExpressed_byCT$celltype)]
degs.filtered$gene_CT = paste0(degs.filtered$geneSym,':',degs.filtered$ct)
geneExpr_nCell_summary$geneCT = paste0(geneExpr_nCell_summary$genesym,':',geneExpr_nCell_summary$celltype)
degs.filtered$isExpressed = geneExpr_nCell_summary$isExpressed[match(degs.filtered$gene_CT,geneExpr_nCell_summary$geneCT)]

nDEG_byCT = degs.filtered %>% filter(isExpressed) %>% group_by(ct,geno,direction,nExpressed) %>% summarise(nDE = n_distinct(geneSym))
nDEG_byCT$fracDE = nDEG_byCT$nDE / nDEG_byCT$nExpressed 

## Add number of cells for each celltype_geno
nCell_byGeno = big.srat@meta.data %>% group_by(finalAnn,Genotype) %>% summarise(nCell = n())
nCell_byGeno$nCellDip = nCell_byGeno$nCell[nCell_byGeno$Genotype == 'diploid'][match(nCell_byGeno$finalAnn,nCell_byGeno$finalAnn[nCell_byGeno$Genotype == 'diploid'])]
nCell_byGeno$totalCell = nCell_byGeno$nCell + nCell_byGeno$nCellDip
nCell_byGeno$id = gsub('_trisomy$','',paste0(nCell_byGeno$finalAnn,':',nCell_byGeno$Genotype))

# Add to nDEG_byCT
nDEG_byCT$id = paste0(nDEG_byCT$ct,':',nDEG_byCT$geno)
nDEG_byCT$totalCell = nCell_byGeno$totalCell[match(nDEG_byCT$id,nCell_byGeno$id)]
nDEG_byCT$AKcell = nCell_byGeno$nCell[match(nDEG_byCT$id,nCell_byGeno$id)]



nDEG_byCT_T21 = nDEG_byCT[nDEG_byCT$geno == 'T21',]
o = nDEG_byCT_T21[nDEG_byCT_T21$direction == 'AK_up',]
o = unique(o$ct[order(o$fracDE,decreasing = F)])
nDEG_byCT_T21$ct = factor(nDEG_byCT_T21$ct,levels = o)
ggplot(nDEG_byCT_T21,aes(ct,fracDE,fill=direction))+
  #geom_boxplot()+
  #geom_point(aes(col=geno))+
  geom_col()+
  scale_color_manual(values = col25[-6])+
  scale_fill_manual(values = col25[-6])+
  facet_wrap(vars(direction),ncol=1)+
  theme_bw() + xlab('') +theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))


nDEG_byCT$direction2 = ifelse(nDEG_byCT$direction == 'AK_up','up regulated','down regulated')
nDEG_byCT$direction2 = factor(nDEG_byCT$direction2,levels = c('up regulated','down regulated'))
colnames(nDEG_byCT)[1] = 'cell type'
colnames(nDEG_byCT)[colnames(nDEG_byCT) == 'AKcell'] = 'Number of AK cells'
ggplot(nDEG_byCT,aes(totalCell,fracDE,size=`Number of AK cells`))+
  #geom_boxplot()+
  geom_smooth(method = 'lm',se = T,col=grey(0.6),size=0.8,alpha=0.2)+
  geom_point(aes(col=`cell type`))+
  scale_color_manual(values = col25[-6])+
  facet_grid(vars(direction2),vars(geno))+
  theme_bw(base_size = 13) + xlab('Total number of cells') + ylab('Fraction of genes being DE')







































# a techinical gene should be changing in the same direction across celltypes
allDEGs.sum = allDEGs.sub %>% group_by(geneSym,direction) %>% summarise(nComp = n_distinct(ct_geno),
                                                                        nGeno = n_distinct(geno),
                                                                        nCT = n_distinct(ct)) %>% 
  group_by(geneSym) %>% mutate(nDirection = n_distinct(direction),
                               total_nComp = sum(nComp),
                               fracComp = nComp/sum(nComp))

allDEGs.sum$gene = allDEGs.sum$geneSym
allDEGs.sum$ensID = geneMap$ensID[match(allDEGs.sum$geneSym,geneMap$geneSym)]
allDEGs.sum = annotateGenes(allDEGs.sum,geneMap = geneMap)
allDEGs.sum$geneSym = allDEGs.sum$gene
allDEGs.sum$id = seq(1:nrow(allDEGs.sum))
allDEGs.sum$chr = geneMap$chr[match(allDEGs.sum$geneSym,geneMap$geneSym)]
# Max 61 CT_geno comparison, 17 celltype, 5 geno
technicalDEGs = allDEGs.sum[allDEGs.sum$nDirection == 1 & allDEGs.sum$nComp >= 15,]
technicalDEGs = rbind(technicalDEGs,
                      allDEGs.sum[allDEGs.sum$nDirection == 2 & allDEGs.sum$geneSym %in% allDEGs.sum$geneSym[allDEGs.sum$fracComp >= 0.7] & 
                                    allDEGs.sum$total_nComp >= 15 & !allDEGs.sum$id %in% technicalDEGs$id,])

View(allDEGs.sum[allDEGs.sum$nDirection == 2 & allDEGs.sum$geneSym %in% allDEGs.sum$geneSym[allDEGs.sum$fracComp >= 0.7] & 
                   allDEGs.sum$total_nComp >= 15 & !allDEGs.sum$id %in% technicalDEGs$id,])

## List of DEGs post filter to remove potential technical genes
degs.filtered = allDEGs.sub[!allDEGs.sub$geneSym %in% technicalDEGs$geneSym,]

## For each cell type, frequency of genes DEGs across different geno
freq_geno = degs.filtered %>% group_by(ct,direction,geneSym) %>% summarise(nGeno = n_distinct(geno)) %>% 
  group_by(ct,direction,nGeno) %>% summarise(nGeno_freq = n())

View(degs.filtered[degs.filtered$geneSym %in% freq_geno.1$geneSym[freq_geno.1$nGeno == 5 & freq_geno.1$direction == 'AK_up'],])
freq_geno.1 = degs.filtered %>% group_by(ct,direction,geneSym) %>% summarise(nGeno = n_distinct(geno))
freq_geno$nGeno = factor(freq_geno$nGeno,levels = c(1,2,3,4,5))
freq_geno = freq_geno %>% group_by(ct) %>% mutate(totalGenes = sum(nGeno_freq),
                                                                  nGeno_freq_frac = nGeno_freq/totalGenes)
freq_geno$direction2 = factor(ifelse(freq_geno$direction == 'AK_up','up regulated','down regulated'),levels = c('up regulated','down regulated'))
colnames(freq_geno)[colnames(freq_geno) == 'ct'] = 'cell type'
ggplot(freq_geno,aes(nGeno,nGeno_freq_frac))+
  #geom_col()+
  #geom_point()+
  geom_boxplot(outlier.size = 0.001)+
  geom_jitter(aes(col=`cell type`),width = 0.3)+
  scale_color_manual(values = col25[-6])+
  facet_wrap(vars(direction2))+
  #facet_grid(vars(direction),vars(ct))+
  theme_bw(base_size = 13) + xlab('Number of genotypes') + ylab('Fraction of all DEGs')




degs.filtered.sub = degs.filtered[degs.filtered$geneSym %in% colnames(tfAct.2),]
View(degs.filtered.sub[degs.filtered.sub$geno == 'T21' & degs.filtered.sub$ct == 'MEMP_MEP' & degs.filtered.sub$direction == 'AK_up',])

allDEGs.sum$nDirection = as.factor(allDEGs.sum$nDirection)
ggplot(allDEGs.sum,aes(nDirection,nComp))+
  geom_boxplot()
View(allDEGs[allDEGs$geneSym == 'TIMD4',])

FeaturePlot(big.srat,'TIMD4',cells = big.srat$cellID[big.srat$finalAnn_broad %in% c('Mono.Mac')])
DimPlot(big.srat,group.by = 'Genotype',cells = big.srat$cellID[big.srat$finalAnn_broad %in% c('MEMP_MEP','MEP')])
Idents(big.srat) = big.srat$finalAnn_broad
DotPlot(big.srat,features = 'NR2F2',group.by = 'Genotype',idents = c('Mono.Mac'))






##--------------------------------------##
##        Enrichment analysis         ####
##--------------------------------------##
# For each cell type - genotype, extract top 200 up/down regulated genes by log2FC --> enrichment

library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021",
         'HDSigDB_Human_2021','KEGG_2021_Human','MSigDB_Hallmark_2020','MSigDB_Oncogenic_Signatures',
         'TF_Perturbations_Followed_by_Expression','WikiPathways_2019_Human',
         'Cancer_Cell_Line_Encyclopedia','CCLE_Proteomics_2020','CellMarker_Augmented_2021',
         'Disease_Perturbations_from_GEO_down','Disease_Perturbations_from_GEO_up',
         'Drug_Perturbations_from_GEO_2014','Drug_Perturbations_from_GEO_down','Drug_Perturbations_from_GEO_up','DrugMatrix','IDG_Drug_Targets_2022',
         'Elsevier_Pathway_Collection','Enrichr_Submissions_TF-Gene_Coocurrence','Gene_Perturbations_from_GEO_down','Gene_Perturbations_from_GEO_up')
up <- enrichr(topUp$geneSym[topUp$nGeno == 4], dbs)
plotEnrich(up[[6]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")

enrichR_output = list()
degs.filtered = allDEGs
for(group in unique(degs.filtered$ct_geno)){
  print(group)
  
  tmp = degs.filtered[degs.filtered$ct_geno == group,] %>% group_by(direction) %>% top_n(200,abs(log2FoldChange))
  print(n_distinct(tmp$geneSym))
  up <- enrichr(tmp$geneSym[tmp$log2FoldChange >0], dbs)
  down <- enrichr(tmp$geneSym[tmp$log2FoldChange <0], dbs)

  enrichR_output[[group]][['up']] = up
  enrichR_output[[group]][['down']] = down
}


plotEnrich(enrichR_output[['MEMP_MEP:T18']][['up']][[6]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(enrichR_output[['MEMP_MEP:T21']][['down']][[6]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")





# ## Let's try to compute DEGs just by celltype for all Genotype at the same time
# 
# ## Loop over cell types
# srat.geno = big.srat
# out = list()
# for(tgtIdx in c(1:(length(unique(srat.geno@meta.data$finalAnn[srat.geno$Genotype!='diploid']))))){
#   #if(tgtIdx %in% c(7)){next}
#   tgtCell = unique(srat.geno@meta.data$finalAnn[srat.geno$Genotype!='diploid'])[tgtIdx]  
#   if(tgtCell %in% c('?','unknown','others')){
#     next
#   }
#   message(sprintf("\n\n------- Consider cell type %s from all geno from tissue %s",tgtCell,tissue))
#   srat = subset(srat.geno,subset = finalAnn == tgtCell)
#   
#   #Check we have a reasonable number of counts from both settings. ie Get number of cells for the relevant cell type and genenotype
#   nCellsGroup = table(factor(srat@meta.data$Genotype))
#   
#   if((!all(nCellsGroup>=50))){
#     message(sprintf('Low number of cells detected'))
#     print(nCellsGroup)
#     next
#   }
#   
#   #Check how many from individual donors
#   nCells = table(srat@meta.data$donorID)
#   if(sum(nCells>50)<3){
#     message(sprintf("Too few effective replicates.  Skipping..."))
#     print(nCells)
#     next
#   }
#   message("Found the following number of high confidence cells")
#   print(nCells)
#   
#   #OK, we're going ahead, create the needed objects
#   toc = srat@assays$RNA@counts
#   toc = toc[rownames(toc) %in% geneMap$geneSym,]
#   m = match(rownames(toc),geneMap$geneSym)
#   sum(is.na(m))
#   rownames(toc) = geneMap$ensID[m]
#   mDat = data.frame(row.names = colnames(toc),
#                     cellID = colnames(toc),
#                     donor = srat@meta.data$donorID)
#   
#   mDat = merge(mDat,sampDat,by='donor')
#   mDat = mDat[match(colnames(toc),mDat$cellID),]
#   rownames(mDat) = colnames(toc)
#   mDat$group = factor(mDat$group,levels = c('diploid','complete_trisomy','T21','T22','T18','MX'))
#   # check that rownames(mDat) is in correct order
#   if(!all(rownames(mDat) == mDat$cellID)){
#     stop(sprintf('Incorrect mDat cellID order for tissue %s cell type %s genotype %s. Please check!',tissue, tgtCell, geno))
#   }
#   # Remove cellID column from mDat
#   mDat = mDat[,colnames(mDat) != 'cellID']
#   
#   coords = gns[rownames(toc)]
#   
#   #formula = ifelse(geno == 'T22','~ %s + sex','~ %s + ageGroup + sex')
#   formula = '~ %s + ageGroup + sex'
#   pdf(file.path(plotDir,paste0(tissue,'___unmerge___',tgtCell,'.pdf')))
#   out[[paste0(tgtCell,'_unmerged')]] = compareCell(toc = toc,
#                                                    mDat = mDat[match(colnames(toc),rownames(mDat)),],
#                                                    coords = gns[rownames(toc)],
#                                                    cellTypeName=tgtCell,
#                                                    formula=formula,
#                                                    donorID='donor',groupID='group')
#   dev.off()
#   
#   pdf(file.path(plotDir,paste0(tissue,'___merge',kMerge,'___',tgtCell,'.pdf')))
#   out[[paste0(tgtCell,'_merge',kMerge)]] = compareCell(toc,mDat,coords,
#                                                        kMerge=kMerge,
#                                                        formula = formula,
#                                                        cellTypeName=tgtCell)
#   dev.off()
# }
# 
# saveRDS(out,paste0(outDir_fp,'/',geno,'_',tissue,'_pseudoBulkDESeq2out.RDS'))
# 
# 
# #############################################
# # Plot median log2FC across celltypes       #
# #############################################
# # Get the median log2FC per chromosome for each cell type --> do a boxplot (x = chromosome) + dots coloured by celltypes
# log2FC = tibble()
# #,paste0('merge',kMerge)
# for(dat in c('unmerged')){
#   #for(dat in c('unmerged')){
#   results = out[grepl(dat,names(out))]
#   celltypes = gsub('_unmerged$','',names(results))
#   for(idx in 1:length(results)){
#     celltype = celltypes[idx]
#     #celltype = sapply(strsplit(names(results)[idx],split='_'),'[',1)
#     print(celltype)
#     log2FC_perChr = results[[idx]][['log2FC_perChr']]
#     
#     #boxplot(log2FC_perChr,
#     #        outline=FALSE,
#     #        xlab='Chromosome',
#     #        ylab='logFC'
#     #)
#     #abline(h=0,col='red')
#     
#     med = sapply(log2FC_perChr,FUN = function(x){median(x)})
#     mu = sapply(log2FC_perChr,FUN = function(x){mean(x)})
#     #plot(x=names(med),y=med,pch=19)
#     #abline(h=0,col='red')
#     
#     tmp = data.frame(celltype = celltype, med_log2FC = med,mu_log2FC=mu, chr = names(med),type=dat)
#     log2FC = rbind(log2FC,tmp)
#     
#   }
# }
# 
# 
# log2FC$chr = factor(log2FC$chr,levels = c(1:22,'X'))
# pdf(file.path(outDir_fp,paste0(geno,'_',tissue,'___',dat,'___medianLog2FCperChr.pdf')),width = 13,height = 10)
# p = ggplot(log2FC,aes(x=chr, y=med_log2FC))+
#   geom_boxplot()+
#   geom_point(aes(col = celltype))+
#   scale_color_manual(values = c(col25,brewer.pal(12,'Paired')))+
#   geom_hline(yintercept = 0)+
#   theme_bw(base_size = 13)+
#   facet_wrap(vars(type),ncol=1,scales = 'free')+
#   xlab('Chromosome')+ylab('median log2FC')+
#   ggtitle(paste0(geno,'_',tissue))
# 
# print(p)
# dev.off()
# 
# write.table(log2FC,file.path(outDir_fp,paste0(geno,'_',tissue,'___',dat,'___medianLog2FCperChr.csv')),sep = ',',row.names = F,col.names = T)
# 
# 
# 
# 
# 
# 
# ######################################################
# # Plot Number of DEGs per Chr across celltypes       #
# ######################################################
# # Get theNumber of DEGs per chromosome for each cell type --> do a boxplot (x = chromosome) + dots coloured by celltypes
# #,paste0('merge',kMerge)
# for(dat in c('unmerged')){
#   celltypes = gsub('_unmerged$','',names(results))
#   results = out[grepl(dat,names(out))]
#   
#   for(idx in 1:length(results)){
#     #celltype = sapply(strsplit(names(results)[idx],split='_'),'[',1)
#     celltype = celltypes[idx]
#     print(celltype)
#     n_deg_perChr = results[[idx]][['deTable']]
#     n_deg_perChr$celltype = celltype
#     n_deg_perChr$geno = geno
#     n_deg_perChr$type = dat
#     
#     n_deg_perChr_allCT = rbind(n_deg_perChr_allCT,n_deg_perChr)
#   }
# }












########## May 10th ###########

#From the results above, some celltypes are extremely homogenous --> very few DEGs
#I'd like to have a look at some of these manually
#eg. Mono.Mac, Kupffer cells


# tgtCell = 'Kupffer.Cell'
# srat = subset(srat.geno,subset = finalAnn == tgtCell)
# srat = standard_clustering(srat)
# 
# DimPlot(srat,group.by = 'gestationalAge')
# degs = as.data.frame(out[['Kupffer.Cell_unmerged']][['de']])
# degs = degs[degs$padj < 0.05,]
# degs$ensID = rownames(degs)
# degs$geneSym = geneMap$geneSym[match(degs$ensID,geneMap$ensID)]
# degs$chr = as.character(seqnames(gns[match(degs$ensID,gns$gene_id)]))
# 
# # Add fraction of cells expressing each gene
# degs$cellPerc = rowSums(srat@assays$RNA@counts[match(degs$geneSym,rownames(srat@assays$RNA@counts)),] > 0.55)
# degs$cellPerc = 100*degs$cellPerc/ncol(srat)
# 
# degs.sub = degs[degs$cellPerc >= 30,]
# degs.sub = degs[degs$cellPerc < 30,]
# 
# VlnPlot(srat,c('HES1','ARL6IP1','DONSON','DYRK1A','UBE2G2'),group.by = 'Genotype')
# 
# FeaturePlot(srat,'DYRK1A')
# 
# 
# # EnrichR
# ## Enrichment analysis ###
# library(enrichR)
# listEnrichrSites()
# setEnrichrSite("Enrichr") # Human genes
# websiteLive <- TRUE
# dbs <- listEnrichrDbs()
# if (is.null(dbs)) websiteLive <- FALSE
# if (websiteLive) head(dbs)
# 
# dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021",
#          'HDSigDB_Human_2021','KEGG_2021_Human','MSigDB_Hallmark_2020','MSigDB_Oncogenic_Signatures',
#          'TF_Perturbations_Followed_by_Expression','WikiPathways_2019_Human',
#          'Cancer_Cell_Line_Encyclopedia','CCLE_Proteomics_2020','CellMarker_Augmented_2021',
#          'Disease_Perturbations_from_GEO_down','Disease_Perturbations_from_GEO_up',
#          'Drug_Perturbations_from_GEO_2014','Drug_Perturbations_from_GEO_down','Drug_Perturbations_from_GEO_up','DrugMatrix','IDG_Drug_Targets_2022',
#          'Elsevier_Pathway_Collection','Enrichr_Submissions_TF-Gene_Coocurrence','Gene_Perturbations_from_GEO_down','Gene_Perturbations_from_GEO_up')
# 
# if (websiteLive) {
#   up <- enrichr(degs.sub$geneSym[degs$log2FoldChange >0], dbs)
#   down <- enrichr(degs.sub$geneSym[degs$log2FoldChange <0], dbs)
# }
# 
# plotEnrich(up[[1]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
# plotEnrich(down[[1]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
# 
# 
# 
# 
# ### Look at an example of cell types which have lots of DEGs, eg. Endo or MEMP
# 
# endo.srat = subset(srat.geno,subset = finalAnn == 'Endo')
# endo.srat = standard_clustering(endo.srat)
# DimPlot(endo.srat,group.by = 'Genotype')
# endo.degs = as.data.frame(out[['Endo_unmerged']][['de']])
# endo.degs = endo.degs[endo.degs$padj < 0.05,]
# endo.degs$ensID = rownames(endo.degs)
# endo.degs$geneSym = geneMap$geneSym[match(endo.degs$ensID,geneMap$ensID)]
# endo.degs$chr = as.character(seqnames(gns[match(endo.degs$ensID,gns$gene_id)]))
# # Add fraction of cells expressing each gene
# endo.degs$cellPerc = rowSums(endo.srat@assays$RNA@counts[match(endo.degs$geneSym,rownames(endo.srat@assays$RNA@counts)),] > 0.55)
# endo.degs$cellPerc = 100*endo.degs$cellPerc/ncol(endo.srat)
# 
# endo.degs.sub = endo.degs[endo.degs$cellPerc >= 10,]
# endo.degs.sub2 = endo.degs[endo.degs$cellPerc < 10,]
# 
# VlnPlot(srat,c('HES1','ARL6IP1','DONSON','DYRK1A','UBE2G2'),group.by = 'Genotype')
# 
# FeaturePlot(endo.srat,'CDKN1A')
# 
# 
# ### Compute frequency of each gene being DEs
# allDEGs2 = do.call(rbind,lapply(seq(1:(length(out)/2)),function(i){
#   de = as.data.frame(out[grepl('_unmerge',names(out))][[i]][['de']])
#   de = de[de$padj < 0.05,]
#   de$ensID = rownames(de)
#   de$geneSym = geneMap$geneSym[match(de$ensID,geneMap$ensID)]
#   de$chr = as.character(seqnames(gns[match(de$ensID,gns$gene_id)]))
#   de$ct = gsub('_unmerged','',names(out[grepl('_unmerge',names(out))])[i])
#   de$geno = 'T21'  
#   # Add fraction of cells expressing each gene
#   de$cellPerc = rowSums(srat.geno@assays$RNA@counts[match(de$geneSym,rownames(srat.geno@assays$RNA@counts)),srat.geno$cellID[srat.geno$finalAnn == unique(de$ct)]] > 0.55)
#   degs$cellPerc = 100*degs$cellPerc/length(srat.geno$cellID[srat.geno$finalAnn == unique(de$ct)])
#   
#   return(de)
# }))
# 
# allDEGs$direction = ifelse(allDEGs$log2FoldChange >0,'up','down')
# allDEGs.sub = allDEGs[allDEGs$cellPerc >= 10 & abs(allDEGs$log2FoldChange) >= log2(1.5),]
# allDEGs.sub = allDEGs.sub %>% group_by(geno,geneSym,direction) %>% mutate(nCT = n_distinct(ct))
# allDEGs.sub$cat = ifelse(allDEGs.sub$nCT >= 10,'>=10',as.character(allDEGs.sub$nCT))
# 
# df = allDEGs.sub %>% group_by(geno,geneSym,direction) %>% summarise(nCT = n_distinct(ct))
# 
# # ~ 7k DEGs across 20 celltypes in T21 vs diploid fLiver comparison
# df$cat = ifelse(df$nCT >= 10,'>=10',as.character(df$nCT))
# df2 = as.data.frame(table(df$cat,df$direction))
# df2$Var1 = factor(df2$Var1,levels = c(as.character(seq(1,9,1)),'>=10'))
# ggplot(df2,aes(x=Var1,y=Freq,fill=Var2))+
#   geom_col(position = 'dodge') + theme_bw(base_size = 15) + xlab('# cell types') + ylab('# DE Genes') + ggtitle('T21 fLiver')
# 
# VlnPlot(srat.geno,group.by = 'finalAnn',split.by = 'Genotype',features = 'A1BG')
# 
# 
# ### Look at genes which are DEGs across multiple cell types
# commonDEGs = allDEGs.sub[allDEGs.sub$cat == '>=10',]
# 
# # Plot normalized average expression level of TF genes 
# # Compute pseudobulk expression by CT:Geno
# mtx = srat.geno@assays$RNA@counts[rownames(srat.geno@assays$RNA@counts) %in% unique(commonDEGs$geneSym),
#                                 colnames(srat.geno@assays$RNA@counts) %in% srat.geno$cellID[srat.geno$finalAnn %in% commonDEGs$ct]]
# mDat = srat.geno@meta.data[rownames(srat.geno@meta.data) %in% colnames(mtx),c('cellID','Genotype','donorID','finalAnn')]
# mDat$ann = paste0(mDat$donorID,':',mDat$finalAnn)
# mDat = mDat %>% group_by(donorID,finalAnn) %>% mutate(nCell = n())
# mDat = mDat[mDat$nCell >= 50,]
# mtx = mtx[,mDat$cellID]
# pb = do.call(cbind,lapply(split(colnames(mtx),mDat[,'ann']),function(e) rowSums(mtx[,e,drop=FALSE])))
# pb = as.matrix(pb)
# 
# # minibulk_normalization: expression of each gene in each minibulk is normalized (divided by) the total count in the minibulk
# pb_norm = sweep(pb,2,colSums(pb),`/`)
# # log transformation
# pb_norm = log(pb_norm * 10^4 + 1)
# 
# #pb_norm = pb_norm[rowSums(pb_norm) != 0,]
# # Scale expression by genotype
# #pb_scaledGeno = do.call(cbind,lapply(split(colnames(pb_norm),mDat[match(gsub(':.*$','',colnames(pb_norm)),mDat$finalAnn_broad),'finalAnn_broad']),function(e) t(scale(t(pb_norm[,e])))))
# 
# library(ComplexHeatmap)
# library(circlize)
# 
# 
# ct_Cols = col25[1:length(unique(mDat$finalAnn))]
# names(ct_Cols) = unique(mDat$finalAnn)
# 
# geno_Cols = col25[1:length(unique(mDat$Genotype))]
# names(geno_Cols) = unique(mDat$Genotype)
# 
# nCell.Cols  = grey(20:0/20)
# pCols = circlize::colorRamp2(seq(min(mDat$nCell),max(mDat$nCell),length.out=length(nCell.Cols)),nCell.Cols)
# 
# 
# 
# mat = pb_norm[!grepl('^RP',rownames(pb_norm)),]
# mat = t(scale(t(mat)))
# 
# q = quantile(mat[!is.na(mat)],seq(0,1,0.15))
# 
# 
# botAnno = HeatmapAnnotation(df = data.frame(Geno = mDat$Genotype[match(colnames(mat),mDat$ann)],
#                                             nCell = mDat$nCell[match(colnames(mat),mDat$ann)]),
#                             annotation_name_side = 'left',
#                             col = list(Geno = geno_Cols,
#                                        nCell = pCols))
# 
# colSplits = mDat$finalAnn[match(colnames(mat),mDat$ann)]
# 
# 
# 
# hm = Heatmap(mat,na_col = 'grey',
#              #col=colorRamp2(q, rev(c('#162668','#1A2D7A','#3A4B8D','#7C87B3','#9CA5C6','#BDC3D9','#DEE1EC'))),
#              #column_split = colSplits,
#              row_split = df$direction[match(rownames(mat),df$geneSym)],
#              show_column_names = T, cluster_columns = T, cluster_row_slices = T,
#              show_row_dend = F,
#              row_names_gp = gpar(fontsize=5),
#              
#              #row_title_gp = gpar(fontsize=12),
#              row_title_rot = 0,
#              column_title_gp = gpar(fontsize=5),
#              column_names_gp = gpar(fontsize=5),
#              column_title_rot = 90,
#              column_gap = unit(2,'mm'),
#              cluster_rows = T,show_row_names = T,
#              bottom_annotation=botAnno)
# draw(hm)
# 
# 
# ### EnrichR on this set of T21_commonDEGs
# 
# up <- enrichr(df$geneSym[df$direction == 'up' & !grepl('^RPL|^RPS',df$geneSym)], dbs)
# down <- enrichr(df$geneSym[df$direction == 'down' & !grepl('^RPL|^RPS',df$geneSym)], dbs)
# 
# 
# plotEnrich(up[[9]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
# plotEnrich(down[[9]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
# 
