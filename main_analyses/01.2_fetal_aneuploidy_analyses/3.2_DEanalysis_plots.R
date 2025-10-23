### Differential Gene Expression analysis between fLiver AK vs diploid ###
### Part 2: analysis ###

# Notes from the analysis:
# amongst the diploid foetuses, only Hsb32 is male, all the other 3 samples are females. this might explain the strange PCA pattern where Hsb32 is always on its own
# Hsb32 markers: EIF1AY (y-linked gene)

outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jul24/without_published2n/genoAssay_withX'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)


#------------------------#
##      Libraries     ####
#------------------------#
library(tidyverse)
library(Seurat)
library(GenomicFeatures)
library(ComplexHeatmap)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
plotDir='~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots'


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


##------------------------------##
##  1. Import the results     ####
##------------------------------##
tissue = 'liver'
## Import fLiver object
# fLiver = filter_sratObj(tissue='liver',ageMatch = T,remove_cyclingCells = F,
#                         srat_in_fp='~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_annotated_0424.RDS')
# fLiver = fLiver[[1]]

fLiver = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS')
fLiver$Sex[fLiver$Sex == 'female'] = 'F'
fLiver$Sex[fLiver$Sex == 'male'] = 'M'
all(fLiver$annot == fLiver$annot_aug24)
fLiver$finalAnn = as.character(fLiver$annot)
fLiver$finalAnn[fLiver$finalAnn %in% c('earlyMK')] = 'MK'
fLiver$finalAnn[fLiver$finalAnn %in% c('promyelocyte','myelocyte')] = 'Myelocyte'
fLiver$finalAnn[fLiver$finalAnn %in% c('proMono')] = 'Monocyte'
fLiver$finalAnn[fLiver$finalAnn %in% c('LMPP_ELP','pro.B.cell','pre.B.cell')] = 'B.cell.prog'
fLiver$finalAnn[fLiver$finalAnn %in% c('ILC.precursor','T.cell','NK_T')] = 'NK.T'
fLiver$ageGroup = fLiver$gestationalAge


fLiver$finalAnn_broad = fLiver$finalAnn

## Import significant DEGs
resultDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jul24/without_published2n/genoAssay_withX/'
allDEGs = import_pbDEGresults(outDir = resultDir,tissue = 'liver')
head(allDEGs)
table(allDEGs$geno,allDEGs$ct)
dim(allDEGs[is.na(allDEGs$geneSym),])




##------------------------##
##  2. enrichR plot     ####
##------------------------##

## Perform enrichR on all DEGs
upDEG_enriched <- enrichr(unique(allDEGs$geneSym[allDEGs$geno == 'T21' & allDEGs$direction == 'AK_up']), dbs)
downDEG_enriched <- enrichr(unique(allDEGs$geneSym[allDEGs$geno == 'T18' & allDEGs$direction == 'AK_down']), dbs)

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





##---------------------------------------##
##  3. Number of DEGs per cell type    ####
##---------------------------------------##
## For each cell type, determine if the gene is expressed in that cell type
if(file.exists('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jul24/without_published2n/genoAssay_withX/fLiver_2n_nCells_perGene_perCelltype.csv')){
  nCell_perGene_perCTGeno = read.csv('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jul24/without_published2n/genoAssay_withX/fLiver_2n_nCells_perGene_perCelltype.csv',row.names = 1)
  #nCell_perGene_perCTGeno$ct[nCell_perGene_perCTGeno$ct == 'myelocyte'] = 'Myelocyte'
  nCell_perGene_perCTGeno = nCell_perGene_perCTGeno[nCell_perGene_perCTGeno$ct != 'Trophoblast',]
  nCell_perGene_perCTGeno$nCell = nCell_perGene_perCTGeno$nCell_0
}else{
  # Calculate number of (diploid only) cells per Geno per CT epxressing each gene
  mtx = fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$Genotype == 'diploid']]
  mtx_list = split(colnames(mtx),fLiver$finalAnn[match(colnames(mtx),fLiver$cellID)])
  nCell_perGene_perCTGeno = data.frame()
  print(length(mtx_list))
  for(i in 1:length(mtx_list)){
    print(i)
    mtx.sub = mtx[,mtx_list[[i]]]
    nCellbyGene = rowSums(mtx.sub >= 0.5)
    nCellbyGene_lowThreshold = rowSums(mtx.sub > 0)
    tmp = data.frame(geneSym = names(nCellbyGene),
                     nCell_0.5 = nCellbyGene,
                     nCell_0 = nCellbyGene_lowThreshold,
                     ct = names(mtx_list)[i])
    nCell_perGene_perCTGeno = rbind(nCell_perGene_perCTGeno,tmp)
  }
  nCell_perCT = as.data.frame(table(fLiver$finalAnn[fLiver$Genotype == 'diploid']))
  colnames(nCell_perCT) = c('ct','totalCell')
  nCell_perGene_perCTGeno = merge(nCell_perGene_perCTGeno,nCell_perCT,by='ct',all=T)
  nCell_perGene_perCTGeno$pct_epxressed_0.5 = 100*nCell_perGene_perCTGeno$nCell_0.5/nCell_perGene_perCTGeno$totalCell
  nCell_perGene_perCTGeno$pct_epxressed_0 = 100*nCell_perGene_perCTGeno$nCell_0/nCell_perGene_perCTGeno$totalCell
  nCell_perGene_perCTGeno$nCell_diff = nCell_perGene_perCTGeno$nCell_0 - nCell_perGene_perCTGeno$nCell_0.5
  
  write.csv(nCell_perGene_perCTGeno,'fLiver_2n_nCells_perGene_perCelltype.csv')  
}



## For each cell type, plot fraction of DEGs being on affected chromosome ##
allDEGs$on_affectedChr = F
allDEGs$on_affectedChr[allDEGs$geno=='T21' & allDEGs$chr == '21'] = T
allDEGs$on_affectedChr[allDEGs$geno=='T18' & allDEGs$chr == '18'] = T
allDEGs$on_affectedChr[allDEGs$geno=='T22' & allDEGs$chr == '22'] = T
allDEGs$on_affectedChr[allDEGs$geno=='MX' & allDEGs$chr == 'X'] = T
allDEGs$on_affectedChr[allDEGs$geno=='3n' & allDEGs$chr != 'X'] = T

## Column plot + up/down regulated
allDEGs.summary = allDEGs %>% group_by(geno,ct,direction,on_affectedChr) %>% summarise(nDEG = n_distinct(geneSym)) %>% 
  group_by(geno,ct,direction) %>% mutate(totalDEG=sum(nDEG),fracDEG = nDEG/totalDEG)
allDEGs.summary$fracDEG[allDEGs.summary$direction == 'AK_down'] = -allDEGs.summary$fracDEG[allDEGs.summary$direction == 'AK_down']
allDEGs.summary$geno = factor(allDEGs.summary$geno,c('MX','T18','T21','T22','3n'))
ggplot(allDEGs.summary,aes(geno,fracDEG,fill=on_affectedChr))+
  geom_col()+
  scale_fill_manual(values = c(grey(0.8),'red'))+
  facet_wrap(vars(ct))+
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + xlab('')


## scatter plot

allDEGs.summary = allDEGs %>% 
  group_by(geno,ct,on_affectedChr) %>% summarise(nDEG_noDirection=n_distinct(geneSym)) %>% 
  group_by(geno,ct) %>% mutate(totalDEG_noDirection=sum(nDEG_noDirection),fracDEG_noDirection = nDEG_noDirection/totalDEG_noDirection)

# Add number of cells going into the comparison
nCell_perCT_perGeno = read.csv(file.path(outDir,'liver_nCell_perCT_perGeno_used_pbDEGs.csv'),row.names = 1)
nCell_perCT_perGeno = nCell_perCT_perGeno[,colnames(nCell_perCT_perGeno) != 'X']
nCell_perCT_perGeno$geno[nCell_perCT_perGeno$geno == 'complete_trisomy'] = 'Triploid'
nCell_perCT_perGeno$ct_geno = paste0(nCell_perCT_perGeno$ct,'_',nCell_perCT_perGeno$geno)
allDEGs.summary$geno[allDEGs.summary$geno == '3n'] = 'Triploid'
allDEGs.summary$ct_geno = paste0(allDEGs.summary$ct,'_',allDEGs.summary$geno)
allDEGs.summary = cbind(allDEGs.summary,nCell_perCT_perGeno[match(allDEGs.summary$ct_geno,nCell_perCT_perGeno$ct_geno),c('nCell_2n','nCell_AK')])
all_nDEG_summary$ct_geno = paste0(all_nDEG_summary$ct,'_',all_nDEG_summary$geno)
all_nDEG_summary$ct_geno[grepl('3n',all_nDEG_summary$ct_geno)] = gsub('3n','Triploid',all_nDEG_summary$ct_geno[grepl('3n',all_nDEG_summary$ct_geno)])
allDEGs.summary$fracDE = all_nDEG_summary$frac_total_nDEG[match(allDEGs.summary$ct_geno,all_nDEG_summary$ct_geno)]

ggplot(allDEGs.summary[allDEGs.summary$ct == 'Endo' & allDEGs.summary$on_affectedChr==T,],aes(nCell_AK,fracDEG_noDirection,col=geno))+
  geom_point()+
  scale_color_manual(values = col25)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + xlab('')

ggplot(allDEGs.summary[allDEGs.summary$on_affectedChr==T,],aes(fracDE,fracDEG_noDirection,col=geno))+
  geom_point(aes(size=nCell_AK),alpha=1)+
  geom_point(aes(size=nCell_AK),shape = 1,colour = "black")+
  #facet_wrap(vars(ct))+
  scale_color_manual(values = geno_cols)+
  scale_x_log10()+
  theme_classic()+
  theme(panel.background = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(colour = 'black')) + 
  xlab('Number of DEG') + ylab('Fraction of DEGs being on Affected chromosome')


ggplot(allDEGs.summary,aes(geno,fracDEG_noDirection,fill=on_affectedChr))+
  geom_col()+
  scale_fill_manual(values = c(grey(0.8),'red'))+
  facet_wrap(vars(ct))+
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)) + xlab('')


## Boxplot
library(ggbeeswarm)
dd = allDEGs.summary[allDEGs.summary$on_affectedChr==T,]
colnames(dd) = c('Genotype','Celltype','on_affectedChr','nDEG_on_TrisomicChrom','total_nDEG','fracDEG_on_TrisomicChrom',
                 'ct_geno','nCell_2n','nCell_AK')
dd = dd[,!colnames(dd) %in% c('on_affectedChr','ct_geno')]

dd = read.delim(file.path(plotDir,'Fig1_2nAKLiver_localImpact_DEGs_rawData.tsv'),sep = '\t',header = T)
plotFun_localImpact = function(noFrame=FALSE,noPlot=FALSE){
  geno_cols = c('diploid' = grey(0.7),
                'T21' = '#93221E',
                'T18' = '#3d5dad',
                'T22' = '#679551',
                'T13' = '#526691',
                'MX' = '#b18db8',
                'Triploid' = '#e07d26')
  p = ggplot(dd,aes(reorder(Genotype,fracDEG_on_TrisomicChrom,FUN=median),fracDEG_on_TrisomicChrom,fill=Genotype))+
    geom_boxplot(width=0.7,color='black',outlier.shape = NA)+
    geom_quasirandom(size=0.8,width = 0.25)+
    scale_fill_manual(values = geno_cols)+
    theme_classic(base_size = 13)+
    theme(panel.background = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
          axis.text = element_text(color = 'black'),
          axis.ticks = element_line(colour = 'black')) + 
    xlab('') + ylab('Local impact')
  
  dd$globalImpact = 1- dd$fracDEG_on_TrisomicChrom
  p2 = ggplot(dd,aes(reorder(Genotype,globalImpact,FUN=median),globalImpact,fill=Genotype))+
    geom_boxplot(width=0.7,color='black',outlier.shape = NA)+
    geom_quasirandom(size=0.5,width = 0.25)+
    scale_fill_manual(values = geno_cols)+
    theme_classic(base_size = 13)+
    theme(panel.background = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
          axis.text = element_text(color = 'black'),
          axis.ticks = element_line(colour = 'black')) + 
    xlab('') + ylab('Global impact')
  
  
  print(p2)
}

saveFig(file.path(plotDir,'Fig1_2nAKLiver_localImpact_DEGs'),plotFun_localImpact,rawData=dd,width = 3.7,height = 3.7,res = 500)
saveFig(file.path(plotDir,'Fig1_2nAKLiver_GlobalImpact_DEGs'),plotFun_localImpact,rawData=dd,width = 3.7,height = 3.7,res = 500)








allDEGs.sub = allDEGs[allDEGs$padj < 0.05 & allDEGs$ct != 'doublets',]
all_nDEG_summary = data.frame()
allDEGs.sub.filtered = data.frame()
ncell_perCT = table(fLiver$finalAnn[fLiver$Genotype == 'diploid'])
ncell_perCT_threshold = ncell_perCT * 0.1

for(celltype in unique(fLiver$finalAnn)){
  if(!celltype %in% unique(allDEGs.sub$ct)){next}
  # Genes which are expressed in 2n cells are,
  # must be expressed in at least 10% of the cells?
  
  #genesExpressed = nCell_perGene_perCTGeno[nCell_perGene_perCTGeno$ct == celltype & nCell_perGene_perCTGeno$nCell >= 5,]
  genesExpressed = nCell_perGene_perCTGeno[nCell_perGene_perCTGeno$ct == celltype & nCell_perGene_perCTGeno$nCell >= ncell_perCT_threshold[celltype],]
  
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

## However, no genes got removed
nrow(allDEGs.sub.filtered) == nrow(allDEGs.sub)
dim(allDEGs.sub.filtered)
dim(allDEGs.sub)
write.csv(allDEGs.sub.filtered,file.path(resultDir,paste0(tissue,'_allDEGs_filtered.csv')))


##----- Calculate number of cells used in pbDESeq2 for each comparison (per Geno per CT)
# nCell used varies because of removal of different individuals depending on the context of the comparison

if(file.exists(file.path(outDir,'liver_nCell_perCT_perGeno_used_pbDEGs.csv'))){
  nCell_perCT_perGeno = read.csv(file.path(outDir,'liver_nCell_perCT_perGeno_used_pbDEGs.csv'),row.names = 1)
  nCell_perCT_perGeno = nCell_perCT_perGeno[,colnames(nCell_perCT_perGeno) != 'X']
  nCell_perCT_perGeno$geno[nCell_perCT_perGeno$geno == 'complete_trisomy'] = 'Triploid'
}
all_nDEG_summary$geno[all_nDEG_summary$geno == '3n'] = 'Triploid'
all_nDEG_summary = merge(all_nDEG_summary,nCell_perCT_perGeno,by=c('geno','ct'),all.x = T)
all_nDEG_summary$nCell = all_nDEG_summary$nCell_2n + all_nDEG_summary$nCell_AK
all_nDEG_summary$nCell_AKto2n = all_nDEG_summary$nCell_AK/all_nDEG_summary$nCell_2n

write.csv(all_nDEG_summary,file.path(resultDir,paste0(tissue,'_nDEGs_perCTperGeno_summary.csv')))


all_nDEG_summary = read.csv(paste0(tissue,'_nDEGs_perCTperGeno_summary.csv'))



fig1f_2nAK_Liv_fracDEG = function(){
  
  geno_cols = c('diploid' = grey(0.7),
                'T21' = '#93221E',
                'T18' = '#3d5dad',
                'T22' = '#679551',
                'T13' = '#526691',
                'MX' = '#b18db8',
                'Triploid' = '#e07d26')
  
  if(file.exists(file.path(plotDir,'Fig2F_2nAKLiver_fractionExpressedGenesDE_log10yScale_withChrX_rawData.tsv'))){
    data = read.delim(file.path(plotDir,'Fig2F_2nAKLiver_fractionExpressedGenesDE_log10yScale_withChrX_rawData.tsv'),sep = '\t')
  }else{
    ## Read in DESeq2 summary
    data = read.csv(file.path(resultDir,paste0(tissue,'_nDEGs_perCTperGeno_summary.csv')))  
    data = data[data$ct !='Cholangiocytes',]
    data = data[,c('geno','ct','total_nDEG','nGeneExpressed','frac_total_nDEG','nCell_2n','nCell_AK')]
    data = unique(data)
    # data$group = ifelse(data$ct %in% c('HSC_MPP'),'HSC_MPP',
    #                     ifelse(data$ct %in% c('MEMP_MEP','MK','Mast.cell','EE','ME','LE'),'Meg/Ery/Mast',
    #                            ifelse(data$ct %in% c('CMP_GMP','Monocyte','Macrophage','Kupffer.cell','myelocyte','Myelocyte','DC2','pDC','myeloid.prog'),'Myeloid',
    #                                   ifelse(data$ct %in% c('B.cell.prog','B.cell','NK.T'),'Lymphoid',
    #                                          ifelse(data$ct %in% c('Hepatocyte','Fibroblast','Endo','Cholangiocytes'),'Stromal','others')))))
    
    data$group = ifelse(data$ct %in% c('HSC_MPP','MEMP_MEP','CMP_GMP','B.cell.prog'),'Progenitors',
                        ifelse(data$ct %in% c('MK','Mast.cell','EE','ME','LE'),'Meg/Ery/Mast',
                               ifelse(data$ct %in% c('Monocyte','Macrophage','Kupffer.cell','myelocyte','Myelocyte','DC2','pDC','myeloid.prog'),'Myeloid',
                                      ifelse(data$ct %in% c('B.cell','NK.T'),'Lymphoid',
                                             ifelse(data$ct %in% c('Hepatocyte','Fibroblast','Endo','Cholangiocytes'),'Stromal','others')))))
    #data = data[data$group != 'Stromal',]
    # p1 = ggplot(data,aes(log10(nCell),frac_total_nDEG))+
    #   geom_smooth(method = 'lm',se = F,col = 'grey',lwd=0.4) +
    #   #geom_point(aes(shape=geno,col=ct),size=3)+
    #   geom_point(aes(col=ct,size=nCell_AKto2n))+
    #   scale_size_continuous(range = c(0.5,5))+
    #   facet_wrap(vars(geno))+
    #   scale_color_manual(values = sample(col25))+
    #   theme_bw() + theme(panel.grid = element_blank())
    # 
    # data$geno = factor(data$geno,c('T21','T18','T22','MX','Triploid'))
    # 
    # cell_cols = c('MEMP_MEP' = '#d10808',
    #               "B.cell.prog" = '#279be3',"B.cell" = '#12436e','MK'='#8d439c','Hepatocyte'=pal37H[c(26)])
    # a = rep(grey(0.8),n_distinct(data$ct)-length(cell_cols))
    # names(a) = unique(data$ct[! data$ct %in% names(cell_cols)])
    # cell_cols = c(cell_cols,a)
    
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
    
    
    # data$ct = factor(as.character(data$ct),c('HSC_MPP',
    #                                          'MEMP_MEP','MK','Mast.cell','EE','ME','LE',
    #                                          'myeloid.prog','CMP_GMP','Monocyte','Macrophage','Kupffer.cell','myelocyte','Myelocyte','DC2','pDC',
    #                                          'Endo','Fibroblast','Hepatocyte','Cholangiocytes',
    #                                          'B.cell.prog','B.cell',"NK.T"
    # ))
    
    data$ct = factor(as.character(data$ct),c('HSC_MPP','MEMP_MEP','CMP_GMP','B.cell.prog',
                                             'MK','Mast.cell','EE','ME','LE',
                                             'Monocyte','Macrophage','Kupffer.cell','myelocyte','Myelocyte','DC2','pDC',
                                             'Endo','Fibroblast','Hepatocyte','Cholangiocytes',
                                             'B.cell',"NK.T"
    ))
    
    #data$group = factor(data$group,c('HSC_MPP','Meg/Ery/Mast','Myeloid','Lymphoid','Stromal'))
    data$group = factor(data$group,c('Progenitors','Meg/Ery/Mast','Myeloid','Lymphoid','Stromal'))
    
    data2 = data %>% group_by(group,ct) %>% summarise(max_frac = max(frac_total_nDEG))
    data$id = 1:nrow(data)
    overlapped_data = data[data$ct %in% c('LE') & data$geno %in% c('Triploid','T21','MX','T18') |
                             data$ct %in% c('Kupffer.cell') & data$geno %in% c('Triploid','T21','T18'),]
    
    
    p3 = ggplot(data,aes(ct,log10(frac_total_nDEG)))+
      geom_segment(data=data2,aes(x = ct,xend=ct,y=min(log10(data$frac_total_nDEG)),yend=log10(max_frac)),col=grey(0.7))+
      facet_grid(.~group,scales = 'free_x',space = 'free_x')+
      geom_quasirandom(data = data[data$id %in% overlapped_data$id,],aes(col=geno,size=log10(nCell_AK)),alpha=1,width = 0.4)+
      geom_quasirandom(data = data[data$id %in% overlapped_data$id,],aes(size=log10(nCell_AK)),shape=1,col='black',alpha=1,width = 0.4)+
      geom_quasirandom(data = data[(!data$id %in% overlapped_data$id) & data$geno =='T21',],aes(col=geno,size=log10(nCell_AK)),alpha=1,width = 0)+
      geom_quasirandom(data = data[(!data$id %in% overlapped_data$id) & data$geno =='T21',],aes(size=log10(nCell_AK)),shape=1,col='black',alpha=1,width = 0)+
      geom_quasirandom(data = data[(!data$id %in% overlapped_data$id) & data$geno !='T21',],aes(col=geno,size=log10(nCell_AK)),alpha=1,width = 0)+
      geom_quasirandom(data = data[(!data$id %in% overlapped_data$id) & data$geno !='T21',],aes(size=log10(nCell_AK)),shape=1,col='black',alpha=1,width = 0)+
      scale_y_continuous(labels = c('0.1','0.01','0.001','0.0001'),breaks = c(-1,-2,-3,-4))+
      scale_color_manual(values = geno_cols)+
      theme_classic(base_size = 15)+
      theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
            axis.text = element_text(colour = 'black'),
            axis.ticks = element_line(colour = 'black'),
            panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.7),
            strip.background = element_blank(),
            axis.line = element_blank()) + 
      xlab('') + ylab('Fraction of DEGs') 
    print(p3)
  }
  
  
  #saveFig(file.path(plotDir,'Fig2F_2nAKLiver_fractionExpressedGenesDE_log10yScale'),plotFun_degFrac_perCT_perGeno,rawData=data,width = 11,height = 5.5,res = 500,useDingbats = F)
  saveFig(file.path(plotDir,'Fig2F_2nAKLiver_fractionExpressedGenesDE_log10yScale_withChrX'),plotFun_degFrac_perCT_perGeno,rawData=data,width = 10,height = 5.5,res = 500,useDingbats = F)
  saveFig(file.path(plotDir,'Fig2F_2nAKLiver_fractionExpressedGenesDE_log10yScale_withChrX_noStromal_v2'),plotFun_degFrac_perCT_perGeno,rawData=data,width = 8.5,height = 5.5,res = 500,useDingbats = F)
  saveFig(file.path(plotDir,'Fig2F_2nAKLiver_fractionExpressedGenesDE_log10yScale_withChrX_wStromal_v2'),plotFun_degFrac_perCT_perGeno,rawData=data,width = 9.5,height = 5.5,res = 500,useDingbats = F)
  #saveFig(file.path(plotDir,'Fig1f_2nAKLiver_fractionExpressedGenesDE_v3'),plotFun_degFrac_perCT_perGeno,rawData=data,width = 11,height = 5.5,res = 500,useDingbats = T)
}



## Frequency of shared genes
nComparison_perGene = allDEGs.sub.filtered %>% group_by(ct,direction,geneSym) %>% mutate(nGeno = n_distinct(geno),
                                                                                         genoShared = paste(unique(geno),collapse = '_'))
table(nComparison_perGene$genoShared[nComparison_perGene$geno == 'T21'],nComparison_perGene$nGeno[nComparison_perGene$geno == 'T21'])

data = allDEGs.sub.filtered %>% group_by(ct,direction,geneSym) %>% mutate(nGeno = n_distinct(geno)) %>% 
  group_by(ct,direction,geno,nGeno) %>% summarise(nGene = n_distinct(geneSym))
#data$nGene[data$direction == 'AK_down'] = -data$nGene[data$direction == 'AK_down']
data$geno[data$geno == '3n'] = 'Triploid'
ggplot(data,aes(geno,y=nGene,fill=direction))+
  geom_boxplot(outlier.shape = NA)+
  #scale_y_log10()+
  geom_point(aes(col=ct))+
  scale_fill_manual(values = geno_cols)+
  scale_color_manual(values = col25)+
  facet_wrap(vars(nGeno),scales = 'free_y')+
  theme_classic()

dd = data[data$geno == 'T21',]
dd$direction = factor(dd$direction,c('AK_up','AK_down'))
plotFun_sharedDEGs_acrossGeno = function(noFrame=FALSE,noPlot=FALSE){
  p = ggplot(dd,aes(as.character(nGeno),y=nGene,fill=direction))+
    geom_boxplot(outlier.shape = NA,width=0.7,color='black')+
    #scale_y_log10()+
    #geom_point(aes(col=ct))+
    geom_quasirandom(width = 0.2,size=0.8)+
    scale_fill_manual(values = col25[c(2,1)])+
    scale_color_manual(values = col25)+
    facet_wrap(vars(direction),scales = 'free_y')+
    theme_classic(base_size = 14)+
    theme(panel.background = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),
          axis.text = element_text(color = 'black'),
          axis.ticks = element_line(colour = 'black'),
          strip.background = element_blank()) + 
    xlab('Number of Genotypes') + ylab('Number of DEG')
  print(p)
}

saveFig(file.path(plotDir,'FigSup2X_nSharedDEGs_acrossGeno'),plotFun_sharedDEGs_acrossGeno,rawData=dd,width = 6,height = 3.5,res = 500)






##-------------------------------------------------------------------##
##    4. Heatmap of log2FC of each gene along genomic position     ####
##-------------------------------------------------------------------##
# Import log2FC output
resultDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jul24/without_published2n/genoAssay_withX'
allGenes = import_pbDEGresults(outDir = resultDir, allGenes = T,tissue='liver')
allGenes$geneSym = geneMap$geneSym[match(allGenes$ensID,geneMap$ensID)]
allGenes$chr = geneMap$chr[match(allGenes$ensID,geneMap$ensID)]
# allGenes$xpos = coords$xpos_linearised[match(allGenes$ensID,coords$gene_id)]
# allGenes$xpos_byChr = coords$TSS[match(allGenes$ensID,coords$gene_id)]
# allGenes$isDE = ifelse(is.na(allGenes$padj) | allGenes$padj > 0.05,F,T) # To extract real DE, I also include criteria of being expressed in >= 10% cells in either group


# a = allGenes %>% filter(!is.na(log2FoldChange)) %>% group_by(ct,geno,chr) %>% summarise(med = median(log2FoldChange))
# ggplot(allGenes[allGenes$ct == 'Endo' & allGenes$chr %in% c('18','21','22','X'),],aes(y=log2FoldChange))+
#   geom_point(aes(x=xpos_byChr,col=isDE),size=0.1,alpha=0.5)+
#   scale_color_manual(values = c('grey','red'))+
#   facet_wrap(vars(geno),ncol = 1)+
#   facet_grid(vars(geno),vars(chr),scales = 'free_x',space='free_x')+
#   #geom_vline(xintercept = chrLens_pos[names(chrLens_pos) %in% c('18','21','22','X')],col='grey',lty=2,size=0.3,alpha=0.5) +
#   theme_classic()+ylim(-2,2)
# 
# df = allGenes[allGenes$ct == 'Hepatocyte',]
# df$geno[df$geno == 'complete'] = 'Triploid'
# df$geno = factor(df$geno,c('T18','T21','T22','MX','Triploid'))
# df$chr2 = ifelse(df$chr %in% c('18','21','22','X'),df$chr,'elsewhere')
# df$chr2 = factor(df$chr2,c('18','21','22','X','elsewhere'))
# 
# df$chr3 = ifelse(as.character(df$chr) %in% as.character(c(1:17)),'1-17',
#                  ifelse(as.character(df$chr) %in% as.character(c(19,20)),'19-20', 
#                         as.character(df$chr)))
# table(df$chr3)
# df$chr3 = factor(df$chr3,c('chr 1 - chr 17','18','chr 19 - chr 20','21','22','X'))
# df$chr3 = factor(df$chr3,c('1-17','18','19-20','21','22','X'))
# df$chr = factor(df$chr,levels = c(as.character(seq(1:22)),'X'))
# 
# 
# ggplot(df,aes(x=chr,y=log2FoldChange))+
#   geom_hline(yintercept = 0,col='grey',lty=2,linewidth=0.3,alpha=1) +
#   #geom_point(aes(col=isDE),size=0.1,alpha=0.5)+
#   geom_jitter(aes(col=isDE),size=0.01,alpha=0.5,width = 0.3)+
#   scale_color_manual(values = c('grey','red'))+
#   facet_wrap(vars(geno),ncol = 1)+
#   facet_grid(vars(geno),vars(chr3),scales = 'free_x',space = 'free_x')+
#   #ylim(-4,4) + 
#   #geom_vline(xintercept = chrLens_pos,col='black',lty=1,size=0.3,alpha=1) +
#   scale_y_continuous(breaks=seq(-4,4,by=4),limits = c(-3,3)) +
#   theme_bw(base_size = 13)+theme(panel.grid = element_blank(),
#                    #axis.text.x = element_blank(),
#                    #axis.ticks.x = element_blank(),
#                    axis.text.x = element_text(size=8),
#                    axis.text.y = element_text(size=8)) + #xlab('Genomic position') + 
#   xlab('Chromosome')+ 
#   ylab('log2 expression fold change (AK / diploid)')
# 



##---- Heatmap of median log2FC per chromosome per genotype
allGenes$geno[allGenes$geno == 'complete'] = 'Triploid'
allGenes$geno = factor(allGenes$geno,c('T18','T21','T22','MX','Triploid','MX.withChrX'))
allGenes$chr = factor(allGenes$chr,levels = c(as.character(seq(1:22)),'X'))
allGenes = allGenes[allGenes$ct != 'doublets',]

## Re-write cell type names
allGenes$ct[allGenes$ct == 'EE'] = 'early Ery'
allGenes$ct[allGenes$ct == 'ME'] = 'mid Ery'
allGenes$ct[allGenes$ct == 'LE'] = 'late Ery'
allGenes$ct[allGenes$ct == 'Endo'] = 'Endothelium'
allGenes$ct[allGenes$ct == 'Mesothelial_cells'] = 'Mesothelial cell'
allGenes$ct[allGenes$ct == 'Cholangiocytes'] = 'Cholangiocyte'
allGenes$ct[allGenes$ct == 'NK.T'] = 'NK / T'

allGenes$ct = gsub('\\.',' ',allGenes$ct)
allGenes$ct = gsub('_',' / ',allGenes$ct)

write.csv(allGenes,file.path(resultDir,paste0(tissue,'_log2FC_allGenes.csv')))

#allGenes = read.csv(file.path(resultDir,paste0(tissue,'_log2FC_allGenes.csv')))


data = allGenes %>% filter(ct != 'doublets') %>% group_by(chr,geno,ct) %>% summarise(median_log2FC = mean(log2FoldChange,na.rm = T))
write.csv(data,file.path(resultDir,paste0(tissue,'_median_log2FC_perCTperGeno.csv')))

tissue = 'liver'
data = read.csv(file.path(resultDir,paste0(tissue,'_median_log2FC_perCTperGeno.csv')))
data = pivot_wider(data,id_cols = c('geno','ct'),names_from = 'chr',values_from = 'median_log2FC')

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.15), c('#255fa8','white','#b8393d'))
#col_fun = colorRamp2(c(-0.1, 0, 1), c('#255fa8','white','#b8393d'))
celltype_order = c('HSC / MPP','MEMP / MEP','early Ery','mid Ery','late Ery','Mast cell','MK',
                   'B cell prog','B cell',
                   'CMP / GMP','Monocyte','Macrophage','Kupffer cell',
                   'DC2','pDC','Myelocyte','NK / T',
                   'Endothelium','Fibroblast','Hepatocyte','Cholangiocyte','Mesothelial cell')


##---   Version 1: 1 heatmap for each Genotype, rows = different celltypes, col = chromosomes  ------##
for(g in unique(data$geno)){
  d = data[data$geno == g,]
  d$ct = factor(d$ct,levels = celltype_order)
  mtx = d[,!colnames(d) %in% c('geno','ct')] %>% as.data.frame()
  rownames(mtx) = paste0(d$ct)
  
  hm = Heatmap(as.matrix(mtx), name = paste0('Average log2FC - ',g),border = F,
               column_split = factor(colnames(mtx),levels = c(1:22,'X')),
               row_split = d$ct,
               show_row_names = T,show_column_names = F,
               column_title_gp = gpar(fontsize=12),
               row_title_gp = gpar(fontsize=0),
               row_names_gp = gpar(fontsize=12),
               row_gap = unit(0.15,'cm'),
               column_gap = unit(0.15,'cm'),
               col = col_fun,
               column_names_rot = 0,row_names_side = 'left',
               cluster_rows = F,cluster_columns = F)
  
  plotFun_log2FC_hm = function(noFrame=FALSE,noPlot=FALSE){
    draw(hm)
  }
  
  if(g == 'T18'){
    h = 5.3
  }else if(g == 'T21'){
    h = 6.0
  }else if(g == 'T22'){
    h = 3.1
  }else if(g == 'MX'){
    h = 5.5
  }else if(g == 'Triploid'){
    h = 4.5
  }else if(g == 'MX.withChrX'){
    h = 5.3
  }
  
  #w = 8.4
  w = 8.55
  #saveFig(file.path(plotDir,paste0('FigSuppXX_2nAK.fLiver_',g,'_meanLog2FC_hm')),plotFun_log2FC_hm,rawData=d,width = w,height = h,res = 500,useDingbats = T)
  saveFig(file.path(plotDir,paste0('FigSuppXX_2nAK.fLiver_',g,'_meanLog2FC_withX_hm')),plotFun_log2FC_hm,rawData=d,width = w,height = h,res = 500,useDingbats = T)
  
  
  
  # pdf(file.path(plotDir,paste0(tissue,'_',g,'_','meanLog2FC_hm.pdf')),width = 8.5,height = 4.85)  
  # draw(hm)
  # dev.off()
  
  # ## Corrected log2FC for Triploid case
  # if(g %in% c('Triploid','3n')){
  #   mtx2 = mtx - mtx[,'X']
  #   
  #   hm = Heatmap(as.matrix(mtx2),name = paste0('Average log2FC - ',g),border = F,
  #                column_split = factor(colnames(mtx),levels = c(1:22,'X')),
  #                row_split = d$ct,
  #                show_row_names = T,show_column_names = F,
  #                column_title_gp = gpar(fontsize=10),
  #                row_title_gp = gpar(fontsize=0),
  #                row_names_gp = gpar(fontsize=10),
  #                row_gap = unit(0.3,'cm'),
  #                col = col_fun,
  #                column_names_rot = 0,row_names_side = 'left',
  #                cluster_rows = F,cluster_columns = F)
  #   
  #   pdf(file.path(plotDir,paste0(tissue,'_',g,'_','corrected_meanLog2FC_hm.pdf')),width = 8,height = 5.2)
  #   draw(hm)
  #   dev.off()
  # }
}















##---   Version 2: Heatmap of a sub-selected set of celltypes, for all Geno  ------##
#ct_toKeep = c('MEMP_MEP','MK','B.cell.prog','Hepatocyte')
ct_toKeep = c('MEMP_MEP','B.cell.prog')

data = data[data$ct %in% ct_toKeep,]
data$ct = factor(data$ct,levels = ct_toKeep)
mtx = data[,!colnames(data) %in% c('geno','ct')] %>% as.data.frame()
rownames(mtx) = paste0(data$ct,':',data$geno)
col_fun = colorRamp2(c(-0.1, 0, 0.2), c('#255fa8','white','#b8393d'))
hm=Heatmap(as.matrix(mtx),name = paste0('Average log2FC'),border = F,
           column_split = factor(colnames(mtx),levels = c(1:22,'X')),
           row_split = data$ct,
           show_row_names = T,show_column_names = F,
           column_title_gp = gpar(fontsize=11),
           row_title_gp = gpar(fontsize=0),
           row_names_gp = gpar(fontsize=10),
           row_gap = unit(0.4,'cm'),
           col = col_fun,
           column_names_rot = 0,row_names_side = 'left',
           cluster_rows = F,cluster_columns = F)

if(length(ct_toKeep) == 4){
  #pdf(file.path(plotDir,paste0(tissue,'_allGeno_subsetCT_meanLog2FC_hm.pdf')),width = 9,height = 5.2)
  pdf(file.path(plotDir,paste0(tissue,'_allGeno_subsetCT_meanLog2FC_hm_withChrX.pdf')),width = 9.3,height = 5.2)  
}else if(length(ct_toKeep) == 2){
  #pdf(file.path(plotDir,paste0(tissue,'_allGeno_subsetCT_meanLog2FC_hm.pdf')),width = 8,height = 2.35)  
  pdf(file.path(plotDir,paste0(tissue,'_allGeno_subsetCT_meanLog2FC_hm_withChrX.pdf')),width = 8.3,height = 2.35)  
}

draw(hm)
dev.off()












##------------------------##
##    T21_HSC DEGs      ####
##------------------------##
fLiver$finalAnn_broad = fLiver$annot_mar24
geno_cols = c('Diploid' = grey(0.7),
              'T21' = '#93221E',
              'T18' = '#3d5dad',
              'T22' = '#679551',
              'T13' = '#526691',
              'MX' = '#b18db8',
              'Triploid' = '#e07d26')
allDEGs.sub.filtered = read.csv(file.path(outDir,paste0(tissue,'_allDEGs_filtered.csv')))
hscDEGs = allDEGs.sub.filtered[allDEGs.sub.filtered$ct == 'HSC_MPP',]
hscDEGs.summary = hscDEGs %>% group_by(geneSym,chr,direction) %>% summarise(nGeno=n_distinct(geno))
hscDEGs.T21 = hscDEGs[hscDEGs$geno == 'T21',]
hscDEGs.T21 = hscDEGs.T21[order(abs(hscDEGs.T21$log2FoldChange),decreasing = T),]


figSupp_fLiver_HSC.MPP_UMAP = function(){
  
  hsc = subset(fLiver,subset = cellID %in% fLiver$cellID[fLiver$finalAnn_broad == 'HSC_MPP'])
  hsc = standard_clustering(hsc)
  hsc$Genotype = as.character(hsc$Genotype)
  hsc$Genotype[hsc$Genotype == 'complete_trisomy'] = 'Triploid'
  hsc$Genotype[hsc$Genotype == 'diploid'] = 'Diploid'
  hsc$Genotype = factor(hsc$Genotype,c('Diploid','T21','Triploid','T18','T22','MX'))
  
  
  hsc_mdat = hsc@meta.data
  hsc_mdat = cbind(hsc_mdat,hsc@reductions$umap@cell.embeddings)
  
  
  plotFun_fLiver_HSC.MPP_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    p1 = ggplot(hsc_mdat[hsc_mdat$Genotype %in% c('Diploid','T21'),],aes(UMAP_1,UMAP_2,col=Genotype))+
      geom_point(size=0.6)+
      geom_point(data = hsc_mdat[!hsc_mdat$Genotype %in% c('Diploid','T21'),],size=0.6)+
      scale_color_manual(values = geno_cols)+
      theme_classic()+theme(panel.border = element_rect(fill=F),
                            axis.line = element_blank(),axis.ticks = element_blank(),
                            axis.text = element_blank()) + ggtitle('HSC_MPP')
    print(p1)  
  }
  
  saveFig(file.path(plotDir,'Supp.Figxx_fLiverHSC_UMAP'),plotFun_fLiver_HSC.MPP_UMAP,rawData=hsc_mdat,width = 3.6,height = 3,res = 500,useDingbats = T)
  
  
  
  
  plotFun_fLiver_HSC.MPP_DEG_dotplot = function(noFrame=FALSE,noPlot=FALSE){
    p = DotPlot(hsc,group.by = 'Genotype',
                cols = c(colAlpha(grey(0.95),0.8),'black'),
                # features = c(hscDEGs$geneSym[hscDEGs$geno == 'T21' & hscDEGs$direction == 'AK_up'],
                #              hscDEGs$geneSym[hscDEGs$geno == 'T21' & hscDEGs$direction == 'AK_down'])
                features = c(hscDEGs.T21$geneSym[hscDEGs.T21$direction == 'AK_up'],
                             hscDEGs.T21$geneSym[hscDEGs.T21$direction == 'AK_down'])
    ) + RotatedAxis()+
      theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                       face=ifelse(hscDEGs.T21$chr == '21','bold','plain'),
                                       colour =ifelse(hscDEGs.T21$isTF == T,'purple',
                                                      ifelse(hscDEGs.T21$isCSM == T,'blue','black'))),
            axis.text.y = element_text(size=11),
            legend.position = 'top',legend.text = element_text(size=9),legend.title = element_text(size=8.5),legend.key.size = unit(0.6,'cm')) + xlab('') + ylab('') 
    
    print(p)
  }
  saveFig(file.path(plotDir,'Supp.Figxx_fLiverHSC_deg_dotPlot'),plotFun_fLiver_HSC.MPP_DEG_dotplot,width = 9,height = 4,res = 500,useDingbats = T)
}







##-------------------------##
##    T21_MEMP DEGs      ####
##-------------------------##
allDEGs.sub.filtered = read.csv(file.path(outDir,paste0(tissue,'_allDEGs_filtered.csv')))
mempDEGs = allDEGs.sub.filtered[allDEGs.sub.filtered$ct == 'MEMP_MEP',] %>% group_by(geneSym,chr,direction) %>% 
  mutate(nGeno=n_distinct(geno),mean_cellFrac = median(max_pct),pct.diff = cellFrac_g1 - cellFrac_g2)

mempDEGs.summary = mempDEGs %>% group_by(geneSym,chr,direction) %>% summarise(nGeno=n_distinct(geno))
mempDEGs.T21 = mempDEGs[mempDEGs$geno == 'T21',]
mempDEGs.T21 = mempDEGs.T21[order(abs(mempDEGs.T21$log2FoldChange),decreasing = T),]


## Perform enrichR on all DEGs
upDEG_enriched <- enrichr(unique(mempDEGs$geneSym[mempDEGs$nGeno >=3 & mempDEGs$direction == 'AK_up']), dbs)
upDEG_enriched <- enrichr(unique(mempDEGs$geneSym[mempDEGs$geno == 'T21' & mempDEGs$direction == 'AK_up']), dbs)
downDEG_enriched <- enrichr(unique(mempDEGs$geneSym[mempDEGs$geno == 'T21' & mempDEGs$direction == 'AK_down']), dbs)

i=5
plotEnrich(upDEG_enriched[[i]][upDEG_enriched[[i]]$Adjusted.P.value <  0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
plotEnrich(downDEG_enriched[[i]][downDEG_enriched[[i]]$Adjusted.P.value <  0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
upDEG_enriched[[i]][upDEG_enriched[[i]]$Adjusted.P.value < 0.05,][1,]
downDEG_enriched[[i]][downDEG_enriched[[i]]$Adjusted.P.value <  0.05,][1,]


fLiver$group = paste0(fLiver$finalAnn,'_',fLiver$Genotype)
Idents(fLiver) = fLiver$Genotype
Idents(fLiver) = fLiver$finalAnn
fLiver$finalAnn = factor(fLiver$finalAnn,levels = c("HSC_MPP","MEMP_MEP","CMP_GMP","LMPP_ELP", # brown
                                                    "earlyMK",'MK', # purple
                                                    'Mast.cell', # slamon
                                                    'EE','ME','LE', # red
                                                    "proMono","Monocyte","Macrophage","Kupffer.cell", # blue
                                                    "DC1","DC2","pDC", # yellow
                                                    "promyelocyte","myelocyte",'Myelocyte', # orange
                                                    'B.cell.prog',"pro.B.cell","pre.B.cell","B.cell", # pink / purple
                                                    "ILC.precursor","T.cell","NK.T", # greys
                                                    "Hepatocyte","Fibroblast","Endo","Mesenchyme","Neuron",'Cholangiocytes','Mesothelial_cells','Trophoblast')) # green
DotPlot(fLiver,group.by = 'group',
        idents = c('diploid','T21'),
        #features = c('IL7','HOXA9','NFKB2','NFKB1','TNFA','IL1B','TRAF3','BCL2L11','BAX')
        features = dd$geneSym[dd$B.cell.prog > 30 & dd$MEMP_MEP < 20 & dd$geneSym %in% allDEGs$geneSym[allDEGs$geno == 'T21' & allDEGs$ct == 'B.cell.prog' & allDEGs$direction == 'AK_down']]
        #features = allDEGs$geneSym[allDEGs$ct == 'B.cell.prog' & allDEGs$geno == 'T21' & allDEGs$direction == 'AK_up'][101:200]
) + RotatedAxis()+
  theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1,size=8),legend.position = 'top')+xlab('')
#idents = c('B.cell','B.cell.prog','HSC_MPP','MEMP_MEP','MK','EE'),
#features = c('PAXBP1','BACH1','IFNAR2','IFNAR1','IFNGR2','IL10RB','PLXNC1','CLEC11A','SLC44A1')) 

tmp = data.frame(nCellExpr_ELP = apply(fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$Genotype == 'diploid' & fLiver$annot_aug24 == 'LMPP_ELP']],1,function(x){sum(x>0)}),
                 nCellExpr_pro = apply(fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$Genotype == 'diploid' & fLiver$annot_aug24 == 'pro.B.cell']],1,function(x){sum(x>0)}),
                 nCellExpr_pre = apply(fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$Genotype == 'diploid' & fLiver$annot_aug24 == 'pre.B.cell']],1,function(x){sum(x>0)}),
                 percCell_ELP = apply(fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$Genotype == 'diploid' & fLiver$annot_aug24 == 'LMPP_ELP']],1,function(x){100*sum(x>0)/length(x)}),
                 percCell_pro = apply(fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$Genotype == 'diploid' & fLiver$annot_aug24 == 'pro.B.cell']],1,function(x){100*sum(x>0)/length(x)}),
                 percCell_pre = apply(fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$Genotype == 'diploid' & fLiver$annot_aug24 == 'pre.B.cell']],1,function(x){100*sum(x>0)/length(x)}))


dd = nCell_perGene_perCTGeno[nCell_perGene_perCTGeno$ct %in% c('HSC_MPP','B.cell.prog','B.cell','MEMP_MEP','CMP_GMP','EE','Monocyte','DC2','MK'),]
dd = pivot_wider(dd,id_cols = c(geneSym,chr),names_from = 'ct',values_from = 'pct_epxressed_0')
dd = cbind(dd,tmp[match(dd$geneSym,rownames(tmp)),])
ggplot(dd,aes(B.cell.prog,MEMP_MEP))+
  geom_point(size=0.1,col=grey(0.8))+
  geom_abline()+
  geom_point(data = dd[dd$geneSym %in% allDEGs$geneSym[allDEGs$geno == 'T21' & 
                                                         allDEGs$ct == 'B.cell.prog' & 
                                                         allDEGs$direction == 'AK_down'],],size=0.7,col='red')+
  theme_classic()

nCell_perGene_perCTGeno$chr = geneMap$chr[match(nCell_perGene_perCTGeno$geneSym,geneMap$geneSym)]
ggplot(nCell_perGene_perCTGeno[nCell_perGene_perCTGeno$chr == '21',],aes(reorder(ct,pct_epxressed_0,median),pct_epxressed_0))+
  geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))


figSupp_fLiver_MEMP.MEP_UMAP = function(){
  
  s = subset(fLiver,subset = cellID %in% fLiver$cellID[fLiver$finalAnn_broad == 'MEMP_MEP'])
  s = standard_clustering(s)
  s$Genotype = as.character(s$Genotype)
  s$Genotype[s$Genotype == 'complete_trisomy'] = 'Triploid'
  s$Genotype[s$Genotype == 'diploid'] = 'Diploid'
  s$Genotype = factor(s$Genotype,c('Diploid','T21','Triploid','T18','T22','MX'))
  
  
  
  s_mdat = s@meta.data
  s_mdat = cbind(s_mdat,s@reductions$umap@cell.embeddings)
  
  plotFun_fLiver_MEMP.MEP_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    if(!noPlot & !noFrame){
      p1 = ggplot(s_mdat,aes(UMAP_1,UMAP_2,col=Genotype))+
        geom_point(size=0.6)+
        scale_color_manual(values = geno_cols)+
        theme_classic(base_size = 9)+theme(panel.border = element_rect(fill=F),
                                           axis.line = element_blank(),axis.ticks = element_blank(),
                                           axis.text = element_blank()) + ggtitle('MEMP.MEP')
      print(p1)    
    }
    
    if(noFrame){
      p1 = ggplot(s_mdat,aes(UMAP_1,UMAP_2,col=Genotype))+
        geom_point(size=0.6)+
        scale_color_manual(values = geno_cols)+
        theme_classic(base_size = 9)+theme(axis.line = element_blank(),axis.ticks = element_blank(),
                                           axis.text = element_blank())+
        xlab('')+ylab('')
      print(p1)    
    }
    
  }
  
  saveFig(file.path(plotDir,'Supp.Figxx_fLiverMEMP_UMAP'),plotFun_fLiver_MEMP.MEP_UMAP,rawData=s_mdat,width = 4,height = 3,res = 500,useDingbats = T)
  
  
  
  plotFun_fLiver_MEMP.MEP_DEG_dotplot = function(noFrame=FALSE,noPlot=FALSE){
    p = DotPlot(s,group.by = 'Genotype',
                cols = c(colAlpha(grey(0.95),0.8),'black'),
                # features = c(mempDEGs$geneSym[mempDEGs$geno == 'T21' & mempDEGs$direction == 'AK_up'][1:30],
                #              mempDEGs$geneSym[mempDEGs$geno == 'T21' & mempDEGs$direction == 'AK_down'])
                features = c(mempDEGs.T21$geneSym[mempDEGs.T21$direction == 'AK_up'],
                             mempDEGs.T21$geneSym[mempDEGs.T21$direction == 'AK_down'])
    ) + RotatedAxis()+
      theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7,
                                       face=ifelse(mempDEGs.T21$chr == '21','bold','plain'),
                                       colour =ifelse(mempDEGs.T21$isTF == T,'purple',
                                                      ifelse(mempDEGs.T21$isCSM == T,'blue','black'))),
            axis.text.y = element_text(size=11),
            legend.position = 'top',legend.text = element_text(size=9),legend.title = element_text(size=8.5),legend.key.size = unit(0.6,'cm')) + xlab('') + ylab('') 
    
    print(p)
  }
  saveFig(file.path(plotDir,'Supp.Figxx_fLiverMEMP_deg_dotPlot'),plotFun_fLiver_MEMP.MEP_DEG_dotplot,width = 9,height = 4,res = 500,useDingbats = T)
}




##-----------------------------------##
##    B.cell progenitors DEGs      ####
##-----------------------------------##
s = subset(fLiver,subset = annot_aug24 %in% c('HSC_MPP','MEMP_MEP','B.cell','LMPP_ELP','pro.B.cell','pre.B.cell'))
s = standard_clustering(s,runHarmony = T,harmonyVar = 'assay')
DimPlot(s,group.by = 'donorID',cols = col25,label = T,label.box = T)
FeaturePlot(s,'MIA2')
DimPlot(s,cells.highlight = s$cellID[s$Genotype == 'T21'])


allDEGs.sub.filtered = read.csv(file.path(outDir,paste0(tissue,'_allDEGs_filtered.csv')))
bProgDEGs = allDEGs.sub.filtered[allDEGs.sub.filtered$ct == 'B.cell.prog',] %>% 
  group_by(geneSym,chr,direction,isTF,isCSM,isCosmic) %>% mutate(nGeno=n_distinct(geno),mean_cellFrac = mean(max_pct))

table(bProgDEGs.summary$nGeno)

fLiver$group = paste0(fLiver$finalAnn,'_',fLiver$Genotype)
Idents(fLiver) = fLiver$finalAnn
DotPlot(fLiver,group.by = 'group',idents = c('B.cell','B.cell.prog','HSC_MPP','Endo'),
        features = c('CD83','CDKN1A','SOD2','KLF6','CD69','CXCR4','ZFP36L1','TCF4','HMGB1','CD79A')) 
## Perform enrichR on all DEGs
upDEG_enriched <- enrichr(unique(bProgDEGs.summary$geneSym[bProgDEGs.summary$nGeno >=3 & bProgDEGs.summary$direction == 'AK_up']), dbs)
upDEG_enriched <- enrichr(unique(bProgDEGs$geneSym[bProgDEGs$geno == 'T21' & bProgDEGs$direction == 'AK_up']), dbs)
downDEG_enriched <- enrichr(unique(bProgDEGs$geneSym[bProgDEGs$geno == 'T21' & bProgDEGs$direction == 'AK_down']), dbs)

i=6
plotEnrich(upDEG_enriched[[i]][upDEG_enriched[[i]]$Adjusted.P.value <  0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
plotEnrich(downDEG_enriched[[i]][downDEG_enriched[[i]]$Adjusted.P.value <  0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
upDEG_enriched[[i]][upDEG_enriched[[i]]$Adjusted.P.value < 0.05,][1,]
downDEG_enriched[[i]][downDEG_enriched[[i]]$Adjusted.P.value <  0.05,][1,]


dd = as.data.frame(table(fLiver$Phase[fLiver$finalAnn == 'B.cell.prog'],fLiver$Genotype[fLiver$finalAnn == 'B.cell.prog']))
colnames(dd) = c('Phase','Geno','Freq')
dd$Phase = factor(dd$Phase,rev(c('G1','S','G2M')))
ggplot(dd,aes(Geno,Freq,fill=Phase))+
  geom_col(position= 'fill')

FeaturePlot(s,'TP53',split.by = 'Genotype',cells = s$cellID[s$Genotype %in% c('diploid','T21')])

##-------------------##
##    MK DEGs      ####
##-------------------##
allDEGs.sub.filtered = read.csv(file.path(outDir,paste0(tissue,'_allDEGs_filtered.csv')))
mkDEGs = allDEGs.sub.filtered[allDEGs.sub.filtered$ct == 'MK',]
mkDEGs = mkDEGs %>% group_by(geneSym,chr,direction) %>% mutate(nGeno=n_distinct(geno),
                                                               mean_cellFrac = median(max_pct))
table(mkDEGs.summary$nGeno)

fLiver$group = paste0(fLiver$finalAnn,'_',fLiver$Genotype)
Idents(fLiver) = fLiver$finalAnn
DotPlot(fLiver,group.by = 'group',idents = c('B.cell','B.cell.prog','HSC_MPP','Endo','MK'),
        features = c('PPP1R15A','NFKBIA','ZFP36','CD69','TNF','HSPA1B','RUNX1')) 

## Perform enrichR on all DEGs
upDEG_enriched <- enrichr(unique(mkDEGs.summary$geneSym[mkDEGs.summary$nGeno >=3 & mkDEGs.summary$direction == 'AK_down']), dbs)
upDEG_enriched <- enrichr(unique(mkDEGs$geneSym[mkDEGs$geno == 'T18' & mkDEGs$direction == 'AK_up']), dbs)
upDEG_enriched <- enrichr(unique(allDEGs$geneSym[allDEGs$geno == 'T18' & allDEGs$ct == 'Kupffer.cell' & allDEGs$direction == 'AK_up']), dbs)
downDEG_enriched <- enrichr(unique(allDEGs$geneSym[allDEGs$geno == 'T18' & allDEGs$direction == 'AK_down']), dbs)

View(upDEG_enriched[[6]][1,])
upDEG_enriched[[6]][1,]
plotEnrich(upDEG_enriched[[9]], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
plotEnrich(downDEG_enriched[[6]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")

