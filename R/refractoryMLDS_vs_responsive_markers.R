##----   Comparing Refractory MLDS vs responsive MLDS blasts    -----##

outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/7_bad.vs.good_MLDS'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)



##---------------##
#### Libraries ####
##---------------##
library(Seurat)
library(GenomicFeatures)
library(tidyverse)
library(UpSetR)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')



##----------------------------##
##   Set Global parameters  ####
##----------------------------##
#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf.gz'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)
## Import Gene Map
geneMap = read.delim('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_45842_SB_Leuk13104278_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz',header = F)
colnames(geneMap) = c('ensID','geneSym','GEX')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))




##------------------------------------------------##
#### 1. Import relevant scRNA datasets: MLDS    ####
##-----------------------------------------------##
## Import MLDS object
# mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
# mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
# mlds$broadLineage = mdat$broadLineage[match(mlds$cellID,mdat$cellID)]
# mlds$timePoint[mlds$orig.ident %in% c('MY.200531.14635833')] = 'D.Relapse'
# mlds$finalAnn_broad = as.character(mlds$annot_mar24)

mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS')
mlds$annot = as.character(mlds$annot_aug24)

# df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_AIRes_240131.csv')
# mlds$AI_res = df$AI_output[match(mlds$cellID,df$cellID)]
# mlds$finalAnn_broad[mlds$finalAnn_broad == 'TumourCNA' & mlds$timePoint == 'Diagnostic'] = 'Tumour'
# mlds$finalAnn_broad[mlds$donorID == 'L038' & mlds$finalAnn_broad == 'Tumour' & mlds$AI_res == 'abbFrac' & mlds$timePoint == 'TP1'] = 'TumourCNA_TP1'
# mlds$finalAnn_broad[mlds$donorID == 'L038' & mlds$finalAnn_broad == 'Tumour' & mlds$AI_res == 'abbFrac' & mlds$timePoint == 'Diagnostic'] = 'TumourCNA_D'
# mlds$finalAnn_broad[mlds$donorID == 'L038' & mlds$finalAnn_broad == 'Tumour' & mlds$AI_res == 'normFrac' & mlds$timePoint == 'Diagnostic'] = 'TumourNorm_D'
# mlds$finalAnn_broad[mlds$donorID == 'L038' & mlds$finalAnn_broad == 'Tumour' & mlds$AI_res == 'normFrac' & mlds$timePoint == 'TP1'] = 'TumourNorm_TP1'



##------------------------------------------------------------------------------------##
##   2a. DEGs between remission MLDS vs relapse/refractory ML-DS BM diagnostic      ####
##------------------------------------------------------------------------------------##
##------- Using pseudobulk -----------
cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in% mlds$cellID[mlds$donorID != 'L041' & mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$tissue =='BM']]]
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv')
mdat$annot = mdat$annot_aug24
cells_toKeep = c(mdat$cellID[mdat$annot_aug24 == 'Tumour' & mdat$disease == 'MLDS' & mdat$timePoint == 'Diagnostic' & mdat$tissue == 'BM' & !mdat$donorID %in% c('L041','CC3','L076','L038')],
                 mdat$cellID[mdat$annot_aug24 == 'Tumour' & mdat$disease == 'MLDS' & mdat$timePoint %in% c('D.Relapse2') & mdat$tissue == 'BM' & mdat$donorID %in% c('L076')],
                 l038$cellID[l038$group %in% c('clone2')])
#cnt_mtx = big.srat@assays$RNA@counts[,cells_toKeep]
cnt_mtx = mlds@assays$RNA@counts[,cells_toKeep]
rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = mdat[mdat$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex','disease','tissue','timePoint')]
mDat$group = ifelse(mDat$donorID %in% c('L076','L038'),'bad','good')
rownames(mDat) = mDat$cellID
mDat$donorID[mDat$donorID == 'L076' & mDat$timePoint == 'D.Relapse2'] = 'L076_R2D'
#mDat$donorID[mDat$donorID == 'L076' & mDat$timePoint == 'D.Relapse'] = 'L076_RD'
mDat$donorID[mDat$donorID == 'L038' & mDat$timePoint == 'TP1'] = 'L038_TP1'
table(mDat$donorID,mDat$group)

bad.vs.good.mlds_res = compareCell_simplified(toc = cnt_mtx,
                                             mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                                             coords = gns[rownames(cnt_mtx)],
                                             cellTypeName='MLDS',
                                             formula='~ %s',tgtChrs = paste0('chr',c(1:22)),
                                             donorID='donorID',groupID='group')

saveRDS(bad.vs.good.mlds_res,'DESeq2_L076L038_R2.clone2TP1.vs.otherGood_dMLDS_res.RDS')
#bad.vs.good.mlds_res = readRDS('DESeq2_L076L038_bad.vs.good_dMLDS_res.RDS')
# 
# bad.vs.good.mlds_res = bad.vs.good.mlds_res[['res']]
# bad.vs.good.mlds_res = as.data.frame(bad.vs.good.mlds_res)
# bad.vs.good.mlds_res = annotateGenes(bad.vs.good.mlds_res,geneMap = geneMap)

bad.vs.good.mlds_deg = bad.vs.good.mlds_res[['mainDE']]
bad.vs.good.mlds_deg$direction = ifelse(bad.vs.good.mlds_deg$log2FoldChange > 0, 'MLDS.good_up','MLDS.good_down')
bad.vs.good.mlds_deg$cellFrac.diff = bad.vs.good.mlds_deg$cellFrac_g1 - bad.vs.good.mlds_deg$cellFrac_g2
bad.vs.good.mlds_deg$cellFrac_max = pmax(bad.vs.good.mlds_deg$cellFrac_g1, bad.vs.good.mlds_deg$cellFrac_g2)

min_l2FC = 0.25

bad.vs.good.mlds_deg = bad.vs.good.mlds_deg[bad.vs.good.mlds_deg$padj < 0.05 & 
                                            abs(bad.vs.good.mlds_deg$log2FoldChange) >= min_l2FC & 
                                              abs(bad.vs.good.mlds_deg$cellFrac.diff) >= 0.05,]


table(bad.vs.good.mlds_deg$direction)
dim(bad.vs.good.mlds_deg)
bad.vs.good.mlds_deg = bad.vs.good.mlds_deg[order(abs(bad.vs.good.mlds_deg$log2FoldChange),decreasing = T),]


Idents(mlds) = mlds$group
DotPlot(mlds,#idents = unique(mlds$group[grepl('MLDS',mlds$group) & !grepl('TP',mlds$group)]),
        features = c(bad.vs.good.mlds_deg$geneSym[bad.vs.good.mlds_deg$direction == 'MLDS.good_up'][1:100]
                     #bad.vs.good.mlds_deg$geneSym[bad.vs.good.mlds_deg$direction == 'MLDS_down']
                     )) + 
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 
FeaturePlot(mlds,c('SLC7A11'))

write.csv(bad.vs.good.mlds_deg,'DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')
bad.vs.good.mlds_deg = read.csv('DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')



p1 = DotPlot(big.srat,
             #idents = unique(big.srat$group_4_dotPlot[grepl('^Tumour.*Diag|Tumour.*D0|Tumour.*L038|Tumour.*L076',big.srat$group_4_dotPlot) & big.srat$donorID !='CC3']),
             idents = unique(big.srat$group_4_dotPlot[!grepl('other',big.srat$group_4_dotPlot)]),
             cols = c(colAlpha(grey(0.98),0.85),'black'),
             features = c(#bad.vs.good.mlds_deg$geneSym[bad.vs.good.mlds_deg$direction == 'MLDS_up'],
               'CD36','EPS8',
                          bad.vs.good.mlds_deg$geneSym[bad.vs.good.mlds_deg$direction == 'MLDS_down'])
             
)+
  
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,face=ifelse(bad.vs.good.mlds_deg$chr == 'chr21','bold','plain'),
                                   colour =ifelse(bad.vs.good.mlds_deg$isTF==T,'purple',ifelse(bad.vs.good.mlds_deg$isCosmic==T,'darkgreen',
                                                                                       ifelse(bad.vs.good.mlds_deg$isCSM == T,'darkblue','black')))),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 


print(p1)  

source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")
upDEG_enriched <- enrichr(unique(bad.vs.good.mlds_deg$geneSym[bad.vs.good.mlds_deg$log2FoldChange < 0]), dbs)

View(upDEG_enriched[[6]])
i=1
plotEnrich(upDEG_enriched[[i]][upDEG_enriched[[i]]$Adjusted.P.value < 0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")




##------------------------------------------------------------------------------------------------------------------------##
#### 1. DEGs between Resistant blasts in L038 vs all other "responsive" MLDS blasts in BM, including Diagnostic L038    ####
##------------------------------------------------------------------------------------------------------------------------##


###---           Resistant clones vs other MLDS       -------------####
L038_dCNA_markers = data.frame()

for(donor in unique(mlds$donorID[!is.na(mlds$finalAnn_broad) & mlds$donorID != 'L041' & mlds$disease == 'MLDS'])){
  print(donor)
  #if(!donor %in% c('L019',"L038",'L039','L040','L042','L076','L091')){next}
  mlds$group_tmp = ifelse(mlds$donorID == 'L038' & 
                            #mlds$Phase == 'G1' & 
                            mlds$finalAnn_broad == 'TumourCNA_TP1' & 
                            mlds$timePoint == 'TP1' & 
                            mlds$tissue == 'BM','L038_TP1_CNA',
                          ifelse(mlds$donorID == 'L038' & 
                                   #mlds$Phase == 'G1' & 
                                   mlds$finalAnn_broad == 'TumourNorm_D' & 
                                   mlds$timePoint == 'Diagnostic' & 
                                   mlds$tissue == 'BM','L038_D_norm',
                                 ifelse(mlds$donorID == 'L038' & 
                                          #mlds$Phase == 'G1' & 
                                          mlds$finalAnn_broad == 'TumourCNA_D' & 
                                          mlds$timePoint == 'Diagnostic' & 
                                          mlds$tissue == 'BM','L038_D_CNA',
                                        ifelse(mlds$donorID == donor & 
                                                 #mlds$Phase == 'G1' & 
                                                 mlds$finalAnn_broad == 'Tumour' & 
                                                 mlds$timePoint == 'Diagnostic' & 
                                                 mlds$tissue == 'BM','other_Diagnostic','others'))))
  Idents(mlds) = mlds$group_tmp
  if(donor != 'L038'){
    markers = FindMarkers(mlds,ident.1 = 'L038_D_CNA',ident.2 = 'other_Diagnostic')  
    markers$gene = rownames(markers)
    markers$donorID = donor
  }else if(donor == 'L038'){
    markers.1 = FindMarkers(mlds,ident.1 = 'L038_D_CNA',ident.2 = 'L038_D_norm')  
    markers.1$gene = rownames(markers.1)
    markers.1$donorID = 'L038_D_norm'
    
    markers.2 = FindMarkers(mlds,ident.1 = 'L038_D_CNA',ident.2 = 'L038_TP1_CNA')
    markers.2$gene = rownames(markers.2)
    markers.2$donorID = 'L038_TP1_CNA'
    
    markers = rbind(markers.1,markers.2)
  }
  
  # Remove rubbish genes
  markers = markers[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',markers$gene),]
  markers = markers[markers$p_val_adj < 0.05,]
  
  L038_dCNA_markers = rbind(L038_dCNA_markers,markers)
  
  if(donor == 'L076'){
    mlds$group_tmp = ifelse(mlds$donorID == donor & 
                              #mlds$Phase == 'G1' & 
                              mlds$finalAnn_broad == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$tissue == 'Blood','other_Diagnostic',
                            ifelse(mlds$donorID == 'L038' & 
                                     mlds$finalAnn_broad == 'TumourCNA_D' & 
                                     mlds$timePoint == 'Diagnostic' & 
                                     mlds$tissue == 'BM','L038_D_CNA','others'))
    Idents(mlds) = mlds$group_tmp
    
    markers = FindMarkers(mlds,ident.1 = 'L038_D_CNA',ident.2 = 'other_Diagnostic')
    markers$gene = rownames(markers)
    # Remove rubbish genes
    markers = markers[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',markers$gene),]
    markers = markers[markers$p_val_adj < 0.05,]
    markers$donorID = 'L076_blood'
    L038_dCNA_markers = rbind(L038_dCNA_markers,markers)
    
    
    
    ##---- Compared to Relapse L076  -------##
    mlds$group_tmp = ifelse(mlds$donorID == donor & 
                              #mlds$Phase == 'G1' & 
                              mlds$finalAnn_broad == 'Tumour' & mlds$timePoint == 'D.Relapse' & mlds$tissue == 'BM','other_Diagnostic',
                            ifelse(mlds$donorID == 'L038' & 
                                     mlds$finalAnn_broad == 'TumourCNA_D' & 
                                     mlds$timePoint == 'Diagnostic' & 
                                     mlds$tissue == 'BM','L038_D_CNA','others'))
    Idents(mlds) = mlds$group_tmp
    
    markers = FindMarkers(mlds,ident.1 = 'L038_D_CNA',ident.2 = 'other_Diagnostic')
    markers$gene = rownames(markers)
    # Remove rubbish genes
    markers = markers[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',markers$gene),]
    markers = markers[markers$p_val_adj < 0.05,]
    markers$donorID = 'L076_BM_relapse'
    L038_dCNA_markers = rbind(L038_dCNA_markers,markers)
  }
}

L038_dCNA_markers$chr = geneMap$chr[match(L038_dCNA_markers$gene,geneMap$geneSym)]

write.csv(L038_dCNA_markers,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/7_DE_Diag.vs.Tum/jan24/L038.dCNA_vs.otherMLDS.markers_FM_240423.csv')
L038_dCNA_markers = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/7_DE_Diag.vs.Tum/jan24/L038.dCNA_vs.otherMLDS.markers_FM_240423.csv')




DimPlot(mlds,group.by = 'donorID',label = T,repel = T,label.size = 3,label.box = T,cols = col25) + NoLegend()
## Extract list of common markers between L038_dCNA and other MLDS samples ##
L038_dCNA_markers$direction = ifelse(L038_dCNA_markers$avg_log2FC > 0, 'L038_dCNA_up','L038_dCNA_down')
L038_dCNA_markers$pct.diff = L038_dCNA_markers$pct.1 - L038_dCNA_markers$pct.2
df = L038_dCNA_markers %>% filter(! donorID %in% c('L076','L076_blood','L076_BM_relapse','CC3')) %>% 
  group_by(direction,gene,chr) %>% summarise(nComp = n_distinct(donorID),
                                             mean_pct.diff = mean(abs(pct.diff)),
                                             mean_log2FC = mean(abs(avg_log2FC)))


extract_markers = function(fcLim = 0.3,pct.threshold = 0,L038_dCNA_markers,geneMap){
  df_woTP1 = L038_dCNA_markers %>% filter(!donorID %in% c('L038_TP1_CNA','L076_blood','L076','L076_BM_relapse','CC3')) %>% filter(abs(avg_log2FC) > fcLim &
                                                                                                                       abs(pct.diff) > pct.threshold) %>% group_by(direction,gene,chr) %>% summarise(nComp = n_distinct(donorID))
  
  
  # Marker genes between L038_dCNA and L038_TP1_CNA
  df_woTP1$in_TP1 = NA
  df_woTP1$in_TP1[df_woTP1$direction == 'L038_dCNA_up'] = (df_woTP1$gene[df_woTP1$direction == 'L038_dCNA_up'] %in% L038_dCNA_markers$gene[L038_dCNA_markers$donorID == 'L038_TP1_CNA' & L038_dCNA_markers$avg_log2FC > fcLim & L038_dCNA_markers$pct.diff > pct.threshold]) 
  df_woTP1$in_TP1[df_woTP1$direction == 'L038_dCNA_down'] = (df_woTP1$gene[df_woTP1$direction == 'L038_dCNA_down'] %in% L038_dCNA_markers$gene[L038_dCNA_markers$donorID == 'L038_TP1_CNA' & L038_dCNA_markers$avg_log2FC < -fcLim & L038_dCNA_markers$pct.diff < -pct.threshold]) 
  
  # Marker genes between L038_dCNA and L038_d
  df_woTP1$in_D = NA
  df_woTP1$in_D[df_woTP1$direction == 'L038_dCNA_up'] = (df_woTP1$gene[df_woTP1$direction == 'L038_dCNA_up'] %in% L038_dCNA_markers$gene[L038_dCNA_markers$donorID == 'L038_D_norm' & L038_dCNA_markers$avg_log2FC > fcLim]) 
  df_woTP1$in_D[df_woTP1$direction == 'L038_dCNA_down'] = (df_woTP1$gene[df_woTP1$direction == 'L038_dCNA_down'] %in% L038_dCNA_markers$gene[L038_dCNA_markers$donorID == 'L038_D_norm' & L038_dCNA_markers$avg_log2FC < -fcLim]) 
  
  df_woTP1$chr = geneMap$chr[match(df_woTP1$gene,geneMap$geneSym)]
  df_woTP1$ensID = geneMap$ensID[match(df_woTP1$gene,geneMap$geneSym)]
  
  df_woTP1_up = as.data.frame(df_woTP1[df_woTP1$direction == 'L038_dCNA_up',])
  rownames(df_woTP1_up) = df_woTP1_up$ensID
  df_woTP1_up = annotateGenes(df_woTP1_up,geneMap = geneMap)
  
  df_woTP1_down = as.data.frame(df_woTP1[df_woTP1$direction == 'L038_dCNA_down',])
  rownames(df_woTP1_down) = df_woTP1_down$ensID
  df_woTP1_down = annotateGenes(df_woTP1_down,geneMap = geneMap)
  
  df_woTP1 = rbind(df_woTP1_up,df_woTP1_down)
  
  return(df_woTP1)
}




##------- Extracting genes with high l2FC but not too different in percentage of cells expression.
##  These genes are expressed by similar amount of cells in both responsive and refractory cells, just differ by expression level.
##  Most of these are genes on Chr5 and chr17 due to the chr17 loss
## Find genes which are DEGs between 2 clones in diagnostic, but NOT DEGs between D vs TP1 refractory clones

fcLim = 0.3
pct.threshold = 0

df_woTP1_l2fc = extract_markers(fcLim = 0.3,pct.threshold = 0,L038_dCNA_markers,geneMap)
l2fc_genes = df_woTP1_l2fc[df_woTP1_l2fc$nComp >= 6 & 
                             df_woTP1_l2fc$in_D == T & df_woTP1_l2fc$in_TP1 == F,] %>% as.data.frame()
table(l2fc_genes$direction)
table(l2fc_genes$chr)
l2fc_genes$mean_pct.diff = df$mean_pct.diff[match(l2fc_genes$gene,df$gene)]
l2fc_genes$mean_l2FC = df$mean_log2FC[match(l2fc_genes$gene,df$gene)]


##------- Extracting genes with high l2FC but also high percentage_cells difference.
##  Looking to narrow down to genes with binary / on-off gene expression 
## Find genes which are DEGs between 2 clones in diagnostic, but NOT DEGs between D vs TP1 refractory clones

df_woTP1_l2fc_binary = extract_markers(fcLim = 0.3,pct.threshold = 0.3,L038_dCNA_markers,geneMap)
l2fc_genes_binary = df_woTP1_l2fc_binary[df_woTP1_l2fc_binary$nComp >= 6 & 
                                           df_woTP1_l2fc_binary$in_D == T & df_woTP1_l2fc_binary$in_TP1 == F,] %>% as.data.frame()
table(l2fc_genes_binary$direction)
table(l2fc_genes_binary$chr)
l2fc_genes_binary$mean_pct.diff = df$mean_pct.diff[match(l2fc_genes_binary$gene,df$gene)]
l2fc_genes_binary$mean_l2FC = df$mean_log2FC[match(l2fc_genes_binary$gene,df$gene)]




##------- Extracting genes with high l2FC but also high percentage_cells difference - Unique to TP1 refractory clones.
fcLim = 0.3
pct.threshold = 0.3
tp1_refractor_markers = L038_dCNA_markers[L038_dCNA_markers$donorID == 'L038_TP1_CNA' & abs(L038_dCNA_markers$avg_log2FC) > fcLim & abs(L038_dCNA_markers$pct.diff) > pct.threshold,]
tp1_refractor_markers$ensID = geneMap$ensID[match(tp1_refractor_markers$gene,geneMap$geneSym)]
rownames(tp1_refractor_markers) = tp1_refractor_markers$ensID
tp1_refractor_markers = annotateGenes(tp1_refractor_markers,geneMap = geneMap)
table(tp1_refractor_markers$direction)
table(tp1_refractor_markers$chr)




##------- Define resistsance module
resistanceModule = df_woTP1[df_woTP1$nComp >= 6 & 
                              df_woTP1$in_D == T & df_woTP1$in_TP1 == F,] %>% as.data.frame()
table(resistanceModule$direction)
table(resistanceModule$chr)
resistanceModule$mean_pct.diff = df$mean_pct.diff[match(resistanceModule$gene,df$gene)]
#resistanceModule = df_woTP1[df_woTP1$nComp >=4,] %>% as.data.frame()
resistanceModule = resistanceModule[order(resistanceModule$mean_pct.diff,decreasing = T),]
#resistanceModule_up = resistanceModule[resistanceModule$direction == 'L038_dCNA_up',]
rownames(resistanceModule) = resistanceModule$ensID
resistanceModule = annotateGenes(resistanceModule,geneMap = geneMap)

dim(resistanceModule)







genes = l2fc_genes_binary[order(abs(l2fc_genes_binary$mean_l2FC),decreasing = T),]
genes = tp1_refractor_markers[order(abs(tp1_refractor_markers$avg_log2FC),decreasing = T),]

##------ DotPlot in MLDS
mlds$broadLineage = mdat$broadLineage[match(mlds$cellID,mdat$cellID)]
mlds$broadLineage[is.na(mlds$broadLineage) & mlds$annot_jan24 == 'Tumour_maybe'] = 'Tumour?'

mlds$group_tmp = ifelse(mlds$finalAnn_broad == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$disease == 'MLDS' & mlds$donorID != 'L041',paste0('Tumour:',mlds$donorID),
                        ifelse(mlds$finalAnn_broad == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$disease == 'TAM',paste0('TAM:',mlds$donorID),
                               ifelse(mlds$donorID == 'L038' & grepl('Tumour.+',mlds$finalAnn_broad),paste0('Tumour:L038:',mlds$finalAnn_broad),
                                      ifelse(mlds$donorID == 'L038' & mlds$finalAnn_broad == 'Tumour',paste0('Tumour:L038:',mlds$timePoint),
                                             ifelse(mlds$donorID == 'L076' & mlds$finalAnn_broad == 'Tumour' & mlds$timePoint == 'D.Relapse',paste0('Tumour:L076:',mlds$timePoint),
                                                    ifelse(mlds$donorID == 'L076' & mlds$finalAnn_broad == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$tissue == 'BM',paste0('Tumour:L076:',mlds$timePoint,'_BM'),
                                             mlds$broadLineage))))))
mlds$group_tmp = factor(mlds$group_tmp,c('Tumour:L038:TumourCNA_TP1','Tumour:L038:TumourCNA_D','Tumour:L038:TumourNorm_D',
                                         unique(mlds$group_tmp[grepl('^Tumour:L038',mlds$group_tmp) & !mlds$group_tmp %in% c('Tumour:L038:TumourCNA_D','Tumour:L038:TumourCNA_TP1','Tumour:L038:TumourNorm_D')]),
                                         unique(mlds$group_tmp[grepl('^Tumour:',mlds$group_tmp) & !grepl('Tumour:L038',mlds$group_tmp)]),
                                         unique(mlds$group_tmp[grepl('^TAM',mlds$group_tmp)]),
                                         unique(mlds$group_tmp[!grepl('^TAM|^Tumour:',mlds$group_tmp)])))


mlds$group_tmp = ifelse(mlds$annot_aug24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & !mlds$donorID %in% c('L041','L076','L038','CC3'),mlds$donorID,
                        ifelse(mlds$annot_aug24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076',paste0('L076_',mlds$tissue),
                               ifelse(mlds$annot_aug24 == 'Tumour' & mlds$timePoint != 'Diagnostic' & mlds$donorID == 'L076',paste0('L076_',mlds$timePoint),
                                      ifelse(mlds$annot_aug24 == 'Tumour' & mlds$donorID == 'L038',paste0('L038_',mlds$timePoint),'others'))))
mlds$group_tmp = as.character(mlds$group_tmp)
mlds$group_tmp[mlds$cellID %in% l038_d$cellID[l038_d$seurat_clusters %in% c(4,18,0,17)]] = 'L038_D_clone2'
mlds$group_tmp[mlds$cellID %in% l038$cellID[l038$group != 'Normal']] = l038$group[match(mlds$cellID[mlds$cellID %in% l038$cellID[l038$group != 'Normal']],l038$cellID)]

mlds$group_tmp = factor(mlds$group_tmp,c('CC1','CC2','CC6','CC7','CC8','L075','L114','L156','L182',
                                         'CC4','CC5','L019','L039','L040','L042','L091','L178',
                                         'L076_Blood','L076_BM','L076_D.Relapse','L076_D.Relapse2',
                                         'clone1','clone2_D','clone2','clone2_ery',
                                         'L038_Diagnostic','L038_D_clone2','L038_TP1',
                                         'others'))

table(is.na(mlds$group_tmp))
Idents(mlds) = mlds$group_tmp

genes = c('PDE4DIP','CREB5','DDIT3','IDI1','CLDND1','IFNGR1',
          'IFI27','GP1BB','ENPP3',
          'H1F0','GJA1',
          'RUNX1','JAK1','CTNNB1',
          
          'SKAP1','RUNX1T1','NCAPG2','FHL2',
          'KLF1','NTRK1','CDKN2A','CD81',
          
          'EPS8',
          'ADAMTS1','CXCL8','NR4A2','ZFP36','MAP3K13','TNFAIP3','NR4A1','ETV3',
          'COL18A1','CD7','CD34','CFD')
DotPlot(mlds,
        idents = unique(mlds$group_tmp[!grepl('others',mlds$group_tmp)]),
        #features = c('EPS8','CD36')
        # features = unique(c(
        #              markers$geneSym[markers$comp == 'RD_vs_D' & markers$pct_diff >= 0.25 & abs(markers$avg_log2FC)>=0.5 & !grepl('^LINC\\d+$|^AC\\d+|^AL\\d+',markers$geneSym) & markers$avg_log2FC > 0],
        #              markers$geneSym[markers$comp == 'R2D_vs_D' & markers$pct_diff >= 0.1 & abs(markers$avg_log2FC)>=0.5 & !grepl('^LINC\\d+$|^AC\\d+|^AL\\d+',markers$geneSym) & markers$avg_log2FC > 0],
        #              'EPS8',
        #              bad.vs.good.mlds_deg$geneSym[!grepl('^LINC\\d+$|^AC\\d+|^AL\\d+',bad.vs.good.mlds_deg$geneSym) & bad.vs.good.mlds_deg$direction == 'MLDS_down' 
        #                                           & abs(bad.vs.good.mlds_deg$cellFrac.diff) > 0.3 & !bad.vs.good.mlds_deg$geneSym %in% c('GABARAPL1')
        #                                           ][1:20],
        #              bad.vs.good.mlds_deg$geneSym[!grepl('^LINC\\d+$|^AC\\d+|^AL\\d+',bad.vs.good.mlds_deg$geneSym) & bad.vs.good.mlds_deg$direction == 'MLDS_up' 
        #                                           & abs(bad.vs.good.mlds_deg$cellFrac.diff) > 0.2
        #                                           ][1:20]))
        features = unique(genes)
        
) + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                   # colour = ifelse(genes$chr[genes$direction == 'L038_dCNA_down'] == '21','purple',
                                   #                 ifelse(genes$isTF[genes$direction == 'L038_dCNA_down'] == T,'blue',
                                   #                        ifelse(genes$isCSM[genes$direction == 'L038_dCNA_down'] == T,'red',
                                   #                               ifelse(genes$isCosmic[genes$direction == 'L038_dCNA_down'][1:50] == T,'darkgreen','black'))))
  ),
  axis.text.y = element_text(size=11)) + xlab('') + ylab('') 



View(bad.vs.good.mlds_deg[bad.vs.good.mlds_deg$cellFrac_max > 0.1,])

DotPlot(mlds,
        idents = unique(mlds$group_tmp[!grepl('others',mlds$group_tmp)]),
        #idents = unique(mlds$group_tmp[!grepl('unsure|maybe|\\?|others|Tumour$|TumourNorm_TP1|Tumour_WT$',mlds$group_tmp)]),
        #idents = unique(mlds$group_tmp[grepl('Tumour:L038',mlds$group_tmp)]),
        #features = c('PRMT9','ANXA2','TUBB6','TAF13','DNM3','MAP3K3','AGPAT4')
        #features = c('TP53','TP53BP1','TSPO','HLA-A','IL2RG')
        #features = df_woTP1$gene[df_woTP1$nComp >=5 & df_woTP1$direction == 'L038_dCNA_up' & df_woTP1$in_D == T & df_woTP1$in_TP1==F]
        #features =  resistanceModule$gene[resistanceModule$direction == 'L038_dCNA_down'][41:84]
        #features =  c(resistanceModule$gene[resistanceModule$direction == 'L038_dCNA_down'][1:50])
        features =  c(genes$gene[genes$direction == 'L038_dCNA_up' & genes$chr == 'chr5'],
                      genes$gene[genes$direction == 'L038_dCNA_up' & genes$chr == 'chr17'],
                      genes$gene[genes$direction == 'L038_dCNA_up' & !genes$chr %in% c('chr5','chr17')],
                      'TP53',
                      genes$gene[genes$direction == 'L038_dCNA_down' & genes$chr == 'chr5'],
                      genes$gene[genes$direction == 'L038_dCNA_down' & genes$chr == 'chr17'],
                      'CA1',
                      genes$gene[genes$direction == 'L038_dCNA_down' & !genes$chr %in% c('chr5','chr17')]#[1:50]
                      )
        # features = c(L038_dCNA_markers$gene[L038_dCNA_markers$direction == 'L038_dCNA_down' & L038_dCNA_markers$donorID == 'L038_TP1_CNA' & abs(L038_dCNA_markers$avg_log2FC) > fcLim & abs(L038_dCNA_markers$pct.diff) > 0.3],
        #              L038_dCNA_markers$gene[L038_dCNA_markers$direction == 'L038_dCNA_up' & L038_dCNA_markers$donorID == 'L038_TP1_CNA' & abs(L038_dCNA_markers$avg_log2FC) > fcLim & abs(L038_dCNA_markers$pct.diff) > 0.3])
        
        #features = df$gene[df$nComp == 8 & df$direction == 'L038_dCNA_up'][1:30]
        #features = df_woTP1$gene[df_woTP1$nComp == 7 & df_woTP1$direction == 'L038_dCNA_down' & df_woTP1$in_TP1 == F][1:30]
) + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                   # colour = ifelse(genes$chr[genes$direction == 'L038_dCNA_down'] == '21','purple',
                                   #                 ifelse(genes$isTF[genes$direction == 'L038_dCNA_down'] == T,'blue',
                                   #                        ifelse(genes$isCSM[genes$direction == 'L038_dCNA_down'] == T,'red',
                                   #                               ifelse(genes$isCosmic[genes$direction == 'L038_dCNA_down'][1:50] == T,'darkgreen','black'))))
                                   ),
        axis.text.y = element_text(size=11)) + xlab('') + ylab('') 





View(genes[genes$direction == 'L038_dCNA_up' & genes$chr == 'chr17',])
View(genes[genes$direction == 'L038_dCNA_up' & !genes$chr %in% c('chr5','chr17'),])


DotPlot(mlds,idents = unique(mlds$group_tmp[!grepl('unsure|maybe|\\?|others|Tumour$|Tumour:L038$|Tumour:L038:TP1$|TumourNorm_TP1|Tumour_WT$',mlds$group_tmp)]),
        #idents = unique(mlds$group_tmp[grepl('Tumour:L038',mlds$group_tmp)]),
        #features = c('PRMT9','ANXA2','TUBB6','TAF13','DNM3','MAP3K3','AGPAT4')
        #features = c('TP53','TP53BP1','TSPO','HLA-A','IL2RG')
        #features = df_woTP1$gene[df_woTP1$nComp >=5 & df_woTP1$direction == 'L038_dCNA_up' & df_woTP1$in_D == T & df_woTP1$in_TP1==F]
        #features =  resistanceModule$gene[resistanceModule$direction == 'L038_dCNA_down'][41:84]
        #features =  c(resistanceModule$gene[resistanceModule$direction == 'L038_dCNA_down'][1:50])
        features =  c('LDLR','FOXO1','LDLRAD3','NRP2'
        )
        # features = c(L038_dCNA_markers$gene[L038_dCNA_markers$direction == 'L038_dCNA_down' & L038_dCNA_markers$donorID == 'L038_TP1_CNA' & abs(L038_dCNA_markers$avg_log2FC) > fcLim & abs(L038_dCNA_markers$pct.diff) > 0.3],
        #              L038_dCNA_markers$gene[L038_dCNA_markers$direction == 'L038_dCNA_up' & L038_dCNA_markers$donorID == 'L038_TP1_CNA' & abs(L038_dCNA_markers$avg_log2FC) > fcLim & abs(L038_dCNA_markers$pct.diff) > 0.3])
        
        #features = df$gene[df$nComp == 8 & df$direction == 'L038_dCNA_up'][1:30]
        #features = df_woTP1$gene[df_woTP1$nComp == 7 & df_woTP1$direction == 'L038_dCNA_down' & df_woTP1$in_TP1 == F][1:30]
) + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                   # colour = ifelse(genes$chr[genes$direction == 'L038_dCNA_down'] == '21','purple',
                                   #                 ifelse(genes$isTF[genes$direction == 'L038_dCNA_down'] == T,'blue',
                                   #                        ifelse(genes$isCSM[genes$direction == 'L038_dCNA_down'] == T,'red',
                                   #                               ifelse(genes$isCosmic[genes$direction == 'L038_dCNA_down'][1:50] == T,'darkgreen','black'))))
  ),
  axis.text.y = element_text(size=11)) + xlab('') + ylab('') 


plotDir='~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')

##------ Plot expression of interesting genes ----  DotPlot in MLDS  -----##
genes_of_interest = c('TP53','PIK3R1','CA2','MAGI1','SKAP1','PRKG1','RUNX1T1','CLCN3','NCAPG2',
                      'CDKN2A','PAG1','ARHGEF6','S1PR1','NTRK1','SELPLG','LDLRAD3')
l2fc_genes = l2fc_genes[order(l2fc_genes$mean_pct.diff,decreasing = T),]
l2fc_genes_binary
tp1_refractor_markers = tp1_refractor_markers[order(tp1_refractor_markers$pct.diff,decreasing = T),]

columns_toKeep = c('direction','chr','ensID','geneSym','isTF','isCSM','isCosmic')
df = rbind(l2fc_genes[l2fc_genes$chr =='chr5' & l2fc_genes$nComp >= 6 & l2fc_genes$direction == 'L038_dCNA_down',columns_toKeep][1:5,],
           l2fc_genes[l2fc_genes$chr =='chr17' & l2fc_genes$ensID %in% gns[end(gns) < 24e6 & seqnames(gns) == 'chr17']$gene_id & l2fc_genes$nComp >= 6 & l2fc_genes$direction == 'L038_dCNA_down',columns_toKeep],
           l2fc_genes[!l2fc_genes$chr %in% c('chr17','chr5') & l2fc_genes$nComp >= 6 & l2fc_genes$direction == 'L038_dCNA_down',columns_toKeep],
           
           #l2fc_genes[l2fc_genes$chr =='chr5' & l2fc_genes$nComp >= 6 & l2fc_genes$direction == 'L038_dCNA_up',columns_toKeep][1:5,],
           l2fc_genes[l2fc_genes$chr =='chr17' & l2fc_genes$nComp >= 6 & l2fc_genes$direction == 'L038_dCNA_up',columns_toKeep],
           l2fc_genes[!l2fc_genes$chr %in% c('chr5') & l2fc_genes$nComp >= 6 & l2fc_genes$direction == 'L038_dCNA_up',columns_toKeep][1:5,],
           
           l2fc_genes_binary[l2fc_genes_binary$direction == 'L038_dCNA_down',columns_toKeep],
           l2fc_genes_binary[l2fc_genes_binary$direction == 'L038_dCNA_up',columns_toKeep],
           tp1_refractor_markers[tp1_refractor_markers$direction == 'L038_dCNA_down' & (tp1_refractor_markers$isCosmic | tp1_refractor_markers$isCSM | tp1_refractor_markers$isTF),columns_toKeep],
           tp1_refractor_markers[tp1_refractor_markers$direction == 'L038_dCNA_down' & !(tp1_refractor_markers$isCosmic | tp1_refractor_markers$isCSM | tp1_refractor_markers$isTF),columns_toKeep][1:15,],
           tp1_refractor_markers[tp1_refractor_markers$direction == 'L038_dCNA_up',columns_toKeep][1:10,]
           )
table(df$geneSym %in% genes_of_interest)
table(genes_of_interest %in% df$geneSym)
genes_of_interest[!genes_of_interest %in% df$geneSym]

df = df[!grepl('AL\\d+|AC\\d+',df$geneSym),]




l038_markers = l038_tp1.vs.clone1
l038_markers = l038_markers[order(abs(l038_markers$avg_log2FC),decreasing = T),]
l038_markers = rbind(l038_markers[l038_markers$avg_log2FC > 0,],
                     l038_markers[l038_markers$avg_log2FC < 0,])
cn_genes = geneMap$geneSym[geneMap$ensID %in% gns[seqnames(gns) == 'chr17' & start(gns) < 24e6]$gene_id]
genes = c(l038_markers$geneSym[l038_markers$geneSym %in% cn_genes & l038_markers$is_DEG == 'DEG' & l038_markers$direction == 'D_down'],
          l038_markers$geneSym[l038_markers$geneSym %in% cn_genes & l038_markers$is_DEG == 'DEG' & l038_markers$direction == 'D_up'],
          l038_markers$geneSym[l038_markers$chr %in% 'chr5' & l038_markers$is_DEG == 'DEG' & l038_markers$direction == 'D_down'][1:5],
          l038_markers$geneSym[l038_markers$chr %in% 'chr5' & l038_markers$is_DEG == 'DEG' & l038_markers$direction == 'D_up'][1:5],
          l038_markers$geneSym[l038_markers$is_DEG == 'DEG' & l038_markers$direction == 'D_up'][1:10],
          l038_markers$geneSym[l038_markers$is_DEG == 'DEG' & l038_markers$direction == 'D_down' & !l038_markers$geneSym %in% c('HBB')][1:10])
genes = bad.vs.good.mlds_deg$geneSym[bad.vs.good.mlds_deg$chr =='chr5']
fig5_resistanceMarkers_L038 = function(){
  
  plotFun_interesting_genes_refractory_L038_vertical = function(noFrame=FALSE,noPlot=FALSE){
    p1 = DotPlot(mlds,idents = unique(mlds$group_tmp[!grepl('unsure|maybe|\\?|others|Tumour$|Tumour:L038$|Tumour:L038:TP1$|TumourNorm_TP1',mlds$group_tmp) & grepl('^Tumour|^TAM',mlds$group_tmp)]),
                 features =  c('TP53',unique(df$gene)),cols = c(grey(0.9),'#5d07a8')
    ) + RotatedAxis()+coord_flip()+
      theme(#legend.position = 'top',
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=11),
        axis.text.y = element_text(size=9,colour = ifelse(df$chr == '21','purple',
                                                          ifelse(df$isTF == T,'blue',
                                                                 ifelse(df$isCSM == T,'red',
                                                                        ifelse(df$isCosmic == T,'darkgreen','black')))))) + xlab('') + ylab('')
    
    print(p1)
  }
  saveFig(file.path(plotDir,'Fig5e_resistanceMarkers_L038'),plotFun_resistanceMarkers_L038,width = 6.5,height = 9.5,res = 500)
  
  
  
  
  plotFun_interesting_genes_refractory_L038_horizontal = function(noFrame=FALSE,noPlot=FALSE){
    # p1 = DotPlot(mlds,idents = unique(mlds$group_tmp[!grepl('unsure|maybe|\\?|others|Tumour$|TumourNorm_TP1|Tumour_WT$|Relapse$',mlds$group_tmp) & grepl('^Tumour|^TAM',mlds$group_tmp)]),
    #              features =  c('TP53',unique(df$gene)),
    #              cols = c(colAlpha(grey(0.95),0.8),'black')
    # ) + RotatedAxis()+
    #   #scale_x_discrete(position = "top")+
    #   theme(#legend.position = 'top',
    #     axis.text.y = element_text(size=11),
    #     axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
    #                                face=c('plain',ifelse(df$chr == 'chr5','bold',ifelse(df$chr %in% c('chr17'),'italic','plain'))),
    #                                colour = c('black',ifelse(df$isTF == T,'purple',ifelse(df$isCSM == T,'darkblue',ifelse(df$isCosmic == T,'darkgreen','black')))))) + 
    #   xlab('') + ylab('')
    
    #print(p1)
    
    
    p2 = DotPlot(mlds,
            idents = unique(mlds$group_tmp[!grepl('others|_ery$',mlds$group_tmp)]),
            #cols = c(grey(0.99),'black'),
            #features = c('EPS8','CD36')
            # features = unique(c(
            #              markers$geneSym[markers$comp == 'RD_vs_D' & markers$pct_diff >= 0.25 & abs(markers$avg_log2FC)>=0.5 & !grepl('^LINC\\d+$|^AC\\d+|^AL\\d+',markers$geneSym) & markers$avg_log2FC > 0],
            #              markers$geneSym[markers$comp == 'R2D_vs_D' & markers$pct_diff >= 0.1 & abs(markers$avg_log2FC)>=0.5 & !grepl('^LINC\\d+$|^AC\\d+|^AL\\d+',markers$geneSym) & markers$avg_log2FC > 0],
            #              'EPS8',
            #              bad.vs.good.mlds_deg$geneSym[!grepl('^LINC\\d+$|^AC\\d+|^AL\\d+',bad.vs.good.mlds_deg$geneSym) & bad.vs.good.mlds_deg$direction == 'MLDS_down' 
            #                                           & abs(bad.vs.good.mlds_deg$cellFrac.diff) > 0.3 & !bad.vs.good.mlds_deg$geneSym %in% c('GABARAPL1')
            #                                           ][1:20],
            #              bad.vs.good.mlds_deg$geneSym[!grepl('^LINC\\d+$|^AC\\d+|^AL\\d+',bad.vs.good.mlds_deg$geneSym) & bad.vs.good.mlds_deg$direction == 'MLDS_up' 
            #                                           & abs(bad.vs.good.mlds_deg$cellFrac.diff) > 0.2
            #                                           ][1:20]))
            features = unique(genes)
            
    ) + RotatedAxis()+
      theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                       # colour = ifelse(genes$chr[genes$direction == 'L038_dCNA_down'] == '21','purple',
                                       #                 ifelse(genes$isTF[genes$direction == 'L038_dCNA_down'] == T,'blue',
                                       #                        ifelse(genes$isCSM[genes$direction == 'L038_dCNA_down'] == T,'red',
                                       #                               ifelse(genes$isCosmic[genes$direction == 'L038_dCNA_down'][1:50] == T,'darkgreen','black'))))
      ),
      axis.text.y = element_text(size=11)) + xlab('') + ylab('') +
      scale_color_gradient2(low = grey(0.98),mid=grey(0.75),high='black')
    
    print(p2)
  }
  
  
  #saveFig(file.path(plotDir,'Fig5e_resistanceMarkers_L038_horizontal'),plotFun_interesting_genes_refractory_L038_horizontal,rawData = df,width = 16,height = 4.8,res = 500)
  saveFig(file.path(plotDir,'SupFig9E_L038_tp1.clone2_vs_clone1.diagnostic_dotPlot_horizontal'),plotFun_interesting_genes_refractory_L038_horizontal,rawData = NULL,width = 11,height = 5.8,res = 500)
}







fig5_resistanceMarkers_L038 = function(){
  # resistanceModule = df_woTP1[df_woTP1$nComp >= 6 & 
  #                               df_woTP1$in_D == T,] %>% as.data.frame()
  # table(resistanceModule$direction)
  # table(resistanceModule$in_D)
  # resistanceModule$mean_pct.diff = df$mean_pct.diff[match(resistanceModule$gene,df$gene)]
  # #resistanceModule = df_woTP1[df_woTP1$nComp >=4,] %>% as.data.frame()
  # resistanceModule = resistanceModule[order(resistanceModule$mean_pct.diff,decreasing = T),]
  # #resistanceModule_up = resistanceModule[resistanceModule$direction == 'L038_dCNA_up',]
  # rownames(resistanceModule) = resistanceModule$ensID
  # resistanceModule = annotateGenes(resistanceModule,geneMap = geneMap)
  # 
  # dim(resistanceModule)
  # 
  # df = rbind(resistanceModule[resistanceModule$direction == 'L038_dCNA_up',][1:20,],
  #            resistanceModule[resistanceModule$direction == 'L038_dCNA_down',][1:10,])
  
  plotFun_resistanceMarkers_L038 = function(noFrame=FALSE,noPlot=FALSE){
    p1 = DotPlot(mlds,idents = unique(mlds$group_tmp[!grepl('unsure|maybe|\\?|others|Tumour$|Tumour:L038$|Tumour:L038:TP1$|TumourNorm_TP1',mlds$group_tmp) & grepl('^Tumour|^TAM',mlds$group_tmp)]),
                 features =  df$gene,cols = c(grey(0.9),'#5d07a8')
    ) + RotatedAxis()+coord_flip()+
      #scale_x_discrete(position = "top")+
      theme(#legend.position = 'top',
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=11),
        axis.text.y = element_text(size=9,colour = ifelse(df$chr == '21','purple',
                                                          ifelse(df$isTF == T,'blue',
                                                                 ifelse(df$isCSM == T,'red',
                                                                        ifelse(df$isCosmic == T,'darkgreen','black')))))) + xlab('') + ylab('')
    
    print(p1)
  }
  saveFig(file.path(plotDir,'Fig5e_resistanceMarkers_L038'),plotFun_resistanceMarkers_L038,width = 6.5,height = 9.5,res = 500)
  
  
  
  df = rbind(resistanceModule[resistanceModule$direction == 'L038_dCNA_up',][1:30,],
             resistanceModule[resistanceModule$direction == 'L038_dCNA_down',][1:30,])
  plotFun_resistanceMarkers_L038_horizontal = function(noFrame=FALSE,noPlot=FALSE){
    p1 = DotPlot(mlds,idents = unique(mlds$group_tmp[!grepl('unsure|maybe|\\?|others|Tumour$|Tumour:L038$|Tumour:L038:TP1$|TumourNorm_TP1',mlds$group_tmp) & grepl('^Tumour|^TAM',mlds$group_tmp)]),
                 features =  df$gene,
                 #cols = c(grey(0.9),'#5d07a8'),
                 cols = c(colAlpha(grey(0.95),0.8),'black')
    ) + RotatedAxis()+
      #scale_x_discrete(position = "top")+
      theme(#legend.position = 'top',
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                   face=ifelse(df$chr == 'chr21','bold',ifelse(df$chr %in% c('chr5','chr17'),'italic','plain')),
                                   colour = ifelse(df$isTF == T,'purple',ifelse(df$isCSM == T,'darkblue',ifelse(df$isCosmic == T,'darkgreen','black'))))) + 
      xlab('') + ylab('')
    
    print(p1)
  }
  saveFig(file.path(plotDir,'Fig5e_resistanceMarkers_L038_horizontal'),plotFun_resistanceMarkers_L038_horizontal,width = 15,height = 4.5,res = 500)
}


