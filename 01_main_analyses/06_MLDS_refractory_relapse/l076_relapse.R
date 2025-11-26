outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/7_bad.vs.good_MLDS'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)


library(Seurat)
library(tidyverse)
source('~/lustre_mt22/generalScripts/utils/misc.R')

#-----------------------------------##
# Transcriptional markers of L076 ####
#-----------------------------------##
l076 = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/L076_sratObj.RDS')
l076$GATA1s_status = l076$GATA1_status
l076$GATA1s_status[l076$GATA1s_status == 'no_GATA1_expr'] = 'No GATA1 expression'
l076$GATA1s_status[l076$GATA1s_status == 'GATA1s_WT'] = 'GATA1 wild type'
l076$GATA1s_status[l076$GATA1s_status == 'GATA1s_mutant'] = 'GATA1s mutation'
l076$GATA1s_status[l076$GATA1s_status %in% c('unsure','noCov','uninformative')] = 'Uninformative'
DimPlot(l076,group.by = 'GATA1s_status',cols = ccs,label = F,repel = T,label.box = T) + NoAxes() + NoLegend()



l076$group = ifelse(l076$annot_aug24 == 'Tumour' & l076$tissue=='Blood','D_Blood',
                    ifelse(l076$annot_aug24 == 'Tumour' & l076$tissue!='Blood',l076$timePoint,'others'))
l076$group[l076$group == 'Diagnostic'] = 'Diagnosis_BM'
l076$group[l076$group == 'D_Blood'] = 'Diagnosis_Blood'
l076$group[l076$group == 'others'] = 'Normal'
DimPlot(l076,group.by = 'group',cols = col25,label = F,repel = T,label.box = T) + NoAxes() + NoLegend()


## --- Plot L076 UMAP
l076_bm = subset(l076,subset = tissue == 'BM')
l076_bm = standard_clustering(l076_bm,clusteringRes = 0.2)
DimPlot(l076_bm,group.by = 'group',cols=col25)



l076_bm$GATA1s_status = l076_bm$GATA1_status
l076_bm$GATA1s_status[l076_bm$GATA1s_status == 'no_GATA1_expr'] = 'No GATA1 expression'
l076_bm$GATA1s_status[l076_bm$GATA1s_status == 'GATA1s_WT'] = 'GATA1 wild type'
l076_bm$GATA1s_status[l076_bm$GATA1s_status == 'GATA1s_mutant'] = 'GATA1s mutation'
l076_bm$GATA1s_status[l076_bm$GATA1s_status %in% c('unsure','noCov','uninformative')] = 'Uninformative'
DimPlot(l076_bm,group.by = 'GATA1s_status',cols = ccs,label = F,repel = T,label.box = T) + NoAxes() + NoLegend()


dd=cbind(l076_bm@meta.data,l076_bm@reductions$umap@cell.embeddings)

dd$umap_1 = dd$UMAP_1
dd$umap_2 = dd$UMAP_2



plotFun_GATA1status = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  ccs = c('No GATA1 expression' = grey(0.9),
          'Uninformative' = grey(0.55),
          'GATA1s mutation' = '#A92821',
          'GATA1 wild type' = '#005579')
  p = ggplot(dd,aes(UMAP_1,UMAP_2))+
    geom_point(size=0.001,aes(col=GATA1s_status),alpha=0.4)+
    scale_color_manual(values = ccs)+
    theme_classic(base_size = 8.3)+ 
    theme(panel.border = element_blank(),axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=3)) + 
    xlab('') + ylab('')
  
  print(p)
}
saveFig(file.path(plotDir,'Fig4C_L076_GATA1s_UMAP'),plotFun_GATA1status,rawData=dd,width = 3.65,height = 2.65,res = 500,useDingbats = F)


dd=cbind(l076.tum@meta.data,l076.tum@reductions$umap@cell.embeddings)

dd$umap_1 = dd$UMAP_1
dd$umap_2 = dd$UMAP_2

plotFun_timePoint = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  ccs = c('D.Relapse' = col25[5],
          'D.Relapse2' = col25[3],
          'Diagnostic' = pal34H[34],
          'Normal' = grey(0.9))
  dd$group = ifelse(dd$annot_aug24 =='Tumour',dd$timePoint,'Normal')
  p = ggplot(dd,aes(UMAP_1,UMAP_2))+
    geom_point(size=0.0001,aes(col=group))+
    scale_color_manual(values = ccs)+
    theme_classic(base_size = 8.3)+ 
    theme(panel.border = element_blank(),axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=3)) + 
    xlab('') + ylab('')
  
  print(p)
}

saveFig(file.path(plotDir,'Fig4C_L076.BM.tum_timePoint_UMAP'),plotFun_timePoint,rawData=dd,width = 4,height = 3.65,res = 500,useDingbats = F)



pp = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/L076/tumourDNA_onlyAnalysis/PD64665a_pp.RDS')

m = match(l076@meta.data$cellID,pp[seqnames(pp) == 'genomeWide',]$cellID)
l076$pp = pp[seqnames(pp) == 'genomeWide',]$probAbberant[m]
l076$pp[l076$annot_aug24 !='Tumour'] = NA
FeaturePlot(l076,'pp',cols = c(grey(0.9),'#c65102')) + NoAxes() + NoLegend()
l076_bm$pp = l076$pp[match(l076_bm$cellID,l076$cellID)]
l076.tum$pp = l076$pp[match(l076.tum$cellID,l076$cellID)]

df = cbind(l076.tum@meta.data,l076.tum@reductions$umap@cell.embeddings)

plotFun_AIresult = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  
  p = ggplot(df,aes(UMAP_1,UMAP_2))+
    geom_point(size=0.0001,aes(col=pp))+
    scale_color_gradient(low='black',high='#e01507',na.value = grey(0.8))+
    theme_classic(base_size = 8.3)+ 
    theme(panel.border = element_blank(),axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key.size = unit(0.45,'cm'),
          legend.text = element_text(size=10)) + 
    xlab('') + ylab('')
  
  print(p)
}

saveFig(file.path(plotDir,'Figure4B_L076.TP1_AIresult_UMAP'),plotFun_AIresult,rawData=df,width = 4,height = 3.8,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'Figure4B_L076.BM.tum_AIresult_UMAP'),plotFun_AIresult,rawData=df,width = 4,height = 3.8,res = 500,useDingbats = T)



saveRDS(l076_bm,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/L076_BM.samples.only_sratObj.RDS')
l076_bm = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/L076_BM.samples.only_sratObj.RDS')
DimPlot(l076_bm,group.by = 'seurat_clusters',label = T)

Idents(l076_bm) = as.character(l076_bm$seurat_clusters)
markers_bm_only = FindMarkers(l076_bm,ident.1 = c('0','1'),ident.2 = c('3','10','5','6'))
markers_bm_only$geneSym = rownames(markers_bm_only)
markers_bm_only$pct_diff = markers_bm_only$pct.1 - markers_bm_only$pct.2
markers_bm_only$comp = 'R2D_vs_RD.D'
markers_bm_only.sub = markers_bm_only[abs(markers_bm_only$pct_diff) > 0.15,]
markers_bm_only.sub = markers_bm_only.sub[order(markers_bm_only.sub$pct_diff,decreasing = T),]
rownames(markers_bm_only.sub) = geneMap$ensID[match(markers_bm_only.sub$geneSym,geneMap$geneSym)]
markers_bm_only.sub = annotateGenes(markers_bm_only.sub,geneMap = geneMap)




l076.tum = subset(l076,subset = cellID %in% l076$cellID[l076$annot_aug24 == 'Tumour' & l076$tissue != 'Blood'])
l076.tum = standard_clustering(l076.tum,clusteringRes = 0.2)

DimPlot(l076.tum,group.by = 'group',cols = c(col25[c(1,2,4)],grey(0.8)),label = F,repel = T,label.box = T) + 
  #ggtitle('L076 all samples all cells') + 
  ggtitle('L076 bone marrow blasts only') + 
  NoAxes()+theme(panel.border = element_rect(fill=F,colour = 'black'))
DimPlot(l076.tum,group.by = 'GATA1_status',cols = col25,label = T,repel = T,label.box = T)
DimPlot(l076,group.by = 'Phase',cols = col25,label = T,repel = T,label.box = T)

l076$GATA1s_status = l076$GATA1_status
l076$GATA1s_status[l076$GATA1s_status == 'no_GATA1_expr'] = 'No GATA1 expression'
l076$GATA1s_status[l076$GATA1s_status == 'GATA1s_WT'] = 'GATA1 wild type'
l076$GATA1s_status[l076$GATA1s_status == 'GATA1s_mutant'] = 'GATA1s mutation'
l076$GATA1s_status[l076$GATA1s_status %in% c('unsure','noCov','uninformative')] = 'Uninformative'
DimPlot(l076,group.by = 'GATA1s_status',cols = ccs,label = F,repel = T,label.box = T) + NoAxes() + NoLegend()


mlds_moduleScore = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/UCell_moduleScore/nov24/MLDS_MLDS_topGenes_2410_cf2_moduleScore.csv')
l076$MLDS_composite_UCell = mlds_moduleScore$MLDS_composite_UCell[match(l076$cellID,mlds_moduleScore$cellID)]
#big.srat$MLDS_composite_UCell = mlds_moduleScore$MLDS_composite_UCell[match(big.srat$cellID,mlds_moduleScore$cellID)]
FeaturePlot(l076,'MLDS_composite_UCell')
FeaturePlot(l076,'nCount_RNA')

Idents(l076) = l076$group
markers_diagnostic = FindMarkers(l076,ident.1 = 'D_Blood',ident.2 = 'Diagnostic')
markers_diagnostic$geneSym = rownames(markers_diagnostic)
markers_diagnostic$pct_diff = markers_diagnostic$pct.1 - markers_diagnostic$pct.2
markers_diagnostic.sub = markers_diagnostic[abs(markers_diagnostic$pct_diff) > 0.6,]
markers_diagnostic.sub = markers_diagnostic.sub[order(markers_diagnostic.sub$pct_diff,decreasing = T),]
FeaturePlot(l076,'PBX4')

DotPlot(l076,group.by = 'group',features = markers_diagnostic.sub$geneSym) + RotatedAxis()

## How is relapse 1 different from diagnostic
markers_rd1_BM = FindMarkers(l076,ident.1 = 'D.Relapse',ident.2 = 'Diagnostic')
markers_rd1_BM$geneSym = rownames(markers_rd1_BM)
markers_rd1_BM$pct_diff = markers_rd1_BM$pct.1 - markers_rd1_BM$pct.2
markers_rd1_BM$comp = 'RD_vs_D'
markers_rd1_BM.sub = markers_rd1_BM[abs(markers_rd1_BM$pct_diff) > 0.5,]
markers_rd1_BM.sub = markers_rd1_BM.sub[order(markers_rd1_BM.sub$pct_diff,decreasing = T),]
rownames(markers_rd1_BM.sub) = geneMap$ensID[match(markers_rd1_BM.sub$geneSym,geneMap$geneSym)]
markers_rd1_BM.sub = annotateGenes(markers_rd1_BM.sub,geneMap = geneMap)
DotPlot(l076,group.by = 'group',features = markers_rd1_BM.sub$geneSym[markers_rd1_BM.sub$avg_log2FC < 0]) + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))

FeaturePlot(l076,'EPS8')

## How is relapse 2 different from diagnostic
markers_rd2_BM = FindMarkers(l076,ident.1 = 'D.Relapse2',ident.2 = 'Diagnostic')
markers_rd2_BM$geneSym = rownames(markers_rd2_BM)
markers_rd2_BM$pct_diff = markers_rd2_BM$pct.1 - markers_rd2_BM$pct.2
markers_rd2_BM$comp = 'R2D_vs_D'
markers_rd2_BM.sub = markers_rd2_BM[abs(markers_rd2_BM$pct_diff) > 0.1,]
markers_rd2_BM.sub = markers_rd2_BM.sub[order(markers_rd2_BM.sub$pct_diff,decreasing = T),]
rownames(markers_rd2_BM.sub) = geneMap$ensID[match(markers_rd2_BM.sub$geneSym,geneMap$geneSym)]
markers_rd2_BM.sub = annotateGenes(markers_rd2_BM.sub,geneMap = geneMap)

markers = rbind(markers_rd1_BM,markers_rd2_BM)
write.csv(markers,'L076_relapse_FindMarkers.csv')

markers = read.csv('L076_relapse_FindMarkers.csv')
markers$direction = ifelse(markers$avg_log2FC > 0 & markers$comp == 'R2D_vs_D','R2D_up',
                           ifelse(markers$avg_log2FC > 0 & markers$comp == 'RD_vs_D','RD_up',
                                  ifelse(markers$avg_log2FC < 0 & markers$comp == 'R2D_vs_D','R2D_down',
                                         ifelse(markers$avg_log2FC < 0 & markers$comp == 'RD_vs_D','RD_down',
                                         'others'))))
#rownames(markers) = geneMap$ensID[match(markers$geneSym,geneMap$geneSym)]
#markers = annotateGenes(markers,geneMap = geneMap)

## Prioritise gene for a module 
markers.sub = markers[abs(markers$pct_diff) >= 0.2 & abs(markers$avg_log2FC) >= 0.5,]
table(markers.sub$comp,markers.sub$direction)




markers_rd2_BM.sub$geneSym[markers_rd2_BM.sub$avg_log2FC < 0]
DotPlot(l076,group.by = 'group',
        features = c(markers.sub$geneSym[markers.sub$avg_log2FC > 0 & markers.sub$comp == 'R2D_vs_D'],
                     markers.sub$geneSym[markers.sub$avg_log2FC < 0 & markers.sub$comp == 'R2D_vs_D'])) + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))

markers = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/findAllMarkers_L076_tumourOnly.csv')
markers$pct_diff = markers$pct.1 - markers$pct.2
markers.sub = markers[markers$pct_diff > 0.4,]
FeaturePlot(l076_bm,'CD99')




## How is relapse 2 different from relapse 1
markers_rd2.vs.rd1_BM = FindMarkers(l076,ident.1 = 'D.Relapse2',ident.2 = 'D.Relapse')
markers_rd2.vs.rd1_BM$geneSym = rownames(markers_rd2.vs.rd1_BM)
markers_rd2.vs.rd1_BM$pct_diff = markers_rd2.vs.rd1_BM$pct.1 - markers_rd2.vs.rd1_BM$pct.2
markers_rd2.vs.rd1_BM$comp = 'R2D_vs_RD'
markers_rd2.vs.rd1_BM$direction = ifelse(markers_rd2.vs.rd1_BM$avg_log2FC > 0,'R2D_up','R2D_down')
rownames(markers_rd2.vs.rd1_BM) = geneMap$ensID[match(markers_rd2.vs.rd1_BM$geneSym,geneMap$geneSym)]
markers_rd2.vs.rd1_BM = annotateGenes(markers_rd2.vs.rd1_BM,geneMap = geneMap)
write.csv(markers_rd2.vs.rd1_BM,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/L076_RD2.vs.RD1_findMarkers.csv')

## Prioritise gene for a module 
markers_rd2.vs.rd1_BM.sub = markers_rd2.vs.rd1_BM[abs(markers_rd2.vs.rd1_BM$pct_diff) >= 0.2 & abs(markers_rd2.vs.rd1_BM$avg_log2FC) >= 0.5,]
markers_rd2.vs.rd1_BM.sub = markers_rd2.vs.rd1_BM.sub[order(markers_rd2.vs.rd1_BM.sub$pct_diff,decreasing = T),]
table(markers_rd2.vs.rd1_BM.sub$direction)

markers.sub = markers.sub[,colnames(markers.sub) !='X']
l076_resistance_geneModule = rbind(markers.sub,markers_rd2.vs.rd1_BM.sub[,match(colnames(markers.sub),colnames(markers_rd2.vs.rd1_BM.sub))])

write.csv(l076_resistance_geneModule,'L076_resistance_geneModules.csv')

l076_resistance_geneModule = read.csv('L076_resistance_geneModules.csv')




## Perform enrichR on all DEGs
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")
upDEG_enriched <- enrichr(unique(markers_rd2_BM.sub$geneSym[markers_rd2_BM.sub$avg_log2FC > 0]), dbs)

View(upDEG_enriched[[6]])
i=6
plotEnrich(upDEG_enriched[[i]][upDEG_enriched[[i]]$Adjusted.P.value < 0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")






FeaturePlot(l076,'GATA1')
library(SoupX)
qm = quickMarkers(l076@assays$RNA@counts[,l076$cellID[l076$annot_aug24 == 'Tumour']],l076$timePoint[l076$annot_aug24 == 'Tumour']) 


##------- Using pseudobulk -----------

## Version 1: only canonical / good TAM
cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in% 
                                                c(mlds$cellID[mlds$donorID %in% c('L038','L076') & mlds$annot == 'Tumour' & mlds$tissue =='BM' & mlds$timePoint == 'Diagnostic'],
                                                  mlds$cellID[mlds$disease == 'MLDS' & mlds$timePoint %in% c('Diagnostic') & mlds$annot == 'Tumour' & grepl('Remission',mlds$clinicalOutcome)])]]

# ## Version 2: using All TAM cases
# cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in% 
#                                                 c(mlds$cellID[mlds$disease == 'TAM' & mlds$annot == 'Tumour'],
#                                                   mlds$cellID[mlds$donorID %in% c('L019','L038','L039','L040') & mlds$annot == 'MEP'])]]


rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex')]

mDat$donorID = as.character(mDat$donorID)
mDat = mDat[match(colnames(cnt_mtx),mDat$cellID),]
mDat$group = ifelse(mDat$donorID %in% c('L038','L076'),'bad','good')

mlds.GoodvsBad_res = compareCell_simplified(toc = cnt_mtx,
                                         mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                                         coords = gns[rownames(cnt_mtx)],
                                         cellTypeName='refractoryMLDS',
                                         formula='~ %s',tgtChrs = paste0('chr',c(1:22)),
                                         donorID='donorID',groupID='group')

saveRDS(tam.vs.gMEP_res,'DESeq2_goodTAM.vs.goshMEP_res.RDS')
tam.vs.gMEP_res = readRDS('DESeq2_goodTAM.vs.goshMEP_res.RDS')

mlds.GoodvsBad_dds = tam.vs.gMEP_res[['dds']]
mlds.GoodvsBad_deg = mlds.GoodvsBad_res[['mainDE']]
mlds.GoodvsBad_deg$direction = ifelse(mlds.GoodvsBad_deg$log2FoldChange > 0, 'TAM_up','TAM_down')
mlds.GoodvsBad_deg$cellFrac.diff = mlds.GoodvsBad_deg$cellFrac_g1 - mlds.GoodvsBad_deg$cellFrac_g2
mlds.GoodvsBad_deg = mlds.GoodvsBad_deg[mlds.GoodvsBad_deg$padj < 0.05 & abs(mlds.GoodvsBad_deg$log2FoldChange) >= min_l2FC & 
                                    abs(mlds.GoodvsBad_deg$cellFrac.diff) >= max_cellfrac,]





##-----
l038_mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_TumourEry_subClustering_AIRes_240131.csv')
l038_d = subset(big.srat,subset = cellID %in% big.srat$cellID[big.srat$donorID == 'L038' & big.srat$timePoint == 'D'])
l038_d = standard_clustering(l038_d)
saveRDS(l038_d,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_diagnostic.RDS')
l038_d = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_diagnostic.RDS')

l038_tp1 = subset(big.srat,subset = cellID %in% big.srat$cellID[big.srat$donorID == 'L038' & big.srat$timePoint == 'TP1'])
l038_tp1 = standard_clustering(l038_tp1)
DimPlot(l038_d,label = T,repel = T,label.box = T)

saveRDS(l038_tp1,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_TP1.RDS')
l038_tp1 = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_TP1.RDS')

l038_tp1$GATA1s_status = l038_tp1$GATA1_status
l038_tp1$GATA1s_status[l038_tp1$GATA1s_status == 'no_GATA1_expr'] = 'No GATA1 expression'
l038_tp1$GATA1s_status[l038_tp1$GATA1s_status == 'GATA1s_WT'] = 'GATA1 wild type'
l038_tp1$GATA1s_status[l038_tp1$GATA1s_status == 'GATA1s_mutant'] = 'GATA1s mutation'
l038_tp1$GATA1s_status[l038_tp1$GATA1s_status %in% c('unsure','noCov','uninformative')] = 'Uninformative'
DimPlot(l038_tp1,group.by = 'GATA1s_status',cols = ccs,label = F,repel = T,label.box = T) + NoAxes() + NoLegend()


ccs = c('No GATA1 expression' = grey(0.8),
        'Uninformative' = grey(0.55),
        'GATA1s mutation' = '#A92821',
        'GATA1 wild type' = '#005579')

l038_d$GATA1s_status = l038_d$GATA1_status
l038_d$GATA1s_status[l038_d$GATA1s_status == 'no_GATA1_expr'] = 'No GATA1 expression'
l038_d$GATA1s_status[l038_d$GATA1s_status == 'GATA1s_WT'] = 'GATA1 wild type'
l038_d$GATA1s_status[l038_d$GATA1s_status == 'GATA1s_mutant'] = 'GATA1s mutation'
l038_d$GATA1s_status[l038_d$GATA1s_status %in% c('unsure','noCov','uninformative')] = 'Uninformative'

DimPlot(l038_d,group.by = 'GATA1s_status',cols = ccs,label = F,repel = T,label.box = T) + NoAxes() + NoLegend()




dd=cbind(l038_d@meta.data,l038_d@reductions$umap@cell.embeddings)
dd=cbind(l038_tp1@meta.data,l038_tp1@reductions$umap@cell.embeddings)

dd$umap_1 = dd$UMAP_1
dd$umap_2 = dd$UMAP_2



plotFun_GATA1status = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  
  ccs = c('No GATA1 expression' = grey(0.9),
          'Uninformative' = grey(0.55),
          'GATA1s mutation' = '#A92821',
          'GATA1 wild type' = '#005579')
  # plot(dd$umap_1,dd$umap_2,
  #      las=1,
  #      type='n',
  #      #xlim=c(-13,17),
  #      #ylim=c(-13,17),
  #      cex.main = 0.85,xaxt='n',yaxt='n',
  #      xlab='',ylab='',
  #      main=ifelse(noFrame,'','L038'),
  #      frame.plot=F)
  # 
  # if(!noPlot){
  #   #Add density contours
  #   #addDensityContours(dd$UMAP_1,dd$UMAP_2,dd$finalAnn,col=colAlpha('black',0.4),nGrid = 2000)
  #   points(dd$umap_1,dd$umap_2,
  #          col = ccs[dd$GATA1s_status],
  #          pch = 19,
  #          cex=0.01)
  #   
  #   
  # }
  # #legend(x=-8.5, y=9,legend=unique(dd$GATA1s_status),fill = ccs[unique(dd$GATA1s_status)],lwd = 0,cex = 0.5,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
  
  
  p = ggplot(dd,aes(UMAP_1,UMAP_2))+
    geom_point(size=0.05,aes(col=GATA1s_status))+
    scale_color_manual(values = ccs)+
    theme_classic(base_size = 8.3)+ 
    theme(panel.border = element_blank(),axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=3)) + 
    xlab('') + ylab('')
  
  print(p)
}

saveFig(file.path(plotDir,'Fig4D_L038.Diag_GATA1s_UMAP'),plotFun_GATA1status,rawData=dd,width = 2.65,height = 1.65,res = 500,useDingbats = F)
saveFig(file.path(plotDir,'Fig4D_L038.TP1_GATA1s_UMAP'),plotFun_GATA1status,rawData=dd,width = 2.65,height = 1.65,res = 500,useDingbats = F)







ai_res = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_TumourOnly_subClustering_mdat_240206.csv')
df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_TumourEry_subClustering_AIRes_240131.csv')

pp = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_pp.RDS')
pp2 = pp

m = match(l038_d@meta.data$cellID,pp2[seqnames(pp2) == 'genomeWide',]$cellID)
l038_d$pp = pp2[seqnames(pp2) == 'genomeWide',]$probAbberant[m]
l038_d$pp_norm = pp2[seqnames(pp2) == 'genomeWide',]$postProb_normFrac[m]
# l038_d$pp_chr5 = srat$AI_output_chr5_pp[match(l038_d$cellID,srat$cellID)]
# l038_d$pp_chr17 = srat$AI_output_chr17_pp[match(l038_d$cellID,srat$cellID)]
l038_d$pp[l038_d$annot_aug24 !='Tumour'] = NA
FeaturePlot(l038_d,'pp',cols = c(grey(0.9),'#c65102')) + NoAxes() + NoLegend()


m = match(l038_tp1@meta.data$cellID,pp2[seqnames(pp2) == 'genomeWide',]$cellID)
l038_tp1$pp = pp2[seqnames(pp2) == 'genomeWide',]$probAbberant[m]
l038_tp1$pp_norm = pp2[seqnames(pp2) == 'genomeWide',]$postProb_normFrac[m]
# l038_tp1$pp_chr5 = srat$AI_output_chr5_pp[match(l038_tp1$cellID,srat$cellID)]
# l038_tp1$pp_chr17 = srat$AI_output_chr17_pp[match(l038_tp1$cellID,srat$cellID)]
l038_tp1$pp[l038_tp1$annot_aug24 !='Tumour'] = NA
FeaturePlot(l038_tp1,'pp',cols = c(grey(0.9),'#c65102')) + NoAxes() + NoLegend()
DimPlot(l038_tp1,cells.highlight = l038_tp1$cellID[l038_tp1$annot_aug24 == 'Tumour'])
# l038_tp1 = FindClusters(l038_tp1,resolution = 1.5)
# table(l038_tp1$GATA1_status[l038_tp1$seurat_clusters %in% c(12,10,9,13,17,1,8,27,29,19) & l038_tp1$annot_aug24 != 'Tumour' & l038_tp1$GATA1_status == 'GATA1s_WT'])
# DimPlot(l038_tp1,group.by = 'seurat_clusters',label = T,repel = F,label.box = T,label.size = 2)
# DimPlot(l038_tp1,cells.highlight = l038_tp1$cellID[l038_tp1$seurat_clusters %in% c(12,10,9,13,17,1,8,27,29,19) & l038_tp1$annot_aug24 != 'Tumour' & l038_tp1$GATA1_status != 'GATA1s_WT'])
# l038_tp1$pp[l038_tp1$seurat_clusters %in% c(12,10,9,13,17,1,8,27,29,19) & l038_tp1$annot_aug24 != 'Tumour' & l038_tp1$GATA1_status != 'GATA1s_WT']


df = cbind(l038_d@meta.data,l038_d@reductions$umap@cell.embeddings)
df = cbind(l038_tp1@meta.data,l038_tp1@reductions$umap@cell.embeddings)
plotFun_AIresult = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  
  p = ggplot(df,aes(UMAP_1,UMAP_2))+
    geom_point(size=0.001,aes(col=pp))+
    scale_color_gradient(low='black',high='red',na.value = grey(0.8))+
    theme_classic(base_size = 8.3)+ 
    theme(panel.border = element_blank(),axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key.size = unit(0.45,'cm')) + 
    xlab('') + ylab('')
  
  print(p)
}

saveFig(file.path(plotDir,'Figure4D_L038.Diag_AIresult_UMAP'),plotFun_AIresult,rawData=df,width = 2.3,height = 1.6,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'Figure4D_L038.TP1_AIresult_UMAP'),plotFun_AIresult,rawData=df,width = 2.3,height = 1.6,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'Figure4B_L076.TP1_AIresult_UMAP'),plotFun_AIresult,rawData=df,width = 2.3,height = 1.6,res = 500,useDingbats = T)



# df = cbind(l038_d@meta.data,l038_d@reductions$umap@cell.embeddings)
# ggplot(df,aes(UMAP_1,UMAP_2))+
#   geom_point(aes(col=pp))+
#   scale_color_gradient2(low=grey(0.99),mid=grey(0.9),high='black')
# l038_d$AI_output = srat$AI_output[match(l038_d$cellID,srat$cellID)]
# table(l038_d$AI_output)
# DimPlot(l038_d,group.by = 'AI_output')

l038_d$AIres = df$AI_output[match(l038_d$cellID,df$cellID)]
l038_d$AIres[l038_d$AIres %in% c('?','Uncalled')] = 'Uninformative'
l038_d$AIres[l038_d$AIres %in% c('normFrac')] = 'without CNA'
l038_d$AIres[l038_d$AIres %in% c('abbFrac')] = 'with CNA'
l038_d$AIres[is.na(l038_d$AIres)] = 'Uninformative'
# l038_d$AIres = ifelse(l038_d$group %in% c('D_noCNA','TP1_noCNA'),'without CNA',
#                   ifelse(l038_d$group %in% c('D_wCNA','TP1_wCNA'),'with CNA','Uninformative'))
l038_d$AIres = factor(l038_d$AIres,c('Uninformative','with CNA','without CNA'))

group_cols = c('Uninformative'=grey(0.8),
               'without CNA' = grey(0.3),
               'with CNA' = '#A92821')

DimPlot(l038_d,group.by = 'AIres',cols = group_cols,label = F,repel = T,label.box = T) + NoAxes() + NoLegend()




l038_d = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_diagnostic.RDS')
l038_tp1 = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_TP1.RDS')
l038 = merge_seurat_objects(l038_d,l038_tp1,keepAllGenes = F,genomeVersions = c('v38','v38'))
l038 = standard_clustering(l038)

l038$group = ifelse(l038$annot_aug24 == 'Tumour' & l038$cellID %in% l038_d$cellID[l038_d$seurat_clusters %in% c(4,18,0,17)],'clone2_D',
                    ifelse(l038$annot_aug24 == 'Tumour' & l038$timePoint == 'D','clone1',
                           ifelse(l038$annot_aug24 == 'Tumour' & l038$timePoint == 'TP1' & l038$seurat_clusters %in% c(4,11,26),'clone2',
                                  ifelse(l038$annot_aug24 == 'Tumour' & l038$timePoint == 'TP1','clone2_ery','Normal'))))

DimPlot(l038,group.by = 'group',label = T)
Idents(l038) = l038$group
l038_tp1.vs.clone1 = FindMarkers(l038,ident.1 = 'clone2',ident.2 = 'clone1')
l038_tp1.vs.clone1$geneSym = rownames(l038_tp1.vs.clone1)
l038_tp1.vs.clone1$pct_diff = l038_tp1.vs.clone1$pct.1 - l038_tp1.vs.clone1$pct.2
l038_tp1.vs.clone1$comp = 'TP1_vs_Dclone1'
rownames(l038_tp1.vs.clone1) = geneMap$ensID[match(l038_tp1.vs.clone1$geneSym,geneMap$geneSym)]
l038_tp1.vs.clone1 = annotateGenes(l038_tp1.vs.clone1,geneMap = geneMap)
l038_tp1.vs.clone1$direction = ifelse(l038_tp1.vs.clone1$avg_log2FC > 0,'clone2_up','clone2_down')

write.csv(l038_tp1.vs.clone1,'L038_tp1.vs.clone1.D_markers.csv')
l038_tp1.vs.clone1.sub = l038_tp1.vs.clone1[abs(l038_tp1.vs.clone1$pct_diff) >= 0.3 & abs(l038_tp1.vs.clone1$avg_log2FC) >= 0.5,]
l038_tp1.vs.clone1.sub = l038_tp1.vs.clone1.sub[order(l038_tp1.vs.clone1.sub$pct_diff,decreasing = T),]

table(l038_tp1.vs.clone1.sub$direction)

genes = c(l038_tp1.vs.clone1.sub$geneSym[l038_tp1.vs.clone1.sub$direction == 'clone2_up'],
          l038_tp1.vs.clone1.sub$geneSym[l038_tp1.vs.clone1.sub$direction == 'clone2_down'])
DotPlot(l038,group.by = 'group',features = genes)+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))


FeaturePlot(l038,'TYROBP')

genes = c(bad.vs.good.mlds_deg$geneSym[bad.vs.good.mlds_deg$direction == 'MLDS.good_up'],
          bad.vs.good.mlds_deg$geneSym[bad.vs.good.mlds_deg$direction == 'MLDS.good_down'])

genes = bad_MLDS_markers$gene_symbol[bad_MLDS_markers$direction=='D_down' & bad_MLDS_markers$is_DEG=='DEG' & bad_MLDS_markers$number_of_comparison == 3]

length(genes)
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
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11)) + xlab('') + ylab('') 



##------------------------------------------------##
##    Quantify expression levels of key genes   ####
##------------------------------------------------##
gene =  c('H1F0','ADAMTS1','FHL2','CXCL8','EPS8',#'IFI27','GJA1',
          #'SMAD3',
          'EPS15','MEIS2',#'ARID1B', # up in R2D_vs_D(BM)
          'CD82',
          'LDLR','CD81',
          'COL18A1')
          #'CD34','EPCAM','NR4A1')#,
          #'CD34','ZFP36','CD7','SAAL1')
          #
tmp = mlds@meta.data
cnt = as.data.frame(t(mlds@assays$RNA@data[gene,]))
cnt$cellID = rownames(cnt)
cnt = pivot_longer(cnt,cols = gene,names_to = 'gene',values_to = 'norm_count')
tmp = cbind(tmp[match(cnt$cellID,tmp$cellID),],cnt[,colnames(cnt)!='cellID'])


tmp$group = ifelse(tmp$group_tmp != 'others' & tmp$disease == 'MLDS' & tmp$timePoint == 'Diagnostic' & !tmp$donorID %in% c('L038','L076'),'responsive_MLDS',
                   ifelse(tmp$group_tmp != 'others' & tmp$disease == 'MLDS' & tmp$donorID %in% c('L038','L076'),as.character(tmp$group_tmp),
                          ifelse(tmp$group_tmp != 'others' & tmp$disease == 'TAM' & tmp$timePoint == 'Diagnostic','TAM',as.character(tmp$group_tmp))))


tmp = tmp %>% dplyr::mutate(group = dplyr::case_when(
  !annot %in% c("Tum_MK?",'Tumour','Tumour_WT','unsure_Tumour') ~ 'normal',
  annot == 'Tumour' & disease == 'TAM' & timePoint == 'Diagnostic' ~ 'TAM',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & !donorID %in% c('L038','L076') ~ 'responsive_MLDS',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & donorID %in% c('L038') & cellID %in% l038_d_mdat$cellID[l038_d_mdat$group == 'clone1_D'] ~ 'L038_D_clone1',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & donorID %in% c('L038') & cellID %in% l038_d_mdat$cellID[l038_d_mdat$group == 'clone2_D'] ~ 'L038_D_clone2',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'TP1' & donorID %in% c('L038') & cellID %in% l038$cellID[l038$group == 'clone2']  ~ 'L038_TP1_clone2',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'TP1' & donorID %in% c('L038') & cellID %in% l038$cellID[l038$group == 'clone2_ery'] ~ 'clone2_ery',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & donorID %in% c('L076') & tissue == 'Blood' ~ 'L076_Blood',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & donorID %in% c('L076') & tissue == 'BM' ~ 'L076_BM',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'D.Relapse' & donorID %in% c('L076') & tissue == 'BM' ~ 'L076_D.Relapse',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'D.Relapse2' & donorID %in% c('L076') & tissue == 'BM' ~ 'L076_D.Relapse2',
  .default = 'others'
  ),
  group = factor(group,c('normal','TAM','responsive_MLDS',
                         'L038_D_clone1','L038_D_clone2','L038_TP1_clone2','clone2_ery',
                         'L076_Blood','L076_BM','L076_D.Relapse','L076_D.Relapse2','others'))) %>% 
  dplyr::filter(group != 'others')
  
table(!tmp$group %in% c('normal','TAM','responsive_MLDS',
                 'L038_D_clone1','L038_D_clone2','L038_TP1_clone2','clone2_ery',
                 'L038_Diagnostic','L038_TP1','L076_Blood','L076_BM','L076_D.Relapse','L076_D.Relapse2'))                  

# tmp$group[tmp$cellID %in% l038$cellID[l038$group == 'clone1']] = 'L038_D_clone1'
# tmp$group[tmp$cellID %in% l038$cellID[l038$group == 'clone2']] = 'L038_TP1_clone2'
# tmp$group[tmp$cellID %in% l038$cellID[l038$group == 'clone2_D']] = 'L038_D_clone2'
# tmp = tmp[tmp$group !='L038_TP1',]

tmp$group2 = ifelse(!tmp$donorID %in% c('L038','L076','L156'),'Remission',
                    ifelse(tmp$donorID == 'L076','Relapse',
                           ifelse(tmp$donorID == 'L156','Recurrent',
                                  ifelse(tmp$donorID == 'L038','Refractory','others'))))
tmp$group2[tmp$group == 'normal'] = 'Normal'
table(tmp$group2)
table(tmp$group)
tmp$group2 = factor(tmp$group2,c('Normal','Remission','Recurrent','Relapse','Refractory'))


tmp = tmp[tmp$gene %in% gene,]
tmp$gene = factor(tmp$gene,gene)
library(ggbeeswarm)
ccs = c('Normal'=grey(0.8),'Remission'='#2D4372','Refractory'=col25[2],'Relapse'=col25[5],'Recurrent'=col25[5])

#tmp = read.delim(file.path(plotDir,'Fig4E_badMLDS_markers_rawData.tsv'),sep = '\t')

plot_normalisedExpression = function(noFrame=FALSE,noPlot=FALSE){
  p = ggplot(tmp[tmp$gene %in% gene,],aes(group,norm_count,fill=group2))+
    geom_boxplot(outliers=F,outlier.shape = NA,colour='black',width=0.7,linewidth=0.3)+
    scale_fill_manual(values = ccs)+
    facet_grid(gene~group2,scales = 'free',space = 'free_x')+
    #geom_quasirandom(size=0.1,alpha=0.1)+
    theme_classic(base_size = 9)+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          strip.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(colour='black'),
          panel.spacing.y = unit(0.1,'cm'),
          axis.ticks = element_line(colour='black'),axis.title = element_text(colour='black'),
          panel.border = element_rect(fill=F,colour='black'))+xlab('')+ylab('Normalized expression')
  print(p)
}

#saveFig(file.path(plotDir,paste0('Fig4E_badMLDS_markers')),plot_normalisedExpression,rawData=tmp,width = 4.5,height = 13,res = 500,useDingbats = F)
saveFig(file.path(plotDir,paste0('Fig4E_badMLDS_markers')),plot_normalisedExpression,rawData=NULL,width = 4,height = 6,res = 500,useDingbats = F)
saveFig(file.path(plotDir,paste0('Fig4E_badMLDS_markers_hor')),plot_normalisedExpression,rawData=tmp,width = 14,height = 4.5,res = 500,useDingbats = F)
saveFig(file.path(plotDir,paste0('Fig4E_badMLDS_markers_hor')),plot_normalisedExpression,rawData=tmp,width = 10,height = 4,res = 500,useDingbats = F)


