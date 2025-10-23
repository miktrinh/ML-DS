## Combine relevant leukemia objects:
#  MLDS + MDS + DS.BALL + DS.BrainLymphoma

outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)



##------------------##
##    Libraries   ####
##------------------##
library(tidyverse)
library(Seurat)
source('~/lustre_mt22/generalScripts/utils/sc_utils.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')





##----    Create other_leuk combined object     -------##
## Import annotated combined seurat metadata
combined_mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/combinedLeuk_2404_mdat.csv',row.names = 1)
combined_mdat$cellID_og = combined_mdat$cellID
combined_mdat$finalAnn_broad = ifelse(combined_mdat$disease %in% c('MLDS','TAM'), as.character(combined_mdat$annot_mar24),as.character(combined_mdat$annot_feb24))
combined_mdat$finalAnn = combined_mdat$finalAnn_broad

##------------------##
##  1. DS B-ALL    ###
##------------------##
## check that new seurat object has all the relevant cells / samples
# dsBALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/BALL_clean_annotated_jan24.RDS')
# write.csv(dsBALL@meta.data,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/BALL_clean_annotated_jan24_mdat.csv')
# dsBALL2 = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/BALL/BALL_clean_noMTCells.RDS')
# table(dsBALL2$orig.ident,dsBALL2$cellID %in% dsBALL$cellID)
# table(dsBALL$annot_feb24,dsBALL$cellID %in% dsBALL2$cellID)
# table(dsBALL$orig.ident %in% dsBALL2$orig.ident)
# table(dsBALL$orig.ident %in% dsBALL$orig.ident)

## all looks good - let's do the merging
dsBALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/BALL/BALL_clean_noMTCells.RDS')
table(dsBALL$orig.ident)
table(combined_mdat$cellID[combined_mdat$disease %in% c('BALL','Lymphoma')] %in% dsBALL$cellID)
dsBALL_mdat = combined_mdat[combined_mdat$cellID %in% dsBALL$cellID,]
# dsBALL_mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/BALL_clean_annotated_jan24_mdat.csv',row.names = 1)
# table(dsBALL_mdat$cellID %in% dsBALL$cellID)
# dsBALL_mdat = dsBALL_mdat[dsBALL_mdat$cellID %in% dsBALL$cellID,]
colnames(dsBALL_mdat)[!colnames(dsBALL_mdat) %in% colnames(dsBALL@meta.data)]
columns_toKeep = colnames(dsBALL_mdat)[!grepl('LR|RNA_snn|UMAP|_v1|_v2|seurat_clusters|X$|GATA1|gata1|may23_ann|cellID_bc|cellID2',colnames(dsBALL_mdat)) & !colnames(dsBALL_mdat) %in% colnames(dsBALL@meta.data)]
dsBALL_mdat = dsBALL_mdat[match(dsBALL$cellID,dsBALL_mdat$cellID),columns_toKeep]

dsBALL@meta.data = cbind(dsBALL@meta.data,dsBALL_mdat)


##------------------------------##
##  2. pAML (including MDS)    ###
##------------------------------##
pAML = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/pAML/pAML_clean_noMTCells.RDS')
table(combined_mdat$cellID[combined_mdat$disease %in% c('MDS','pAML')] %in% pAML$cellID)
pAML_mdat = combined_mdat[combined_mdat$cellID %in% pAML$cellID,]
columns_toKeep = colnames(pAML_mdat)[!grepl('LR|RNA_snn|UMAP|_v1|_v2|seurat_clusters|X$|GATA1|gata1|may23_ann|cellID_bc|cellID2',colnames(pAML_mdat)) & !colnames(pAML_mdat) %in% colnames(pAML@meta.data)]
pAML_mdat = pAML_mdat[match(pAML$cellID,pAML_mdat$cellID),columns_toKeep]

pAML@meta.data = cbind(pAML@meta.data,pAML_mdat)


# ##---  old MDS  ----###
# mds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MDS/MDS_clean_annotated_tmp.RDS')
# mds$tissue = 'BM'
# mds$dataset = 'MDS'
# mds$donorID = 'L067'
# mds$disease = 'MDS'
# mds$age_yrs = '?'
# mds$sex = 'M'
# mds$mutation = '?'
# mds$timePoint = 'Diagnostic'
# mds$blastPerc = '?'
# mds$clinicalOutcome = '?'
# mds$assay = 'GEX5p'
# mds$Genotype = 'diploid'
# mds$finalAnn_broad = mds$finalAnn
# 
# 
# ##---  old pAML  ----###
# # pAML2 = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML/AML_clean_noMTcells_annotated_0923.RDS')
# pAML2 = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML/AML_clean_noMTcells_annotated_0923.RDS')
# pAML2$dataset = 'pAML'
# pAML2$disease = 'pAML'
# pAML2$Genotype = 'diploid'
# pAML2$annot_feb24 = pAML2$finalAnn_broad
# pAML2$finalAnn = pAML2$finalAnn_broad



##------- Combine DS-BALL and pAML
print('Merging BALL and pAML')
big.srat = merge_seurat_objects(dsBALL,pAML,keepAllGenes = F,genomeVersions = c('v38','v38'))
big.srat@misc$geneMap = dsBALL@misc$geneMap
big.srat@misc$geneMapCommon = dsBALL@misc$geneMapCommon
rm(dsBALL)
rm(pAML)
big.srat = standard_clustering(big.srat)

# DimPlot(big.srat,group.by = 'annot_feb24',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
# 
# 



##----------------------##
##    3.  AMKL        ####
##----------------------##

## Ellie's infant ALL paper - 8 diploid ALL + 1 AMKL
infantALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/ek12_infantALL/ek12_infantALL_clean_annotated.RDS')
infantALL$finalAnn_broad = infantALL$annot
infantALL$annot_feb24= infantALL$finalAnn_broad
infantALL$tissue = 'BM'
infantALL$dataset = 'infantALL'
infantALL$disease = '?'
infantALL$disease[infantALL$donorID == 'P9_iAML'] = 'AMKL'
infantALL$disease[infantALL$donorID != 'P9_iAML'] = 'iALL'
infantALL$finalAnn_broad[infantALL$finalAnn_broad=='Cancer'] = 'Tumour'
infantALL$age_yrs = '?'
infantALL$sex = ifelse(grepl('P8_|P7_|P1_',infantALL$donorID),'M','F')
infantALL$mutation = '?'
infantALL$timePoint = ifelse(infantALL$donorID == 'P9_iAML_TP1','TP1','Diagnostic')
infantALL$blastPerc = '?'
infantALL$clinicalOutcome = '?'
infantALL$assay = '?'
infantALL$Genotype = 'diploid'
infantALL$donorID[infantALL$donorID == 'P9_iAML_TP1'] = 'P9_iAML'


table(combined_mdat$cellID %in% infantALL$cellID,combined_mdat$disease)
table(infantALL$cellID %in% combined_mdat$cellID)
iALL_mdat = combined_mdat[combined_mdat$cellID %in% infantALL$cellID,]
columns_toKeep = colnames(iALL_mdat)[!grepl('LR|RNA_snn|UMAP|_v1|_v2|seurat_clusters|X$|GATA1|gata1|may23_ann|cellID_bc|cellID2',colnames(iALL_mdat)) & !colnames(iALL_mdat) %in% colnames(infantALL@meta.data)]
iALL_mdat = iALL_mdat[match(infantALL$cellID,iALL_mdat$cellID),columns_toKeep]

##------- Combine DS-BALL + pAML and iALL
print('Merging iALL')
colnames(big.srat@meta.data)[!colnames(big.srat@meta.data) %in% colnames(infantALL@meta.data)]
big.srat = merge_seurat_objects(big.srat,infantALL,keepAllGenes = F,genomeVersions = c('v38','v38'))
rm(infantALL)

big.srat = standard_clustering(big.srat)



##----  Save the objects  ------##
print('Merging completed - Saving the object')
mdat = cbind(big.srat@meta.data,big.srat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/combinedLeuk_mdat_2404.csv')
saveRDS(big.srat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/combinedLeuk_2404.RDS')


big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/combinedLeuk_2404.RDS')
# 
DimPlot(big.srat,group.by = 'finalAnn_broad',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()




##--------     Annotate    ---------
annotDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_scProcessing_leukaemia'
if(!dir.exists(annotDir)){
  dir.create(annotDir,recursive = T)
}


big.srat$annot_apr24 = big.srat$annot_feb24
## Annotate each cluster with NA cells
i=0
for(clust in unique(as.character(big.srat$seurat_clusters[is.na(big.srat$annot_feb24)]))){
  i=i+1
  print(sprintf('Cluster %d out of %d',i,length(unique(as.character(big.srat$seurat_clusters[is.na(big.srat$annot_feb24)])))))
  nNewCells = length(big.srat$cellID[is.na(big.srat$annot_feb24) & as.character(big.srat$seurat_clusters) == clust])
  #nNewCells = length(big.srat$cellID[big.srat$annot_feb24 == 'NA' & as.character(big.srat$seurat_clusters) == clust])
  nTot = length(big.srat$cellID[as.character(big.srat$seurat_clusters) == clust])
  p = DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$seurat_clusters == clust]) + ggtitle(clust) + NoLegend()
  p1 = DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$seurat_clusters == clust & is.na(big.srat$annot_feb24)]) + ggtitle(clust,subtitle = 'new cells')+ NoLegend()
  #p1 = DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$seurat_clusters == clust & big.srat$annot_feb24 == 'NA']) + ggtitle(clust,subtitle = 'new cells')+ NoLegend()
  print(p+p1)

  current_annot = unique(big.srat$annot_feb24[!is.na(big.srat$annot_feb24) & as.character(big.srat$seurat_clusters) == clust])
  #current_annot = unique(big.srat$annot_feb24[big.srat$annot_feb24 != 'NA' & as.character(big.srat$seurat_clusters) == clust])
  print(sprintf('Cluster %s: %d new cells out of %d total cell count, representing a fraction of %f',clust,nNewCells,nTot,nNewCells/nTot))
  if(length(current_annot) == 1){
    assign = readline(sprintf('Cluster %s: current annot is %s. Assigned cluster annot to this now? ', clust, current_annot))
  }else{
    print(sprintf('Cluster %s current annotation:',clust))
    print(table(big.srat$annot_feb24[!is.na(big.srat$annot_feb24) & as.character(big.srat$seurat_clusters) == clust]))
    #print(table(big.srat$annot_feb24[big.srat$annot_feb24 != 'NA' & as.character(big.srat$seurat_clusters) == clust]))
    assign = readline(sprintf('Cluster %s: Assigned cluster annot now? ', clust))
  }

  if(assign == 'n'){
    big.srat$annot_apr24[is.na(big.srat$annot_apr24) & big.srat$seurat_clusters == clust] = 'NA'
    #big.srat$annot_feb24[big.srat$annot_feb24 =='NA' & big.srat$seurat_clusters == clust] = 'NA'
  }else if(assign == 'y'){
    big.srat$annot_apr24[is.na(big.srat$annot_apr24) & big.srat$seurat_clusters == clust] = current_annot
    #big.srat$annot_feb24[big.srat$annot_feb24 =='NA' & big.srat$seurat_clusters == clust] = current_annot
  }else{
    big.srat$annot_apr24[is.na(big.srat$annot_apr24) & big.srat$seurat_clusters == clust] = assign
    #big.srat$annot_feb24[big.srat$annot_feb24 =='NA' & big.srat$seurat_clusters == clust] = assign
  }
}


complete_mdat = big.srat@meta.data[!is.na(big.srat$disease),]
big.srat$donorID[is.na(big.srat$disease)] = complete_mdat$donorID[match(big.srat$orig.ident[is.na(big.srat$disease)],complete_mdat$orig.ident)]
big.srat$disease[is.na(big.srat$disease)] = complete_mdat$disease[match(big.srat$orig.ident[is.na(big.srat$disease)],complete_mdat$orig.ident)]
big.srat$dataset[is.na(big.srat$dataset)] = complete_mdat$dataset[match(big.srat$orig.ident[is.na(big.srat$dataset)],complete_mdat$orig.ident)]
big.srat$timePoint[is.na(big.srat$timePoint)] = complete_mdat$timePoint[match(big.srat$orig.ident[is.na(big.srat$timePoint)],complete_mdat$orig.ident)]
DimPlot(big.srat,group.by = 'timePoint',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$annot_apr24 == 'Tumour' & big.srat$disease == 'MDS'])
DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$timePoint == 'RelapseD29' & big.srat$donorID == 'L001'])
DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$annot_apr24 == '0'])

write.csv(big.srat@meta.data,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/combinedLeuk_mdat_2405.csv')

##---- Subcluster T / NK cells -------##
s = subset(big.srat,subset = annot_apr24 %in% c(unique(big.srat$annot_apr24)[29:53],'T_CD8','T_CD4','T_gd','NK_T','maybe_NK','NK','T_MAIT','T.cell','?'))
s = standard_clustering(s,clusteringRes = 5)
DimPlot(s,group.by = 'seurat_clusters',label.size = 3,label = T,repel = T,label.box = T)+NoLegend()
DimPlot(s,group.by = 'annot_apr24',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
DimPlot(s,cells.highlight = s$cellID[s$percent.mt > 30])
library(SoupX)
qm = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)


FeaturePlot(s,c('CD8A','CD8B','CD3G','CD3D','CD3E','CD4'))
FeaturePlot(s,c('TRDV2','TRGV9'))
FeaturePlot(s,c('CD3D','CD3E'))
FeaturePlot(s,c('KLRD1','NKG7','CD3D','CD3E','TRAV1-2','SLC4A10'))

DimPlot(s,cells.highlight = s$cellID[s$seurat_clusters %in% as.character(c('42','84','41','77','29'))])
s$annot_apr24 = '?'

non_NKT = c(62,45,56,33,58,72,34,83,40,44,84,75,70,87,89,65,68,42,"74","41","30","47","78","77","86","18","29","79","82","80","88","71","85","69","64","14","43","59","4","21")
doublets = c(32,52,76,88)
NKT_clusters = as.character(unique(s$seurat_clusters[!as.character(s$seurat_clusters) %in% c(non_NKT,doublets,'39')]))

s$annot_apr24[s$seurat_clusters %in% NKT_clusters] = 'T.cell'
tCD8 = WhichCells(s,expression = (CD8A > 1.5 | CD8B > 1.5))
s$annot_apr24[s$seurat_clusters %in% NKT_clusters & s$cellID %in% tCD8] = 'T_CD8'
s$annot_apr24[s$seurat_clusters %in% c('11','51','24','48','61')] = 'NK'
s$annot_apr24[s$seurat_clusters %in% c('63','10')] = 'T_gd'
s$annot_apr24[s$seurat_clusters %in% c('67')] = 'T_MAIT'

other_clusters = non_NKT[!non_NKT %in% c(NKT_clusters,doublets,'39')]
s$annot_apr24[s$seurat_clusters %in% c('83')] = 'Mono.Mac'
s$annot_apr24[s$seurat_clusters %in% c('44')] = 'Glial.cell'
s$annot_apr24[s$seurat_clusters %in% c('42','84','41','77','29')] = 'B.cell'
s$annot_apr24[s$seurat_clusters %in% c('45','56')] = 'Cancer'
s$annot_apr24[s$seurat_clusters %in% doublets] = 'doublets'

newLab = s@meta.data[s$annot_apr24 != '?',]
## Add back to big.srat
big.srat$annot_apr24[big.srat$cellID %in% newLab$cellID] = newLab$annot_apr24[match(big.srat$cellID[big.srat$cellID %in% newLab$cellID],newLab$cellID)]
table(big.srat$annot_apr24[!big.srat$cellID %in% newLab$cellID])
DimPlot(s,cells.highlight = big.srat$cellID[!big.srat$cellID %in% newLab$cellID & big.srat$annot_apr24 %in% c('T/NK')])
big.srat$annot_apr24[big.srat$annot_apr24 %in% c('maybe_NK.T','NK_T') & !big.srat$cellID %in% newLab$cellID ] = 'NK'
big.srat$annot_apr24[big.srat$annot_apr24 %in% c('T_CD4','T_CD8','T.cell','T/NK') & !big.srat$cellID %in% newLab$cellID] = 'doublets'


table(big.srat$annot_apr24)
##---- Subcluster T / NK cells -------##
s = subset(big.srat,subset = cellID %in% 
             c(big.srat$cellID[big.srat$annot_apr24 %in% 
                               c(unique(big.srat$annot_apr24)[29:52],'naive.B','HSC_MPP','Plasma.cell','pDC','B.cell','Mono_CD14','?')],
               big.srat$cellID[big.srat$disease == 'iALL']))
s = standard_clustering(s,clusteringRes = 2)
DimPlot(s,group.by = 'seurat_clusters',label.size = 3,label = T,repel = T,label.box = T)+NoLegend()
DimPlot(s,group.by = 'disease',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
qm = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)
DimPlot(big.srat,cells.highlight = s$cellID[s$seurat_clusters %in% c()])

s$annot_may24 = '?'
s$annot_may24[s$seurat_clusters %in% c(1,2,42,30,39,4)] = 'B.cell'
s$annot_may24[s$seurat_clusters %in% c(55)] = 'doublets'
s$annot_may24[s$seurat_clusters %in% c(59,58,37)] = 'Plasma.cell'
s$annot_may24[s$seurat_clusters %in% c(32,28,
                                       0,7,24,50,
                                       48,
                                       3,8,5,
                                       11
                                       )] = 'Cancer'
s$annot_may24[s$seurat_clusters %in% c(53)] = 'unknown'

table(big.srat$annot_apr24[big.srat$cellID %in% s$cellID[s$seurat_clusters %in% c(61)]])
DimPlot(s,cells.highlight = s$cellID[s$annot_apr24 == 'B.cell'])
DimPlot(s,cells.highlight = s$cellID[s$seurat_clusters %in% c(40)])
DimPlot(s,cells.highlight = s$cellID[as.character(s$seurat_clusters) %in% other_clusters])
FeaturePlot(s,'percent.mt')

## For each seurat cluster, calculate fraction of cells CD8A|CD8B|CD4 positive
nCell_expressed_perCluster = do.call(rbind,lapply(split(s$cellID,s$seurat_clusters),function(e){rowSums(s@assays$RNA@counts[c('CD8A','CD8B','CD4'),e] >= 0.5)})) %>% as.data.frame()
nCell_expressed_perCluster$seurat_clusters = rownames(nCell_expressed_perCluster)
nCell_perCluster = table(s$seurat_clusters)
nCell_expressed_perCluster$totalCell = nCell_perCluster[match(nCell_expressed_perCluster$seurat_clusters,names(nCell_perCluster))]
nCell_expressed_perCluster$frac_CD8A = nCell_expressed_perCluster$CD8A / nCell_expressed_perCluster$totalCell
nCell_expressed_perCluster$frac_CD8B = nCell_expressed_perCluster$CD8B / nCell_expressed_perCluster$totalCell
nCell_expressed_perCluster$frac_CD4 = nCell_expressed_perCluster$CD4 / nCell_expressed_perCluster$totalCell
DimPlot(s,cells.highlight = s$cellID[s$seurat_clusters %in% nCell_expressed_perCluster$seurat_clusters[nCell_expressed_perCluster$frac_CD8A >= 0.5]])

big.srat$annot_feb24[big.srat$annot_feb24 == 'T_CD4' & !is.na(big.srat$finalAnn) & big.srat$finalAnn %in% c('T_CD8','?') ] = 'T_CD8'
big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$finalAnn) & big.srat$finalAnn %in% c('T_CD4') ] = 'T_CD8'
big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$finalAnn) & big.srat$finalAnn %in% c('MonoMac') ] = 'Mono_CD14'
big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$finalAnn) & big.srat$finalAnn %in% c('NK') & big.srat$seurat_clusters == 37 ] = 'MK'
big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$finalAnn) & big.srat$finalAnn %in% c('NK') & big.srat$seurat_clusters == 19 ] = 'NK'
big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$finalAnn) & big.srat$finalAnn %in% c('Tumour') ] = 'Cancer'
big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$finalAnn) & big.srat$finalAnn %in% c('naive.B') ] = 'naive.B'
big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$finalAnn) & big.srat$finalAnn %in% c('T_CD8') ] = 'T_CD8'

View(as.data.frame(table(big.srat$annot_feb24[big.srat$donorID == 'L067'],big.srat$finalAnn[big.srat$donorID == 'L067'])))
DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$annot_feb24 == 'naive.B' & !is.na(big.srat$finalAnn) & big.srat$finalAnn =='?'])
DimPlot(mds,cells.highlight = big.srat$cellID[big.srat$annot_feb24 == 'T_CD4' & !is.na(big.srat$finalAnn) & big.srat$finalAnn =='?'])

big.srat$group_tmp = paste0(big.srat$annot_feb24,':',big.srat$finalAnn)

Idents(big.srat) = big.srat$group_tmp
keyGenes = unique(c('HLA-DRA','CLEC9A','CD34','CD38','HLF','SPINK2','MLLT3','PRSS57', # HSC_MPP
                    'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1',#MEMP
                    'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
                    #'KIT',
                    'CTSG',	'PRTN3', # CMP
                    'AZU1','MPO','FLT3','PTPRC', # GMP
                    'ZBTB16','LTB', 'CD52',# ILC precursor
                    'IL7R','CD3D','GZMA','CD4','CD8A', #Early.lymphoid_T.lymphocyte
                    'TRDV2','TRGV9', # gamma delta T-cell
                    'SLC4A10','TRAV1-2', #MAIT t-cell
                    'PRF1', # effector T cells
                    'FOXP3',	'CDH1', # regulatory T cells
                    'NKG7','KLRD1', #NK
                    'IGLL1','CD99', # Pre-pro
                    'DNTT','CD79B','VPREB1','EBF1','CD19','RAG1',# pro-b-cells
                    'MME','CD79A',# pre-b-cells
                    'TCL1A','MME','RAG1','MS4A1',  # B-cells
                    'ITGB2','SELL','ITGAM','CD14','CCR2','FCGR3A','CX3CR1',# Monocytes
                    'CD68','MSR1',#Macrophages
                    'IRF8',	'CLEC10A', #DC.precursor
                    'ITGAX','CD1C','CLEC9A','THBD',#DC1
                    'CLEC10A', # DC2 
                    'IL3RA', 'CLEC4C',#pDC
                    'CD27',#plasma cell?
                    'CSF2RB','HDC','SERPINB1','TPSAB1','KIT', # Mast.cell
                    'PF4','ITGA2B', #Megakaryocyte
                    'PPBP','TUBB1', # Platelets
                    'GATA1','KLF1', # early Erythroid
                    'ALAS2', # mid.erythroid
                    'HBA1','BPGM', # late.erythroid
                    'CSF3R','FPR1','FCGR3B','NAMPT','MNDA', # NEUTROPHILS
                    'DEFA3','DEFA4','CAMP', 'LCN2',
                    'ESAM','PECAM1', 'KDR', 'PTPRB', 'PLVAP', # Endothelium
                    'DCN','SERPINF1','COL3A1','BGN','COL1A1','VIM','ECM1', # endosteal fibroblast / fibroblast
                    'APOA1','SCD','ALB','TTR' # Hepatocyte
))
DotPlot(big.srat,idents = c('NA:NK'),group.by = 'seurat_clusters',
        features = keyGenes)+RotatedAxis()




##------ Annotate neuronal clusters in L062 Lymphoma  -------
library(SoupX)
m = quickMarkers(big.srat@assays$RNA@counts,big.srat$seurat_clusters)

DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$seurat_clusters == 47])

big.srat$annot_feb24[big.srat$seurat_clusters %in% c(64,65,61,57)] = 'Neuronal.cells'
big.srat$annot_feb24[big.srat$seurat_clusters %in% c(47)] = 'Astrocytes'




##-------------------------##
##    Save the object    ####
##-------------------------##
## Finalizing annotation
table(big.srat$annot_feb24)
big.srat$annot_feb24[big.srat$annot_feb24 == 'Cancer'] = 'Tumour'
big.srat$annot_feb24[is.na(big.srat$annot_feb24)] = 'NA'
big.srat$finalAnn_broad = as.character(big.srat$annot_feb24)
big.srat$finalAnn_broad[big.srat$finalAnn_broad == 'MEMP_MEP'] = 'MEP'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('naive.B','pre.B.cell','pro.B.cell','transitional.B.cell','Plasma.cell')] = 'B lineage'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('DC2','DC1','pDC')] = 'Dendritic cells'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('EE','ME','LE')] = 'Erythroblasts'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('HSC_MPP','CMP_GMP','MEP','MEMP_MEP')] = 'HSC & prog.'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('MK','MK_WT')] = 'Megakaryocytes'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('Mono_CD14','Mono_CD16')] = 'Monocyte/Macrophage'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('Neutrophil','Myelocytes','activated_neutrophil')] = 'Neutrophil'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('T_CD4','T_CD8','NK','T_gd','T_reg','NK_T','T_MAIT')] = 'T/NK lineage'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('Tumour')] = 'Tumour'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('Neuronal.cells','Astrocytes')] = 'Neuronal'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('doublets','lowQual','NA')] = 'lowQual'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('Tumour_maybe','Normal','P9_iAML_TP1')] = 'others'

big.srat$Genotype[big.srat$donorID == 'L156'] = 'T21'

## Add alleleIntegrator results
df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_AIRes_240131.csv')
big.srat$AI_res = df$AI_output[match(big.srat$cellID,df$cellID)]
big.srat$finalAnn_broad[big.srat$donorID == 'L038' & big.srat$finalAnn_broad == 'Tumour' & big.srat$AI_res == 'abbFrac' & big.srat$timePoint == 'TP1'] = 'TumourCNA_TP1'
big.srat$finalAnn_broad[big.srat$donorID == 'L038' & big.srat$finalAnn_broad == 'Tumour' & big.srat$AI_res == 'abbFrac' & big.srat$timePoint == 'Diagnostic'] = 'TumourCNA_D'

big.srat@meta.data = big.srat@meta.data[,colnames(big.srat@meta.data) != 'group_tmp']
umap = big.srat@reductions$umap@cell.embeddings %>% as.data.frame()
umap$cellID = big.srat$cellID
mdat = merge(big.srat@meta.data,umap,by='cellID',all=T)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24_mdat.csv')
saveRDS(big.srat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24.RDS')


big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24.RDS')


















##-------------##
##  1. MLDS   ###
##-------------##
# #mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_clean_annotated_noUnknowns_jan24.RDS')
# #mlds$dataset[mlds$donorID == 'L156'] = 'GOSH'
# #mlds$finalAnn_broad = as.character(mlds$annot_jan24)
# mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
# mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
# mlds$broadLineage = mdat$broadLineage[match(mlds$cellID,mdat$cellID)]
# mlds$timePoint[mlds$orig.ident %in% c('MY.200531.14635833')] = 'D.Relapse'
# mlds$finalAnn_broad = as.character(mlds$annot_mar24)




# ## 5. Combine DS-BALL and MLDS
# mlds@meta.data = mlds@meta.data[,!colnames(mlds@meta.data) %in% c('group','RNA_snn_res.2','RNA_snn_res.10')]
# colnames(mlds@meta.data)[!colnames(mlds@meta.data) %in% colnames(dsBALL@meta.data)]
# 
# big.srat = merge_seurat_objects(mlds,dsBALL,keepAllGenes = F,genomeVersions = c('v38','v38'))
# big.srat$annot_feb24 = mdat$annot_feb24[match(big.srat$cellID,mdat$cellID)]
# colnames(mdat)[!colnames(mdat) %in% colnames(big.srat@meta.data)]
# 
# rm(mlds)
# # rm(tgt.srat)
# rm(dsBALL)
# 
# 
# 
# 
# colnames(big.srat@meta.data)[!colnames(big.srat@meta.data) %in% colnames(mds@meta.data)]
# 
# 
# ## 6. Combine DS-BALL and MLDS and MDS
# colnames(big.srat@meta.data)[!colnames(big.srat@meta.data) %in% colnames(mds@meta.data)]
# big.srat = merge_seurat_objects(big.srat,mds,keepAllGenes = F,genomeVersions = c('v38','v38'))
# big.srat = standard_clustering(big.srat)






# DimPlot(big.srat,group.by = 'annot_feb24',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
# 
# 




# ##-------------------------------------------##
# ##    March 2024: remove i-BALL, add AML   ####
# ##-------------------------------------------##
#big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24.RDS')

## Combine big.srat and new MLDS
# table(big.srat$cellID %in% mlds$cellID)
# table(mlds$cellID %in% big.srat$cellID)
# 
# mlds = subset(mlds,subset = cellID %in% mlds$cellID[!mlds$cellID %in% big.srat$cellID])
# print(table(mlds$cellID %in% big.srat$cellID))
# 
# big.srat = merge_seurat_objects(big.srat,mlds,keepAllGenes = F,genomeVersions = c('v38','v38'))
# 
# 
# rm(mlds)

## Combine big.srat and pAML
# colnames(big.srat@meta.data)[!colnames(big.srat@meta.data) %in% colnames(pAML@meta.data)]
# big.srat = merge_seurat_objects(big.srat,pAML,keepAllGenes = F,genomeVersions = c('v38','v38'))
# 
# rm(pAML)
# 
# 
# big.srat = standard_clustering(big.srat)
# 
# ## Import final annotation to add
# mdat_latest = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/combinedLeuk_2404_mdat.csv')
# big.srat$annot_feb24 = mdat_latest$annot_feb24[match(big.srat$cellID,mdat_latest$cellID)]
# big.srat$broadLineage = mdat_latest$broadLineage[match(big.srat$cellID,mdat_latest$cellID)]
# big.srat$annot_apr24 = big.srat$annot_feb24
# big.srat$annot_apr24[big.srat$disease %in% c('MLDS','TAM')] = big.srat$annot_mar24[big.srat$disease %in% c('MLDS','TAM')]  
# 
# mdat = cbind(big.srat@meta.data,big.srat@reductions$umap@cell.embeddings)
# write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/combinedLeuk_2404_v2_mdat.csv')
# saveRDS(big.srat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/combinedLeuk_2404_v2.RDS')
# 
# 
# 
# # DimPlot(big.srat,group.by = 'disease',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
# # big.srat$disease[big.srat$dataset == 'pAML'] = 'pAML'
# # write.csv(big.srat@meta.data,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_mar24_mdat.csv')
# # saveRDS(big.srat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24.RDS')
# # mergedSrat = big.srat
# # 
# # 
# # 
# # big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24.RDS')
# 
