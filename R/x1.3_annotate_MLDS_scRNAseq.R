##----   Processing leukaemia scRNAseq datasets    -----##

outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)

##------------------##
##    Libraries   ####
##------------------##
library(Seurat)
library(tidyverse)
source('~/lustre_mt22/generalScripts/utils/sc_utils.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')


##----------------------------------##
##       Annotation - ML-DS       ####
##----------------------------------##
skipIfExists=T

## 1. Import object
## Cluster new seurat object
cleanSrat_fp = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_noMTCells.RDS'
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)
}
## for some reason, RNA14831998 and RNA14831999 - most cells got flagged as doublets by scrublets... might be due to the fact that majority of cells are blasts --> only 1 cell types...
cleanSrat@misc = cleanSrat@misc[names(cleanSrat@misc) != 'preQC']
cleanSrat = standard_clustering(cleanSrat)
DimPlot(cleanSrat,label = T,repel = T,label.size = 3)+NoLegend() + NoAxes()


## 2. Add metadata
# Import metadata
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2406_mdat.csv',row.names = 1)
n_distinct(cleanSrat$orig.ident)
table(cleanSrat$orig.ident %in% mdat$orig.ident)
table(mdat$orig.ident %in% cleanSrat$orig.ident)
table(mdat$orig.ident[!mdat$orig.ident %in% cleanSrat$orig.ident])
table(mdat$cellID %in% cleanSrat$cellID)

# Add to cleanSrat seurat object
# Subset mdat to only keep cells present in cleanSrat
table(cleanSrat$orig.ident,cleanSrat$cellID %in% mdat$cellID)
table(mdat$annot_jun24,mdat$cellID %in% cleanSrat$cellID)
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[!cleanSrat$cellID %in% mdat$cellID])
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$nCount_RNA > 1000])
FeaturePlot(cleanSrat,c('nCount_RNA','nFeature_RNA','percent.mt','GATA1'))

mdat = mdat[mdat$cellID %in% cleanSrat$cellID,]

## Select columns to be added to cleanSrat
colnames(mdat)[!colnames(mdat) %in% colnames(cleanSrat@meta.data)]
columns_toKeep = colnames(mdat)[!grepl('LR|RNA_snn|UMAP|_v1|_v2',colnames(mdat)) & !colnames(mdat) %in% colnames(cleanSrat@meta.data)]
columns_toKeep

mdat = mdat[match(cleanSrat$cellID,mdat$cellID),columns_toKeep]
cleanSrat@meta.data = cbind(cleanSrat@meta.data,mdat)
cleanSrat@meta.data = cleanSrat@meta.data[,colnames(cleanSrat@meta.data) != 'X']
table(cleanSrat$cellID == cleanSrat$cellID_og)


## Add missing data
columns_na = unlist(sapply(1:ncol(cleanSrat@meta.data),function(i){if(sum(is.na(cleanSrat@meta.data[,i])) > 0){return(colnames(cleanSrat@meta.data)[i])}}))

projMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'MLDS_GOSH')
projMani = projMani[!is.na(projMani$assay) & projMani$assay == "5' V2 Dual Index",]
table(cleanSrat$orig.ident %in% gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID)))
table(cleanSrat$orig.ident[!cleanSrat$orig.ident %in% gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID))])

## donorID
cleanSrat$donorID = projMani$donorID[match(cleanSrat$orig.ident, gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID)))]

## age
cleanSrat$age_yrs = projMani$`age (months)`[match(cleanSrat$orig.ident, gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID)))]
# convert to years
cleanSrat$age_yrs[cleanSrat$age_yrs == '12'] = 1
cleanSrat$age_yrs[cleanSrat$age_yrs == '14'] = round(14/12,digits = 1)
cleanSrat$age_yrs[cleanSrat$age_yrs == '19'] = round(19/12,digits = 1)
cleanSrat$age_yrs[cleanSrat$age_yrs == '20'] = round(20/12,digits = 1)
cleanSrat$age_yrs[cleanSrat$age_yrs == '24'] = round(24/12,digits = 1)
cleanSrat$age_yrs[cleanSrat$age_yrs == '8'] = round(8/12,digits = 1)
cleanSrat$age_yrs[cleanSrat$age_yrs == '2y, 10m'] = 2+round(10/12,digits = 1)
cleanSrat$age_yrs[cleanSrat$age_yrs == '3y, 7m'] = 3+round(7/12,digits = 1)
cleanSrat$age_yrs[cleanSrat$age_yrs == "1y,1m, 21d"] = 1+round(1.5/12,digits = 1)
cleanSrat$age_yrs[cleanSrat$age_yrs == '1d'] = 0
cleanSrat$age_yrs[cleanSrat$age_yrs == '20'] = 0
cleanSrat$age_yrs[is.na(cleanSrat$age_yrs)] = '?'

## Sex
avgExpr = AverageExpression(cleanSrat,group.by = 'donorID',features = c('XIST','RPS4Y1'))
avgExpr = as.data.frame(t(avgExpr$RNA))
avgExpr$sex = ifelse(avgExpr$XIST > avgExpr$RPS4Y1,'F','M')
avgExpr$sex_current = cleanSrat$sex[match(rownames(avgExpr),cleanSrat$donorID)]
cleanSrat$sex = avgExpr$sex[match(cleanSrat$donorID,rownames(avgExpr))]

## tissue
tissue = projMani$Tissue[match(cleanSrat$orig.ident, gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID)))]
tissue[tissue == 'Bone Marrow'] = 'BM'
tissue[tissue == 'PBMC'] = 'Blood'
cleanSrat$tissue = tissue

## TimePoint
tp = projMani$`Point in treatment`[match(cleanSrat$orig.ident, gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID)))]
tp = gsub('postChemo ','',tp)
tp[tp=='Relapse'] = 'D.Relapse'
tp[tp=='Relapse2_D0'] = 'D.Relapse2'
cleanSrat$timePoint = tp

#clinicalOutcome
cleanSrat$clinicalOutcome = projMani$`Clinical outcome`[match(cleanSrat$orig.ident, gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID)))]
cleanSrat$blastPerc = projMani$`Blast Percentage`[match(cleanSrat$orig.ident, gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID)))]

# Assay
cleanSrat$assay = 'GEX5p'

#dataset
cleanSrat$dataset[!grepl('^CC',cleanSrat$donorID)] = 'GOSH'
cleanSrat$dataset[grepl('^CC',cleanSrat$donorID)] = 'Hennings'
cleanSrat$Genotype = 'T21'
cleanSrat$disease = projMani$Disease[match(cleanSrat$orig.ident, gsub('_','.',gsub('SB_|ALeuk_','',projMani$sangerSampleID)))]



DimPlot(cleanSrat,group.by = 'donorID',cols = c(col25,pal34H), label = T,repel = T,label.size = 3,label.box = T)+NoLegend() + NoAxes()
DimPlot(cleanSrat,group.by = 'seurat_clusters',label.box = T,label = T,repel = T)  +NoLegend()
DimPlot(cleanSrat,group.by = 'Phase',label.box = T,label = F,repel = T,cols = col25)  +NoAxes()


FeaturePlot(cleanSrat,c('GATA1','KLF1','XIST','RPS4Y1'))

mdat = cbind(cleanSrat@meta.data,cleanSrat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_mdat.csv')




##  3. Annotate new cells
cleanSrat$annot_aug24 = cleanSrat$annot_jun24
nCluster = n_distinct(as.character(cleanSrat$seurat_clusters[is.na(cleanSrat$annot_aug24)]))
i=0
for(clust in unique(as.character(cleanSrat$seurat_clusters[is.na(cleanSrat$annot_aug24)]))){
  i=i+1
  print(sprintf('%s / %s',i,nCluster))
  nNewCells = length(cleanSrat$cellID[is.na(cleanSrat$annot_jun24) & as.character(cleanSrat$seurat_clusters) == clust])
  nTot = length(cleanSrat$cellID[as.character(cleanSrat$seurat_clusters) == clust]) 
  p = DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$seurat_clusters == clust]) + ggtitle(clust)
  p2 = DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$seurat_clusters == clust & is.na(cleanSrat$annot_jun24)]) + ggtitle(clust,subtitle = 'New cells')
  print(p+p2)
  current_annot = unique(cleanSrat$annot_jun24[!is.na(cleanSrat$annot_jun24) & as.character(cleanSrat$seurat_clusters) == clust])
  print(sprintf('Cluster %s: %d new cells out of %d total cell count, representing a fraction of %f',clust,nNewCells,nTot,nNewCells/nTot))
  if(length(current_annot) == 1){
    assign = readline(sprintf('Cluster %s: current annot is %s. Assigned cluster annot to this now? ', clust, current_annot))
  }else{
    print(sprintf('Cluster %s current annotation:',clust))
    print(table(cleanSrat$annot_jun24[!is.na(cleanSrat$annot_jun24) & as.character(cleanSrat$seurat_clusters) == clust]))
    assign = readline(sprintf('Cluster %s: Assigned cluster annot now? ', clust))
  }
  
  if(assign == 'n'){
    cleanSrat$annot_aug24[is.na(cleanSrat$annot_aug24) & cleanSrat$seurat_clusters == clust] = 'NA'
  }else if(assign == 'y'){
    cleanSrat$annot_aug24[is.na(cleanSrat$annot_aug24) & cleanSrat$seurat_clusters == clust] = current_annot
  }else{
    cleanSrat$annot_aug24[is.na(cleanSrat$annot_aug24) & cleanSrat$seurat_clusters == clust] = assign
  }
}
#cleanSrat$annot_jun24[cleanSrat$donorID == 'CC1' & cleanSrat$annot_jun24 == ''] = 'NA'

DimPlot(cleanSrat,group.by = 'annot_aug24',cols = c(col25,pal34H),label = T,repel = T,label.box = T,label.size = 3)+NoLegend()

##------    Save tmp metadata annotation  -------##
mdat = cbind(cleanSrat@meta.data,cleanSrat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_mdat.csv')
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_mdat.csv',row.names = 1) 


##--------------------------------##
##  Add celltypist annotation   ####
##--------------------------------##
celltypist = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/annotation_wCT/Immune_All_Low.REF_LRwCT_MLDS.tgt_output_results.csv')
cleanSrat$LR_predicted_label_ImmuneAllLow = celltypist$LR_predicted_label_ImmuneAllLow[match(cleanSrat$cellID,celltypist$cellID)]
cleanSrat$LR_softmax_predicted_label_ImmuneAllLow = celltypist$LR_softmax_predicted_label_ImmuneAllLow[match(cleanSrat$cellID,celltypist$cellID)]

##--------------------------------##
##  Subclustering MLDS object   ####
##--------------------------------##
library(SoupX)
keyMarkers = unique(c('CA1','CLEC9A','CD34','CD38','HLF','SPINK2','MLLT3','PRSS57', # HSC_MPP
                      'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1',#MEMP
                      'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
                      #'KIT',
                      'CTSG',	'PRTN3', # CMP
                      'AZU1','MPO','FLT3','PTPRC', # GMP
                      'ZBTB16','LTB', 'CD52',# ILC precursor
                      'IL7R','CD3D','GZMA','CD4',
                      'FHIT', # naive cd4
                      'CD8A', #Early.lymphoid_T.lymphocyte
                      'CRTAM',
                      'TRDV2','TRGV9', # gamma delta T-cell
                      'SLC4A10','TRAV1-2', #MAIT t-cell
                      'TRAV8-2',
                      'PRF1', # effector T cells
                      'FOXP3',	'CDH1', # regulatory T cells
                      'NKG7','KLRD1', #NK
                      'IGLL1','CD99', # Pre-pro
                      'DNTT','CD79B','VPREB1','EBF1','CD19','RAG1',# pro-b-cells
                      'MME','CD79A',# pre-b-cells
                      'TCL1A','MME','RAG1','MS4A1',  # B-cells
                      'ITGB2','SELL','ITGAM','CD14','CCR2','FCGR3A','CX3CR1',# Monocytes (SELL is a good monocyte vs macro marker)
                      'S100A8', 'CD52', 'MS4A6A', 'CD14', 'CXCR4', 'CCR2', 'IL1B',
                      'CD68','MSR1',#Macrophages
                      'LYVE1', 'C1QA', 'CD163', 'SPP1', 'FCGR3A', 'FCGR3B', 'FCGR1A', 'FCGR1B', #MACROPHAGE
                      'FABP3',#Kupffer
                      'IRF8',	'CLEC10A','ACY3', 'TIFAB', 'KIF17', #DC.precursor
                      'ITGAX','CLEC9A','CD1C','THBD','BATF3', 'ANPEP',#DC1
                      'CLEC10A', # DC2
                      'CCR7', 'LAMP3', 'IL7R', #Migratory DC
                      'IL3RA', 'CLEC4C','LILRA4',#pDC
                      'CD27',#plasma cell?
                      'CSF2RB','HDC','SERPINB1','TPSAB1','KIT', # Mast.cell
                      'PF4','ITGA2B', #Megakaryocyte
                      'PPBP','TUBB1', # Platelets
                      'GATA1','KLF1','APOC1', # early Erythroid
                      'ALAS2', # mid.erythroid
                      'HBA1','BPGM', # late.erythroid
                      'CSF3R','FPR1','FCGR3B','NAMPT','MNDA', # NEUTROPHILS
                      'DEFA3','DEFA4', # pro-myelocytes
                      'CAMP','LCN2', #myelocytes
                      'CXCR2','CSF3R','FCGR3B','FUT4', # Neutrophil
                      'ESAM','PECAM1', 'KDR', 'PTPRB', 'PLVAP', # Endothelium
                      'DCN','SERPINF1','COL3A1','BGN','COL1A1','VIM','ECM1' # endosteal fibroblast / fibroblast
                      #'APOA1','SCD','ALB','TTR' # Hepatocyte
))


##------    NK/T lineage  ------##
s = subset(cleanSrat,subset = annot_aug24 %in% c('-','?','HSC_MPP','NK','NK_T','T_CD4','T_CD8','T_gd','T_reg'))
s = standard_clustering(s)
qm = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)
DimPlot(s,group.by = 'annot_aug24_new',cols=c(col25,pal34H),label = T,repel = T,label.box = T,label.size = 2) + NoLegend()
FeaturePlot(s,c('CD3D','CD4','CD8A','CD3E'))

DimPlot(cleanSrat,cells.highlight = t_cells)
DotPlot(s,group.by = 'seurat_clusters',features = keyMarkers) + 
  RotatedAxis() + 
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'right') + xlab('') + ylab('')

s$annot_aug24 = cleanSrat$annot_aug24[match(s$cellID,cleanSrat$cellID)]
s$annot_aug24_new = '?'
s$annot_aug24_new[s$annot_aug24 %in% c('?','-','HSC_MPP')] = s$annot_aug24[s$annot_aug24 %in% c('?','-','HSC_MPP')]
s$annot_aug24_new[s$seurat_clusters %in% c(9,32)] = 'T_gd'
cd8 = WhichCells(s,cells = s$cellID[s$seurat_clusters %in% c(10,13,15,21,24:27,30)],expression = (CD3D>0 | CD3E >0) & (CD8A>0|CD8B>0))
s$annot_aug24_new[s$cellID %in% cd8] = 'T_CD8'
s$annot_aug24_new[s$seurat_clusters %in% c(10,13,15,21,24,25,27,30)] = 'T_CD8'
s$annot_aug24_new[s$seurat_clusters %in% c(1,2,6,20,23,31)] = 'T_CD4'
s$annot_aug24_new[s$seurat_clusters %in% c(8,14,28)] = 'NK'
s$annot_aug24_new[s$seurat_clusters %in% c(21,30,25)] = 'T_CD8'
t_cells = WhichCells(s,cells = s$cellID[s$annot_aug24_new == '?'],expression = (CD3D>0 | CD3E >0))
s$annot_aug24_new[s$cellID %in% t_cells] = 'T_cells'

## Add to big.srat
cleanSrat$annot_aug24_new = cleanSrat$annot_aug24
cleanSrat$annot_aug24_new[cleanSrat$cellID %in% s$cellID] = s$annot_aug24_new[match(cleanSrat$cellID[cleanSrat$cellID %in% s$cellID],s$cellID)]





##------    Myeloid lineage  ------##
s = subset(cleanSrat,subset = annot_aug24 %in% c('-','?','HSC_MPP','activated_neutrophil','CMP_GMP',
                                                 'DC1','DC2','Mono_CD14','Mono_CD16','Myelocytes','Neutrophil','pDC','MEP'))
s = standard_clustering(s,clusteringRes = 0.5)
qm = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)
DimPlot(s,group.by = 'annot_aug24_new',cols=c(col25,pal34H),label = T,repel = T,label.box = T,label.size = 4) + NoLegend()
s$group_tmp = ifelse(s$seurat_clusters %in% c(17),paste0(s$annot_aug24,'_test'),s$annot_aug24)
DotPlot(s,group.by = 'group_tmp',features = keyMarkers) + 
  RotatedAxis() + 
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'right') + xlab('') + ylab('')

DimPlot(cleanSrat,cells.highlight = s$cellID[s$seurat_clusters %in% c(12) & s$annot_aug24 == '-'])
table(s$timePoint[s$seurat_clusters==23])
table(s$annot_aug24[s$seurat_clusters %in% c(17)])

s$annot_aug24_new = '?'
s$annot_aug24_new[s$seurat_clusters %in% c(20,19,23,22)] = 'doublets'
s$annot_aug24_new[s$seurat_clusters == 8] = 'pDC'
s$annot_aug24_new[s$seurat_clusters == 8 & s$annot_aug24 %in% c('Mono_CD16','MEP','HSC_MPP','DC1','DC2')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(18,9,15,13)] = s$annot_aug24[s$seurat_clusters  %in% c(18,9,15,13)]
s$annot_aug24_new[s$seurat_clusters == 18 & s$annot_aug24_new %in% c('?','Mono_CD14')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(9,15,13) & s$annot_aug24 == 'HSC_MPP' ] = 'CMP_GMP'
s$annot_aug24_new[s$seurat_clusters %in% c(9,15,13) & s$annot_aug24 %in% c('DC2','Mono_CD14','pDC','?','-') ] = 'doublets'
s$annot_aug24_new[s$seurat_clusters == 21 & s$annot_aug24 !='DC2'] = 'DC1'
s$annot_aug24_new[s$seurat_clusters == 21 & s$annot_aug24 =='DC2'] = 'DC2'
s$annot_aug24_new[s$seurat_clusters == 5 & s$annot_aug24 %in% c('pDC','Mono_CD16')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters == 5 & s$annot_aug24 %in% c('DC2','Mono_CD16')] = 'DC2'
s$annot_aug24_new[s$seurat_clusters == 10 & s$annot_aug24 %in% c('pDC','DC1','Mono_CD16')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters == 10 & s$annot_aug24 %in% c('Mono_CD14','DC2','?')] = 'DC2'
s$annot_aug24_new[s$seurat_clusters == 16 & s$annot_aug24 %in% c('Mono_CD14','-','?','Mono_CD16')] = 'Mono_CD16'
s$annot_aug24_new[s$seurat_clusters %in% c(2,6,1) & s$annot_aug24 %in% c('Myelocytes','Mono_CD16','Mono_CD14','activated_neutrophil','CMP_GMP','pDC')] = s$annot_aug24[s$seurat_clusters %in% c(2,6,1) & s$annot_aug24 %in% c('Myelocytes','Mono_CD16','Mono_CD14','activated_neutrophil','CMP_GMP','pDC')]
s$annot_aug24_new[s$seurat_clusters %in% c(2,6,1) & s$annot_aug24 %in% c('HSC_MPP')] = 'CMP_GMP'
s$annot_aug24_new[s$seurat_clusters %in% c(2,6,1) & s$annot_aug24 %in% c('?','DC2')] = 'Mono_CD14'
s$annot_aug24_new[s$seurat_clusters %in% c(14,7,11) & s$annot_aug24 %in% c('-','activated_neutrophil','Neutrophil','?')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(14,7,11) & s$annot_aug24 %in% c('DC2')] = 'Mono_CD14'
s$annot_aug24_new[s$seurat_clusters %in% c(14,7,11) & s$annot_aug24 %in% c('CMP_GMP','Mono_CD16','Mono_CD14','pDC')] = s$annot_aug24[s$seurat_clusters %in% c(14,7,11) & s$annot_aug24 %in% c('CMP_GMP','Mono_CD16','Mono_CD14','pDC')]
s$annot_aug24_new[s$seurat_clusters %in% c(4,3) & s$annot_aug24 %in% c('pDC')] = 'B?'
s$annot_aug24_new[s$seurat_clusters %in% c(4,3) & s$annot_aug24 %in% c('Myelocytes')] = 'CMP_GMP'
s$annot_aug24_new[s$seurat_clusters %in% c(4,3) & s$annot_aug24 %in% c('MEP','DC2','DC1','activated_neutrophil','HSC_MPP')] = s$annot_aug24[s$seurat_clusters %in% c(4,3) & s$annot_aug24 %in% c('MEP','DC2','DC1','activated_neutrophil','HSC_MPP')]
gmp = WhichCells(s,cells = s$cellID[s$seurat_clusters %in% c(4,3) & s$annot_aug24 %in% c('CMP_GMP','Mono_CD14')],expression = (MPO > 0 | AZU1 > 0))
hsc = WhichCells(s,cells = s$cellID[s$seurat_clusters %in% c(4,3) & s$annot_aug24 %in% c('CMP_GMP')],expression = (SPINK2 > 0))
s$annot_aug24_new[s$cellID %in% gmp & s$annot_aug24_new =='?'] = 'CMP_GMP'
s$annot_aug24_new[s$cellID %in% hsc & s$annot_aug24_new =='?'] = 'HSC_MPP'
s$annot_aug24_new[s$seurat_clusters==17 & s$annot_aug24 %in% c('pDC','Neutrophil','Mono_CD14','DC2','CMP_GMP')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters==17 & s$annot_aug24 %in% c('HSC_MPP','MEP')] = 'MEP'

## Add to big.srat
cleanSrat$annot_aug24_new[cleanSrat$cellID %in% s$cellID] = s$annot_aug24_new[match(cleanSrat$cellID[cleanSrat$cellID %in% s$cellID],s$cellID)]






##------    Ery / B lineage  ------##
s = subset(cleanSrat,subset = cellID %in% cleanSrat$cellID[cleanSrat$annot_aug24 %in% c('-','?','HSC_MPP','EE','ME','LE','MEP','MEP_EE','MK','naive.B','Plasma.cell',
                                                 'pre.B.cell','pro.B.cell','unsure_EE','unsure_LE','unsure_ME','unsure_MEP','unsure_Tum_MK?') |
                                                   cleanSrat$annot_aug24_new %in% c('-','?')])

s = standard_clustering(s,clusteringRes = 0.5)
s = FindClusters(s, resolution = 1)
DimPlot(s,group.by = 'annot_aug24_new',cols=c(col25,pal34H),label = T,repel = T,label.box = T,label.size = 4) + NoLegend()
qm = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)
s$group_tmp = ifelse(s$seurat_clusters %in% c(20,18),paste0(s$annot_aug24,'_test'),s$annot_aug24)
DotPlot(s,group.by = 'group_tmp',features = keyMarkers) + 
  RotatedAxis() + 
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'right') + xlab('') + ylab('')

DimPlot(s,cells.highlight = s$cellID[s$seurat_clusters %in% c(20,18) & s$annot_aug24 == 'MEP_EE'],pt.size = 3)
DimPlot(s,cells.highlight = s$cellID[s$seurat_clusters %in% c(3) & s$LR_predicted_label_ImmuneAllLow == 'Tem/Trm cytotoxic T cells'])
table(cleanSrat$annot_aug24[cleanSrat$cellID %in% s$cellID[s$seurat_clusters == 17]],
      cleanSrat$annot_aug24_new[cleanSrat$cellID %in% s$cellID[s$seurat_clusters == 17]])
table(s$annot_aug24[s$seurat_clusters %in% c(20,18)],s$annot_aug24_new[s$seurat_clusters %in% c(20,18)])
table(s$LR_predicted_label_ImmuneAllLow[s$seurat_clusters %in% c(8)],s$annot_aug24_new[s$seurat_clusters %in% c(8)])

s$annot_aug24_new = '?'
s$annot_aug24_new[s$seurat_clusters %in% c(3,4) & s$annot_aug24 %in% c('ME','unsure_LE','LE')] = 'LE'
s$annot_aug24_new[s$seurat_clusters %in% c(17,9) & s$annot_aug24 %in% c('unsure_LE','LE')] = 'LE'
s$annot_aug24_new[s$seurat_clusters %in% c(17,9) & s$annot_aug24 %in% c('unsure_ME','ME')] = 'ME'
s$annot_aug24_new[s$seurat_clusters %in% c(9) & s$annot_aug24 %in% c('MEP_EE','MEP')] = 'MEP'
s$annot_aug24_new[s$seurat_clusters == 3 & s$annot_aug24 %in% c('pro.B.cell')] = 'pro.B.cell'
s$annot_aug24_new[s$seurat_clusters == 17 & s$annot_aug24 %in% c('T_CD8','T_CD4')] = 'T_cells'
s$annot_aug24_new[s$seurat_clusters == 17 & s$annot_aug24 %in% c('NK')] = 'NK'
s$annot_aug24_new[s$seurat_clusters == 17 & s$annot_aug24 %in% c('pro.B.cell','naive.B','MEP')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters == 17 & s$annot_aug24 %in% c('Mono_CD14')] = 'Mono_CD14'
s$annot_aug24_new[s$seurat_clusters %in% c(15,11) & s$annot_aug24 %in% c('?')] = 'pro.B.cell'
s$annot_aug24_new[s$seurat_clusters %in% c(15,11) & s$annot_aug24 %in% c('pro.B.cell')] = 'pro.B.cell'
s$annot_aug24_new[s$seurat_clusters %in% c(15,11) & s$annot_aug24 %in% c('pre.B.cell','naive.B')] = 'pre.B.cell'
s$annot_aug24_new[s$seurat_clusters %in% c(15,11) & s$annot_aug24 %in% c('HSC_MPP','CMP_GMP','Mono_CD14','NK')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(14,6,12,5,18) & s$annot_aug24 %in% c('T_CD8','Mono_CD14','-','?')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(14,6,12,5,18) & s$annot_aug24 %in% c('T_CD4')] = 'T_CD4'
b.cells = WhichCells(s,cells = s$cellID[s$seurat_clusters %in% c(14,6,12,5,18) & s$annot_aug24 %in% c('pro.B.cell','pre.B.cell')],expression = (TCL1A > 0 | MS4A1 > 0))
s$annot_aug24_new[s$seurat_clusters %in% c(14,6,12,5,18) & s$cellID %in% b.cells] = 'pre.B.cell'
s$annot_aug24_new[s$seurat_clusters %in% c(14,6,12,5,18) & !s$cellID %in% b.cells & s$annot_aug24 %in% c('pro.B.cell')] = 'pro.B.cell'
s$annot_aug24_new[s$seurat_clusters %in% c(14,6,12,5,18) & !s$cellID %in% b.cells & s$annot_aug24 %in% c('pre.B.cell')] = 'pre.B.cell'
s$annot_aug24_new[s$seurat_clusters %in% c(14,6,12,5,18) & s$annot_aug24 %in% c('NK')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(14,6,12,5,18) & s$annot_aug24 %in% c('naive.B')] = 'naive.B'
s$annot_aug24_new[s$seurat_clusters %in% c(14,6,12,5,18) & s$annot_aug24 %in% c('Plasma.cell')] = 'Plasma.cell'
s$annot_aug24_new[s$seurat_clusters %in% c(19)] = 'Plasma.cell'
s$annot_aug24_new[s$seurat_clusters %in% c(22) & s$annot_aug24 %in% c('unsure_Tum_MK?','MK')] = 'MK'
s$annot_aug24_new[s$seurat_clusters %in% c(22) & s$annot_aug24 %in% c('naive.B')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(23,21)] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(7) & s$annot_aug24 %in% c('EE')] = 'MEP'
s$annot_aug24_new[s$seurat_clusters %in% c(7) & s$annot_aug24 %in% c('LE','Mono_CD14','naive.B')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(7) & s$annot_aug24 == 'MEP'] = 'MEP'
s$annot_aug24_new[s$seurat_clusters %in% c(7) & s$annot_aug24 == 'CMP_GMP'] = 'CMP_GMP'
s$annot_aug24_new[s$seurat_clusters %in% c(7) & s$annot_aug24 == 'HSC_MPP'] = 'HSC_MPP'
gmp = WhichCells(s,cells = s$cellID[s$seurat_clusters %in% c(7) & s$annot_aug24 %in% c('HSC_MPP','?')],expression = (MPO > 0 | AZU1 > 0))
s$annot_aug24_new[s$seurat_clusters %in% c(7) & s$cellID %in% gmp] = 'CMP_GMP'
s$annot_aug24_new[s$seurat_clusters %in% c(7) & s$annot_aug24_new == '?'] = 'HSC_MPP'
s$annot_aug24_new[s$seurat_clusters %in% c(7) & s$annot_aug24 %in% c('pro.B.cell','pre.B.cell')] = 'pro.B.cell'


s = FindClusters(s, resolution = 1)
s$annot_aug24_new[s$seurat_clusters %in% c(8) & s$annot_aug24 %in% c('HSC_MPP')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(8) & s$annot_aug24 %in% c('EE','unsure_EE','ME','unsure_ME','MEP_EE')] = 'EE'
s$annot_aug24_new[s$seurat_clusters %in% c(8) & s$annot_aug24 %in% c('unsure_MEP','MEP')] = 'MEP'
s$annot_aug24_new[s$seurat_clusters %in% c(8) & s$annot_aug24 %in% c('LE')] = 'ME'
s$annot_aug24_new[s$seurat_clusters %in% c(5) & s$annot_aug24 %in% c('EE','unsure_EE','unsure_ME','unsure_MEP','-')] = 'EE'
s$annot_aug24_new[s$seurat_clusters %in% c(5) & s$annot_aug24 %in% c('HSC_MPP','pre.B.cell')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(5) & s$annot_aug24 %in% c('ME','LE')] = 'ME'
s$annot_aug24_new[s$seurat_clusters %in% c(5) & s$annot_aug24 %in% c('MEP','MEP_EE')] = 'MEP'
s$annot_aug24_new[s$seurat_clusters %in% c(20,18) & s$annot_aug24 %in% c('-','activated_neutrophil','LE','Mono_CD14','Neutrophil','NK','pro.B.cell','T_CD4')] = 'doublets'
s$annot_aug24_new[s$seurat_clusters %in% c(20,18) & s$annot_aug24 %in% c('?','EE','HSC_MPP','MEP','unsure_MEP','MEP_EE')] = 'MEP'
s$annot_aug24_new[s$seurat_clusters %in% c(20,18) & s$annot_aug24 %in% c('EE')] = 'EE'

## Add to big.srat
cleanSrat$annot_aug24_new[cleanSrat$cellID %in% s$cellID & cleanSrat$annot_aug24_new == '?'] = s$annot_aug24_new[match(cleanSrat$cellID[cleanSrat$cellID %in% s$cellID & cleanSrat$annot_aug24_new == '?'],s$cellID)]





##------    L178 myeloid  ------##
s = subset(cleanSrat,subset = cellID %in% cleanSrat$cellID[cleanSrat$donorID == 'L178'])
s = standard_clustering(s,clusteringRes = 0.5)
qm = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)
DimPlot(s,group.by = 'annot_aug24_new',cols = c(col25,pal34H),label = T,repel = T,label.box = T,label.size = 3)+ NoLegend()
table(s$LR_predicted_label_ImmuneAllLow[s$se])
View(table(s$LR_predicted_label_ImmuneAllLow[s$annot_aug24_new == '?']))

DimPlot(cleanSrat,cells.highlight = s$cellID[s$LR_predicted_label_ImmuneAllLow == 'Neutrophil-myeloid progenitor' & s$annot_aug24_new =='?' & s$seurat_clusters %in% c(1,5,13) ])
DimPlot(s,cells.highlight = s$cellID[s$annot_aug24_new=='?' & s$LR_predicted_label_ImmuneAllLow == 'Intermediate macrophages'])
FeaturePlot(s,c('CD14','FCGR3A','CLEC9A','CLEC10A'))
FeaturePlot(s,c('GATA1','KIT'))
FeaturePlot(cleanSrat,c('ACSL1'))
View(table(s$LR_predicted_label_ImmuneAllLow))


s$group_tmp = ifelse(s$seurat_clusters %in% c(1,5,13) & s$LR_predicted_label_ImmuneAllLow == 'DC2','DC2_test',s$annot_aug24_new)
DotPlot(s,group.by = 'group_tmp',features = keyMarkers) + 
  RotatedAxis() + 
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'right') + xlab('') + ylab('')


s$annot_aug24_new[s$seurat_clusters %in% c(1,5,13) & s$LR_predicted_label_ImmuneAllLow %in% c('Classical monocytes','Macrophages','Alveolar macrophages','Neutrophil-myeloid progenitor','Intermediate macrophages',
                                                                                              'Erythrophagocytic macrophages','Monocyte precursor','Intestinal macrophages','Kupffer cells')] = 'Mono_CD14'
s$annot_aug24_new[s$seurat_clusters %in% c(1,5,13) & s$LR_predicted_label_ImmuneAllLow %in% c('DC2','DC precursor','DC','Transitional DC')] = 'DC2'
s$annot_aug24_new[s$seurat_clusters %in% c(1,5,13) & s$LR_predicted_label_ImmuneAllLow %in% c('Non-classical monocytes')] = 'Mono_CD16'
s$annot_aug24_new[s$seurat_clusters %in% c(1,5,13) & s$LR_predicted_label_ImmuneAllLow %in% c('GMP')] = 'CMP_GMP'

s$annot_aug24_new[s$seurat_clusters %in% c(1,5,13) & s$LR_predicted_label_ImmuneAllLow %in% c('DC1')] = 'DC1'
s$annot_aug24_new[s$seurat_clusters %in% c(1,5,13) & s$LR_predicted_label_ImmuneAllLow %in% c('CD16- NK cells','CD16+ NK cells','Naive B cells','Tcm/Naive helper T cells','Memory B cells','Type 1 helper T cells','Double-positive thymocytes','Age-associated B cells',
                                                                                              'Myelocytes','Promyelocytes','ILC3','pDC',
                                                                                              'Plasmablasts','Proliferative germinal center B cells','Plasma cells','Tem/Effector helper T cells','Tem/Trm cytotoxic T cells','NK cells','Regulatory T cells','Type 17 helper T cells','Pro-B cells',
                                                                                              'B cells','Tcm/Naive cytotoxic T cells','Trm cytotoxic T cells','CD8a/b(entry)','Cycling NK cells','Epithelial cells','Follicular helper T cells','Transitional B cells','CRTAM+ gamma-delta T cells','Double-negative thymocytes')] = 'doublets'
s$annot_aug24_new[s$LR_predicted_label_ImmuneAllLow == 'pDC' & s$seurat_clusters %in% c(1,5,13) & s$cellID %in% cleanSrat$cellID[cleanSrat$seurat_clusters==48]] = 'pDC'
s$annot_aug24_new[s$LR_predicted_label_ImmuneAllLow == 'Classical monocytes' & s$seurat_clusters %in% c(10) & s$cellID %in% cleanSrat$cellID[cleanSrat$seurat_clusters==56]] = 'activated_neutrophil'
s$annot_aug24_new[s$LR_predicted_label_ImmuneAllLow == 'pDC' & s$seurat_clusters %in% c(10)] = 'pDC'
s$annot_aug24_new[s$LR_predicted_label_ImmuneAllLow %in% c('Macrophages','Monocytes') & s$seurat_clusters %in% c(10)] = 'Mono_CD14'
s$annot_aug24_new[s$LR_predicted_label_ImmuneAllLow %in% c('Promyelocytes','Neutrophil-myeloid progenitor','Intermediate macrophages','Alveolar macrophages') & s$annot_aug24_new == '?'] = 'Mono_CD14'
s$annot_aug24_new[s$annot_aug24_new == '?'] = 'doublets'

## Add to big.srat
cleanSrat$annot_aug24_new[cleanSrat$cellID %in% s$cellID] = s$annot_aug24_new[match(cleanSrat$cellID[cleanSrat$cellID %in% s$cellID],s$cellID)]



##---   Refined annotation with Celltypist prediction -----##
DimPlot(cleanSrat,group.by = 'annot_aug24_new',cols=c(col25,pal34H,col25),label = T,repel = T,label.box = T,label.size = 4) + NoLegend()
View(table(cleanSrat$LR_predicted_label_ImmuneAllLow[cleanSrat$annot_aug24_new == 'MEP_EE']))
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$annot_aug24_new == 'MEP_EE'])
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$annot_aug24_new == 'pre.B.cell' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Pro-B cells')])


table(cleanSrat$seurat_clusters[cleanSrat$annot_aug24_new == 'MEP'])
table(cleanSrat$annot_aug24_new)
cleanSrat$group_tmp = ifelse(cleanSrat$annot_aug24_new=='pre.B.cell' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Pro-B cells'),'pre_test',cleanSrat$annot_aug24_new)
DotPlot(cleanSrat,group.by = 'group_tmp',features = keyMarkers) + 
  RotatedAxis() + 
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'right') + xlab('') + ylab('')


##--- T/NK lineage
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(23,25,40,35) & cleanSrat$annot_aug24_new == '?' &  cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Classical monocytes','Monocytes')] = 'Mono_CD14'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(23,25) & cleanSrat$annot_aug24_new == '?' &  cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('DC2','DC precursor')] = 'DC2'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(23,25) & cleanSrat$annot_aug24_new == '?' &  cleanSrat$LR_predicted_label_ImmuneAllLow == 'Non-classical monocytes'] = 'Mono_CD16'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(23,25) & cleanSrat$annot_aug24_new == '?' &  cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Intermediate macrophages','Macrophages')] = 'Macrophage'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(10,3,37,3,4,44) & cleanSrat$annot_aug24_new == '?' &  cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Tcm/Naive helper T cells','Regulatory T cells','Tem/Trm cytotoxic T cells','Tcm/Naive helper T cells',
                                                                                                                                                                  'Tcm/Naive cytotoxic T cells','	Tcm/Naive helper T cells','Tem/Trm cytotoxic T cells','Tem/Temra cytotoxic T cells')] = 'T_cells'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(10,3,37,3,4,44) & cleanSrat$annot_aug24_new == '?' &  cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('CD16- NK cells','CD16+ NK cells','NK cells')] = 'NK'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(40,35,25,23) & cleanSrat$annot_aug24_new == '?' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('CD16- NK cells','Naive B cells')] = 'doublets'

cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(0,10,4,37,33,44,50,3) & cleanSrat$annot_aug24_new == '?' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Memory B cells','Naive B cells','Classical monocytes','Alveolar macrophages','Monocyte precursor',
                                                                                                                                                                       'DC','pDC','DC3','Epithelial cells','Late erythroid','Germinal center B cells','Late erythroid',
                                                                                                                                                                       'Age-associated B cells','ETP','Plasma cells','B cells','Plasmablasts','Monocytes','Macrophages','Myelocytes')] = 'doublets'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(0,10,4,37,33,44,50,3) & cleanSrat$annot_aug24_new == '?' & cleanSrat$LR_predicted_label_ImmuneAllLow == 'CD16+ NK cells'] = 'NK'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(0,10,4,37,33,44,50,3) & cleanSrat$annot_aug24_new == '?'] = 'T_cells'

##---- B-cells

cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(35,36,14,22,57,61) & cleanSrat$annot_aug24_new == '?'] = 'doublets'

##---- Myeloid
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(56,46,39,48) & cleanSrat$annot_aug24_new == '?'] = 'doublets'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(35,25,60,18,40,16,23,32) & cleanSrat$annot_aug24_new == '?' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Classical monocytes','Macrophages','Intermediate macrophages','Monocyte precursor','Monocytes')] = 'Mono_CD14'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(35,25,60,18,40,16,23,32) & cleanSrat$annot_aug24_new == '?' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('DC precursor','DC')] = 'DC2'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(35,25,60,18,40,16,23,32) & cleanSrat$annot_aug24_new == '?' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('DC1')] = 'DC1'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(35,25,60,18,40,16,23,32) & cleanSrat$annot_aug24_new == '?' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Non-classical monocytes')] = 'Mono_CD16'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(35,25,60,18,40,16,23,32) & cleanSrat$annot_aug24_new == '?' ] = 'doublets'

cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(11,13,15,17,19,20,21,26,42,45,1,12) & cleanSrat$annot_aug24_new == '?'] = 'Tumour'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(59) & cleanSrat$annot_aug24_new == '?'] = 'NK'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(53) & cleanSrat$annot_aug24_new == '?'] = 'T_cells'
cleanSrat$annot_aug24_new[cleanSrat$seurat_clusters %in% c(30,38,28,46,39,56,8,56) & cleanSrat$annot_aug24_new == '?'] = 'doublets'

cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'unsure_LE'] = 'LE'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'unsure_ME'] = 'ME'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'unsure_EE'] = 'EE'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'unsure_?'] = 'EE'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'B?'] = 'doublets'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'unsure_Tum_MK?'] = 'Tum_MK?'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'unsure_MEP'] = 'EE'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'MEP_EE' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('MEMP','Early erythroid','Megakaryocyte-erythroid-mast cell progenitor')] = 'MEP'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'MEP_EE' & cleanSrat$LR_predicted_label_ImmuneAllLow == 'Mid erythroid'] = 'EE'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'MEP_EE' & cleanSrat$LR_predicted_label_ImmuneAllLow == 'Classical monocytes'] = 'doublets'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'MEP_EE' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Tcm/Naive helper T cells','CD16- NK cells','Double-positive thymocytes','Memory B cells','NK cells','Naive B cells','ILC3')] = 'MEP'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'MEP_EE'] = 'EE'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == '?'] = 'doublets'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new == 'MEP' & !cleanSrat$seurat_clusters %in% c(8,35,25,40,16,23,32,36)] = 'Tumour'

dd = as.data.frame(table(cleanSrat$annot_aug24_new,cleanSrat$LR_predicted_label_ImmuneAllLow))
dd = dd[dd$Freq >0,]
View(dd[as.character(dd$Var1) != as.character(dd$Var2),])
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$annot_aug24_new == 'DC2' & cleanSrat$LR_predicted_label_ImmuneAllLow=='Classical monocytes'])
View(table(cleanSrat$LR_predicted_label_ImmuneAllLow[cleanSrat$annot_aug24 == 'Mono_CD14' & cleanSrat$annot_aug24_new=='CMP_GMP']))
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24 == 'Mono_CD14' & cleanSrat$annot_aug24_new=='CMP_GMP' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Classical monocytes','Intermediate macrophages')] = 'Mono_CD14'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new=='naive.B' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Large pre-B cells','Small pre-B cells')] = 'pre.B.cell'
cleanSrat$annot_aug24_new[cleanSrat$annot_aug24_new=='pre.B.cell' & cleanSrat$LR_predicted_label_ImmuneAllLow %in% c('Pro-B cells')] = 'pro.B.cell'

cleanSrat$annot_aug24 = as.character(cleanSrat$annot_aug24_new)

##------------------------------------------##
##  Remove cells with nCount_RNA < 1000   ####
##------------------------------------------##
## I noticed that cells from Hennings batch often have high proportion with %MT > 10%
## I wont filter these out because cells have been through more wet lab processes (FACS --> freezed --> sequencing with more treatment?)
##  I think they are more stressed because they have been handled a lot but ok enough to keep them for analysis?

# DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$nCount_RNA < 1000])
# table(cleanSrat$percent.mt > 10,is.na(cleanSrat$annot_mar24))
# cleanSrat$lowQual = rowSums(data.frame(lowCount = cleanSrat$nCount_RNA < 1000,
#                                        lowFeature= cleanSrat$nCount_RNA < 500,
#                                        highMT = cleanSrat$percent.mt > 10))
# tmp = data.frame(cellID = cleanSrat$cellID,
#                  annot_jun24 = cleanSrat$annot_jun24,
#                  donorID = cleanSrat$donorID,
#                  lowCount = cleanSrat$nCount_RNA < 1000,
#                  lowFeature= cleanSrat$nCount_RNA < 500,
#                  highMT = cleanSrat$percent.mt > 10)
# table(cleanSrat$donorID[cleanSrat$lowQual == 1],cleanSrat$annot_jun24[cleanSrat$lowQual == 1])
# table(cleanSrat$annot_jun24[cleanSrat$nCount_RNA < 1000])
# 
# 
# 
# ## Subset to remove cells with nCount < 1000
# cleanSrat = subset(cleanSrat,subset = nCount_RNA >= 1000)
# cleanSrat = standard_clustering(cleanSrat)
# cleanSrat = RunUMAP(cleanSrat, dims=1:75,seed.use = 3456)

##-----------------------------------------##
##  Save final object + metadata + UMAP  ####
##-----------------------------------------##
## Define broad lineages
cleanSrat$broadLineage = '?'
table(cleanSrat$annot_aug24,cleanSrat$broadLineage)
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('naive.B','Plasma.cell','pre.B.cell','pro.B.cell')] = 'B lineage'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('DC1','DC2','pDC')] = 'Dendritic cells'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('EE','ME','LE')] = 'Erythroblasts'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('CMP_GMP','HSC_MPP','MEP')] = 'HSC & prog.'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('MK')] = 'Megakaryocytes'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('Mono_CD14','Mono_CD16','Macrophage')] = 'Monocyte/Macrophage'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('Myelocytes','activated_neutrophil','Neutrophil')] = 'Neutrophil'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('NK','NK_T','T_CD4','T_CD8','T_gd','T_reg','T_cells')] = 'T/NK lineage'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('Tumour','Tum_MK?')] = 'Tumour'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('doublets','NA','?','unsure_EE','unsure_?','unsure_LE','unsure_ME','unsure_MEP','unsure_Tum_MK?')] = 'lowQual'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('Tumour_WT','unsure_Tumour')] = 'Tumour_unsure'



cleanSrat$finalAnn_broad = as.character(cleanSrat$annot_aug24)
mdat = cbind(cleanSrat@meta.data,cleanSrat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_mdat.csv')






##---------------------------------------------------------##
##  GATA1s mutation - refine Tumour / normal annotation  ####
##---------------------------------------------------------##
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_mdat.csv',row.names = 1)
## Remove doublets
table(is.na(mdat$annot_aug24[match(cleanSrat$cellID,mdat$cellID)]))
table(cleanSrat$annot_aug24 == mdat$annot_aug24[match(cleanSrat$cellID,mdat$cellID)])
cleanSrat$annot_aug24 = mdat$annot_aug24[match(cleanSrat$cellID,mdat$cellID)]
cleanSrat = subset(cleanSrat,subset = annot_aug24 != 'doublets')
cleanSrat = standard_clustering(cleanSrat)
colnames(mdat)[!colnames(mdat) %in% colnames(cleanSrat@meta.data)]
colnames(cleanSrat@meta.data)[!colnames(cleanSrat@meta.data) %in% colnames(mdat)]
cleanSrat@meta.data = mdat[match(cleanSrat$cellID,mdat$cellID),!colnames(mdat) %in% c('group_tmp','UMAP_1','UMAP_2'),]
mdat = cbind(cleanSrat@meta.data,cleanSrat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv')
saveRDS(cleanSrat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS')






mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2404.RDS')
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2404_mdat_tmp.csv',row.names = 1)
mlds = FindClusters(mlds, resolution = 0.4)
DimPlot(mlds,group.by = 'orig.ident',label = T,label.box = T,cols = c(col25,pal34H),label.size = 3,repel = T) + NoLegend() + NoAxes()

library(SoupX)
qm = quickMarkers(mlds@assays$RNA@counts,mlds$seurat_clusters,N = 20)
write.csv(qm,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_quickMarkers_0724.csv')

mlds$GATA1_expression = (mlds@assays$RNA@counts['GATA1',]>0)
mlds$group_tmp = ifelse(mlds$annot_jun24 == 'Tumour' & mlds$GATA1_expression == T,'Tumour_GATA1pos',
                        ifelse(mlds$annot_jun24 == 'Tumour' & mlds$GATA1_expression == F,'Tumour_GATA1neg','others'))
table(mlds$group_tmp)
qm = quickMarkers(mlds@assays$RNA@counts,mlds$group_tmp,N = 20)
write.csv(qm,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_cancer_GATA1pos.neg_quickMarkers_0724.csv')

FeaturePlot(mlds,'SLC25A21')
FeaturePlot(mlds,'percent.mt') + NoAxes()




## Define broad lineages
mlds$broadLineage = '?'
table(mlds$annot_jun24,mlds$broadLineage)
mlds$broadLineage[mlds$annot_jun24 %in% c('naive.B','Plasma.cell','pre.B.cell','pro.B.cell')] = 'B lineage'
mlds$broadLineage[mlds$annot_jun24 %in% c('DC1','DC2','pDC')] = 'Dendritic cells'
mlds$broadLineage[mlds$annot_jun24 %in% c('EE','ME','LE')] = 'Erythroblasts'
mlds$broadLineage[mlds$annot_jun24 %in% c('CMP_GMP','HSC_MPP','MEP')] = 'HSC & prog.'
mlds$broadLineage[mlds$annot_jun24 %in% c('MK')] = 'Megakaryocytes'
mlds$broadLineage[mlds$annot_jun24 %in% c('Mono_CD14','Mono_CD16')] = 'Monocyte/Macrophage'
mlds$broadLineage[mlds$annot_jun24 %in% c('Myelocytes','activated_neutrophil','Neutrophil')] = 'Neutrophil'
mlds$broadLineage[mlds$annot_jun24 %in% c('NK','NK_T','T_CD4','T_CD8','T_gd','T_reg')] = 'T/NK lineage'
mlds$broadLineage[mlds$annot_jun24 %in% c('Tumour','Tum_MK?')] = 'Tumour'
mlds$broadLineage[mlds$annot_jun24 %in% c('doublets','NA','?','unsure_EE','unsure_?','unsure_LE','unsure_ME','unsure_MEP','unsure_Tum_MK?')] = 'lowQual'
mlds$broadLineage[mlds$annot_jun24 %in% c('Tumour_WT','unsure_Tumour')] = 'Tumour_unsure'

mdat$broadLineage = '?'
##------------------------------------------------##
##  Sub-clustering cells from individual donor  ####
##------------------------------------------------##


cleanSrat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_withMTCells.RDS')




DimPlot(cleanSrat,group.by = 'seurat_clusters',cols = c(col25,pal34H,col25),label = T,repel = T,label.box = T,label.size = 3)+NoLegend()

###----------------------- Cancer - TAM and MLDS ####
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$seurat_clusters %in% c(1,3,10,9,12,13,14,46,25,30,34)])

s = subset(cleanSrat,subset = cellID %in% cleanSrat$cellID[cleanSrat$annot_jun24 %in% c('-','NK','HSC_MPP','MEP','EE','Tumour','Tumour_WT','unsure_Tumour','unsure_MEP','Tum_MK?')])
s = standard_clustering(s,clusteringRes = 0.5)

library(SoupX)
markers = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)

DimPlot(s,label = T,label.box = T) + theme(legend.position = 'none')
DimPlot(s,group.by = 'annot_jun24',label = T,label.box = T,cols = col25) + theme(legend.position = 'none')
DimPlot(s,group.by = 'seurat_clusters',label = T,label.box = T) + theme(legend.position = 'none')
DimPlot(s,group.by = 'finalAnn',label = T,label.box = T,repel = T) + theme(legend.position = 'none')
FeaturePlot(s,features = c('HBB','HBA2','ALB'))

DotPlot(s,group.by = 'finalAnn',#idents = unique(s$ann_tmp[grepl('B.cell',s$ann_tmp)]),
        features = keyMarkers)+RotatedAxis()+
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,size=8))


View(table(s$seurat_clusters,s$LR_predicted_label_ImmuneAllLow))
DimPlot(cleanSrat,cells.highlight = s$cellID[s$seurat_clusters == 39])

s$finalAnn = s$annot_jun24
s$finalAnn[s$finalAnn == '-'] = 'Tumour'
s$finalAnn[s$seurat_clusters == 39] = paste0('39-',s$annot_jun24[s$seurat_clusters == 39])
s$finalAnn[s$finalAnn == '39-Tumour'] = 'doublets'
s$finalAnn[grepl('39-',s$finalAnn)] = gsub('39-','',s$finalAnn[grepl('39-',s$finalAnn)])

cleanSrat$annot_jun24[cleanSrat$cellID %in% s$cellID] = s$finalAnn[match(cleanSrat$cellID[cleanSrat$cellID %in% s$cellID],s$cellID)]

###----------------------- HSC / GMP / MEP ####
s = subset(cleanSrat,subset = cellID %in% cleanSrat$cellID[cleanSrat$annot_jun24 %in% c('-','NK','HSC_MPP','MEP','EE','CMP_GMP','DC1','DC2','Mono_CD14','Mono_CD16','unsure_MEP','unsure_EE')])
s = standard_clustering(s,clusteringRes = 0.5)
DimPlot(s,group.by = 'seurat_clusters',label = T,label.box = T,cols = col25) + theme(legend.position = 'none')

s$finalAnn = s$annot_jun24
s$finalAnn[s$seurat_clusters %in% c(5,1) & s$finalAnn == 'DC2'] = paste0('DC2-',s$seurat_clusters[s$seurat_clusters %in% c(5,1) & s$finalAnn == 'DC2'])

###----------------------- T/NK ####
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$seurat_clusters %in% c(1,3,10,9,12,13,14,46,25,30,34)])

s = subset(cleanSrat,subset = cellID %in% cleanSrat$cellID[cleanSrat$seurat_clusters %in% c(28,27,0,7,35,8,9) | cleanSrat$annot_jun24 %in% c('-','T_CD4','T_CD8','T_gd','T_reg','NK','NK_T')])
s = standard_clustering(s,clusteringRes = 0.5)

library(SoupX)
markers = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)

DimPlot(s,label = T,label.box = T) + theme(legend.position = 'none')
DimPlot(s,group.by = 'may23_ann',label = T,label.box = T) + theme(legend.position = 'none')
DimPlot(s,group.by = 'seurat_clusters',label = T,label.box = T) + theme(legend.position = 'none')
DimPlot(s,group.by = 'finalAnn',label = T,label.box = T,repel = T) + theme(legend.position = 'none')
FeaturePlot(s,features = c('HBB','HBA2','ALB'))

DotPlot(s,group.by = 'finalAnn',#idents = unique(s$ann_tmp[grepl('B.cell',s$ann_tmp)]),
        features = keyMarkers)+RotatedAxis() +
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,size=8))

View(table(s$seurat_clusters,s$LR_predicted_label_ImmuneAllLow))
DimPlot(s,cells.highlight = tgd_cells)

s$finalAnn[s$seurat_clusters %in% c(4,5,13,20)] = 'NK'

