## Run Logistic Regression v1, with REF = fLiver 2n, tgt = ML-DS ##

outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/goodTAM_vs_badTAM/LRv1_tam.mlds.ref'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}


setwd(outDir)



##------------------##
##    Libraries   ####
##------------------##
library(tidyverse)
library(readxl)
library(Seurat)
library(RColorBrewer)

source("~/lustre_mt22/generalScripts/utils/logisticRegression.R")
source("~/lustre_mt22/generalScripts/utils/runLR.R")
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")


message('\n--------  Hello! --------\n')


##----------##
##  Params  ##
##----------##


seed = 2397
keepMTCells = T
numPCs = 75

library(GenomicFeatures)
#Define genomic coordinates
# gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf.gz'
# txdb = makeTxDbFromGFF(gtf)
# gns = genes(txdb)
# 
# geneMap = read.delim('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_45842_SB_Leuk13104278_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz',header = F)
# colnames(geneMap) = c('ensID','geneSym','GEX')
# geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
# geneMap$geneSym = gsub('_','-',geneMap$geneSym)
# geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))


##--------------------------------------------------------##
####   1. LR_orig with REF = canonical TAM vs all ML-DS ####
##--------------------------------------------------------##
ref_dataset = c('publishedFL','2n_FL','T21_FL','MLDS')
ref_dataset = 'MLDS'

ref.srat.fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS'


if(file.exists(ref.srat.fp)){
  REF.srat = readRDS(ref.srat.fp)
  REF.srat$timePoint[REF.srat$orig.ident %in% c('MY.200531.14635833')] = 'D.Relapse'
  
  mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
  mdat$broadLineage[mdat$broadLineage %in% c('Tumour_unsure','lowQual')] = 'others'
  mdat$broadLineage[mdat$broadLineage =='Tumour' & mdat$donorID == 'L041'] = 'others'
  #mdat$broadLineage[mdat$broadLineage == 'Tumour'] = paste0('Leuk:',mdat$disease[mdat$broadLineage == 'Tumour'],':',mdat$donorID[mdat$broadLineage == 'Tumour'])
  
  mdat$group = ifelse(mdat$broadLineage == 'Tumour' & mdat$disease == 'TAM' & mdat$donorID %in% c('L075','CC1','CC2'),'good_TAM',
                      ifelse(mdat$broadLineage == 'Tumour' & mdat$disease == 'TAM' & !mdat$donorID %in% c('L075','CC1','CC2'),paste0(mdat$donorID,'_TAM'),
                             ifelse(mdat$broadLineage == 'Tumour' & mdat$disease == 'MLDS' & !mdat$donorID %in% c('L041','CC3') & mdat$timePoint == 'Diagnostic','dMLDS','others')))
  
  
  mdat$group = ifelse(mdat$broadLineage == 'Tumour' & mdat$disease == 'TAM' & mdat$donorID %in% c('L075','CC1','CC2'),'good_TAM',
                      ifelse(mdat$broadLineage == 'Tumour' & mdat$disease == 'TAM' & !mdat$donorID %in% c('L075','CC1','CC2'),paste0(mdat$donorID,'_TAM'),
                             ifelse(mdat$broadLineage == 'Tumour' & mdat$disease == 'MLDS' & !mdat$donorID %in% c('L041','CC3') & mdat$timePoint == 'Diagnostic',paste0(mdat$donorID,'_MLDS'),'others')))
  
  REF.srat$group = mdat$group[match(REF.srat$cellID,mdat$cellID)]
  
  
  ## tgt.srat = 'L156'
  tgt.srat = subset(REF.srat,subset = cellID %in% c(mdat$cellID[mdat$group %in% c('L156_TAM','CC6_TAM','CC7_TAM','CC8_TAM') & mdat$annot_mar24 == 'Tumour']))

  REF.srat = subset(REF.srat,subset = group %in% c('good_TAM',unique(mdat$group[grepl('_MLDS',mdat$group)])))
  REF.srat$finalAnn[REF.srat$cellID %in% mdat$cellID] = REF.srat$group[REF.srat$cellID %in% mdat$cellID]
}




## Import fAdr
fAdr = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/adrenal/oct22/adrenal_clean_filtered_annotated.RDS')
fAdr@meta.data$cellID = rownames(fAdr@meta.data)
cells_toKeep = fAdr@meta.data$cellID[fAdr$Phase == 'G1' & fAdr@meta.data$finalAnn %in% c('SCPs')]
fAdr = subset(fAdr,subset = cellID %in% cells_toKeep)


REF.srat = merge_seurat_objects(REF.srat,fAdr,keepAllGenes = F,genomeVersions = c('v38','v38'))
message('1. REF.srat loaded')







message(('2. tgt.srat loaded'))



## Configure REF and tgt sratObj for LR
# Only keep the same genes between REF and tgt.srat
REF.srat$type = 'ref'
tgt.srat$type = 'tgt'

genesToKeep = intersect(rownames(REF.srat),rownames(tgt.srat))
# remove rubbish genes
genesToKeep = genesToKeep[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS|GATA1$',genesToKeep)]
REF.srat@assays$RNA@counts = REF.srat@assays$RNA@counts[genesToKeep,]
tgt.srat@assays$RNA@counts = tgt.srat@assays$RNA@counts[genesToKeep,]


##----------------------------##
##   Run Logistic Regression  ##
##----------------------------##
skipIfExists=F
ref_annot='finalAnn'
maxCells=4000
#outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_MLDSref/')
model_fp = file.path(outDir,'gTAM.MLDS_REF_trainModel_4kmaxcells_70perc.RDS')  # REF = gTAM + ML-DS + SCPs
model_fp = file.path(outDir,'gTAM.dMLDS_REF_trainModel_4kmaxcells_70perc.RDS')  # REF = gTAM + diagnostic_ML-DS + SCPs  
model_fp = file.path(outDir,'gTAM.dMLDS_filteredGenes_REF_trainModel_4kmaxcells_70perc.RDS')   # REF = gTAM + diagnostic_ML-DS + SCPs (removed rubbished genes before training the model)
model_fp = file.path(outDir,'gTAM.indivMLDS_REF_trainModel_4kmaxcells_70perc.RDS')  # REF = gTAM + individual donor diagnostic_ML-DS + SCPs  

LR_level='both'
srat_annot='group'
minGeneMatch = 0.99
maxCells=4000
tissue = 'MLDS'


out_prefix = 'L156_gTAM.MLDSref_maxCells_70perc_'
plot_prefix = 'L156_gTAM.MLDSref_maxCells_70perc_'

out_prefix = 'otherTAM_gTAM.indivMLDSref_maxCells_70perc_'
plot_prefix = 'otherTAM_gTAM.indivMLDSref_maxCells_70perc_'

out_prefix = 'otherTAM_gTAM.dMLDSref_maxCells_70perc_'
plot_prefix = 'otherTAM_gTAM.dMLDSref_maxCells_70perc_'

out_prefix = 'otherTAM_gTAM.dMLDS.filteredGenes.ref_maxCells_70perc_'
plot_prefix = 'otherTAM_gTAM.dMLDS.filteredGenes.ref_maxCells_70perc_'

plot_prefix = NULL

outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
                model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
                minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
                scLR_TGTtype='',scLR_REFtype='')

message(sprintf('3. LR completed for tissue %s',tissue))


##---------------------------------##
##   Do some LR_similarity plots   ##
##---------------------------------##

ref_order = c('SCPs','good_TAM','MLDS')
ref_order = c('SCPs','good_TAM','dMLDS')
ref_order = c('SCPs','good_TAM',colnames(output)[grepl('_MLDS',colnames(output))])





outputs = readRDS(file.path(outDir,'otherTAM_gTAM.dMLDSref_maxCells_70perc_raw_LR_outputs.RDS'))
## annotated Cluster level #
if(length(outputs) == 2){
  output = outputs[[2]][[1]]  
}else{
  output = outputs[[1]]
}

type = ifelse(grepl('ref_',rownames(output)),'REF','TGT')
show_row_names = T

table(colnames(output) %in% ref_order)
colnames(output)[!colnames(output) %in% ref_order]
ref_order = ref_order[ref_order %in% colnames(output)]


row_order = c(rownames(output)[grepl('ref_',rownames(output))],
              rownames(output)[grepl('Tumour',rownames(output))],
              rownames(output)[!grepl('Tumour|ref_',rownames(output))])

row_order = unique(c(rownames(output)[grepl('ref_',rownames(output))],
              rownames(output)[grepl('TAM',rownames(output))],
              rownames(output)[!grepl('TAM|ref_',rownames(output))]))

table(rownames(output) %in% row_order)




pdf(file.path(outDir,paste0(plot_prefix,'clusterLR.pdf')),width = 10,height = 55)

hm = similarityHeatmap(output,
                       row_order=row_order,
                       column_order = ref_order,
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize=10),row_names_gp = gpar(fontsize=10),row_names_max_width = unit(6,'cm'),
                       column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                       split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = F)
draw(hm)
dev.off()





## annotated single-cell level #
if(length(outputs) == 2){
  output = outputs[[2]][[2]]  
}else{
  output = outputs[[2]]
}
output = output[!rownames(output) %in% tgt.srat$cellID[tgt.srat$annot_mar24 == 'Tumour' & tgt.srat$donorID %in% c('CC6','CC7','CC8') ],]
type = ifelse(grepl('ref_',rownames(output)),REF.srat$donorID[match(gsub('ref_','',rownames(output)),REF.srat$cellID)],
       tgt.srat$donorID[match(rownames(output),tgt.srat$cellID)])
type[grepl('_adrenal',type)] = 'SCPs'

hm = similarityHeatmap(output,
                       #row_order=row_order,
                       column_order = ref_order,
                       row_title_rot = 0,
                       row_title_gp = gpar(fontsize=8),row_names_gp = gpar(fontsize=7),row_names_max_width = unit(6,'cm'),
                       column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
                       split = type, gap = unit(2,'mm'), show_row_names = F, cluster_rows = T)
ht = draw(hm)
cell_order = data.frame(cellID = rownames(output),group='?')
cell_order$group[cell_order$cellID %in% rownames(output)[row_order(ht)[['1,Tumour']]]] = 'pale_both'
cell_order$group[cell_order$cellID %in% rownames(output)[row_order(ht)[['2,Tumour']]]] = 'pale_gTAM'
cell_order$group[cell_order$cellID %in% rownames(output)[row_order(ht)[['3,Tumour']]]] = 'mid_gTAM'
cell_order$group[cell_order$cellID %in% rownames(output)[row_order(ht)[['4,Tumour']]]] = 'strong_gTAM'

DimPlot(tgt.srat,group.by = 'donorID')
DimPlot(tgt.srat,cells.highlight = cell_order$cellID[cell_order$group == 'pale_both'])

tgt.srat@assays$RNA@data = tgt.srat@assays$RNA@counts
tgt.srat = standard_clustering(tgt.srat)
FeaturePlot(tgt.srat,'CD74')


df = output[rownames(output) %in% tgt.srat$cellID[tgt.srat$donorID == 'L156'],]
df$nCount = tgt.srat$nCount_RNA[match(rownames(df),tgt.srat$cellID)]
df$nGene = tgt.srat$nFeature_RNA[match(rownames(df),tgt.srat$cellID)]
df$BST2 = tgt.srat@assays$RNA@data['BST2',match(rownames(df),colnames(tgt.srat@assays$RNA@data))]
ggplot(df,aes(BST2,dMLDS))+
  geom_point(size=0.1)

##-----------------------------------------------------------------##
##     Plot expression of model markers in scRNAseq dataset      ####
##-----------------------------------------------------------------##

## Import MLDS object
mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
mlds$broadLineage = mdat$broadLineage[match(mlds$cellID,mdat$cellID)]
mlds$timePoint[mlds$orig.ident %in% c('MY.200531.14635833')] = 'D.Relapse'
mlds$finalAnn_broad = as.character(mlds$annot_mar24)

## Import model
model_fp = file.path(outDir,'gTAM.MLDS_REF_trainModel_4kmaxcells_70perc.RDS')  
model_fp = file.path(outDir,'gTAM.dMLDS_REF_trainModel_4kmaxcells_70perc.RDS')  
model_fp = file.path(outDir,'gTAM.indivMLDS_REF_trainModel_4kmaxcells_70perc.RDS')  

model = readRDS(model_fp)
model_markers = getMarkers(model)


## DotPlot
mlds$group = ifelse(mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID != 'L041',paste0(mlds$disease,':',mlds$donorID),mlds$annot_mar24)
mlds$group[mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076'] = paste0('MLDS:L076:',mlds$tissue[mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076'])
mlds$group[mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'D.Relapse' & mlds$donorID == 'L076'] = 'MLDS:L076:DR'
mlds$group[mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'TP1' & mlds$donorID == 'L038'] = 'MLDS:L038:TP1'
mlds$group = gsub('MLDS:|TAM:','',mlds$group)
mlds$group = factor(mlds$group,c(unique(mlds$group[grepl('CC',mlds$group) & !grepl('L076|L038',mlds$group) & mlds$disease == 'MLDS']),
                                 unique(mlds$group[grepl('L\\d+',mlds$group) & !grepl('L076|L038',mlds$group) & mlds$disease == 'MLDS']),
                                 unique(mlds$group[grepl('L038',mlds$group)]),
                                 unique(mlds$group[grepl('L076:DR',mlds$group)]),
                                 unique(mlds$group[grepl('L076:BM',mlds$group)]),
                                 unique(mlds$group[grepl('L076:Blood',mlds$group)]),
                                 unique(mlds$group[grepl('L075',mlds$group)]),
                                 unique(mlds$group[grepl('CC2',mlds$group)]),
                                 unique(mlds$group[grepl('CC1',mlds$group)]),
                                 unique(mlds$group[grepl('CC6',mlds$group)]),
                                 unique(mlds$group[grepl('CC7',mlds$group)]),
                                 unique(mlds$group[grepl('CC8',mlds$group)]),
                                 unique(mlds$group[grepl('L156',mlds$group)]),
                                 unique(mlds$group[mlds$annot_mar24 != 'Tumour']),
                                 'Tumour'
                                 
))

Idents(mlds) = mlds$group
model_markers = model_markers[order(abs(model_markers$coef),decreasing = T),]

DotPlot(mlds,
        idents = unique(Idents(mlds)[grepl('CC\\d+|L\\d+',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT|Tumour',Idents(mlds))]),
        scale = T,
        features = c(model_markers$gene[model_markers$coef > 0 & model_markers$class == 'good_TAM'][1:20],
                     model_markers$gene[model_markers$coef > 0 & model_markers$class == 'dMLDS'][1:20])
)+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('')

