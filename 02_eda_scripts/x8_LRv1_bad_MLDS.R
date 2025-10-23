## Logistic regression between all bad ML-DS clusters

outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/x8_LRv1_bad_MLDS'
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


##--------------------------------------------------------##
####   1. LR_orig with REF = canonical TAM vs all ML-DS ####
##--------------------------------------------------------##
mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS')
mlds$annot = as.character(mlds$annot_aug24)

## Subset to just bad MLDS
REF.srat = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID %in% c('L076','L038') & mlds$annot == 'Tumour' & mlds$tissue == 'BM'])
REF.srat$finalAnn = REF.srat$annot

## Import fAdr
fAdr = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/adrenal/oct22/adrenal_clean_filtered_annotated.RDS')
fAdr@meta.data$cellID = rownames(fAdr@meta.data)
cells_toKeep = fAdr@meta.data$cellID[fAdr$Phase == 'G1' & fAdr@meta.data$finalAnn %in% c('SCPs')]
fAdr = subset(fAdr,subset = cellID %in% cells_toKeep)

## Combine to make REF.srat
REF.srat = merge_seurat_objects(REF.srat,fAdr,keepAllGenes = F,genomeVersions = c('v38','v38'))
message('1. REF.srat loaded')

df = REF.srat@meta.data[REF.srat$finalAnn == 'Tumour',]
df$group = paste0(df$timePoint,'_',df$tissue,'_',df$donorID)
df = df %>% group_by(group) %>% mutate(idx = 1:n(),
                                       toKeep = (idx %in% sample(1:n(),n()/5)))
table(df$group,df$toKeep)
tgt.srat = subset(REF.srat,subset=cellID %in% df$cellID[df$toKeep])
REF.srat = subset(REF.srat,subset=cellID %in% REF.srat$cellID[!REF.srat$cellID %in% df$cellID[df$toKeep]])
REF.srat$type = 'ref'
tgt.srat$type = 'tgt'

REF.srat$finalAnn = ifelse(REF.srat$finalAnn == 'SCPs','SCPs',paste0(REF.srat$donorID,'_',REF.srat$timePoint,'_',REF.srat$tissue))
tgt.srat$finalAnn = paste0(tgt.srat$donorID,'_',tgt.srat$timePoint,'_',tgt.srat$tissue)




##----------------------------##
##   Run Logistic Regression  ##
##----------------------------##
skipIfExists=F
ref_annot='finalAnn'
maxCells=4000
#outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_MLDSref/')
model_fp = file.path(outDir,'badMLDS_BM_REF_trainModel_4kmaxcells_70perc.RDS')  # REF = gTAM + ML-DS + SCPs

LR_level='both'
srat_annot='finalAnn'
minGeneMatch = 0.99
maxCells=4000
tissue = 'MLDS'


out_prefix = 'badMLDSref_maxCells_70perc_'
plot_prefix = 'badMLDSref_maxCells_70perc_'


plot_prefix = NULL

outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
                model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
                minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
                scLR_TGTtype='',scLR_REFtype='')

message(sprintf('3. LR completed for tissue %s',tissue))


##---------------------------------##
##   Do some LR_similarity plots   ##
##---------------------------------##
ref_order = c('SCPs','L076_Diagnostic_Blood','L076_Diagnostic_BM','L076_D.Relapse_BM','L076_D.Relapse2_BM','L038_Diagnostic_BM','L038_TP1_BM')





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
              ref_order[ref_order!='SCPs'])
row_order[!row_order %in% rownames(output)]
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
