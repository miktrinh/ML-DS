## Run Logistic Regression v1, with REF = fLiver 2n, tgt = ML-DS ##

outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/LRv1_published_FLref'
outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/2_LRv1_FLref'
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
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf.gz'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

geneMap = read.delim('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_45842_SB_Leuk13104278_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz',header = F)
colnames(geneMap) = c('ensID','geneSym','GEX')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))


##-------------------------------------------------------##
####   1. LR_orig with REF = published FL ; tgt = MLDS ####
##-------------------------------------------------------##
ref_dataset = c('publishedFL','2n_FL','T21_FL')
ref_dataset = '2n_FL'

if(ref_dataset == 'publishedFL'){
  # Import REF sratObj = scRNA MLDS
  REF.srat = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver_2n/jan24/liver_liverREFmerged_clean_processed_annotated_noUnknowns_0124.RDS')
  ## Only keep published cells
  REF.srat = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$published_ann_2 != 'NA'])
  published_mdat = read.csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/Muz_fLiver_SoupXed_sratObj_metadata.csv')
  REF.srat$cellID = gsub('-.*$','',REF.srat$cellID)
  message(sprintf('Of the %d published cells, %d cells are retained',nrow(published_mdat),sum(REF.srat$cellID %in% published_mdat$cellID)))
  
  REF.srat$finalAnn = published_mdat$cell.labels[match(REF.srat$cellID,published_mdat$cellID)]
  table(REF.srat$finalAnn)
}else if(ref_dataset == '2n_FL'){
  akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'
  akLiv_mdat_fp = gsub('_0824.RDS','_0824_mdat.csv',akLiv_srat_fp)
  
  REF.srat = readRDS(akLiv_srat_fp)
  REF.srat$finalAnn = as.character(REF.srat$annot_aug24)
  REF.srat = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$finalAnn != 'doublets' & REF.srat$Genotype == 'diploid'])
  table(REF.srat$Genotype,REF.srat$finalAnn)
  #REF.srat = subset(REF.srat,subset = Genotype == 'diploid')
}else if(ref_dataset == 'T21_FL'){
  
}



## Import fAdr
fAdr = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/adrenal/oct22/adrenal_clean_filtered_annotated.RDS')
fAdr@meta.data$cellID = rownames(fAdr@meta.data)
cells_toKeep = fAdr@meta.data$cellID[fAdr$Phase == 'G1' & fAdr@meta.data$finalAnn %in% c('SCPs')]
fAdr = subset(fAdr,subset = cellID %in% cells_toKeep)


REF.srat = merge_seurat_objects(REF.srat,fAdr,keepAllGenes = F,genomeVersions = c('v38','v38'))
message('1. REF.srat loaded')




## Import tgt sratObj = MLDS
tgt.srat.fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS'
mlds_mdat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv'

if(file.exists(tgt.srat.fp)){
  tgt.srat = readRDS(tgt.srat.fp)
  
  mdat = read.csv(mlds_mdat_fp)
  mdat$broadLineage[mdat$broadLineage != 'Tumour' & mdat$annot_aug24 == 'Tumour']='Tumour'
  mdat$broadLineage[mdat$broadLineage %in% c('Tumour_unsure','lowQual')] = 'others'
  mdat$broadLineage[mdat$broadLineage =='Tumour' & mdat$donorID == 'L041'] = 'others'
  #mdat$broadLineage[mdat$broadLineage == 'Tumour'] = paste0('Leuk:',mdat$disease[mdat$broadLineage == 'Tumour'],':',mdat$donorID[mdat$broadLineage == 'Tumour'])
  
  mdat$group = ifelse(mdat$broadLineage == 'Tumour',paste0(mdat$disease,':',mdat$donorID,':',mdat$timePoint,':',mdat$tissue),
                      ifelse(mdat$broadLineage == 'others','others',mdat$annot_aug24))
  mdat$group[mdat$group %in% c('T_gd','T_reg','T_CD4','T_CD8','T_cells')] = 'T cell'
  mdat$group[mdat$group %in% c('activated_neutrophil','Neutrophil')] = 'Neutrophil'
  mdat$group[mdat$group %in% c('Mono_CD14','Mono_CD16')] = 'Monocyte'
  mdat$group[mdat$group %in% c('naive.B')] = 'B cell'
  mdat$group[grepl('\\.',mdat$group) & !grepl('Leuk',mdat$group)] = gsub('\\.',' ',mdat$group[grepl('\\.',mdat$group) & !grepl('Leuk',mdat$group)])
  
  tgt.srat$group = mdat$group[match(tgt.srat$cellID,mdat$cellID)]
}

message(('2. tgt.srat loaded'))



## Configure REF and tgt sratObj for LR
# Only keep the same genes between REF and tgt.srat
REF.srat$type = 'fLiver'
tgt.srat$type = 'MLDS'

genesToKeep = intersect(rownames(REF.srat),rownames(tgt.srat))
REF.srat@assays$RNA@counts = REF.srat@assays$RNA@counts[genesToKeep,]
tgt.srat@assays$RNA@counts = tgt.srat@assays$RNA@counts[genesToKeep,]
#merged.srat = merge_seurat_objects(srat1 = REF.srat,srat2 = tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))

#REF.srat = subset(merged.srat,subset = type == 'MLDS')
#tgt.srat = subset(merged.srat,subset = type == 'fliver')

##----------------------------##
##   Run Logistic Regression  ##
##----------------------------##
skipIfExists=T
ref_annot='finalAnn'
maxCells=4000
#outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_MLDSref/')
#model_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/AdrKidLiv_REF_trainModel_4kmaxcells_70perc.RDS'
if(ref_dataset == 'publishedFL'){
  model_fp = file.path(outDir,'publishedFL_REF_trainModel_4kmaxcells_70perc.RDS')  
}else if(ref_dataset == '2n_FL'){
  model_fp = file.path(outDir,'diploiFL_REF_trainModel_4kmaxcells_70perc.RDS')  
}else{
  model_fp = file.path(outDir,'diploiFL_REF_trainModel_4kmaxcells_70perc.RDS')
}



LR_level='both'
srat_annot='group'
minGeneMatch = 0.99
maxCells=4000
tissue = 'MLDS'

if(ref_dataset == 'publishedFL'){
  out_prefix = 'MLDS_publishedFLref_maxCells_70perc_'
  plot_prefix = 'MLDS_publishedFLref_maxCells_70perc_'
}else if(ref_dataset == '2n_FL'){
  out_prefix = 'MLDS_2nFLref_maxCells_70perc_'
  plot_prefix = 'MLDS_2nFLref_maxCells_70perc_'
  
}else{
  out_prefix = 'MLDS_T21FLref_maxCells_70perc_'
  plot_prefix = 'MLDS_T21FLref_maxCells_70perc_'
}


outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
                model_fp = model_fp,outDir=outDir,plot_prefix=NULL,out_prefix=out_prefix,
                minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
                scLR_TGTtype='',scLR_REFtype='')

message(sprintf('3. LR completed for tissue %s',tissue))


##---------------------------------##
##   Do some LR_similarity plots   ##
##---------------------------------##

if(ref_dataset == 'publishedFL'){
  ref_order = c("HSC_MPP",
                "MEMP","Megakaryocyte","Mast cell",
                "Early Erythroid", "Mid Erythroid","Late Erythroid",
                "Neutrophil-myeloid progenitor",
                "Monocyte precursor","Monocyte","Mono-Mac","Kupffer Cell",'VCAM1..EI.macrophage',
                "DC precursor","DC1","DC2","pDC precursor",
                "Pre pro B cell","pro-B cell","pre-B cell","B cell",
                "ILC precursor", "Early lymphoid_T lymphocyte", "NK",
                "Hepatocyte","Endothelial cell","Fibroblast",'SCPs' 
  )
}else if(ref_dataset == '2n_FL'){
  ref_order = c("HSC_MPP","MEMP_MEP",'CMP_GMP','LMPP_ELP',
                'earlyMK',"MK","Mast.cell",
                "EE", "ME","LE",
                'promyelocyte','myelocyte','proMono','Monocyte','Macrophage','Kupffer.cell',
                "DC precursor","DC1","DC2","pDC",
                "pro.B.cell","pre.B.cell","B.cell",
                "ILC.precursor", "T.cell", "NK_T",
                "Hepatocyte","Endo",'Mesenchyme','Fibroblast',
                'Cholangiocytes','Mesothelial_cells','Neuron','SCPs' 
  )
}






outputs = readRDS(file.path(outDir,paste0(out_prefix,'raw_LR_outputs.RDS')))
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
              rownames(output)[grepl('TAM',rownames(output))],
              rownames(output)[grepl('MLDS',rownames(output))],
              rownames(output)[!grepl('ref_|TAM|MLDS',rownames(output))])

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












# ##-------------------------------------------------------##
# ####   2. LR_orig with REF = 2n fLiver ; tgt = MLDS   ####
# ##------------------------------------------------------##
# 
# ## Import REF sratObj = scRNA fLiver (all genotypes)
# REF.srat.fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/sept23/liver_liverREFmerged_clean_processed_annotated_noUnknowns_0923.RDS'
# #REF.srat.fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/may23/liver_liverREFmerged_clean_processed_annotated_v2.RDS'
# 
# if(file.exists(REF.srat.fp)){
#   REF.srat = readRDS(REF.srat.fp)
#   REF.srat$published_ann_2[is.na(REF.srat$published_ann_2)] == '?'
#   REF.srat$ann = paste0(REF.srat$Genotype,'_',REF.srat$finalAnn_broad)
#   # Remove unpublished diploid cells
#   REF.srat$cells_toKeep = ifelse(REF.srat$Genotype == 'diploid' & REF.srat$published_ann_2 == '?',F,T)
#   print(table(REF.srat$cells_toKeep,REF.srat$Genotype))
#   #REF.srat = subset(REF.srat,subset = cells_toKeep == TRUE)
#   print(table(REF.srat$published_ann_2 == '?',REF.srat$Genotype))
#   
#   REF.srat = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$Genotype %in% c('T21')])
#   REF.srat = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$annot_sept23 %in% c('HSC_MPP','MEMP_MEP','MK','EE','Mast.cell','NK_T','Monocyte','B.cell')])
#   print(table(REF.srat$finalAnn_broad,REF.srat$annot_sept23))
#   print(table(REF.srat$Genotype))
# }
# 
# 
# message(('1. REF.srat loaded'))
# 
# ## Import tgt sratObj = scRNA MLDS
# #tgt.srat = readRDS(('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov22/MLDSonly_clean_LRwCTannotated_v3.RDS'))
# #tgt.srat = readRDS(('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov22/MLDS_clean_LRwCTannotated_jan23.RDS'))
# #tgt.srat = readRDS(('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/may23/MLDS_clean_LRwCTannotated_may23.RDS'))
# tgt.srat = readRDS('~/lustre_mt22/Aneuploidy/Results/6c_cancerModuleScore_along_2n_T21_TAM_MLDS/combined_srat.RDS')
# tgt.srat = subset(tgt.srat,subset = cellID %in% tgt.srat$cellID[!tgt.srat$cellID %in% REF.srat$cellID])
# message('2. tgt.srat loaded')
# 
# 
# ## Configure REF and tgt sratObj for LR
# # Only keep the same genes between REF and tgt.srat
# REF.srat$type = 'fLiver'
# tgt.srat$type = 'MLDS'
# 
# genesToKeep = intersect(rownames(REF.srat),rownames(tgt.srat))
# # Remove chr21 genes
# genesToKeep = genesToKeep[!genesToKeep %in% geneMap$geneSym[geneMap$chr == 'chr21']]
# REF.srat_full = REF.srat
# tgt.srat_full = tgt.srat
# REF.srat@assays$RNA@counts = REF.srat@assays$RNA@counts[genesToKeep,]
# tgt.srat@assays$RNA@counts = tgt.srat@assays$RNA@counts[genesToKeep,]
# 
# print(dim(REF.srat@assays$RNA@counts))
# print(dim(tgt.srat@assays$RNA@counts))
# 
# # merged.srat = merge_seurat_objects(srat1 = REF.srat,srat2 = tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # 
# # tgt.srat = subset(merged.srat,subset = type == 'MLDS')
# # REF.srat = subset(merged.srat,subset = type == 'fLiver')
# 
# tgt.srat$ann[grepl('Tumour:',tgt.srat$ann)] = paste0(tgt.srat$donorID[grepl('Tumour:',tgt.srat$ann)],'_',tgt.srat$ann[grepl('Tumour:',tgt.srat$ann)])
# 
# # If we want to train LR model on specific cell types only - subset REF.srat to include these CTs and remove others
# #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Hepatocyte','Early.Erythroid','Megakaryocyte','Fibroblast','CMP','GMP','Mono.Mac'))
# #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Early.Erythroid','Megakaryocyte'))
# #REF.srat = subset(REF.srat,subset = Genotype %in% c('T21','diploid'))
# 
# ##----------------------------##
# ##   Run Logistic Regression  ##
# ##----------------------------##
# skipIfExists=F
# ref_annot='ann'
# maxCells=4000
# outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/')
# #model_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/AdrKidLiv_REF_trainModel_4kmaxcells_70perc.RDS'
# #model_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/fLiver_MKlineage_dipT21_REF_trainModel_4kmaxcells_95perc.RDS'
# #model_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/fLiver_REF_trainModel_4kmaxcells_70perc_jan23.RDS'
# model_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/fLiverT21_noChr21_REF_trainModel_4kmaxcells_70perc_oct23_MKlineage.RDS'
# 
# LR_level='both'
# srat_annot='ann'
# minGeneMatch = 0.99
# maxCells=4000
# tissue = 'liver'
# 
# out_prefix = 'MLDS1023_fLiverT21.noChr21REF1023_maxCells_70perc_MKlineage_'
# plot_prefix = 'MLDS1023_fLiverT21.noChr21REF1023_maxCells_70perc_MKlineagge_'
# plot_prefix = NULL
# 
# write.csv(tgt.srat@meta.data,file.path(outDir,paste0(out_prefix,'TGTsrat_mdat.csv')))
# write.csv(REF.srat@meta.data,file.path(outDir,paste0(out_prefix,'REFsrat_mdat.csv')))
# 
# 
# outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
#                 model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
#                 minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
#                 scLR_TGTtype='',scLR_REFtype='')
# 
# message(sprintf('3. LR completed for tissue %s',tissue))
# 
# 
# #---------------------------------##
# #   Do some LR_similarity plots   ##
# #---------------------------------##
# model_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/fLiverT21_noChr21_REF_trainModel_4kmaxcells_70perc_oct23_MKlineage.RDS'
# outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/')
# out_prefix = 'MLDS1023_fLiverT21.noChr21REF1023_maxCells_70perc_MKlineage_'
# tgt.srat = read.csv(file.path(outDir,paste0(out_prefix,'TGTsrat_mdat.csv')))
# REF.srat = read.csv(file.path(outDir,paste0(out_prefix,'TGTsrat_mdat.csv')))
# 
# 
# outputs = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/MLDS1023_fLiverT21.noChr21REF1023_maxCells_70perc_MKlineage_raw_LR_outputs.RDS')
# # ## annotated Cluster level #
# # if(length(outputs) > 2){
# #   output = outputs[[1]]
# # }else{
# #   output = outputs[[2]][[1]]
# # }
# # 
# # type = ifelse(grepl('ref_',rownames(output)),'REF','TGT')
# # show_row_names = T
# # column_order = colnames(output)[order(colnames(output))]
# # column_order = column_order[order(gsub('diploid_|T21_','',column_order))]
# # row_order = rownames(output)[!grepl('Tumour',rownames(output))]
# # row_order = row_order[order(row_order)]
# # row_order = c(row_order,rownames(output)[grepl('Tumour',rownames(output))])
# # 
# # 
# # plot_prefix = 'MLDS1023_fLiver2nREF1023_maxCells_70perc_MKlineagge_'
# # 
# # pdf(file.path(outDir,paste0(plot_prefix,'clusterLR.pdf')),width = 25,height = 35)
# # 
# # hm = similarityHeatmap(output,
# #                        row_order=row_order,
# #                        column_order = column_order,
# #                        row_title_rot = 0,
# #                        row_title_gp = gpar(fontsize=10),row_names_gp = gpar(fontsize=10),row_names_max_width = unit(6,'cm'),
# #                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
# #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = F)
# # draw(hm)
# # dev.off()
# # 
# # 
# # 
# #------- Single-cell level
# if(length(outputs)>2){
#   output = outputs[['scLR_all']]
# }else{
#   output = outputs[[2]][['scLR_all']]
# }
# 
# # 
# # # Only plots those of interest
# # ref_toKeep = c('diploid_HSC_MPP','T21_HSC_MPP','diploid_MEMP_MEP','T21_MEMP_MEP','diploid_MK','T21_MK','T21_EE','diploid_EE',"diploid_Mast.cell", "T21_Mast.cell",'T21_NK_T','diploid_NK_T')
# # tgt_cellstoKeep = tgt.srat$cellID[grepl('Tumour|HSC|MEP|EE|NK|MK',tgt.srat$finalAnn)]
# # ref_cellstoKeep = REF.srat$cellID[REF.srat$ann %in% ref_toKeep]
# # 
# # in_mtx = output[rownames(output) %in% c(tgt_cellstoKeep,paste0('ref_',ref_cellstoKeep)),colnames(output) %in% ref_toKeep]
# # type = tgt.srat$finalAnn[match(rownames(in_mtx),tgt.srat$cellID)]
# # type[is.na(type)] = REF.srat$ann[match(rownames(in_mtx)[is.na(type)],paste0('ref_',REF.srat$cellID))]
# # 
# # 
# # pdf('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/MLDS1023_fLiver2nT21REF1023_maxCells_70perc_scLR_MKlineage.pdf',width = 5,height = 35)
# # show_row_names=F
# # hm = similarityHeatmap(in_mtx,
# #                        column_order = ref_toKeep[grepl('diploid',ref_toKeep)],
# #                        row_title_rot = 0,
# #                        #row_title_gp = gpar(fontsize=5),row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
# #                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
# #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
# # draw(hm)
# # 
# # dev.off()
# # 
# # 
# # in_mtx = output
# # type = tgt.srat$finalAnn[match(rownames(output),tgt.srat$cellID)]
# # type[is.na(type)] = REF.srat$ann[match(rownames(output)[is.na(type)],paste0('ref_',REF.srat$cellID))]
# # 
# # pdf('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/MLDS1023_fLiver2nT21REF1023_maxCells_70perc_scLR_all_Heatmap.pdf',width = 50,height = 200)
# # show_row_names=F
# # hm = similarityHeatmap(in_mtx,
# #                        #column_order = ref_toKeep,
# #                        row_title_rot = 0,
# #                        #row_title_gp = gpar(fontsize=5),row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
# #                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
# #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
# # draw(hm)
# # 
# # dev.off()
# # 
# # 
# # 
# # 
# # 
# # 
# # ## For every tumour cell, assign it to the reference cell type of best match....
# # mtx = output[rownames(output) %in% tgt.srat$cellID[grepl('_Tumour$',tgt.srat$finalAnn)],grepl('diploid',colnames(output))]
# # mtx = output[rownames(output) %in% tgt.srat$cellID[grepl('_Tumour$',tgt.srat$finalAnn)],ref_toKeep[grepl('diploid',ref_toKeep)]]
# # mtx = mtx[,!is.na(colSums(mtx))]
# # bestMatch = do.call(c,lapply(1:nrow(mtx),function(i){
# #   n = colnames(mtx)[which(mtx[i,] == max(mtx[i,]))]
# #   if(length(n) > 1){
# #     n = paste(n,collapse = ':')
# #   }
# #   return(n)
# # }))
# # 
# # 
# # df = data.frame(cellID = rownames(mtx),
# #                 bestMatch = bestMatch)
# # df$tumourAnn = paste0(tgt.srat$donorID[match(df$cellID,tgt.srat$cellID)],':Tumour')
# # df$isMLDS = ifelse(is.na(tgt.srat$finalAnn_broad[match(df$cellID,tgt.srat$cellID)]),F,T)
# # df$bestMatch2 = df$bestMatch
# # df$bestMatch2[grepl(':',df$bestMatch2)] = 'multiple'
# # df$bestMatch2[grepl('B.cell',df$bestMatch2)] = 'diploid_Blineage'
# # df$bestMatch2[grepl('NK|T|ILC',df$bestMatch2)] = 'diploid_NK_T'
# # df$bestMatch2[grepl('Mono|DC|MOP|Mac|myelocyte',df$bestMatch2)] = 'diploid_myeloid'
# # 
# # 
# # ggplot(df[df$isMLDS & df$tumourAnn != 'L041:Tumour',],aes(tumourAnn,fill=bestMatch2))+
# #   geom_bar(position="fill")+
# #   scale_fill_manual(values = c(col25))+
# #   theme_classic() + theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))+ylab('Fraction of Blasts') + xlab('')
# # 
# # 
# # 
# # ## subset to MK/EE/Mast lineage only
# # mtx = output[rownames(output) %in% tgt.srat$cellID[grepl('_Tumour$',tgt.srat$finalAnn)],c(ref_toKeep[grepl('diploid',ref_toKeep)],'diploid_proMono','diploid_CMP_GMP','diploid_Monocyte')]
# # mtx = mtx[,!is.na(colSums(mtx))]
# # bestMatch = do.call(c,lapply(1:nrow(mtx),function(i){
# #   n = colnames(mtx)[which(mtx[i,] == max(mtx[i,]))]
# #   if(length(n) > 1){
# #     n = paste(n,collapse = ':')
# #   }
# #   return(n)
# # }))
# # 
# # 
# # df = data.frame(cellID = rownames(mtx),
# #                 bestMatch = bestMatch)
# # df$tumourAnn = paste0(tgt.srat$donorID[match(df$cellID,tgt.srat$cellID)],':Tumour')
# # df$isMLDS = ifelse(is.na(tgt.srat$finalAnn_broad[match(df$cellID,tgt.srat$cellID)]),F,T)
# # df$bestMatch = factor(df$bestMatch,levels = rev(c('diploid_HSC_MPP','diploid_MEMP_MEP','diploid_MK','diploid_EE','diploid_Mast.cell','diploid_NK_T',unique(bestMatch[!bestMatch %in% c('diploid_HSC_MPP','diploid_MEMP_MEP','diploid_MK','diploid_EE','diploid_Mast.cell','diploid_NK_T')]))))
# # 
# # df$tumourAnn[df$tumourAnn == 'L075:Tumour'] = 'TAM:L075'
# # ggplot(df[df$isMLDS & df$tumourAnn != 'L041:Tumour',],aes(tumourAnn,fill=bestMatch))+
# #   geom_bar(position="fill")+
# #   scale_fill_manual(values = c(col25))+
# #   theme_classic(base_size = 14) + theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))+ylab('Fraction of Blasts') + xlab('')
# # 
# # ## Heatmap
# # logitCols = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
# # cols = circlize::colorRamp2(seq(-5,5,length.out=length(logitCols)),logitCols)
# # 
# # mtx = output[rownames(output) %in% tgt.srat$cellID[grepl('_Tumour$',tgt.srat$finalAnn)],c(ref_toKeep[grepl('diploid',ref_toKeep)],'diploid_proMono','diploid_CMP_GMP','diploid_Monocyte')]
# # pdf(file.path(outDir,'MLDS1023_fLiver2nT21REF1023_maxCells_70perc_scLR_MKlineage_subset.pdf'),width = 20,height = 50)
# # hm = Heatmap(as.matrix(mtx),show_row_dend = F,show_column_dend = F,
# #         col = cols,
# #         show_column_names = T,show_row_names = F,row_title_rot = 0,
# #         split = df$tumourAnn[match(rownames(mtx),df$cellID)])
# # draw(hm)
# # dev.off()
# # 
# # 
# # 
# # 
# # 
# ##----------------------------       Ellie's style plot     -------------------------------------#####
# ## Plot the range of similarity score for each Tumour_patientID against relevant Reference cell type #
# model = readRDS(model_fp)
# 
# mtx = predictSimilarity(model, tgt.srat@assays$RNA@counts,logits = T)
# pp = predictSimilarity(model, tgt.srat@assays$RNA@counts,logits = F)
# 
# ## subset to MK/EE/Mast lineage only
# ref_toKeep = c('diploid_HSC_MPP','T21_HSC_MPP','diploid_MEMP_MEP','T21_MEMP_MEP','diploid_MK','T21_MK','T21_EE','diploid_EE',"diploid_Mast.cell", "T21_Mast.cell",'T21_NK_T','diploid_NK_T')
# mtx = pp[rownames(pp) %in% tgt.srat$cellID[grepl('_Tumour:TAM|_Tumour:MLDS|Tumour:pAML|Tumour:infantALL|MEP:fLiver:diploid|MEP:fLiver:T21|MK:fLiver:diploid|MK:fLiver:T21',tgt.srat$ann)],
#          c(ref_toKeep[grepl('T21',ref_toKeep)],'T21_B.cell','T21_Monocyte')]
# mtx = mtx[,!is.na(colSums(mtx))]
# source('~/lustre_mt22/generalScripts/sharedCode-main/logisticRegressionCellTypist.R')
# mtx = as.data.frame(mtx)
# mtx$cellID = rownames(mtx)
# df = pivot_longer(mtx,cols = colnames(mtx)[colnames(mtx) != 'cellID'],names_to = 'REF_celltype',values_to = 'LRpp_score')
# df$group = tgt.srat$ann[match(df$cellID,tgt.srat$cellID)]
# df$timePoint = tgt.srat$timePoint[match(df$cellID,tgt.srat$cellID)]
# df$group[df$group == 'L038_Tumour:MLDS:MLDS:T21' & df$timePoint == 'postChemo TP1'] = 'L038:TP1'
# 
# #df = df[df$group == 'L019:Tumour',]
# df$LR_score_group = ifelse(df$LRpp_score > 0.99, 'high',
#                            ifelse(df$LRpp_score > 0.9, 'mid2',
#                                   ifelse(df$LRpp_score < 0.2,'low','mid')))
# df$LR_score_group = factor(df$LR_score_group,levels = rev(c('high','mid2','mid','low')))
# df$REF_celltype = factor(df$REF_celltype,levels = c(ref_toKeep[grepl('T21',ref_toKeep)],'T21_B.cell','T21_Monocyte'))
# df$group2 = ifelse(grepl('L075',df$group),'TAM',
#                    ifelse(grepl('L019|L039|L040|L042|L76|L091|L038',df$group),'MLDS',
#                           ifelse(grepl('infantALL',df$group),'infantALL',
#                                  ifelse(grepl('MEP|MK',df$group),'normalCtrl','pAML'))))
# df$timePoint = tgt.srat$timePoint[match(df$cellID,tgt.srat$cellID)]
# 
# #df$group2 = factor(df$group2,c('TAM','MLDS','pAML'))
# df$toKeep = ifelse(grepl('Tumour',df$group) & !grepl('infantALL',df$group) & !df$timePoint %in% c('D0','Diagnostic')  ,F,T)
# ggplot(df[df$toKeep == T,],aes(y=group,fill = LR_score_group))+
#   geom_bar(position = 'fill',col='black')+
#   facet_grid(vars(group2),vars(REF_celltype),scales = 'free_y',space = 'free')+
#   scale_fill_manual(values = rev(c('#EA2335','#F06571','white','#2D3F90'))) + theme_classic(base_size = 13) + xlab('Fraction of cells') + ylab('')
# 
# # plotByGroup(fit=mtx,cellGroupings=df$group,colSpaceTruncate=T,colSpaceBreaks = c(-5,0,5),cluster_columns = F)
# # 
# # plotByGroup = function(fit,cellGroupings,
# #                        colSpaceBreaks=c(0,0.5,1),colSpaceCols =c('#2D3F90','white','#EA2335'),colSpaceBins=20,colSpaceTruncate=FALSE,bgColor='grey',...){
# #   #Group data
# #   collapsedDat = collapseToClusters(fit,cellGroupings,outputForFancyHeatmap=TRUE)
# #   #Define the colour scheme
# #   col_fun = circlize::colorRamp2(colSpaceBreaks,colSpaceCols)
# #   colSpaceRange = range(colSpaceBreaks)
# #   #Fix colour space if required
# #   if(colSpaceTruncate){
# #     collapsedDat$mat[collapsedDat$mat < min(colSpaceBreaks)] = min(colSpaceBreaks)
# #     collapsedDat$mat[collapsedDat$mat > max(colSpaceBreaks)] = max(colSpaceBreaks)
# #     for(i in seq_along(collapsedDat$x)){
# #       for(j in seq_along(collapsedDat$x[[i]])){
# #         collapsedDat$x[[i]][[j]][collapsedDat$x[[i]][[j]] < min(colSpaceBreaks)] = min(colSpaceBreaks)
# #         collapsedDat$x[[i]][[j]][collapsedDat$x[[i]][[j]] > max(colSpaceBreaks)] = max(colSpaceBreaks)
# #       }
# #     }
# #   }
# #   #Define the cell plotting function
# #   cell_fun = function(j, i, x, y, width, height, fill) {
# #     #Define number of boxes to split into
# #     N=colSpaceBins
# #     #Do tons of teny tiny rectangles
# #     dat = collapsedDat$x[[i]][[j]]
# #     #dat = (1+exp(-sortDat[[i]][[j]]))**-1
# #     #White out the area
# #     grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = NA, fill = bgColor))
# #     #Draw N boxes per thingy
# #     #Have a bit less width available for this than the whole region
# #     subWidth = width*0.95
# #     subHeight = height*1.0
# #     sDat = split(dat,cut(dat,seq(colSpaceRange[1],colSpaceRange[2],length.out=N+1),include.lowest=TRUE))
# #     boxWidths = lengths(sDat)/length(dat)
# #     for(ii in seq_along(sDat)){
# #       grid.rect(x = (x-subWidth/2) + (sum(boxWidths[which(seq_along(boxWidths)<=(ii-1))]) + boxWidths[ii]/2)*subWidth,
# #                 y = y,
# #                 width = boxWidths[ii]*subWidth,
# #                 height = subHeight,
# #                 gp = gpar(fill=col_fun(mean(sDat[[ii]])),col=NA)
# #       )
# #     }
# #     #Draw the bounding box for the whole thing
# #     grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = bgColor, fill = NA))
# #   }
# #   hmParams = list(matrix = collapsedDat$mat,
# #                   name='Similarity',
# #                   col=col_fun,
# #                   cell_fun = cell_fun,
# #                   show_row_dend=FALSE,
# #                   show_column_dend=FALSE)
# #   theDots = list(...)
# #   for(nom in names(theDots))
# #     hmParams[[nom]] = theDots[[nom]]
# #   hm = do.call(Heatmap,hmParams)
# #   return(hm)
# # }
# 
# 
# 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # srat = standard_clustering(srat)
# # ### Plot distribution of logit scores of Cancer_cell matched to REF categories of interest
# # cancer_cell_LRscore = output[rownames(output) %in% tgt.srat$cellID[tgt.srat$cluster_ann == 'Cancer_cell'],]
# # cancer_cell_LRscore = cancer_cell_LRscore[,!is.na(colSums(cancer_cell_LRscore))]
# # cancer_cell_LRscore = as.data.frame(cancer_cell_LRscore)
# # cancer_cell_LRscore$cellID = rownames(cancer_cell_LRscore)
# # cancer_cell_LRscore = pivot_longer(cancer_cell_LRscore,cols = c(1:ncol(cancer_cell_LRscore)-1),names_to = 'REF_type',values_to = 'LR_score')
# # sum_LR_score = cancer_cell_LRscore[,colnames(cancer_cell_LRscore) != "cellID"] %>% group_by(REF_type) %>% summarise(med_LR_score = median(LR_score))
# # sum_LR_score = sum_LR_score[order(sum_LR_score$med_LR_score,decreasing = T),]
# # cancer_cell_LRscore$REF_type = factor(cancer_cell_LRscore$REF_type,levels = sum_LR_score$REF_type)
# # pdf('2_LRorig_fLiverREF_onMLDS/mldsCancerCell_LRscore_allREFcelltype_1.pdf',width = 30,height = 25)
# # p=ggplot(cancer_cell_LRscore,aes(LR_score,fill=REF_type))+
# #   geom_density()+
# #   geom_vline(xintercept = 0)+
# #   theme_bw()+
# #   theme(legend.position = 'none')+
# #   facet_wrap(vars(REF_type),scales = 'free_y')
# # 
# # print(p)
# # dev.off()
# # pdf('2_LRorig_fLiverREF_onMLDS/mldsCancerCell_LRscore_allREFcelltype_2_sub.pdf',width = 7,height = 5)
# # p=ggplot(cancer_cell_LRscore[cancer_cell_LRscore$REF_type %in% ref_toKeep,],aes(LR_score,fill=REF_type))+
# #   geom_density(alpha=0.7)+
# #   geom_vline(xintercept = 0)+
# #   scale_fill_brewer(palette = 'Dark2')+
# #   theme_bw()
# #   #theme(legend.position = 'none')
# #   #facet_wrap(vars(REF_type))
# # 
# # print(p)
# # dev.off()
# # 
# # 
# # #### April 2023:
# # cancer_cell_LRscore = output[rownames(output) %in% mlds$cellID[mlds$finalAnn_broad=='Tumour'],]
# # cancer_cell_LRscore = cancer_cell_LRscore[,grepl('HSC|MEMP|MK|Mega|EE|Early.Er',colnames(cancer_cell_LRscore))]
# # cancer_cell_LRscore = as.data.frame(cancer_cell_LRscore)
# # cancer_cell_LRscore$cellID = rownames(cancer_cell_LRscore)
# # cancer_cell_LRscore = pivot_longer(cancer_cell_LRscore,cols = c(1:ncol(cancer_cell_LRscore)-1),names_to = 'REF_type',values_to = 'LR_score')
# # cancer_cell_LRscore$donorID = mlds$donorID[match(cancer_cell_LRscore$cellID,mlds$cellID)]
# # cancer_cell_LRscore$REF_ct = gsub('T21_|T18_|T22_|MX_|complete_trisomy_|diploid_','',cancer_cell_LRscore$REF_type)
# # cancer_cell_LRscore$REF_geno = sapply(strsplit(cancer_cell_LRscore$REF_type,split='_'),'[',1)
# # ggplot(cancer_cell_LRscore,aes(REF_geno,LR_score,fill=donorID))+
# #   geom_boxplot()+
# #   geom_hline(yintercept = 0)+
# #   theme_bw()+ylim(-10,10)+theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
# #   facet_wrap(vars(REF_ct),scales = 'free_x')
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ##------------------------------------------------##
# # ####   3. LR_orig with REF = fBM ; tgt = MLDS   ####
# # ##------------------------------------------------##
# # 
# # ## Import REF sratObj = scRNA fBM (diploid + T21)
# # 
# # # fBM diploid
# # fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_diploid_Laura_sratObj.RDS'
# # if(file.exists(fBM_fp)){
# #   fBM = readRDS(fBM_fp)
# # }else{
# #   stop(sprintf('Cannot find file(s) below, please check!\n%s',fBM_fp))
# # }
# # 
# # # fBM T21
# # t21_fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_T21_Laura_sratObj.RDS'
# # if(file.exists(t21_fBM_fp)){
# #   t21_fBM = readRDS(t21_fBM_fp)
# # }else{
# #   stop(sprintf('Cannot find file(s) below, please check!\n%s',t21_fBM_fp))
# # }
# # 
# # REF.srat = merge_seurat_objects(fBM,t21_fBM,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # rm(fBM)
# # rm(t21_fBM)
# # 
# # 
# # REF.srat$finalAnn2 = REF.srat$finalAnn
# # REF.srat$finalAnn2[grepl('MEMP|MEP',REF.srat$finalAnn2)] = 'MEMP_MEP'
# # REF.srat$finalAnn2[grepl('MK',REF.srat$finalAnn2)] = 'MK'
# # REF.srat$finalAnn2[grepl('early erythroid',REF.srat$finalAnn2)] = 'EE'
# # 
# # REF.srat$ann = paste0(REF.srat$Genotype,'_',REF.srat$finalAnn2)
# # 
# # # Remove categories with < 10 cells
# # #ct_to_remove = names(table(REF.srat$ann)[table(REF.srat$ann) < 10])
# # #REF.srat = subset(REF.srat,subset = ann %in% unique(REF.srat$ann[!REF.srat$ann %in% ct_to_remove]))
# # 
# # 
# # message(('1. REF.srat loaded'))
# # 
# # ## Import tgt sratObj = scRNA MLDS
# # #tgt.srat = readRDS(('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov22/MLDSonly_clean_LRwCTannotated_v3.RDS'))
# # tgt.srat = readRDS(('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov22/MLDS_clean_LRwCTannotated_jan23.RDS'))
# # message('2. tgt.srat loaded')
# # 
# # 
# # ## Configure REF and tgt sratObj for LR
# # # Only keep the same genes between REF and tgt.srat
# # REF.srat$type = 'fBM'
# # tgt.srat$type = 'MLDS'
# # merged.srat = merge_seurat_objects(srat1 = REF.srat,srat2 = tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # 
# # tgt.srat = subset(merged.srat,subset = type == 'MLDS')
# # REF.srat = subset(merged.srat,subset = type == 'fBM')
# # 
# # tgt.srat$final_broadAnn[tgt.srat$final_broadAnn == 'Tumour'] = paste0(tgt.srat$donorID[tgt.srat$final_broadAnn == 'Tumour'],'_',tgt.srat$final_broadAnn[tgt.srat$final_broadAnn == 'Tumour'])
# # 
# # # If we want to train LR model on specific cell types only - subset REF.srat to include these CTs and remove others
# # #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Hepatocyte','Early.Erythroid','Megakaryocyte','Fibroblast','CMP','GMP','Mono.Mac'))
# # #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Early.Erythroid','Megakaryocyte'))
# # #REF.srat = subset(REF.srat,subset = Genotype %in% c('T21','diploid'))
# # 
# # ##----------------------------##
# # ##   Run Logistic Regression  ##
# # ##----------------------------##
# # skipIfExists=F
# # ref_annot='ann'
# # maxCells=4000
# # outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fBoneMarrowREF_onMLDS/')
# # model_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fBoneMarrowREF_onMLDS/fBM_diploid_T21_REF_trainModel_4kmaxcells_70perc_jan23.RDS'
# # 
# # LR_level='both'
# # srat_annot='final_broadAnn'
# # minGeneMatch = 0.99
# # maxCells=4000
# # tissue = 'BoneMarrow'
# # 
# # out_prefix = 'MLDS0123_fBoneMarrowREF_maxCells_70perc_'
# # plot_prefix = 'MLDS0123_fBoneMarrowREF_maxCells_70perc_'
# # plot_prefix = NULL
# # 
# # outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
# #                 model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
# #                 minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
# #                 scLR_TGTtype='',scLR_REFtype='')
# # 
# # message(sprintf('3. LR completed for tissue %s',tissue))
# # 
# # 
# # ##---------------------------------##
# # ##   Do some LR_similarity plots   ##
# # ##---------------------------------##
# # 
# # ## annotated Cluster level #
# # output = outputs[[2]][[1]]
# # type = ifelse(grepl('ref_',rownames(output)),'REF','TGT')
# # show_row_names = T
# # column_order = colnames(output)[order(colnames(output))]
# # column_order = column_order[order(gsub('diploid_|T21_','',column_order))]
# # row_order = rownames(output)[!grepl('Tumour',rownames(output))]
# # row_order = row_order[order(row_order)]
# # row_order = c(row_order,rownames(output)[grepl('Tumour',rownames(output))])
# # 
# # 
# # plot_prefix = 'MLDS0123_fBoneMarrowREF_maxCells_70perc_'
# # 
# # pdf(file.path(outDir,paste0(plot_prefix,'_clusterLR.pdf')),width = 20,height = 25)
# # 
# # hm = similarityHeatmap(output,
# #                        row_order=row_order,
# #                        column_order = column_order,
# #                        row_title_rot = 0,
# #                        row_title_gp = gpar(fontsize=10),row_names_gp = gpar(fontsize=10),row_names_max_width = unit(6,'cm'),
# #                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
# #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = F)
# # draw(hm)
# # dev.off()
# # 
# # 
# # 
# # output = outputs[[2]][['scLR_all']]
# # # Only plots those of interest
# # #ref_toKeep = c('diploid_HSC_MPP','T21_HSC_MPP','diploid_MEMP_MEP','T21_MEMP_MEP','diploid_Megakaryocyte','T21_Megakaryocyte','T21_Early.Erythroid','diploid_Early.Erythroid')
# # #ref_toKeep = unique(REF.srat$ann[grepl('erythroid','EE','')])[grepl()]
# # tgt_cellstoKeep = tgt.srat$cellID[grepl('Cancer',tgt.srat$cluster_ann)]
# # ref_cellstoKeep = REF.srat$cellID[REF.srat$ann %in% ref_toKeep]
# # 
# # in_mtx = output[rownames(output) %in% c(tgt_cellstoKeep,paste0('ref_',ref_cellstoKeep)),colnames(output) %in% ref_toKeep]
# # type = tgt.srat$cluster_ann[match(rownames(in_mtx),tgt.srat$cellID)]
# # type[is.na(type)] = REF.srat$ann[match(rownames(in_mtx)[is.na(type)],paste0('ref_',REF.srat$cellID))]
# # 
# # 
# # type = tgt.srat$cluster_ann[match(rownames(output),tgt.srat$cellID)]
# # type[is.na(type)] = REF.srat$ann[match(rownames(output)[is.na(type)],paste0('ref_',REF.srat$cellID))]
# # 
# # pdf('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/MLDS_fLiverREF_maxCells_70perc_scLR_all_Heatmap.pdf',width = 50,height = 200)
# # pdf('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fLiverREF_onMLDS/MLDS_fLiverREF_maxCells_70perc_scLR_MKlineage.pdf',width = 5,height = 35)
# # hm = similarityHeatmap(in_mtx,
# #                        column_order = ref_toKeep,
# #                        #show_row_names=F,
# #                        row_title_rot = 0,
# #                        #row_title_gp = gpar(fontsize=5),row_names_gp = gpar(fontsize=5),row_names_max_width = unit(6,'cm'),
# #                        column_names_gp = gpar(fontsize=10),column_names_max_height = unit(6,'cm'),
# #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
# # draw(hm)
# # 
# # dev.off()
# # 
# # srat = standard_clustering(srat)
# # ### Plot distribution of logit scores of Cancer_cell matched to REF categories of interest
# # cancer_cell_LRscore = output[rownames(output) %in% tgt.srat$cellID[tgt.srat$cluster_ann == 'Cancer_cell'],]
# # cancer_cell_LRscore = cancer_cell_LRscore[,!is.na(colSums(cancer_cell_LRscore))]
# # cancer_cell_LRscore = as.data.frame(cancer_cell_LRscore)
# # cancer_cell_LRscore$cellID = rownames(cancer_cell_LRscore)
# # cancer_cell_LRscore = pivot_longer(cancer_cell_LRscore,cols = c(1:ncol(cancer_cell_LRscore)-1),names_to = 'REF_type',values_to = 'LR_score')
# # sum_LR_score = cancer_cell_LRscore[,colnames(cancer_cell_LRscore) != "cellID"] %>% group_by(REF_type) %>% summarise(med_LR_score = median(LR_score))
# # sum_LR_score = sum_LR_score[order(sum_LR_score$med_LR_score,decreasing = T),]
# # cancer_cell_LRscore$REF_type = factor(cancer_cell_LRscore$REF_type,levels = sum_LR_score$REF_type)
# # pdf('2_LRorig_fLiverREF_onMLDS/mldsCancerCell_LRscore_allREFcelltype_1.pdf',width = 30,height = 25)
# # p=ggplot(cancer_cell_LRscore,aes(LR_score,fill=REF_type))+
# #   geom_density()+
# #   geom_vline(xintercept = 0)+
# #   theme_bw()+
# #   theme(legend.position = 'none')+
# #   facet_wrap(vars(REF_type),scales = 'free_y')
# # 
# # print(p)
# # dev.off()
# # pdf('2_LRorig_fLiverREF_onMLDS/mldsCancerCell_LRscore_allREFcelltype_2_sub.pdf',width = 7,height = 5)
# # p=ggplot(cancer_cell_LRscore[cancer_cell_LRscore$REF_type %in% ref_toKeep,],aes(LR_score,fill=REF_type))+
# #   geom_density(alpha=0.7)+
# #   geom_vline(xintercept = 0)+
# #   scale_fill_brewer(palette = 'Dark2')+
# #   theme_bw()
# # #theme(legend.position = 'none')
# # #facet_wrap(vars(REF_type))
# # 
# # print(p)
# # dev.off()
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ##---------------------------------------------------------##
# # ####   4. LR_orig with REF = fBM + fLiver ; tgt = MLDS   ####
# # ##---------------------------------------------------------##
# # 
# # ## Import REF sratObj = scRNA fBM (diploid + T21) + fLiver (all genotypes)
# # 
# # # fBM diploid
# # fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_diploid_Laura_sratObj.RDS'
# # if(file.exists(fBM_fp)){
# #   fBM = readRDS(fBM_fp)
# # }else{
# #   stop(sprintf('Cannot find file(s) below, please check!\n%s',fBM_fp))
# # }
# # 
# # # fBM T21
# # t21_fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_T21_Laura_sratObj.RDS'
# # if(file.exists(t21_fBM_fp)){
# #   t21_fBM = readRDS(t21_fBM_fp)
# # }else{
# #   stop(sprintf('Cannot find file(s) below, please check!\n%s',t21_fBM_fp))
# # }
# # 
# # REF.srat = merge_seurat_objects(fBM,t21_fBM,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # rm(fBM)
# # rm(t21_fBM)
# # 
# # 
# # REF.srat$finalAnn2 = REF.srat$finalAnn
# # REF.srat$finalAnn2[grepl('MEMP|MEP',REF.srat$finalAnn2)] = 'MEMP_MEP'
# # REF.srat$finalAnn2[grepl('MK',REF.srat$finalAnn2)] = 'MK'
# # REF.srat$finalAnn2[grepl('early erythroid',REF.srat$finalAnn2)] = 'EE'
# # 
# # REF.srat$ann = paste0(REF.srat$Genotype,'_',REF.srat$finalAnn2)
# # REF.srat$type = 'fBM'
# # 
# # 
# # ## Import fLiver
# # fLiver_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/oct22/liver_liverREFmerged_clean_processed_annotated_v2.RDS'
# # 
# # if(file.exists(fLiver_fp)){
# #   fLiver = readRDS(fLiver_fp)
# #   fLiver$type = 'fLiver'
# #   fLiver$ann = paste0(fLiver$Genotype,'_',fLiver$finalAnn_broad)
# # }
# # 
# # 
# # REF.srat = merge_seurat_objects(REF.srat,fLiver,keepAllGenes = FALSE,genomeVersions = c('v38','v38'))
# # rm(fLiver)
# # REF.srat$ann = paste0(REF.srat$ann,'_',REF.srat$type)
# # 
# # message(('1. REF.srat loaded'))
# # 
# # # Remove categories with < 10 cells
# # ct_to_remove = names(table(REF.srat$ann)[table(REF.srat$ann) < 10])
# # REF.srat = subset(REF.srat,subset = ann %in% unique(REF.srat$ann[!REF.srat$ann %in% ct_to_remove]))
# # 
# # 
# # message(('1. REF.srat loaded'))
# # 
# # ## Import tgt sratObj = scRNA MLDS
# # tgt.srat = readRDS(('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov22/MLDSonly_clean_LRwCTannotated_v3.RDS'))
# # 
# # message('2. tgt.srat loaded')
# # 
# # 
# # ## Configure REF and tgt sratObj for LR
# # # Only keep the same genes between REF and tgt.srat
# # REF.srat$type = 'fBM'
# # tgt.srat$type = 'MLDS'
# # merged.srat = merge_seurat_objects(srat1 = REF.srat,srat2 = tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # 
# # tgt.srat = subset(merged.srat,subset = type == 'MLDS')
# # REF.srat = subset(merged.srat,subset = type == 'fBM')
# # 
# # # If we want to train LR model on specific cell types only - subset REF.srat to include these CTs and remove others
# # #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Hepatocyte','Early.Erythroid','Megakaryocyte','Fibroblast','CMP','GMP','Mono.Mac'))
# # #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Early.Erythroid','Megakaryocyte'))
# # #REF.srat = subset(REF.srat,subset = Genotype %in% c('T21','diploid'))
# # 
# # ##----------------------------##
# # ##   Run Logistic Regression  ##
# # ##----------------------------##
# # skipIfExists=F
# # ref_annot='ann'
# # maxCells=4000
# # outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fBoneMarrowREF_onMLDS/')
# # #model_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/AdrKidLiv_REF_trainModel_4kmaxcells_70perc.RDS'
# # model_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fBoneMarrowREF_onMLDS/fBM_fLiver_allGeno_REF_trainModel_4kmaxcells_70perc.RDS'
# # 
# # LR_level='both'
# # srat_annot='cluster_ann'
# # minGeneMatch = 0.99
# # maxCells=4000
# # tissue = 'fLiver_BoneMarrow'
# # 
# # out_prefix = 'MLDS_fLiver_fBoneMarrowREF_maxCells_70perc_'
# # plot_prefix = 'MLDS_fLiver_fBoneMarrowREF_maxCells_70perc_'
# # plot_prefix = NULL
# # 
# # outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
# #                 model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
# #                 minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
# #                 scLR_TGTtype='',scLR_REFtype='')
# # 
# # message(sprintf('3. LR completed for tissue %s',tissue))
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ##------------------------------------------------------------##
# # ####   3. LR_orig with REF = fBM ; tgt = MLDS_tumour only   ####
# # ##-----------------------------------------------------------##
# # 
# # ## Import REF sratObj = scRNA fBM (diploid + T21)
# # 
# # # fBM diploid
# # fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_diploid_Laura_sratObj.RDS'
# # if(file.exists(fBM_fp)){
# #   fBM = readRDS(fBM_fp)
# # }else{
# #   stop(sprintf('Cannot find file(s) below, please check!\n%s',fBM_fp))
# # }
# # 
# # # fBM T21
# # t21_fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_T21_Laura_sratObj.RDS'
# # if(file.exists(t21_fBM_fp)){
# #   t21_fBM = readRDS(t21_fBM_fp)
# # }else{
# #   stop(sprintf('Cannot find file(s) below, please check!\n%s',t21_fBM_fp))
# # }
# # 
# # REF.srat = merge_seurat_objects(fBM,t21_fBM,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # rm(fBM)
# # rm(t21_fBM)
# # 
# # 
# # REF.srat$finalAnn2 = REF.srat$finalAnn
# # REF.srat$finalAnn2[grepl('MEMP|MEP',REF.srat$finalAnn2)] = 'MEMP_MEP'
# # REF.srat$finalAnn2[grepl('MK',REF.srat$finalAnn2)] = 'MK'
# # REF.srat$finalAnn2[grepl('early erythroid',REF.srat$finalAnn2)] = 'EE'
# # 
# # REF.srat$ann = paste0(REF.srat$Genotype,'_',REF.srat$finalAnn2)
# # 
# # # Remove categories with < 10 cells
# # #ct_to_remove = names(table(REF.srat$ann)[table(REF.srat$ann) < 10])
# # #REF.srat = subset(REF.srat,subset = ann %in% unique(REF.srat$ann[!REF.srat$ann %in% ct_to_remove]))
# # 
# # 
# # message(('1. REF.srat loaded'))
# # 
# # ## Import tgt sratObj = scRNA MLDS
# # #tgt.srat = readRDS(('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov22/MLDSonly_clean_LRwCTannotated_v3.RDS'))
# # tgt.srat = readRDS(('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov22/MLDS_cancerOnly_jan23.RDS'))
# # message('2. tgt.srat loaded')
# # 
# # 
# # ## Configure REF and tgt sratObj for LR
# # # Only keep the same genes between REF and tgt.srat
# # REF.srat$type = 'fBM'
# # tgt.srat$type = 'MLDS'
# # merged.srat = merge_seurat_objects(srat1 = REF.srat,srat2 = tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # 
# # tgt.srat = subset(merged.srat,subset = type == 'MLDS')
# # REF.srat = subset(merged.srat,subset = type == 'fBM')
# # 
# # tgt.srat$final_broadAnn[tgt.srat$final_broadAnn == 'Tumour'] = paste0(tgt.srat$donorID[tgt.srat$final_broadAnn == 'Tumour'],'_',tgt.srat$final_broadAnn[tgt.srat$final_broadAnn == 'Tumour'])
# # 
# # # If we want to train LR model on specific cell types only - subset REF.srat to include these CTs and remove others
# # #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Hepatocyte','Early.Erythroid','Megakaryocyte','Fibroblast','CMP','GMP','Mono.Mac'))
# # #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Early.Erythroid','Megakaryocyte'))
# # #REF.srat = subset(REF.srat,subset = Genotype %in% c('T21','diploid'))
# # 
# # ##----------------------------##
# # ##   Run Logistic Regression  ##
# # ##----------------------------##
# # skipIfExists=F
# # ref_annot='ann'
# # maxCells=4000
# # outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fBoneMarrowREF_onMLDS/')
# # model_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fBoneMarrowREF_onMLDS/fBM_diploid_T21_REF_trainModel_4kmaxcells_70perc_jan23.RDS'
# # 
# # LR_level='both'
# # srat_annot=''
# # minGeneMatch = 0.99
# # maxCells=4000
# # tissue = 'BoneMarrow'
# # 
# # out_prefix = 'MLDS0123_fBoneMarrowREF_maxCells_70perc_'
# # plot_prefix = 'MLDS0123_fBoneMarrowREF_maxCells_70perc_'
# # plot_prefix = NULL
# # 
# # outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
# #                 model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
# #                 minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
# #                 scLR_TGTtype='',scLR_REFtype='')
# # 
# # message(sprintf('3. LR completed for tissue %s',tissue))
# # 
# # 
# # 
# # 
# # 
# # 
# # ##------------------------------------------------------------##
# # ####   3. LR_orig with REF = fBM ; tgt = fLiver   ####
# # ##-----------------------------------------------------------##
# # 
# # ## Import REF sratObj = scRNA fBM (diploid + T21)
# # 
# # # fBM diploid
# # fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_diploid_Laura_sratObj.RDS'
# # if(file.exists(fBM_fp)){
# #   fBM = readRDS(fBM_fp)
# # }else{
# #   stop(sprintf('Cannot find file(s) below, please check!\n%s',fBM_fp))
# # }
# # 
# # # fBM T21
# # t21_fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_T21_Laura_sratObj.RDS'
# # if(file.exists(t21_fBM_fp)){
# #   t21_fBM = readRDS(t21_fBM_fp)
# # }else{
# #   stop(sprintf('Cannot find file(s) below, please check!\n%s',t21_fBM_fp))
# # }
# # 
# # REF.srat = merge_seurat_objects(fBM,t21_fBM,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # rm(fBM)
# # rm(t21_fBM)
# # 
# # 
# # REF.srat$finalAnn2 = REF.srat$finalAnn
# # REF.srat$finalAnn2[grepl('MEMP|MEP',REF.srat$finalAnn2)] = 'MEMP_MEP'
# # REF.srat$finalAnn2[grepl('MK',REF.srat$finalAnn2)] = 'MK'
# # REF.srat$finalAnn2[grepl('early erythroid',REF.srat$finalAnn2)] = 'EE'
# # 
# # REF.srat$ann = paste0(REF.srat$Genotype,'_',REF.srat$finalAnn2)
# # 
# # ## Cluster fBM ref.srat
# # REF.srat = standard_clustering(REF.srat)
# # DimPlot(REF.srat,group.by = 'Genotype')
# # DimPlot(REF.srat,group.by = 'finalAnn2',label = T,label.size = 3,repel = T) + NoLegend()
# # DimPlot(REF.srat,cells.highlight = rownames(REF.srat@meta.data[REF.srat$finalAnn2 == 'HSC/MPP' &
# #                                                                  REF.srat$Genotype != 'T21',]))
# # # Remove categories with < 10 cells
# # #ct_to_remove = names(table(REF.srat$ann)[table(REF.srat$ann) < 10])
# # #REF.srat = subset(REF.srat,subset = ann %in% unique(REF.srat$ann[!REF.srat$ann %in% ct_to_remove]))
# # 
# # 
# # message(('1. REF.srat loaded'))
# # 
# # ## Import tgt sratObj = scRNA fLiver (all genotypes)
# # tgt.srat.fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/oct22/liver_liverREFmerged_clean_processed_annotated_v2.RDS'
# # 
# # if(file.exists(tgt.srat.fp)){
# #   tgt.srat = readRDS(tgt.srat.fp)
# #   tgt.srat$ann = paste0(tgt.srat$Genotype,'_',tgt.srat$finalAnn_broad_2)
# # }
# # 
# # message('2. tgt.srat loaded')
# # 
# # 
# # ## Configure REF and tgt sratObj for LR
# # # Only keep the same genes between REF and tgt.srat
# # REF.srat$type = 'fBM'
# # tgt.srat$type = 'fLiver'
# # merged.srat = merge_seurat_objects(srat1 = REF.srat,srat2 = tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# # 
# # tgt.srat = subset(merged.srat,subset = type == 'fLiver')
# # REF.srat = subset(merged.srat,subset = type == 'fBM')
# # 
# # merged.srat = standard_clustering(merged.srat)
# # DimPlot(merged.srat, group.by = 'type')
# # DimPlot(merged.srat, group.by = 'Genotype')
# # DimPlot(merged.srat, group.by = 'ann',label = T,label.size = 3,repel = T) + NoLegend()
# # DimPlot(merged.srat, cells.highlight = rownames(merged.srat@meta.data[merged.srat$ann=='T21_HSC_MPP' & merged.srat$type == 'fLiver',]))
# # Idents(merged.srat) = 'ann'
# # m = FindMarkers(merged.srat,ident.1 = 'T21_EE',ident.2 = 'T21_Early.Erythroid')
# # # If we want to train LR model on specific cell types only - subset REF.srat to include these CTs and remove others
# # #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Hepatocyte','Early.Erythroid','Megakaryocyte','Fibroblast','CMP','GMP','Mono.Mac'))
# # #REF.srat = subset(REF.srat,subset = finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','Early.Erythroid','Megakaryocyte'))
# # #REF.srat = subset(REF.srat,subset = Genotype %in% c('T21','diploid'))
# # 
# # ##----------------------------##
# # ##   Run Logistic Regression  ##
# # ##----------------------------##
# # skipIfExists=F
# # ref_annot='ann'
# # maxCells=4000
# # outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fBoneMarrowREF_onMLDS/')
# # model_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/2_LRorig_fBoneMarrowREF_onMLDS/fBM_diploid_T21_REF_trainModel_4kmaxcells_70perc_jan23.RDS'
# # 
# # LR_level='both'
# # srat_annot='ann'
# # minGeneMatch = 0.99
# # maxCells=4000
# # tissue = 'BoneMarrow'
# # 
# # out_prefix = 'fLiver_fBoneMarrowREF_maxCells_70perc_'
# # plot_prefix = 'fLiver_fBoneMarrowREF_maxCells_70perc_'
# # plot_prefix = NULL
# # 
# # outputs = runLR(REF.srat,ref_annot=ref_annot,srat=tgt.srat,LR_level=LR_level,srat_annot=srat_annot,
# #                 model_fp = model_fp,outDir=outDir,plot_prefix=plot_prefix,out_prefix=out_prefix,
# #                 minGeneMatch=minGeneMatch,maxCells = maxCells,tissue=tissue,
# #                 scLR_TGTtype='',scLR_REFtype='')
# # 
# # message(sprintf('3. LR completed for tissue %s',tissue))
# # 
# # ### Plotting ClusterLR heatmap #
# # output = outputs[[2]][[1]]
# # 
# # type = ifelse(grepl('ref_',rownames(output)),'REF','TGT')
# # show_row_names = T
# # plot_prefix = 'fLiver_fBoneMarrowREF_maxCells_70perc_'
# # 
# # pdf(file.path(outDir,paste0(plot_prefix,'_clusterLR','.pdf')),width = 25,height = 70)
# # 
# # hm = similarityHeatmap(output,
# #                        row_title_rot = 0,
# #                        row_title_gp = gpar(fontsize=12),row_names_gp = gpar(fontsize=10),row_names_max_width = unit(6,'cm'),
# #                        column_names_gp = gpar(fontsize=12),column_names_max_height = unit(6,'cm'),
# #                        split = type, gap = unit(2,'mm'), show_row_names = show_row_names, cluster_rows = T)
# # draw(hm)
# # dev.off()
# # 
# # 
