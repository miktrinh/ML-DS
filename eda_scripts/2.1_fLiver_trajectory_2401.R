# AIM: to construct MK-trajectory for 2n and T21, then use G2G to compare the 2 trajectories.

##------------------##
##    Libraries   ####
##------------------##
library(Seurat)
library(tidyverse)
source('~/lustre_mt22/generalScripts/utils/sc_utils.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')
source('~/lustre_mt22/generalScripts/utils/pseudobulk.R')
source('~/lustre_mt22/na15/trajectoryProjection_utils.R')

outDir = '~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output'

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}



setwd(outDir)





##-----------------------------------##
##    Convert RDS to h5ad object   ####
##-----------------------------------##
ANNOTATION_KEY = 'annot_mar24'
tissue = 'liver'

## Import fLiver object
# REF.srat = filter_sratObj(tissue = tissue,ageMatch=ageMatch,remove_cyclingCells=remove_cyclingCells,
#                                      srat_in_fp='~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_annotated_0424.RDS')
# REF.srat = REF.srat[[1]]
# mdat = fLiver@meta.data

fLiver = readRDS(ref.srat_fp)
mdat = read.csv('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_annotated_0424_mdat.csv')
fLiver$broadLineage = mdat$broadLineage[match(fLiver$cellID,mdat$cellID)] 
fLiver$finalAnn_broad = as.character(fLiver@meta.data[[ANNOTATION_KEY]])
fLiver$cellID = rownames(fLiver@meta.data)



## Define cellID to keep
ct_toKeep = c('HSC_MPP',
              "MEMP_MEP","earlyMK","MK","Mast.cell",
              "EE", "ME","LE",
              "CMP_GMP","promyelocyte","myelocyte",
              "proMono","Monocyte","Macrophage",
              "LMPP_ELP","pro.B.cell","pre.B.cell","B.cell")
mdat$cellID_toKeep = ifelse(as.character(mdat$finalAnn_broad %in% ct_toKeep),T,F)

## Down-sample some cell types
ct_toDS = table(as.character(mdat$finalAnn_broad[mdat$finalAnn_broad %in% ct_toKeep]))
ct_toDS = names(ct_toDS[ct_toDS>10000])

for(ct in ct_toDS){
  set.seed(1234)
  cellID_toRemove = mdat[mdat$finalAnn_broad == ct,] %>% group_by(donorID) %>% mutate(id=1:n(),nCell = n())
  cellID_toRemove = cellID_toRemove[cellID_toRemove$nCell >= 2000,]
  cellID_toRemove = cellID_toRemove %>% group_by(donorID) %>% mutate(selected = ifelse(id %in% sample(1:n(),2000),T,F))
  cellID_toRemove = cellID_toRemove$cellID[cellID_toRemove$selected == F]
  mdat$cellID_toKeep[mdat$cellID %in% cellID_toRemove] = F
}

table(mdat$finalAnn_broad[mdat$cellID_toKeep==T],mdat$donorID[mdat$cellID_toKeep==T])
write.csv(mdat,file.path(outDir,'liver_liverREFmerged_clean_processed_annotated_noUnknowns_0124_mdat_forPalantir.csv'))

colnames(mdat)[!colnames(mdat) %in% colnames(fLiver@meta.data)]
fLiver$cellID_toKeep = mdat$cellID_toKeep[match(fLiver$cellID,mdat$cellID)]
fLiver$UMAP_1 = fLiver@reductions$umap@cell.embeddings[,1]
fLiver$UMAP_2 = fLiver@reductions$umap@cell.embeddings[,2]

## Convert to h5ad for Palantir (in python)
#devtools::install_github("cellgeni/sceasy")
library(sceasy)
use_condaenv('base',required = T)
# Converting Raw counts to adata
convertFormat(fLiver, from="seurat", to="anndata",assay = "RNA", main_layer = "counts",
              outFile='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_annotated_0424.h5ad')



##--------------------------------##
##  Investigate 2n trajectory   ####
##--------------------------------##
nonHaem_celltypes = c('Hepatocyte','Endo','Kupffer.cell','Fibroblast','NPC','Mesenchyme')

celltypes_toKeep = c('HSC_MPP','MEMP_MEP','EE','ME','Mast.cell','earlyMK','MK')
                     #'LMPP_ELP','pro.B.cell','pre.B.cell','B.cell',
                     #'CMP_GMP','proMono','MOP','Monocyte')
diploid = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$Genotype == 'diploid' & REF.srat$annot_jan24 %in% celltypes_toKeep & REF.srat$published_ann_2 == 'NA'])

#diploid = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$Phase == 'G1' & REF.srat$Genotype == 'diploid' & !REF.srat$annot_sept23 %in% nonHaem_celltypes])
diploid = standard_clustering(diploid)

DimPlot(diploid,group.by = 'annot_jan24',label = T,repel = T,cols = col25,label.box = T) + NoLegend()

# Define start cell
start_cell = names(diploid@assays$RNA@counts['CD34',])[diploid@assays$RNA@counts['CD34',] == max(diploid@assays$RNA@counts['CD34',])]


## Import palantir trajectory
diploid_trajectory = get_palantir_ref_traj(palantir_output_fp='~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/inhouse2n_haemTraj.csv',
                                           nTerm=5,split_ref=FALSE,
                                           test_ratio=0.2,ann_key='annot_jan24')

table(diploid_trajectory$branch2,diploid_trajectory$annot_jan24)

## Define cells to include along MK trajectory
# 1. HSC_MPP or MEMP_MEP or MK with entropy > 1 and belongs to MK/Mast/EE
# 2. HSC_MPP with the difference between branchProb_Monocyte and branchProb_MK/ME/Mast < 0.2
DimPlot(diploid,cells.highlight = diploid_trajectory$cellID[diploid_trajectory$annot_jan24 %in%  c('MK') &
                                                              #diploid_trajectory$palantir_entropy > 1 &
                                                              diploid_trajectory$branch2 %in% c('EE')])
diploid_trajectory$branch3 = diploid_trajectory$branch2
diploid_trajectory$branch3[diploid_trajectory$annot_jan24 %in%  c('HSC_MPP') &
                             diploid_trajectory$palantir_entropy > 1 &
                             diploid_trajectory$branch2 %in% c('EE','Mast.cell')] = 'MK'

diploid_trajectory$branch3[diploid_trajectory$annot_jan24 %in%  c('earlyMK','MK')] = 'MK'
diploid_trajectory$branch3[diploid_trajectory$annot_jan24 %in%  c('Mast.cell')] = 'Mast.cell'
diploid_trajectory$branch3[diploid_trajectory$annot_jan24 %in%  c('EE') & diploid_trajectory$branch3 == 'MK'] = 'EE'
diploid_trajectory$branch3[diploid_trajectory$annot_jan24 %in%  c('pre.B.cell') & diploid_trajectory$branch3 == 'MK'] = 'B.cell'

table(diploid_trajectory$branch3,diploid_trajectory$annot_jan24)


write.csv(diploid_trajectory,'~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/inhouse2n_haemTraj_assigned.csv')
diploid_trajectory = read.csv('~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/inhouse2n_haemTraj_assigned.csv')






##--------------------------------##
##  Investigate T21 trajectory   ####
##--------------------------------##
t21 = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$Genotype == 'T21' & REF.srat$annot_jan24 %in% celltypes_toKeep])
t21 = standard_clustering(t21)
DimPlot(t21,group.by = 'annot_jan24',label = T,repel = T,label.box = T) + NoLegend()

## Import palantir trajectory
t21_trajectory = get_palantir_ref_traj(palantir_output_fp='~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/T21_haemTraj_noB.csv',
                                       nTerm=4,split_ref=FALSE,
                                       test_ratio=0.2,ann_key='annot_jan24')

ggplot(t21_trajectory,aes(palantir_pseudotime,palantir_entropy,col=annot_jan24))+
  geom_point(size=0.8,alpha=0.9)+
  scale_colour_manual(values = col25[-6])+
  theme_bw(base_size = 15) + ggtitle('Palantir trajectory')



table(t21_trajectory$branch2,t21_trajectory$annot_jan24)

## Define cells to include along MK trajectory
DimPlot(t21,cells.highlight = t21_trajectory$cellID[t21_trajectory$annot_jan24 %in%  c('HSC_MPP') &
                                                      t21_trajectory$diff_branchProb < 0.2 &
                                                      t21_trajectory$branch2 %in% c('Monocyte')])

t21$branch2 = t21_trajectory$branch2[match(t21$cellID,t21_trajectory$cellID)]
t21$branch2[is.na(t21$branch2)] = 'NA'
Idents(t21) = t21$annot_jan24
DotPlot(t21,idents = 'MEMP_MEP',group.by = 'branch2',features = unique(c(
  'CD34','SPINK2',
  'CD38',#'MLLT3','PRSS57', # HSC_MPP
  #'CTSG',	'PRTN3', # CMP
  
  
  #'SERPINB1', 
  'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1',#MEMP
  'HBD','KLF1','FLI1','FCER1A', 'ITGA2B', 'PLEK', # MEP
  
  'PF4','ITGA2B', #Megakaryocyte
  'PPBP','TUBB1', # Platelets
  'CSF2RB','HDC','TPSAB1','KIT', # Mast.cell
  
  'GATA1','KLF1', # early Erythroid
  'ALAS2', # mid.erythroid
  'HBB','BPGM', # late.erythroid
  
  
  'CD14',#'FCGR3A',# Monocytes
  'C1QA','CD68',#'MSR1',#Macrophages
  'FABP3',#Kupffer
  #'CCR2','CX3CR1',
  
  'ITGAX','IRF8',	 #DC.precursor
  'CLEC10A','CD1C','CLEC9A',#DC1
  #'CLEC10A', # DC2 
  'CLEC4C','IL3RA', #pDC
  
  'AZU1','MPO',
  #'CSF3R','FPR1','FCGR3B',
  'MNDA', # NEUTROPHILS
  'DEFA3','DEFA4', # pro-myelocytes
  'CAMP','LCN2', #myelocytes
  #'CXCR2',
  'CSF3R','FCGR3B',#'FUT4', # Neutrophil
  
  #'IGLL1','CD99', # Pre-pro
  'DNTT',
  #'EBF1',
  'CD19',#'RAG1',# pro-b-cells
  'MME','VPREB1','CD79A','CD79B',# pre-b-cells
  'TCL1A','MME','RAG1',
  'MS4A1',  # B-cells
  #'CD27',#plasma cell?
  
  'CD52','IL7R',# ILC precursor
  'ZBTB16',#'LTB', 
  'CD3D',#'GZMA',
  #'CD4','CD8A', #Early.lymphoid_T.lymphocyte
  #'TRDV2','TRGV9', # gamma delta T-cell
  #'SLC4A10','TRAV1-2', #MAIT t-cell
  #'PRF1', # effector T cells
  #'FOXP3',	'CDH1', # regulatory T cells
  'NKG7','KLRD1', #NK
  
  
  #'ESAM','PECAM1', 
  'KDR', #'PTPRB', 
  #'PLVAP', # Endothelium
  #'DCN','SERPINF1',
  'COL3A1',#'COL1A1',
  #'BGN','VIM','ECM1', # endosteal fibroblast / fibroblast
  #'APOA1',
  #'SCD',
  'ALB','TTR', # Hepatocyte
  'PTPRC'
  
))) + RotatedAxis()



t21_trajectory$branch3 = t21_trajectory$branch2
t21_trajectory$max_branchProb_MK.EE.Mast = sapply(seq(1:nrow(t21_trajectory)),FUN = function(i){
  return(max(t21_trajectory[i,c('branchProb_Mast.cell','branchProb_EE','branchProb_MK')]))})
t21_trajectory$diff_branchProb = t21_trajectory$branchProb_Monocyte - t21_trajectory$max_branchProb_MK.EE.Mast


t21_trajectory$branch3[t21_trajectory$annot_jan24 %in%  c('HSC_MPP') &
                             t21_trajectory$palantir_entropy > 1 &
                             t21_trajectory$branch2 %in% c('EE','Mast.cell')] = 'MK'
t21_trajectory$branch3[t21_trajectory$annot_jan24 %in%  c('HSC_MPP') &
                         t21_trajectory$palantir_entropy > 1 &
                         t21_trajectory$diff_branchProb < 0.2 &
                         t21_trajectory$branch2 %in% c('Monocyte')] = 'MK'


t21_trajectory$branch3[t21_trajectory$annot_jan24 %in%  c('earlyMK')] = 'MK'
t21_trajectory$branch3[t21_trajectory$annot_jan24 %in%  c('Mast.cell')] = 'Mast.cell'
t21_trajectory$branch3[t21_trajectory$annot_jan24 %in%  c('EE') & t21_trajectory$branch3 == 'MK'] = 'EE'
t21_trajectory$branch3[t21_trajectory$annot_jan24 %in%  c('CMP_GMP','Monocyte') & t21_trajectory$branch3 == 'MK'] = 'Monocyte'
t21_trajectory$branch3[t21_trajectory$annot_jan24 %in%  c('pro.B.cell','B.cell','LMPP_ELP') & t21_trajectory$branch3 == 'MK'] = 'B.cell'

table(t21_trajectory$branch3,t21_trajectory$annot_jan24)


write.csv(t21_trajectory,'~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/T21_haemTraj_noB_assigned.csv')
t21_trajectory = read.csv('~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/T21_haemTraj_noB_assigned.csv')




##-------------------------------------##
##  Investigate allGeno trajectory   ####
##-------------------------------------##
nonHaem_celltypes = c('Hepatocyte','Endo','Kupffer.cell','Fibroblast','NPC','Mesenchyme')

celltypes_toKeep = c('HSC_MPP','MEMP_MEP','EE','ME','Mast.cell','earlyMK','MK')
#'LMPP_ELP','pro.B.cell','pre.B.cell','B.cell',
#'CMP_GMP','proMono','MOP','Monocyte')
diploid = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$Genotype == 'diploid' & REF.srat$annot_jan24 %in% celltypes_toKeep & REF.srat$published_ann_2 == 'NA'])

#diploid = subset(REF.srat,subset = cellID %in% REF.srat$cellID[REF.srat$Phase == 'G1' & REF.srat$Genotype == 'diploid' & !REF.srat$annot_sept23 %in% nonHaem_celltypes])
diploid = standard_clustering(diploid)

DimPlot(diploid,group.by = 'annot_jan24',label = T,repel = T,cols = col25,label.box = T) + NoLegend()

# Define start cell
start_cell = names(diploid@assays$RNA@counts['CD34',])[diploid@assays$RNA@counts['CD34',] == max(diploid@assays$RNA@counts['CD34',])]


## Import palantir trajectory
allGeno_trajectory = get_palantir_ref_traj(palantir_output_fp='~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/allGeno_haemTraj_noPublished2n_mar24.csv',
                                           nTerm=5,split_ref=FALSE,
                                           test_ratio=0.2,ann_key='annot_jan24')

table(allGeno_trajectory$branch2,allGeno_trajectory$annot_jan24)
ggplot(allGeno_trajectory[allGeno_trajectory$finalAnn_broad == 'HSC_MPP',],aes(palantir_pseudotime,palantir_entropy,col=finalAnn_broad))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = col25)

# ## Define cells to include along MK trajectory
# # 1. HSC_MPP or MEMP_MEP or MK with entropy > 1 and belongs to MK/Mast/EE
# # 2. HSC_MPP with the difference between branchProb_Monocyte and branchProb_MK/ME/Mast < 0.2
# DimPlot(diploid,cells.highlight = allGeno_trajectory$cellID[allGeno_trajectory$annot_jan24 %in%  c('MK') &
#                                                               #allGeno_trajectory$palantir_entropy > 1 &
#                                                               allGeno_trajectory$branch2 %in% c('EE')])
# allGeno_trajectory$branch3 = allGeno_trajectory$branch2
# allGeno_trajectory$branch3[allGeno_trajectory$annot_jan24 %in%  c('HSC_MPP') &
#                              allGeno_trajectory$palantir_entropy > 1 &
#                              allGeno_trajectory$branch2 %in% c('EE','Mast.cell')] = 'MK'
# 
# allGeno_trajectory$branch3[allGeno_trajectory$annot_jan24 %in%  c('earlyMK','MK')] = 'MK'
# allGeno_trajectory$branch3[allGeno_trajectory$annot_jan24 %in%  c('Mast.cell')] = 'Mast.cell'
# allGeno_trajectory$branch3[allGeno_trajectory$annot_jan24 %in%  c('EE') & allGeno_trajectory$branch3 == 'MK'] = 'EE'
# allGeno_trajectory$branch3[allGeno_trajectory$annot_jan24 %in%  c('pre.B.cell') & allGeno_trajectory$branch3 == 'MK'] = 'B.cell'
# 
# table(allGeno_trajectory$branch3,allGeno_trajectory$annot_jan24)
# 
# 
# write.csv(allGeno_trajectory,'~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/inhouse2n_haemTraj_assigned.csv')
# allGeno_trajectory = read.csv('~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/inhouse2n_haemTraj_assigned.csv')



## Cell Rank fate probability
fateProb = read.csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/fLiver_2nAK_nonDS_palantired_CR_fateProbs.csv')

