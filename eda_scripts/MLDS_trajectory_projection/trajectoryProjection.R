### Projecting query cells onto a Reference trajectory  ###
###         TAM/MLDS onto fLiver MK/Ery/Mast           ###



outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/jul24/wholeTranscriptome'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

# Also make 3 sub_folders for 1_palantir_output; 2_tradeSeq_output; 3_projection_output
subDirs = file.path(outDir,c('1_palantir_output','2_tradeSeq_output','3_projection_output'))
for(subDir in subDirs){
  if(!dir.exists(subDir)){
    dir.create(subDir,recursive = T)  
  }
}

setwd(outDir)



#----------------------------##
#          Libraries       ####
#----------------------------##
library(tidyverse)
library(Seurat)
source('~/lustre_mt22/Aneuploidy/scripts/projection_method/trajectoryProjection_utils.R')


##----------------------------##
##   Set Global parameters  ####
##----------------------------##


# Path to reference seurat object
ref.srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'


# Path to save h5ad reference object for palantir (must be full path, no "~")
ref.srat_h5ad_fp = '/home/jovyan/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0724.h5ad'


# Path to Palantir output csv file
palantir_output_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/6_MK_lineage/jan24/fLiver2n_allcells/1_palantir_output/inhouse2n_haemTraj.csv'
palantir_output_fp = '~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/inhouse2n_MK.EE.Mast.Traj.csv'
palantir_output_fp = '~/lustre_mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/inhouse2n_MK.LE.Mast.Traj_2407.csv'
palantir_output_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/jul24/wholeTranscriptome/1_palantir_output/inhouse2n_haemTraj_2407.csv'
palantir_output_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/jul24/wholeTranscriptome/1_palantir_output/inhouse2n_MK.LE.Mast.Traj_2407.csv'
# Path to seurat object containing query cells
tgt.srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS'

# Number of terminal states
nTerm = 3
ANNOTATION_KEY = 'annot_aug24'
#tradeSeq_outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/6_MK_lineage/jan24/fLiver2n_allcells/2_tradeSeq_output'
tradeSeq_outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/6_MK_lineage/sept23/fLiver2n_allCells/2_tradeSeq_output'

skipIfExist = T
parallel = T




#------------------------------------------------------------##
##    Step 1: Build REFERENCE trajectory using Palantir    ####
#------------------------------------------------------------##

## Import REF dataset
fLiver = readRDS(ref.srat_fp)
ANNOTATION_KEY = 'annot_aug24'
fLiver$finalAnn_broad = as.character(fLiver@meta.data[[ANNOTATION_KEY]])
fLiver$cellID = rownames(fLiver@meta.data)

# Subset REF.srat to include cell of interest only
REF.srat = subset(fLiver, subset = cellID %in% fLiver$cellID[fLiver$Genotype == 'diploid' & 
                                                               fLiver$finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','earlyMK','MK',
                                                                                            'EE','ME','LE','Mast.cell'#,'Monocyte','CMP_GMP','proMono'
                                                                                            #'B.cell','Monocyte','CMP_GMP','pre.B.cell'
)])




## Convert to h5ad for Palantir (in python)
#devtools::install_github("cellgeni/sceasy")
# library(sceasy)
# use_condaenv('base',required = F,conda = '/opt/conda/etc/profile.d/conda.sh')
# # Converting Raw counts to adata
# convertFormat(REF.srat, from="seurat", to="anndata",assay = "RNA", main_layer = "counts",
#              outFile=ref.srat_h5ad_fp)

## Run Palantir in python
## Down-sampling EE/ME/LE for diploid and T21 samples
set.seed(4321)
fLiver$cells_toKeep = T
for(d in unique(fLiver$donorID)){
  nCells = table(fLiver$finalAnn_broad[fLiver$finalAnn_broad %in% c('EE','ME','LE') & fLiver$donorID == d])
  if(any(nCells > 5000)){
    print(unique(fLiver$Genotype[fLiver$donorID == d]))
    for(ct in names(nCells[nCells > 5000])){
      cells_toKeep = sample(fLiver$cellID[fLiver$finalAnn_broad == ct & fLiver$donorID == d],5000)
      fLiver$cells_toKeep[fLiver$finalAnn_broad == ct & fLiver$donorID == d & !fLiver$cellID %in% cells_toKeep] = F
    }
  }
}

write.csv(fLiver@meta.data,file.path(outDir,'1_palantir_output','fLiver_cells_downsampled.csv'))



## Import Palantir output
print('Importing Palantir results')
ANNOTATION_KEY = 'annot_jun24'
ref_trajectory = get_palantir_ref_traj(palantir_output_fp=palantir_output_fp,
                                       nTerm=nTerm,split_ref=TRUE,
                                       test_ratio=0.2,ann_key=ANNOTATION_KEY,cell_key = 'cellID')
ref_trajectory$annot = fLiver$annot[match(ref_trajectory$cellID,fLiver$cellID)]
table(ref_trajectory$branch2,ref_trajectory$annot)
# dd = ref_trajectory
# dd$tmp = ifelse(dd$annot == 'EE' & dd$branch2 == 'MK','1',dd$annot)
# ggplot(dd,aes(palantir_pseudotime,palantir_entropy,col=tmp))+
#   geom_point(size=0.2)+
#   scale_color_manual(values = c(col25[1],rep(grey(0.8),10)))+
#   theme_classic()

ref_trajectory$branch2[ref_trajectory$annot %in% c('EE') & ref_trajectory$branch2 == 'MK'] = 'LE'
ref_trajectory$branch2[ref_trajectory$annot %in% c('earlyMK','MK')] = 'MK'
ref_trajectory$branch2[ref_trajectory$annot %in% c('Mast.cell')] = 'Mast.cell'
ref_trajectory$branch2[ref_trajectory$annot %in% c('B.cell','LMPP_ELP','pre.B.cell','pro.B.cell')] = 'B.cell'
ref_trajectory$branch2[ref_trajectory$annot %in% c('proMono','Monocyte','CMP_GMP')] = 'Monocyte'
ref_trajectory$branch2[ref_trajectory$annot %in% c('LE','ME')] = 'LE'

table(ref_trajectory$branch2,ref_trajectory$annot)
table(ref_trajectory$cellID %in% REF.srat$cellID)
ref_trajectory = ref_trajectory[ref_trajectory$cellID %in% REF.srat$cellID,]
##----- Annotate branch structure
table(ref_trajectory$finalAnn_broad[ref_trajectory$palantir_pseudotime <= 0.4])
# cells before 0.5 are common progenitor to all 3 branches
ref_trajectory$differentiation_state = '?'
ref_trajectory$differentiation_state[ref_trajectory$palantir_pseudotime <= 0.4] = 'prog'
ref_trajectory$differentiation_state[ref_trajectory$palantir_pseudotime > 0.4] = ref_trajectory$branch2[ref_trajectory$palantir_pseudotime > 0.4]

ggplot(ref_trajectory,aes(palantir_pseudotime,palantir_entropy,col=annot))+
  geom_point(size=0.01)+
  #facet_wrap(vars(annot))+
  geom_vline(xintercept = 0.4)+
  scale_color_manual(values = col25)+
  theme_classic(base_size = 13)+
  theme(panel.border = element_rect(fill=F,color='black',linewidth = 1),axis.line = element_blank()) 


ggplot(ref_trajectory,aes(palantir_pseudotime,nFeature_RNA,col=differentiation_state))+
  geom_point(size=0.01)+
  theme_classic()+
  theme(panel.border = element_rect(fill=F,color='black',linewidth = 1),axis.line = element_blank()) 



# ref_trajectory_sept23 = get_palantir_ref_traj(palantir_output_fp=palantir_output_fp,
#                                        nTerm=nTerm,split_ref=TRUE,
#                                        test_ratio=0.2,ann_key=ANNOTATION_KEY,cell_key = 'cellID')
# ref_trajectory_sept23$Phase = mdat$Phase[match(ref_trajectory_sept23$cellID,mdat$cellID)]
# table(ref_trajectory_sept23$Phase)
# ref_trajectory_sept23$pseudotime_nov23 = ref_trajectory$palantir_pseudotime[match(ref_trajectory_sept23$cellID,ref_trajectory$cellID)]
# ggplot(ref_trajectory_sept23,aes(palantir_pseudotime,pseudotime_nov23,col=finalAnn_broad))+
#   geom_point()+
#   geom_abline()+
#   theme_bw()


# Scaling pseudotime
if(max(ref_trajectory$palantir_pseudotime) > 1){
  ref_trajectory$palantir_pseudotime = (ref_trajectory$palantir_pseudotime - min(ref_trajectory$palantir_pseudotime))/(max(ref_trajectory$palantir_pseudotime) - min(ref_trajectory$palantir_pseudotime))
}

ggplot(ref_trajectory,aes(palantir_pseudotime,palantir_entropy,col=differentiation_state))+
  geom_point()

# Scaling individual branches to stretch between branch beginning to the end
for(b in unique(ref_trajectory$branch2)){
  d = ref_trajectory[ref_trajectory$branch2 == b & ref_trajectory$differentiation_state != 'prog',]
  # ggplot(d,aes(palantir_pseudotime,palantir_entropy,col=differentiation_state))+
  #   geom_point()
  
  d$palantir_pseudotime = (d$palantir_pseudotime - min(d$palantir_pseudotime))/(max(d$palantir_pseudotime) - min(d$palantir_pseudotime))*(1-min(d$palantir_pseudotime)) + min(d$palantir_pseudotime)
  ref_trajectory$palantir_pseudotime[match(d$cellID,ref_trajectory$cellID)] = d$palantir_pseudotime
}

ggplot(ref_trajectory,aes(palantir_pseudotime,palantir_entropy,col=differentiation_state))+
  geom_point(size=0.01,alpha=0.2)+
  scale_color_manual(values = col25)+
  facet_wrap(vars(differentiation_state))+
  ggtitle('REF trajectory post branch-wise scaling')+
  theme_classic()



# ref_trajectory2 = binning_trajectory(ref_trajectory[ref_trajectory$branch2=='MK',],mergeBin = F)
# ref_trajectory$finalAnn_broad2 = as.character(ref_trajectory$finalAnn_broad)
# ref_trajectory$finalAnn_broad2 = factor(ref_trajectory$finalAnn_broad2,levels = c('HSC_MPP','MEMP_MEP','Mast.cell','MK','EE'))
# ref_trajectory$palantir_pseudotime_round = round(ref_trajectory$palantir_pseudotime,1.5)
# df = ref_trajectory %>% group_by(branch2,palantir_pseudotime_round,finalAnn_broad) %>% summarise(nCell = n())
# df$finalAnn_broad2 = factor(df$finalAnn_broad,levels = c('HSC_MPP','MEMP_MEP','Mast.cell','MK','EE'))
# table(ref_trajectory$finalAnn_broad2,ref_trajectory$branch2)
# 
# ggplot(ref_trajectory[ref_trajectory$branch2 == 'MK',],aes(palantir_pseudotime))+
#   #geom_jitter(height = 3,size=0.2,aes(y=3,col=finalAnn_broad2),alpha=0.7)+
#   #geom_density(aes(palantir_pseudotime,fill=finalAnn_broad2),alpha=0.7)+
#   geom_col(data = df[df$branch2 == 'MK' & df$finalAnn_broad2 != 'Mast.cell',],aes(palantir_pseudotime_round,log10(nCell),fill = finalAnn_broad2),alpha=0.6)+
#   
#   #facet_wrap(vars(branch2),ncol=1,scales = 'free')+
#   scale_color_manual(values = c(colAlpha('#EF6548',0.4),'#990000','#A0D7E2',colAlpha('#D360A3',0.6),colAlpha('#8C96C6',0.6)))+
#   scale_fill_manual(values = c(colAlpha('#EF6548',0.9),'#990000','#A0D7E2',colAlpha('#D360A3',0.9),colAlpha('#8C96C6',0.9)))+
#   theme_classic(base_size = 15) + xlab('pseudotime') + ylab('log10 cell count')
#   
# 
# 
# 
# print(p)
# 
# df$finalAnn_broad2 = factor(df$finalAnn_broad2,levels = c('MK','MEMP_MEP','HSC_MPP'))
# plotFun = function(noFrame=FALSE,noPlot=FALSE){
#   
#   p = ggplot(df[df$branch2 == 'MK' & df$finalAnn_broad2 != 'Mast.cell',],aes(palantir_pseudotime_round,log10(nCell),fill = finalAnn_broad2))+
#     geom_col(alpha=0.6)+
#     scale_color_manual(values = c(colAlpha('#EF6548',0.4),'#990000','#A0D7E2',colAlpha('#D360A3',0.6),colAlpha('#8C96C6',0.6)))+
#     scale_fill_manual(values = c(colAlpha('#D360A3',0.9),'#990000',colAlpha('#EF6548',0.9),colAlpha('#8C96C6',0.9)),name='cell type')+
#     theme_classic(base_size = 11) + xlab('pseudotime') + ylab('log10 cell count')
#   
#   print(p)
# }
# 
# saveFig(file.path(plotDir,'example_trajectory'),plotFun,rawData=df,width = 4.8,height = 1.8,res = 500)


# Define training set
ref_trajectory.train = ref_trajectory[ref_trajectory$set == 'train',]


# REF.srat.sub = subset(REF.srat,subset = cellID %in% ref_trajectory.train$cellID)
# REF.srat.sub = standard_clustering(REF.srat.sub)
# REF.srat.sub$pseudotime = ref_trajectory.train$palantir_pseudotime[match(REF.srat.sub$cellID,ref_trajectory.train$cellID)]
# REF.srat.sub$entropy = ref_trajectory.train$palantir_entropy[match(REF.srat.sub$cellID,ref_trajectory.train$cellID)]
# REF.srat.sub$UMAP_1 = REF.srat.sub@reductions$umap@cell.embeddings[,'UMAP_1']
# REF.srat.sub$UMAP_2 = REF.srat.sub@reductions$umap@cell.embeddings[,'UMAP_2']
# 
# FeaturePlot(REF.srat.sub,'pseudotime')
# DimPlot(REF.srat.sub,group.by = 'finalAnn_broad',label = T,label.box = T,repel = T)+ NoLegend()
# assoRes.sub = assoRes.sub[order(assoRes.sub$waldStat,decreasing = T),]
# ## Expression by group
# res = avgExpr_cellFrac_byGroup(srat = REF.srat.sub,group = 'finalAnn_broad',genes = assoRes.sub$gene,doPlot = F)
# cellfrac_df = res[["cellFrac_df"]]
# avgExpr = res[["avgExpr_df"]]
# avgExpr_long = pivot_wider(avgExpr,id_cols = 1,names_from = 'group',values_from = 'avgExpr')
# avgExpr_long$maxBranch = sapply(1:nrow(avgExpr_long),function(i){
#   n = colnames(avgExpr_long)[c(which(avgExpr_long[i,-1] == max(avgExpr_long[i,-1])) + 1)]
#   return(n)
# })
# 
# DotPlot(l076,group.by = 'tissue',features = avgExpr_long$gene[avgExpr_long$maxBranch == 'Mast.cell'][211:300]) +RotatedAxis()
# 
# cellfrac_df = merge(cellfrac_df,avgExpr,by=c('gene','group'))
# mast.markers = cellfrac_df[cellfrac_df$group == 'Mast.cell' & cellfrac_df$frac > 0.8,]
# DotPlot(REF.srat.sub,group.by = 'finalAnn_broad',features = assoRes.sub$gene[])
# umap = as.data.frame(REF.srat@reductions$umap@cell.embeddings)
# ref_trajectory$umap1 = umap$UMAP_1[match(ref_trajectory$cellID,rownames(umap))]
# ref_trajectory$umap2 = umap$UMAP_2[match(ref_trajectory$cellID,rownames(umap))]
# 
# ggplot(ref_trajectory,aes(umap1,umap2,col=palantir_pseudotime))+
#   geom_point(size=0.2)+
#   theme_bw()






#--------------------------------##
##    Step 2: Run TradeSeq     ####
#--------------------------------##
# print('Running TradeSeq')
# # TradeSeq identify genes associated with the trajectory 
# # i.e. genes with significant changes in expression pattern along the trajectory
# #ref_trajectory.train$cellID = gsub('-','.',ref_trajectory.train$cellID)
# #REF.srat$cellID = gsub('\\.','-',REF.srat$cellID)
# 
# assoRes = run_tradeSeq(REF.srat=REF.srat,ref_trajectory=ref_trajectory.train,
#                        outDir=tradeSeq_outDir,skipIfExist=T,
#                        term_states=NULL,nTerm=NULL,parallel = parallel)
# 
# ## Remove genes with low waldStat - keep top 75% only
# assoRes$gene = assoRes$X
# class(assoRes)
# ## Remove rubbish genes from list of associated genes
# assoRes = assoRes[!is.na(assoRes$pvalue) & !grepl('^RPL|^RPS|^MT-',assoRes$gene),]
# assoRes$padj = p.adjust(assoRes$pvalue,method = 'BH')
# assoRes = assoRes[assoRes$padj <0.01,]
# assoRes = assoRes[!assoRes$gene %in% c('MALAT1','NEAT1'),]
# assoRes.sub = assoRes[assoRes$waldStat > quantile(assoRes$waldStat,0.5),]
# 
# REF.srat$group = paste0(REF.srat$annot,"_",REF.srat$donorID)
# avgExpr = AverageExpression(REF.srat,group.by = 'group',features = assoRes$gene)
# avgExpr = avgExpr[['RNA']]
# 
# hm = Heatmap(t(scale(t(avgExpr))),show_row_dend = F,show_column_dend = F,
#              column_split = REF.srat$annot[match(colnames(avgExpr),REF.srat$group)],km=4,
#              cluster_rows = T,cluster_columns = T,show_row_names = F,show_column_names = T)
# ht = draw(hm)
# 
# geneOrder = data.frame(geneSym = rownames(avgExpr)[do.call(c,row_order(ht))],
#                        group = rep(c(names(row_order(ht))),sapply(row_order(ht),length)))
# geneOrder$celltype = ifelse(geneOrder$group == 1,'HSC',
#                             ifelse(geneOrder$group == 2,'Mast',
#                                    ifelse(geneOrder$group == 3,'MK',
#                                           ifelse(geneOrder$group == 4,'EE','others'))))
# 
# View(assoRes[assoRes$gene %in% geneOrder$geneSym[geneOrder$celltype == 'Mast'],])
# mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24_mdat.csv')
# mlds$broadLineage = mdat$broadLineage[match(mlds$cellID,mdat$cellID)]
# mlds$group = ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic',paste0('Tumour_D_',mlds$donorID),
#                     ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'TP1' & mlds$donorID == 'L038','Tumour_TP1_L038',
#                            ifelse(mlds$annot == 'Mast.cell','Mast',mlds$broadLineage)))
# avgExpr_mlds = AverageExpression(mlds,features = assoRes$gene,group.by = 'group')
# avgExpr_mlds = avgExpr_mlds[['RNA']]
# avgExpr_mlds = avgExpr_mlds[!is.na(rownames(avgExpr_mlds)),]
# 
# 
# hm2 = Heatmap(t(scale(t(avgExpr_mlds[,grepl('Tumour_',colnames(avgExpr_mlds))]))),show_row_dend = F,show_column_dend = F,
#               column_split = mlds$broadLineage[match(colnames(avgExpr_mlds)[grepl('Tumour_',colnames(avgExpr_mlds))],mlds$group)],#km=4,
#               column_title_rot = 90,
#               split = geneOrder$celltype[match(rownames(avgExpr_mlds),geneOrder$geneSym)],
#               cluster_rows = T,cluster_columns = T,show_row_names = F,show_column_names = T)
# ht2 = draw(hm2)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
##-------------------------------------------------##
##  Step 3. Prepare query cell count matrix      ####
##-------------------------------------------------##

#### Test set is made of:
# 1. fAdrenal SCPs + Megakaryocyte + .cells
# 2. fLiver 2n test set
# 3. fLiver AK (T21) MK.lineage
# 5. fBM 2n + T21 MK.lineage
# 4. Ellie's infant ALL - 8 diploid ALL + 1 AMKL
# 5. T21 ALL
# 6. MLDS

test_datasets = c('fAdr','fLiver2n','fLiverAK','fBM_2n','fBM_T21','MLDS','infantALL','pAML','BALL','MDS')
test_datasets = c('fLiver2n','MLDS','MDS')
cols_toKeep = c('cellID','donorID','Phase','nFeature_RNA','nCount_RNA','percent.mt','finalAnn_broad','Genotype','dataset','timePoint','tissue')

#test_datasets = c('fAdr','fLiver2n','MLDS')
inputData = import_datasets(test_datasets = test_datasets,
                            cols_toKeep = cols_toKeep,
                            ref_trajectory=ref_trajectory)

# tum = inputData[[2]][['MLDS']]
# inputData[[1]][['TAM']] = inputData[[1]][['MLDS']][,colnames(inputData[[1]][['MLDS']]) %in% tum$cellID[tum$finalAnn_broad == 'Tumour' & tum$dataset == 'TAM']]
# inputData[[1]][['MLDS_tum']] = inputData[[1]][['MLDS']][,colnames(inputData[[1]][['MLDS']]) %in% tum$cellID[tum$finalAnn_broad == 'Tumour' & tum$dataset == 'MLDS']]
# inputData[[1]][['MLDS_norm']] = inputData[[1]][['MLDS']][,colnames(inputData[[1]][['MLDS']]) %in% tum$cellID[tum$finalAnn_broad != 'Tumour']]
# 
# inputData[[2]][['TAM']] = tum[tum$finalAnn_broad == 'Tumour' & tum$dataset == 'TAM',]
# inputData[[2]][['MLDS_tum']] = tum[tum$finalAnn_broad == 'Tumour' & tum$dataset == 'MLDS',]
# inputData[[2]][['MLDS_norm']] = tum[tum$finalAnn_broad != 'Tumour',]
# 
# inputData[[1]] = inputData[[1]][names(inputData[[1]]) != 'MLDS']
# inputData[[2]] = inputData[[2]][names(inputData[[2]]) != 'MLDS']
# names(inputData[[1]])



##---------------------------------------------------------------##
##    Step 4. Project query cells onto reference trajectory    ####
##---------------------------------------------------------------##
nBin_toUse=50 # Number of bin to split the pseudotime by
refTraj='fLiver2n_MKlin_traj1200_240904_combinedBranches' 

geneSelectionMethod = 'allGenes'

useTradeseqGenes = F # Do we want to use tradeSeq_associated genes? if FALSE, all genes in test_toc will be used
tradeSeq_outDir=tradeSeq_outDir # Directory containing tradeSeq output

useHVG=F # Use only top_n highly variable genes to compute PCA space
n_HVG = 5000 # Number of top HVG to use


group = 'celltype' # name of column in tgt_mdat to inform how to split query cell output for plotting
nPCs = 20 # Number of PCs to be used in calculating correlation
branch_ofInterest = NULL # if only interested in a subset of terminal states, specify a vector of these terminal states (eg. c('B.cell'))
skipIfExist = T
doPlot=T
mergeBin=F
minCellperBin = 300
projection_outDir = file.path(outDir,'3_projection_output')
branchMethod = 'combined_branches'
dimRedMethod = NULL #'PCA'
binMethod = 'binByPseudotime' #'binByPseudotime.noDS'
compute_ED = F
print('Starting projection....')
# ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/6_MK_lineage/sept23/fLiver2n_allCells//fLiver2n_MKlin_traj1200_nov23.REF_wholeTrans_PCA_30PCs_50binByPseudotime.noDS_projection_output.RDS
# tgt_mdat = data.frame()
# test_toc = data.frame()
# REF.srat=c()
# save test_toc for future use

# #saveRDS(test_toc,file.path(outDir,'3_projection_output/fLiver2nAK_MLDS_AMKL_MDS_pAML.test_toc.RDS'))
# saveRDS(test_toc,file.path(outDir,'3_projection_output/fLiver2nAK_MLDS_MDS_fAdr.test_toc.RDS'))
# write.csv(tgt_mdat,file.path(outDir,'3_projection_outputfLiver2nAK_MLDS_MDS_fAdr.tgt_mdat_240406.csv'))
# results is a list of raw_output and processed_output
results = project_cells_v2(REF.srat=REF.srat,
                        ann_key='finalAnn_broad',
                        
                        ref_trajectory=ref_trajectory.train,
                        branch_ofInterest=branch_ofInterest,
                        
                        inputData=inputData,
                        # tgt_toc=test_toc,
                        # tgt_mdat=tgt_mdat,
                        group = group,
                        
                        
                        
                        geneSelectionMethod = geneSelectionMethod,
                        
                        useTradeseqGenes=useTradeseqGenes,
                        tradeSeq_outDir=tradeSeq_outDir,
                        assoRes = assoRes.sub,
                        
                        useHVG = useHVG,
                        n_HVG = n_HVG,
                        
                        
                        dimRedMethod=dimRedMethod,#"PCA",
                        nBin_toUse=nBin_toUse,
                        mergeBin=mergeBin,
                        branchMethod=branchMethod,
                        refTraj=refTraj,
                        minCellperBin = minCellperBin,
                        skipIfExist = skipIfExist,
                        nPCs = nPCs,
                        outDir=projection_outDir,
                        doPlot=doPlot,
                        parallel=F,
                        binMethod = binMethod,
                        compute_ED=compute_ED)

print('Cell projection completed!')



# out_list = results[['raw_output']]
# pre = paste0('wholeTrans_',ifelse(dimRedMethod == 'PCA',paste0(nPCs,'PCs_'),''))
# parallel = T
# data = summarise_output_mtx(out_list=out_list,nBin=nBin_toUse,mergeBin=mergeBin,
#                             maxScoreVariable = 'ed',doPlot=doPlot,
#                             title = paste0(refTraj,' ',gsub('_$','',pre),' projection results'),
#                             parallel = parallel)
# data_ed = data
# data_ed$branchAssign = gsub('^.*/|:all','',data_ed$binName)
# table(data_ed$celltype,data_ed$branchAssign)
# data_cor = processed_output  
# data_ed

# 
# df = results[[2]]
# ggplot(df,aes(avg_pseudotime,avg_entropy,col=max_correlation))+
#   geom_point()+
#   facet_wrap(vars(group))+
#   theme_classic(base_size = 13)+
#   theme(panel.border = element_rect(fill=F,color='black',linewidth = 1),axis.line = element_blank()) 
# a = df[grepl('Tumour$',df$group) & !grepl('unsure',df$group),] %>% group_by(group,dataset) %>% summarise(med=mean(avg_pseudotime))
# ggplot(a,aes(dataset,med,fill=dataset))+
#   geom_boxplot(outlier.size = 0.001)+
#   #facet_wrap(vars)+
#   theme_classic(base_size = 13)+
#   theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
#         panel.border = element_rect(fill=F,color='black',linewidth = 1),axis.line = element_blank()) 
# 
# ggplot(,aes(group,avg_pseudotime,fill=dataset))+
#   geom_boxplot(outlier.size = 0.001)+
#   #facet_wrap(vars)+
#   theme_classic(base_size = 13)+
#   theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
#         panel.border = element_rect(fill=F,color='black',linewidth = 1),axis.line = element_blank()) 
# 
# 
# 
# results = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/wholeTranscriptome/3_projection_output/fLiver2n_MKlin_traj1200_240708_combinedBranches.REF_wholeTrans_PCA_20PCs_50binByPseudotime.noDS_projection_output.RDS')
# df = res[[2]]




