##-------- Prepare Reference foetal objects (fLiver, fAdr, fKidney)  --------##
## 2. Train LR_v2 (celltypist) on remapped foetal references
## [230324]: Only completed for fLiver, using published annotation


setwd('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/')

#------------------------#
##      Libraries     ####
#------------------------#
library(tidyverse)
library(Seurat)
library(RColorBrewer)
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/sharedCode-main/logisticRegressionCellTypist.R")
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/misc.R")
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/sc_utils.R")


#---------------------------------------------------------#
##    1. Train LRv2 model on remapped foetal REF       ####
##       Train a different model for each tissue         ##  
##       Using published annotations                     ##
#---------------------------------------------------------#
use_published_annot = T
for(tissue in c('liver')){
  message(sprintf('Training new model for tissue %s',tissue))  
  new_model_fp = paste0('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/2_annotation/remapped_fetal_',tissue,'_REF_publishedAnnoLRwCT_model.pkl')
  ## Import reference seurat object 
  #srat = readRDS(paste0('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/',tissue,'/oct22/',tissue,'_clean_filtered_annotated.RDS'))
  
  ## Sub-sampling the object if it's too large...
  if(tissue == 'liver'){
    if(use_published_annot){
      ##--- Import QC_preFiltered remapped ref object
      srat = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/liver/oct22/liver_clean_withMTCells.RDS')
      ##--- Only keep published cells
      published_mdat = read.csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/Muz_fLiver_SoupXed_sratObj_metadata.csv')
      srat$cellID = gsub('-.*$','',srat$cellID)
      message(sprintf('Of the %d published cells, %d cells are retained',nrow(published_mdat),sum(srat$cellID %in% published_mdat$cellID)))
      srat = subset(srat,subset = cellID %in% published_mdat$cellID)
      srat$finalAnn = published_mdat$cell.labels[match(srat$cellID,published_mdat$cellID)]
    }
    
    ### Remove 50% of any cell types with >10000 cells
    nCell_per_ct = table(srat$finalAnn)
    ct_to_50perc = names(nCell_per_ct[nCell_per_ct>10000])

    srat$rownames = rownames(srat@meta.data)
    set.seed(123)
    # cellID_tokeep = srat@meta.data %>% group_by(finalAnn) %>%
    #   mutate(id = seq(1,n()),
    #          selected = ifelse(finalAnn %in% ct_to_25perc, id %in% sample(c(1:max(id)),round(0.25*max(id))),
    #                            ifelse(finalAnn %in% ct_to_50perc,id %in% sample(c(1:max(id)),round(0.5*max(id))),T)))
    cellID_tokeep = srat@meta.data %>% group_by(finalAnn) %>%
      mutate(id = seq(1,n()),
             selected = ifelse(finalAnn %in% ct_to_50perc, id %in% sample(c(1:max(id)),round(0.5*max(id))),T))
    cellID_tokeep = cellID_tokeep$rownames[cellID_tokeep$selected]

    srat = subset(srat, subset = rownames %in% cellID_tokeep)
  }else if(tissue == 'adrenal'){
    ### Remove 50% of cortex cells
    nCell_per_ct = table(srat$finalAnn)
    ct_to_25perc = names(nCell_per_ct[nCell_per_ct>10000])
    ct_to_50perc = c()
    #ct_to_50perc = names(nCell_per_ct[nCell_per_ct>5000 & !names(nCell_per_ct) %in% ct_to_25perc])

    srat$rownames = rownames(srat@meta.data)
    set.seed(123)
    cellID_tokeep = srat@meta.data %>% group_by(finalAnn) %>%
      mutate(id = seq(1,n()),
             selected = ifelse(finalAnn %in% ct_to_25perc, id %in% sample(c(1:max(id)),round(0.25*max(id))),
                               ifelse(finalAnn %in% ct_to_50perc,id %in% sample(c(1:max(id)),round(0.5*max(id))),T)))
    cellID_tokeep = cellID_tokeep$rownames[cellID_tokeep$selected]

    srat = subset(srat, subset = rownames %in% cellID_tokeep)
  }

  message(sprintf('Training new model on srat object of %d cells.',ncol(srat)))

  if(!file.exists(new_model_fp)){
    trainCelltypistModel(cnts = srat@assays$RNA@counts,
                         outPath = new_model_fp,
                         labels = srat@meta.data$finalAnn,
                         n_jobs = 24)

  }else{
    message('Existing model found. Skipping...')
  }

}






#--------------------------------------------------------------------#
##    2. Train LRv2 model on ALL remapped foetal REF              ####
##       Train 1 model on combined foetal REF with all 3 tissue     ##  
##       Using published annotations?                               ##
#--------------------------------------------------------------------#

#### Make a massive REF seurat object, combining Adrenal Liver and Kidney ####
## big.srat object is 17.2 Gb (322,949 cells)
# big.srat_fp = paste0('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/allREF_clean_filtered_annotated.RDS')
# if(!file.exists(big.srat_fp)){
#   big.srat = NULL
#   for(tissue in c('adrenal','kidney','liver')){
#     srat_fp = paste0('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/',tissue,'/oct22/',tissue,'_clean_filtered_annotated.RDS')
#     srat = readRDS(srat_fp)
#     if(is.null(big.srat)){
#       big.srat = srat
#     }else{
#       big.srat = merge_seurat_objects(srat1 = big.srat,srat2 = srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
#     }
#   }
#   
#   saveRDS(big.srat,big.srat_fp)
# }else{
#   big.srat = readRDS(big.srat_fp)
# }
# 
# big.srat$tissue = gsub('^.*_','',big.srat$donorID)
# big.srat$finalAnn[big.srat$finalAnn %in% c('early T cell','Early.lymphoid_T.lymphocyte')] = 'early.T.cell'
# big.srat$finalAnn[big.srat$finalAnn %in% c('Endothelial.cell')] = 'Endothelium'
# big.srat$finalAnn[big.srat$finalAnn %in% c('Mast cell')] = 'Mast.cell'
# big.srat$finalAnn_tissue = paste0(big.srat$finalAnn,'_',big.srat$tissue)
# 
# new_model_fp = paste0('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/2_annotation/remapped_fetal_allREF_LRwCT_model.pkl')



# ## Remove 50% of any cell types with >5000 cells, and 25% of those with >10k cells
# nCell_per_ct = table(big.srat$finalAnn_tissue)
# ct_to_25perc = names(nCell_per_ct[nCell_per_ct>10000])
# ct_to_50perc = names(nCell_per_ct[nCell_per_ct>5000 & !names(nCell_per_ct) %in% ct_to_25perc])
# 
# big.srat$rownames = rownames(big.srat@meta.data)
# set.seed(123)
# cellID_tokeep = big.srat@meta.data %>% group_by(finalAnn_tissue) %>%
#   mutate(id = seq(1,n()),
#          selected = ifelse(finalAnn_tissue %in% ct_to_25perc, id %in% sample(c(1:max(id)),round(0.25*max(id))),
#                            ifelse(finalAnn_tissue %in% ct_to_50perc,id %in% sample(c(1:max(id)),round(0.5*max(id))),T)))
# cellID_tokeep = cellID_tokeep$rownames[cellID_tokeep$selected]
# 
# big.srat = subset(big.srat, subset = rownames %in% cellID_tokeep)
# 
# 
# message(sprintf('Training new model on srat object of %d cells.',ncol(big.srat)))
# 
# 
# if(!file.exists(new_model_fp)){
#   trainCelltypistModel(cnts = big.srat@assays$RNA@counts,
#                        outPath = new_model_fp,
#                        labels = big.srat@meta.data$finalAnn_tissue,
#                        n_jobs = 23)
# 
# }else{
#   message('Existing model found. Skipping...')
# }













