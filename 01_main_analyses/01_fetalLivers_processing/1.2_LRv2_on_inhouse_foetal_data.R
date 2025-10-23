## Run LRv2_w_celltypist on inhouse fLiver 2n+AK object ##
## REF = foetal liver (2n) atlas (2019)

##------------------##
##    Libraries   ####
##------------------##
library(Seurat)
library(tidyverse)
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/sharedCode-main/logisticRegressionCellTypist.R")

outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/2_annotation/liver/jun24'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}


setwd(outDir)



##-------------------------------------------##
##    1. Import inhouse fLiver dataset     ####
##-------------------------------------------##

big.srat = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0724.RDS')



##----------------------------------------------------##
####    2. Run LR_v2 on inhouse fLiver dataset      ####
##----------------------------------------------------##
tissue = 'liver'
skipIfExists = T

mtx = big.srat@assays$RNA@counts
mdat = big.srat@meta.data
rm(big.srat)
for(ref_tissue in c('liver')){
  #### Logistic Regression v2 (with CT) to annotate the dataset
  message(sprintf('\nTraining model with reference tissue %s ....',ref_tissue))
  
  model_fp = paste0('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/2_annotation/remapped_fetal_',ref_tissue,'_REF_publishedAnnoLRwCT_model.pkl')
  
  if(!file.exists(model_fp)){
    stop(sprintf('LR_v2 model for %s REF does not seem to exist.... Please check and/or re-run script 0.2_train_LRv2_model_on_new_fetalREF.R\n'))
  }else{
    message('Existing model found. Skipping...')
  }
  
  
  
  #### LR prediction with fetal **tissue** ####
  message(sprintf('Predicting label against fetal %s',ref_tissue))
  
  
  lr_output_fp = file.path(outDir,paste0(tissue,'_fetalREF',ref_tissue,'_LRwCT_output.RDS'))
  if(file.exists(lr_output_fp) & skipIfExists){
    print('Reading existing output')
    lr_output = readRDS(lr_output_fp)
  }else{
    if(!dir.exists(dirname(lr_output_fp))){
      dir.create(dirname(lr_output_fp))
    }
    lr_output = runCelltypist(cnts=mtx,
                              model=model_fp)
    saveRDS(lr_output,lr_output_fp)
    
    message('Doing plots....')
    # 
    # pdf(gsub('output.RDS','logit.pdf',lr_output_fp),width = 30,height = 30)
    # library(ComplexHeatmap)
    # p = plotByGroup(fit=lr_output[['logitMat']],cellGroupings=lr_output[['labMat']][,2],colSpaceTruncate=T,colSpaceBreaks = c(-5,0,5))
    # print(p)
    # #p = plotByGroup(fit=lr_output[['logitMat']],cellGroupings=AK.srat$seurat_clusters[match(rownames(lr_output[[1]]),rownames(AK.srat@meta.data))],colSpaceTruncate=T,colSpaceBreaks = c(-5,0,5))
    # #print(p)
    # dev.off()
  }  
  
  
  # message('Calculating softmax')
  # ## Softmax to assign low_conf cells
  # softmax_p = lr_output[['logitMat']]
  # softmax_p = t(apply(softmax_p,1,function(x){sapply(x,function(xi){xi = exp(xi)/sum(exp(x))})}))
  # 
  # ## Assign labels to each cell
  # #  for each cell, find the name with the best match - ie. highest probability
  # #  assign the celltype if prob >=0.95
  # #  assign as low_conf if prob >=0.8
  # #  else - unknown
  # softmax_label = apply(softmax_p, 1, function(x){
  #   max_score = max(x)
  #   second_max_score = max(x[x!=max_score])
  #   diff = max_score - second_max_score
  #   
  #   
  #   lab = ifelse(sum(x == max(x)) > 1,
  #                paste(c('ambiguous',names(x[x==max(x)])),sep = ':'),
  #                ifelse(max(x) >= 0.95,
  #                       ifelse(diff >= 0.1, names(x[x==max_score]), 
  #                              paste(c('ambiguous',names(x[x==max_score]),names(x[x==second_max_score])),sep = ':')),
  #                       ifelse(max(x) >= 0.85,
  #                              ifelse(diff >= 0.1, paste0('lowConf_',names(x[x==max_score])), 
  #                                     paste(c('lowConf_ambiguous',names(x[x==max_score]),names(x[x==second_max_score])),sep = ':')),
  #                              'unknown')))
  #   return(lab)
  # })
  # 
  # # pdf(gsub('output.RDS','softmax.pdf',lr_output_fp),width = 30,height = 30)
  # # library(ComplexHeatmap)
  # # p = plotByGroup(fit=softmax_p,cellGroupings=AK.srat$seurat_clusters[match(rownames(softmax_p),rownames(AK.srat@meta.data))],colSpaceTruncate=F,colSpaceBreaks = c(0,0.5,1))
  # # print(p)
  # # dev.off()
  # 
  # #### Add annotation to the new object
  mdat[[paste0('LR_predicted_label_',ref_tissue)]] = lr_output[[2]]$predicted_labels[match(rownames(mdat),lr_output[[2]]$X)]
  #mdat[[paste0('LR_softmax_predicted_label_',ref_tissue)]] = softmax_label[match(rownames(mdat),names(softmax_label))]

  
  ##---- SAVE LR output
  write.csv(mdat,file.path(outDir,'inhouseFoetal.REF_LRwCT_fLiver.tgt_output_results.csv'))
  
}

## Adjust annotation accordingly ##
big.srat = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0724.RDS')
mdat2 = read.csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/2_annotation/liver/jun24/inhouseFoetal.REF_LRwCT_fLiver.tgt_output_results.csv')
dd = as.data.frame(table(mdat2$annot_jun24,mdat2$LR_predicted_label_liver))
dd = dd[dd$Freq > 0,]

s = subset(big.srat,subset=annot_jun24 %in% c('MEMP_MEP','earlyMK','EE','Mast.cell'))
s = standard_clustering(s,runHarmony = T,harmonyVar = 'assay',clusteringRes = 1.3)
DimPlot(s,group.by = 'annot_aug24',cols = c(col25,pal34H),label = T,repel = T,label.box = T) 

table(s$seurat_clusters[s$annot_jun24=='EE'])
DimPlot(big.srat,cells.highlight = s$cellID[s$seurat_clusters %in% c(3) & s$annot_aug24 == 'EE'])
DimPlot(big.srat,cells.highlight = s$cellID[s$seurat_clusters == 4 & s$annot_aug24 == 'Mast.cell'])
DimPlot(s,cells.highlight = rownames(lr_output)[lr_output$MEMP < 0 & lr_output$Early.Erythroid > 0 & rownames(lr_output) %in% s$cellID[s$seurat_clusters == 17]])

s$group_tmp = ifelse(s$seurat_clusters == 1, s$annot_aug24,s$seurat_clusters)
DotPlot(s,group.by = 'group_tmp',
        features = fLiver_markers)+RotatedAxis()+
  theme(axis.text.x = element_text(size=10,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'top') + xlab('') + ylab('')

table(s$annot_aug24[s$seurat_clusters %in% c(18,19,20,29,30,22)])

s$annot_aug24 = s$annot_jun24
s$annot_aug24[s$seurat_clusters %in% c(0,1,2,4,5,6,7,8,9,10,11,14,15,17,21,23,24,25,27,31) & s$annot_aug24 == 'MEMP_MEP'] = 'EE'
s$annot_aug24[s$seurat_clusters %in% c(0,1,2,4,5,6,7,8,9,10,11,14,15,17,21,23,24,25,27,31) & s$annot_aug24 == 'earlyMK'] = 'MEMP_MEP'
s$annot_aug24[s$seurat_clusters %in% c(3)] = 'MEMP_MEP'
s$annot_aug24[s$seurat_clusters %in% c(5) & s$annot_aug24 == 'Mast.cell'] = 'MEMP_MEP'
s$annot_aug24[s$seurat_clusters %in% c(12) & s$annot_aug24 == 'MEMP_MEP'] = 'earlyMK'



## Add back to big.srat
big.srat$annot_aug24 = as.character(big.srat$annot_jun24)
big.srat$annot_aug24[big.srat$cellID %in% s$cellID] = s$annot_aug24[match(big.srat$cellID[big.srat$cellID %in% s$cellID],s$cellID)]
DotPlot(big.srat,group.by = 'annot_aug24',
        features = fLiver_markers)+RotatedAxis()+
  theme(axis.text.x = element_text(size=10,angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(size=10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'top') + xlab('') + ylab('')

DimPlot(big.srat,group.by = 'annot_aug24',cols = c(col25,pal34H),label = T,repel = T,label.box = T) 

# earlyMK_klf1 = WhichCells(big.srat,cells = big.srat$cellID[big.srat$annot_aug24 == 'earlyMK'],expression = (IL1B>0.5 & KLF1 > 0.5) )
# DimPlot(big.srat,cells.highlight = earlyMK_klf1)
# big.srat$annot_aug24[big.srat$cellID %in% earlyMK_klf1] = 'earlyMK_1'
# big.srat$annot_aug24[big.srat$annot_aug24 == 'B.cell' & big.srat$LR_predicted_label_liver == 'B.cell'] = 'pre.B.cell_2'
# big.srat$annot_aug24[big.srat$annot_aug24 == 'pre.B.cell_2'] = 'B.cell'
DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$annot_aug24 == 'earlyMK'])
DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$annot_aug24 == 'B.cell' & big.srat$LR_predicted_label_liver == 'B.cell'])


big.srat$finalAnn_broad = as.character(big.srat$annot_aug24)
big.srat$annot = as.character(big.srat$annot_aug24)
big.srat$finalAnn = as.character(big.srat$annot_aug24)
big.srat$broadLineage = as.character(big.srat$finalAnn_broad)
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('pro.B.cell','pre.B.cell','B.cell')] = 'B lineage'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('DC1','DC2','pDC')] = 'Dendritic cells'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('EE','ME','LE')] = 'Erythroblasts'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','CMP_GMP','LMPP_ELP')] = 'HSC & prog.'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('Mast.cell')] = 'Mast.cell'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('MK','earlyMK')] = 'Megakaryocytes'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('MOP','proMono','Monocyte','Macrophage','Kupffer.cell')] = 'Monocyte/Macrophage'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('promyelocyte','myelocyte')] = 'Myelocytes'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('Endo','NPC','Fibroblast','Mesenchyme','Hepatocyte')] = 'Stromal'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('ILC.precursor','T.cell','NK_T')] = 'T/NK lineage'
big.srat$broadLineage[big.srat$finalAnn_broad %in% c('doublets','HB.cluster','others','unknown')] = 'others'


## Without doublets
saveRDS(big.srat,'~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS')
mdat = cbind(big.srat@meta.data,big.srat@reductions$umap@cell.embeddings)
mdat$broadLineage = as.character(mdat$finalAnn_broad)
mdat$broadLineage[mdat$finalAnn_broad %in% c('pro.B.cell','pre.B.cell','B.cell')] = 'B lineage'
mdat$broadLineage[mdat$finalAnn_broad %in% c('DC1','DC2','pDC')] = 'Dendritic cells'
mdat$broadLineage[mdat$finalAnn_broad %in% c('EE','ME','LE')] = 'Erythroblasts'
mdat$broadLineage[mdat$finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','CMP_GMP','LMPP_ELP')] = 'HSC & prog.'
mdat$broadLineage[mdat$finalAnn_broad %in% c('Mast.cell')] = 'Mast.cell'
mdat$broadLineage[mdat$finalAnn_broad %in% c('MK','earlyMK')] = 'Megakaryocytes'
mdat$broadLineage[mdat$finalAnn_broad %in% c('MOP','proMono','Monocyte','Macrophage','Kupffer.cell')] = 'Monocyte/Macrophage'
mdat$broadLineage[mdat$finalAnn_broad %in% c('promyelocyte','myelocyte')] = 'Myelocytes'
mdat$broadLineage[mdat$finalAnn_broad %in% c('Endo','NPC','Fibroblast','Mesenchyme','Hepatocyte','Neuron','Cholangiocytes','Mesothelial_cells')] = 'Stromal'
mdat$broadLineage[mdat$finalAnn_broad %in% c('ILC.precursor','T.cell','NK_T')] = 'T/NK lineage'
mdat$broadLineage[mdat$finalAnn_broad %in% c('doublets','HB.cluster','others','unknown','Trophoblast')] = 'others'

write.csv(mdat,'~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824_mdat.csv')


##---- SAVE tmp big.srat
#write.csv(big.srat@meta.data,'~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_clean_LRv2_mdat.csv')

