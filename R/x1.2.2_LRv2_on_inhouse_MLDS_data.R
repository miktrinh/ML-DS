## Automated celltype annotation using CellTypist 
## Run LRv2_w_celltypist on inhouse TAM - MLDS object ##
## REF = foetal liver (2n) atlas (2019)
## REF = Immune_All_Low.pkl	

##------------------##
##    Libraries   ####
##------------------##
library(Seurat)
library(tidyverse)
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/sharedCode-main/logisticRegressionCellTypist.R")

outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/annotation_wCT'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}


setwd(outDir)


##-----------------------------------------##
##    1. Import inhouse MLDS dataset     ####
##-----------------------------------------##
## Import MLDS object
mlds = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_noMTCells.RDS')

##--------------------------------------------------##
####    2. Run LR_v2 on inhouse MLDS dataset      ####
##--------------------------------------------------##
skipIfExists = F

model_fp = '/nfs/users/nfs_m/mt22/.celltypist/data/models/Immune_All_Low.pkl'
tissue = 'MLDS'
mtx = mlds@assays$RNA@counts
mdat = mlds@meta.data
rm(mlds)
lr_output_fp = file.path(outDir,paste0('Immune_All_Low.REF_LRwCT_MLDS.tgt_output.RDS'))

if(file.exists(lr_output_fp) & skipIfExists){
  lr_output = readRDS(lr_output_fp)
}else{
  if(!dir.exists(dirname(lr_output_fp))){
    dir.create(dirname(lr_output_fp))
  }
  lr_output = runCelltypist(cnts=mtx,
                            model=model_fp)
  saveRDS(lr_output,lr_output_fp)
  
  message('Doing plots....')
}
  
  
  
message('Calculating softmax')
## Softmax to assign low_conf cells
softmax_p = lr_output[['logitMat']]
softmax_p = t(apply(softmax_p,1,function(x){sapply(x,function(xi){xi = exp(xi)/sum(exp(x))})}))

## Assign labels to each cell
#  for each cell, find the name with the best match - ie. highest probability
#  assign the celltype if prob >=0.95
#  assign as low_conf if prob >=0.8
#  else - unknown
softmax_label = apply(softmax_p, 1, function(x){
  max_score = max(x)
  second_max_score = max(x[x!=max_score])
  diff = max_score - second_max_score


  lab = ifelse(sum(x == max(x)) > 1,
               paste(c('ambiguous',names(x[x==max(x)])),sep = ':'),
               ifelse(max(x) >= 0.95,
                      ifelse(diff >= 0.1, names(x[x==max_score]),
                             paste(c('ambiguous',names(x[x==max_score]),names(x[x==second_max_score])),sep = ':')),
                      ifelse(max(x) >= 0.85,
                             ifelse(diff >= 0.1, paste0('lowConf_',names(x[x==max_score])),
                                    paste(c('lowConf_ambiguous',names(x[x==max_score]),names(x[x==second_max_score])),sep = ':')),
                             'unknown')))
  return(lab)
})

# pdf(gsub('output.RDS','softmax.pdf',lr_output_fp),width = 30,height = 30)
# library(ComplexHeatmap)
# p = plotByGroup(fit=softmax_p,cellGroupings=mdat$group[match(rownames(softmax_p),rownames(mdat))],colSpaceTruncate=F,colSpaceBreaks = c(0,0.5,1))
# print(p)
# dev.off()

#### Add annotation to the new object
mdat[[paste0('LR_predicted_label_ImmuneAllLow')]] = lr_output[[2]]$predicted_labels[match(rownames(mdat),lr_output[[2]]$X)]
mdat[[paste0('LR_softmax_predicted_label_ImmuneAllLow')]] = softmax_label[match(rownames(mdat),names(softmax_label))]


##---- SAVE LR output
write.csv(mdat,file.path(outDir,'Immune_All_Low.REF_LRwCT_MLDS.tgt_output_results.csv'))






