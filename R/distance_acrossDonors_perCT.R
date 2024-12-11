##----   An attempt to quantify the "distance" between donorID across different cell types    -----##
## Motivation: We want to show that TAM cells between patients are more transcriptionally similar to each other than MLDS (full blown leukaemia - highly patient specific)
## Expectation: Some metrics showing that:
  # 1. absolute normal celltypes: high degree of overlaps between donorIDs
  # 2. TAM: donor specific, but TAM cells between donors are NOT too different from each other
  # 3. MLDS: full-blown cancer --> highly specific to patient --> distinct clusters by patient ID

## Method: I thought of a few methods to try:
# 1. Literally calculate the correlation between every single cell --> compare correlation distribution within patient vs between patients
# 2. pseudobulk by donorID per celltype --> Heatmap of top HVG --> hierarchical clustering
# 3. pseudobulk by donorID per celltype --> PCA --> calculate distance between donorID per celltype
# 4. mahalanobis distance of some sort??
# 5. average silhouette width
# 6. manhattan distance of each cell towards patients centroid
# 7. correlation / eucledian distance between each cell towards its own patient-level pseudobulk, and other patient's pseudobulk

##------------------------##
outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/distance_acrossDonors_perCT/'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)


##---------------##
#### Libraries ####
##---------------##
library(Seurat)
library(tidyverse)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
#source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")




##------------------------------------------------##
#### 1. Import relevant scRNA datasets: MLDS    ####
##-----------------------------------------------##
## Import MLDS object
mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
mlds$annot_mar24 = mdat$annot_mar24[match(mlds$cellID,mdat$cellID)]
mlds$broadLineage = mdat$broadLineage[match(mlds$cellID,mdat$cellID)]

#DimPlot(mlds,group.by = 'annot_mar24',label = T,repel = T,label.box = T) + NoLegend()
mlds$group = ifelse(mlds$annot_mar24 == 'Tumour',mlds$disease,mlds$annot_mar24)

##-------------------------------------------------------------------------------------##
#### METHOD 1. Raw correlation between every single cell per celltype per donorID    ####
##-------------------------------------------------------------------------------------##





##------------------------------------------------------##
#### METHOD 3. pb per celltype per donorID   --> PCA  ####
##------------------------------------------------------##

# mlds$group4pb = paste0(mlds$group,'_',mlds$donorID)
# 
# pb = do.call(cbind,lapply(split(colnames(mlds),mlds$group4pb[match(colnames(mlds),mlds$cellID)]),function(e) rowSums(mlds@assays$RNA@counts[,e,drop=FALSE])))
# pb = pb[,!grepl('unsure|NA|\\?|Tumour_WT|doublets',colnames(pb))]
# ## Only use group4pb with >= 100 cells
# nCell_perGroup = as.data.frame(table(mlds$group4pb))
# nCell_perGroup = nCell_perGroup[nCell_perGroup$Freq >= 100,]
# 
# 
# message('Computing reference PCA')
# s = CreateSeuratObject(counts = pb[,colnames(pb) %in% nCell_perGroup$Var1])
# s = NormalizeData(s,verbose = F)
# print(sprintf('Number of pseudobulk groups: %d',ncol(s)))
# s = FindVariableFeatures(s,nfeatures = 5000)
# s = ScaleData(s,verbose = F)
# s = RunPCA(s, npcs = 70,verbose = F)
# # Elbow plot
# p = ElbowPlot(s)
# print(p)
# 
# pca_df = as.data.frame(s@reductions$pca@cell.embeddings)
# pca_df$celltype = mlds$group[match(rownames(pca_df),mlds$group4pb)]
# pca_df$donorID = mlds$donorID[match(rownames(pca_df),mlds$group4pb)]
# pca_df$celltype_group = ifelse(pca_df$celltype %in% c('TAM','MLDS'),pca_df$celltype,'Normal')
# # pca_df.sub$pseudoTime = traj$pseudotime[match(rownames(pca_df.sub),as.character(traj$cellID))]
# # pca_df.sub$mainRefCT = traj$finalAnn_broad[match(rownames(pca_df.sub),as.character(traj$cellID))]
# # #pca_df.sub$mainRefCT_ct = sapply(strsplit(pca_df.sub$mainRefCT,split=':'),'[',1)
# 
# p1 = ggplot(pca_df[!pca_df$celltype %in% c('TAM','MLDS'),],aes(PC_1,PC_2,col=donorID,shape=celltype_group))+
#   geom_point() +
#   geom_point(data = pca_df[pca_df$celltype %in% c('TAM','MLDS'),],size=3) +
#   scale_color_manual(values = c(col25,pal34H))+
#   scale_shape_manual(values = c(15,19,17))+
#   theme_classic(base_size = 15) + 
#   theme(panel.border = element_rect(fill=F),axis.line = element_blank())+
#   ggtitle('pb by celltype by donorID')
# 
# print(p1)
# 
# 
# ## Eucledian distance between donor, using the top 20PC
# pca_mtx = pca_df[,c(colnames(pca_df)[1:20],'celltype','donorID')]
# out_df = data.frame()
# for(ct in unique(pca_df$celltype)){
#   df = pca_mtx[pca_mtx$celltype == ct,colnames(pca_mtx)[!colnames(pca_mtx) %in% c('celltype','donorID')]]
#   if(nrow(df) <= 1){next}
#   dist_df = as.matrix(dist(df, method = "euclidean", diag = F))
#   dist_df = dist_df[upper.tri(dist_df,diag = F)]
#   tmp = data.frame(celltype = ct,
#                    dist_betweenDonor = dist_df)
#   
#   if(nrow(out_df) == 0){
#     out_df = tmp
#   }else{
#     out_df = rbind(out_df,tmp)  
#   }
#   
# }
# 
# library(ggbeeswarm)
# ggplot(out_df,aes(x=reorder(celltype,dist_betweenDonor,FUN = median),dist_betweenDonor))+
#   geom_quasirandom(width = 0.2,size=0.5)+
#   stat_summary(
#     aes(group = celltype), fun = median, fun.min = median, fun.max = median,
#     geom = "crossbar", color = "black", width = 0.7, lwd = 0.2,
#     
#     # add this bit here to your stat_summary function
#     position=position_dodge(width=0.75))+
#   #geom_boxplot(alpha=0.4,width=0.3)+
#   theme_classic(base_size = 14)+
#   theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
#         axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))+
#   xlab('cell type') + ylab('Eucledian distance between donors')



##------------------------------------------------------------------------------------------------------##
#### METHOD 6. Manhattan distance of each cell towards its own centroid and other donor_CT's centroid ####
##------------------------------------------------------------------------------------------------------##
# celltype = 'pDC'
# celltype = 'TAM'
# 
# manDist_matrix_list=list()
# manDist_matrix_long_list = list()
# for(ct in unique(mlds$annot_mar24[!mlds$broadLineage %in% c('lowQual','Tumour_unsure')])){
#   print(sprintf('Working on %s',ct))
#   
#   mtx = mlds@assays$RNA@data[,colnames(mlds@assays$RNA@data) %in% mlds$cellID[mlds$annot_mar24 == ct]]
#   # Compute donor centroid
#   centroid_byDonor = do.call(cbind,lapply(split(colnames(mtx),mlds$donorID[match(colnames(mtx),mlds$cellID)]),function(e) rowMeans(mtx[,e,drop=FALSE])))
#   # Calculate Manhattan distance of each cell towards each centroid
#   manDist_matrix = do.call(rbind,lapply(1:ncol(mtx),function(i){
#     tmp_mtx = cbind(mtx[,i],centroid_byDonor)
#     manDist = as.matrix(dist(t(tmp_mtx),method = 'manhattan'))
#     manDist = manDist[1,-1]
#     return(manDist)
#   }))
#   rownames(manDist_matrix) = colnames(mtx)
#   manDist_matrix = as.data.frame(manDist_matrix)
#   manDist_matrix$donorID = mlds$donorID[match(rownames(manDist_matrix),mlds$cellID)]
#   manDist_matrix$cellID = rownames(manDist_matrix)
#   manDist_matrix_long = pivot_longer(manDist_matrix,cols = colnames(manDist_matrix)[!colnames(manDist_matrix) %in% c('donorID','cellID')],names_to = 'partner_donor',values_to = 'distance')
#   manDist_matrix_long$type = ifelse(manDist_matrix_long$donorID == manDist_matrix_long$partner_donor,'within','between')
#   manDist_matrix_long$celltype = ct
#   manDist_matrix_list[[ct]] = manDist_matrix
#   manDist_matrix_long_list[[ct]] = manDist_matrix_long
# }
# 
# 
# 
# 
# ggplot(manDist_matrix_long,aes(type,distance,col=type))+
#   geom_point()+
#   facet_wrap(vars(donorID),nrow=1)
# # Correlation between cells for each donor
# for(d in unique(mlds$donorID[mlds$cellID %in% colnames(mtx)])){
#   withinD_mtx = mtx[,colnames(mtx) %in% mlds$cellID[mlds$annot_mar24 == ct & mlds$donorID == d]]
#   withinD_cor = cor(as.matrix(withinD_mtx))
#   
#   mtx
# }
# 
# 
# 
# ## Extract normalized count matrix for relevant cells
# s = subset(mlds,subset = cellID %in% mlds$cellID[mlds$group == celltype])
# mtx = s@assays$RNA@data
# ## calculate correlation of any one cell to all other cells
# cor_mtx = cor(as.matrix(mtx))
# 
# ## Split correlation matrix by donorID
# cor_mtx_byDonor =  lapply(split(colnames(cor_mtx),mlds$donorID[match(colnames(cor_mtx),mlds$cellID)]),function(e) cor_mtx[,e,drop=FALSE])
# output = data.frame()
# 
# for(d in unique(s$donorID)){
#   print(d)
#   # within donorID correlation
#   withinDonor_mtx = cor_mtx[rownames(cor_mtx) %in% s$cellID[s$donorID == d],colnames(cor_mtx) %in% s$cellID[s$donorID == d]]
#   withinDonor_mtx = withinDonor_mtx[upper.tri(withinDonor_mtx,diag = F)]
#   
#   betweenDonor_mtx = cor_mtx[rownames(cor_mtx) %in% s$cellID[s$donorID == d],!colnames(cor_mtx) %in% s$cellID[s$donorID == d]]
#   # Not a diagonal matrix
#   table(colnames(betweenDonor_mtx) %in% rownames(betweenDonor_mtx))
#   betweenDonor_mtx = as.vector(betweenDonor_mtx)
#   
#   tmp = data.frame(celltype = celltype,
#                    donorID = d,
#                    withinDonor_meanCor = mean(withinDonor_mtx),
#                    withinDonor_medCor = median(withinDonor_mtx),
#                    betweenDonor_meanCor = mean(betweenDonor_mtx),
#                    betweenDonor_medCor = median(betweenDonor_mtx))
#   if(nrow(output) == 0){
#     output = tmp
#   }else{
#     output = rbind(output,tmp)  
#   }
#   
# }
# 
# 
# ggplot(output,aes(withinDonor_meanCor,betweenDonor_meanCor,col=celltype))+
#   geom_point()+
#   theme_classic()+
#   xlim(0,1)+ylim(0,1)+
#   theme(panel.border = element_rect(fill=F))
# 


##---------------------------------------------------------------------------------------------##
#### METHOD 7a. Correlation between each cell and its patient-level pb vs other patients pb  ####
##---------------------------------------------------------------------------------------------##
metaCell_method = c('pb','centroid')
metaCell_method = 'centroid'
metric = 'manhattan'
## Calculate pb of the celltype for each donor
mlds$group4pb = paste0(mlds$group,'_',mlds$donorID,'_',mlds$timePoint)
## Only use group4pb with >= 100 cells
nCell_perGroup = as.data.frame(table(mlds$group4pb))
nCell_perGroup = nCell_perGroup[nCell_perGroup$Freq >= 100,]



if(metaCell_method == 'pb'){
  pb = do.call(cbind,lapply(split(colnames(mlds),mlds$group4pb[match(colnames(mlds),mlds$cellID)]),function(e) rowSums(mlds@assays$RNA@counts[,e,drop=FALSE])))
  pb = pb[,!grepl('unsure|NA|\\?|Tumour_WT|doublets',colnames(pb))]
  pb = pb[,colnames(pb) %in% nCell_perGroup$Var1]
  
}else if(metaCell_method == 'centroid'){
  # Compute donor centroid
  centroid_byGroup = do.call(cbind,lapply(split(colnames(mlds)[colnames(mlds) %in% mlds$cellID[mlds$group4pb %in% nCell_perGroup$Var1]],
                                                mlds$group4pb[match(colnames(mlds)[colnames(mlds) %in% mlds$cellID[mlds$group4pb %in% nCell_perGroup$Var1]],mlds$cellID)]),function(e) rowMeans(mlds@assays$RNA@counts[,e,drop=FALSE])))
  centroid_byGroup = centroid_byGroup[,!grepl('unsure|NA|\\?|Tumour_WT|doublets',colnames(centroid_byGroup))]
  centroid_byGroup = centroid_byGroup[,colnames(centroid_byGroup) %in% nCell_perGroup$Var1]
}


compute_distance_metric = function(celltype,metaCell,mdat, metric = c('cor','euclidean','manhattan'),nCores=1){
  if(nCores > 1){
    library(doParallel)
    
    #Setup backend to use many processors
    totalCores = detectCores()
    
    #Leave some core to avoid overload your computer
    nCore_touse = min(totalCores[1]-5,nCores)
    cluster <- makeCluster(nCore_touse)
    registerDoParallel(cluster)  
  }
  
  out_df = data.frame()
  metaCell.sub = metaCell[,colnames(metaCell)[grepl(celltype,colnames(metaCell))]]
  
  ## Define query cell matrix
  mtx = mlds@assays$RNA@counts[,colnames(mlds@assays$RNA@counts) %in% mlds$cellID[mlds$group == celltype & !is.na(mlds$group)]]
  ## Split the matrix into multiple batches to run in parallel
  mtx_list = lapply(split(1:ncol(mtx),c(1:ncol(mtx)) %/% 1000), function(e){mtx[,e]})
  
  if(nCores > 1){
    start_time <- Sys.time()
    
    dist_df <- foreach(i = 1:length(mtx_list), .combine=rbind) %dopar% {
      print(i)
      x = mtx_list[[i]]
      out = lapply(1:ncol(metaCell.sub),function(j){
        y = metaCell.sub[,j]
        df = cbind(as.matrix(x),y[match(rownames(x),names(y))])
        colnames(df)[ncol(df)] = colnames(metaCell.sub)[j]
        dist_matrix = as.matrix(dist(t(df),method = 'manhattan'))
        dist_value = dist_matrix[,colnames(metaCell.sub)[j]]
        tmp = data.frame(cellID = colnames(x),
                         refID = colnames(metaCell.sub)[j],
                         dist_value = dist_value[match(colnames(x),names(dist_value))])
        return(tmp)
      })
      
      return(out)
    }
    
    #Stop cluster
    stopCluster(cluster)
    end_time <- Sys.time()
    print(end_time - start_time)
  }else{
    dist_df <- do.call(rbind,lapply(1:length(mtx_list),function(i){
      print(i)
      x = mtx_list[[i]]
      out = lapply(1:ncol(metaCell.sub),function(j){
        y = metaCell.sub[,j]
        df = cbind(as.matrix(x),y[match(rownames(x),names(y))])
        colnames(df)[ncol(df)] = colnames(metaCell.sub)[j]
        dist_matrix = as.matrix(dist(t(df),method = 'manhattan'))
        dist_value = dist_matrix[,colnames(metaCell.sub)[j]]
        tmp = data.frame(cellID = colnames(x),
                         refID = colnames(metaCell.sub)[j],
                         dist_value = dist_value[match(colnames(x),names(dist_value))])
        return(tmp)
      })
      
      return(out)
      
    }))
    
  }
  
  
  dist_df$annot = mdat$annot_mar24[match(dist_df$cellID,mdat$cellID)]
  dist_df$donorID = mdat$donorID[match(dist_df$cellID,mdat$cellID)]
  dist_df$timePoint = mdat$timePoint[match(dist_df$cellID,mdat$cellID)]

  return(dist_df)
}






# for(d in unique(mlds$donorID[mlds$group == celltype & !is.na(mlds$group)])){
#   if(length(mlds$cellID[mlds$group == celltype & !is.na(mlds$group) & mlds$donorID == d])<=1){next}
#   mtx = mlds@assays$RNA@counts[,colnames(mlds@assays$RNA@counts) %in% mlds$cellID[mlds$group == celltype & !is.na(mlds$group) & mlds$donorID == d]]
#   
#   cor_mtx = apply(mtx, 2,function(x){
#     cor_val = cor(x,metaCell.sub)
#     names(cor_val) = colnames(metaCell.sub)
#     return(cor_val)
#   })
#   
#   cor_mtx = t(cor_mtx)
#   
#   # Calculate Manhattan distance of each cell towards each centroid
#   manDist_matrix = do.call(rbind,lapply(1:ncol(metaCell.sub),function(i){
#     tmp_mtx = cbind(as.matrix(mtx),metaCell.sub[,i])
#     manDist = as.matrix(dist(t(tmp_mtx),method = 'manhattan'))
#     manDist = manDist[1,-1]
#     return(manDist)
#   }))
#   
#   
#   out_df = rbind(out_df,as.data.frame(cor_mtx))
# }
# 
# out_df$cellID = rownames(out_df)
# out_df = pivot_longer(out_df,cols = colnames(out_df)[colnames(out_df) != 'cellID'])




## Pseudobulk correlation  
# tam_cor = compute_distance_metric(celltype = 'TAM',metaCell = pb,mdat = mdat)
# mlds_cor = compute_cor(celltype = 'MLDS',metaCell = pb,mdat = mdat)
# b_cor = compute_cor(celltype = 'naive.B',metaCell = pb,mdat = mdat)

## Centroid Manhattan distance
tam_cor_centroid = compute_distance_metric(celltype = 'TAM',metaCell = centroid_byGroup,mdat = mdat,metric='manhattan',nCores=10)
print('completed TAM')
mlds_cor_centroid = compute_distance_metric(celltype = 'MLDS',metaCell = centroid_byGroup,mdat = mdat,metric='manhattan',nCores=10)
print('completed MLDS')

data = do.call(rbind,list(tam_cor_centroid,mlds_cor_centroid))
data$group = mlds$group[match(data$cellID,mlds$cellID)]


write.csv(data,paste0(metaCell_method,'_',metric,'_distance_2406.csv'))
# ggplot(data[#!grepl('CC',mlds_cor$donorID) & 
#   data$timePoint == 'Diagnostic',],aes(name,value,fill=donorID))+
#   geom_boxplot(outlier.shape = NA)+
#   #geom_jitter(size=0.1,alpha=0.3)+
#   scale_fill_manual(values = col25)+
#   facet_grid(timePoint~group,scales = 'free_x')+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
# 
# 
# 
# 
# 
# ## Subclustering cancer cells ####
# tam = subset(mlds,subset = donorID %in% mlds$donorID[mlds$disease == 'TAM'])
# tam = standard_clustering(tam)
# 
# DimPlot(tam,group.by = 'annot_mar24',label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend() + ggtitle('TAM') + NoAxes()
# 
# mlds_sub = subset(mlds,subset = donorID %in% mlds$donorID[mlds$disease == 'MLDS' & mlds$timePoint %in% c('Diagnostic')])
# mlds_sub = standard_clustering(mlds_sub)
# 
# DimPlot(mlds_sub,group.by = 'donorID',label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend() + ggtitle('ML-DS') + NoAxes()
# 
# 
# 
# 
