##----   Processing leukaemia scRNAseq datasets    -----##

outDir = "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia"

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)



##---------------##
#### Libraries ####
##---------------##
library(tidyverse)
library(Seurat)
library(SoupX)
library(RColorBrewer)
library(readxl)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/sc_basicQC.R")


##----------------------------##
##   Set Global parameters  ####
##----------------------------##

maxMT = 30
minGenes = 300
minUMIs = 500  
maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup = F
scrubScoreMax = 0.5
scrubPath='./scrubletScores.tsv'
scPath="./strainedCounts"
doPlot=T
verbose = T
skipIfExists=T
keepMTCells=T



##----------------------##
##     TAM / ML-DS    ####
##----------------------##

sub_outDir = file.path(outDir,'MLDS')
plotDir = file.path(sub_outDir,'MLDS_')
outPath = file.path(sub_outDir,'MLDS')
# sub_outDir = file.path(outDir,'test')
# plotDir = file.path(sub_outDir,'test_')
# outPath = file.path(sub_outDir,'test')
metadata = NULL
matchBy = NULL
cleanCountDir = file.path(sub_outDir,'cleanCount')



cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/MLDS_clean_withMTCells.RDS'),paste0(sub_outDir,'/MLDS_clean_noMTCells.RDS'))
# if(file.exists(cleanSrat_fp) & skipIfExists){
#   cleanSrat = readRDS(cleanSrat_fp)
# }else{
  if(!dir.exists(sub_outDir)){
    message(sprintf('Creating output directory'))
    dir.create(sub_outDir,recursive = T)
  }
  
  setwd(sub_outDir)
  
  # List location of cellranger outputs
  dataDirs = c(list.dirs('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX',full.names = T,recursive = T),
               list.dirs('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/',recursive = T,full.names = T))
  dataDirs = dataDirs[!grepl('tic-3108',dataDirs)]
  dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
  dataDirs = dataDirs[!grepl('GRCh38-1_2_0',dataDirs)]
  
  
  names(dataDirs) = basename(dirname(dataDirs))
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_ALeuk_','',names(dataDirs)))
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',names(dataDirs)))
  names(dataDirs) = gsub('_','.',gsub('_GRCh38-2020-A.*$','',gsub('^.*_MY_','MY_',names(dataDirs))))
  
  ## Remove L062 - LPD samples
  dataDirs = dataDirs[!names(dataDirs) %in% c("Leuk13234200", "Leuk13234201", "Leuk13234202", "NB14406184", "NB14406185")]
  length(names(dataDirs))
  #table(names(dataDirs) %in% mlds2$orig.ident)
  #names(dataDirs)[!names(dataDirs) %in% mlds2$orig.ident]
  dataDirs=dataDirs[file.exists(dataDirs)]
  #dataDirs = dataDirs[grepl('RNA14831998|RNA14831999',dataDirs)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))
  
  
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose)
  print(n_distinct(QC.output[[1]]$orig.ident))
  

  # cleanSrat = QC.output[[1]]
  #
  # df.out = QC.output[[2]]
  #
  # write.csv(df.out,paste0(sub_outDir,'/MLDS_240422_qc_summary.csv'))


#}

## Below are the code I used before July 2024
# numPCs_soupX=75
# metadata=NULL
# matchBy=NULL
# rho_max_limit=NULL
# excludeGenes = c()
# geneSets=list()
# mtPattern = '^MT-'
# is10X=TRUE
# sPhaseGenes=NULL
# g2mPhaseGenes=NULL
# doPlot=TRUE
# #plotDir=NULL
# verbose=TRUE
# 
# require(tidyverse)
# require(Seurat)
# require(ggVennDiagram)
# require(RColorBrewer)
# 
# # Check if plotDir is provided
# if(doPlot & is.null(plotDir)){
#   warning('Please provide plotDir or set doPlot=F')
#   return()
# }
# if(is.null(dataDirs)){
#   warning('Please provide dataDirs')
#   return()
# }
# 
# if(!is.null(outPath) & file.exists(outPath) & skipIfExists==F){
#   # Read in existing output
#   message(sprintf('Existing clean output detected:\noutPath = %s \nTo overwite existing file, please specify skipIfExists = F',outPath))
#   srat = read.RDS(outPath)
# 
#   df = srat@misc$preQC@meta.data
#   df.out = df %>% group_by(orig.ident) %>% summarise(nCells_prefilter=n(),medMT=median(percent.mt),muMT=mean(percent.mt),nPASS = sum(PASS),nDoublet = sum(PASS_doublet==FALSE),medCnt = median(nCount_RNA),medGenes=median(nFeature_RNA),nBadClust = sum(PASS_cluster==FALSE),
#                                                      lowCnt = sum(PASS_nCounts==F),highMT = sum(PASS_MT==F),lowGenes = sum(PASS_nGenes==F))
#   return(list(srat,df.out))
# }
# 
# if(is.null(names(dataDirs))){
#   warning("Data directories not named, generating automatic names.")
#   names(dataDirs) = paste0('DataSource',seq_along(dataDirs))
# }
# 
# if(verbose)
#   message('Loading data.')
# if(!is.null(dataDirs)){
#   srat = Read10X(dataDirs)
#   if(class(srat) == 'list' & 'Gene Expression' %in% names(srat)){
#     srat = srat[['Gene Expression']]
#   }
# }else{
#   if(!is.null(mtx)){
#     srat = mtx
#   }else{
#     stop('No data provided. Please check!')
#   }
# }
# 
# #if(ncol(srat) > 50000){
# #  srat = CreateSeuratObject(srat,min.cells = 3,min.features = minGenes)
# #}else{
# #  srat = CreateSeuratObject(srat)
# #}
# srat = CreateSeuratObject(srat)
# 
# 
# #### =================================================== #
# # Load and store the gene objects
# if(verbose)
#   message("Loading gene meta-data.")
# srat@misc$geneMap = lapply(dataDirs,function(e) {
#   if(file.exists(file.path(e,'features.tsv.gz'))){
#     x = read.table(file.path(e,'features.tsv.gz'),sep='\t',header=FALSE)
#   }else{
#     x = read.table(file.path(e,'genes.tsv'),sep='\t',header=FALSE)
#   }
#   #Do the Seurat gene name conversion.  Assume column 2.
#   x[,2] = make.unique(gsub('_','-',x[,2]))
#   return(x)
# })
# 
# #Check if they're all the same.
# srat@misc$geneMapCommon=TRUE
# for(i in seq_along(dataDirs)){
#   if(i==1){
#     base = srat@misc$geneMap[[i]]
#   }else{
#     tst = srat@misc$geneMap[[i]]
#     #Check they have the same dimensions
#     if(!all(dim(base)==dim(tst))){
#       srat@misc$geneMapCommon=FALSE
#       next
#     }
#     #Check they're all the same
#     for(j in seq(ncol(base)))
#       if(!all(base[,j]==tst[,j]))
#         srat@misc$geneMapCommon=FALSE
#   }
# }
# #If they're all the same, store the same one
# if(srat@misc$geneMapCommon){
#   srat@misc$geneMapCommon = base
#   rownames(srat@misc$geneMapCommon) = base[,2]
# }else{
#   srat@misc$geneMapCommon = NULL
# }
# 
# 
# 
# #### =================================================== #
# srat = NormalizeData(srat,verbose=FALSE)
# # Get cell cycle scores
# data('cc.genes.updated.2019',package='Seurat')
# sPhaseGenes = cc.genes.updated.2019$s.genes
# g2mPhaseGenes = cc.genes.updated.2019$g2m.genes
# #The seed parameter magic is needed because the Seurat authors are dicks
# srat = CellCycleScoring(srat, s.features=sPhaseGenes,g2m.features=g2mPhaseGenes,seed=sample(1e9,1))
# 
# # Get MT content
# if(!is.null(mtPattern)){
#   srat[['percent.mt']] = PercentageFeatureSet(srat, pattern = mtPattern)
# }
# 
# # Get content of any other gene sets
# if(is.null(names(geneSets))){
#   if(length(geneSets)==0){
#     message('No extra gene sets provided.')
#   }else if (length(geneSets) > 0){
#     warning("No names given for gene sets.  Setting to uninformative auto-generated values.")
#     names(geneSets) = paste0('geneSet',seq_along(geneSets))
#   }
# 
# }
# if(verbose)
#   message("Generating extra meta-data.")
# for(nom in names(geneSets)){
#   tmp = geneSets[[nom]]
#   tmp = tmp[tmp %in% rownames(srat@assays$RNA@data)]
#   if(length(tmp)==0)
#     next
#   srat@meta.data[,nom] = colSums(srat@assays$RNA@counts[tmp,,drop=FALSE])/srat@meta.data$nCount_RNA
# }
# 
# 
# #### =================================================== #
# # Call doublets with Scrublet
# if(verbose)
#   message("Loading or running scrublet.")
# srat@meta.data$scrubScore=NA
# srat@meta.data$scrubCall=NA
# srat@meta.data$scrubSim1=NA
# srat@meta.data$scrubSim2=NA
# 
# # Determine where to save scrublet output (with 10X data or locally somewhere)
# if(is.null(cleanCountDir)){
#   scrubDirs = dataDirs
# }else{
#   scrubDirs = file.path(cleanCountDir,names(dataDirs))
# }
# 
# if(skipScrub){
#   scrubRun = rep(FALSE,length(scrubDirs))
# }else{
#   if(!is.null(scrubPath)){
#     scrubRun = !file.exists(file.path(scrubDirs,scrubPath))
#   }else{
#     scrubRun = rep(TRUE,length(scrubDirs))
#   }
# 
#   for(i in seq_along(dataDirs)){
#     if(scrubRun[i]){
#       if(verbose)
#         message(sprintf('Running for channel %s',names(dataDirs)[i]))
# 
#       if(!file.exists(file.path(scrubDirs[i],ifelse(grepl('\\.tsv$',scrubPath),dirname(scrubPath),scrubPath)))){
#         if(verbose)
#           message(sprintf('Creating Clean Count dir for channel %s',names(dataDirs)[i]))
#         dir.create(file.path(scrubDirs[i],ifelse(grepl('\\.tsv$',scrubPath),dirname(scrubPath),scrubPath)),recursive = T)
#       }
# 
#       scrubScores = runScrublet(dataDirs[i],nPCs=numPCs)
#       write_delim(scrubScores,file.path(scrubDirs[i],scrubPath),delim = '\t')
#     }else{
#       scrubScores = read.table(file.path(scrubDirs[i],scrubPath),sep='\t',header=TRUE)
#     }
# 
#     #Save the relevant information in meta.data
#     m = match(paste0(names(dataDirs)[i],'_',scrubScores$barcode),rownames(srat@meta.data))
#     o = !is.na(m)
#     if(sum(o) == 0){
#       m = match(scrubScores$barcode,rownames(srat@meta.data))
#       o = !is.na(m)
#     }
#     if(sum(o) == 0){
#       warning(sprintf('Skipping Scrublet for channel %s as scrublet cell barcodes do not match with seurat object. Please check!',names(dataDirs)[i]))
#     }else{
#       srat@meta.data$scrubScore[m[o]] = scrubScores$score[o]
#       srat@meta.data$scrubCall[m[o]] = scrubScores$call[o]
#       srat@meta.data$scrubSim1[m[o]] = scrubScores$simScore1[o]
#       srat@meta.data$scrubSim2[m[o]] = scrubScores$simScore2[o]
#     }
#   }
# }
# 
# 
# #### =================================================== #
# #Do basic processing and clustering
# if(verbose)
#   message("Basic processing to complete QC.")
# srat = FindVariableFeatures(srat,verbose=FALSE)
# srat = ScaleData(srat,verbose=FALSE)
# srat = RunPCA(srat,npcs=numPCs,approx=FALSE,verbose=FALSE)
# srat = FindNeighbors(srat,dims=seq(numPCs),verbose=FALSE)
# srat = FindClusters(srat,res=clusteringRes,verbose=FALSE)
# srat = RunUMAP(srat,dims=seq(numPCs),verbose=FALSE)
# 
# #### =================================================== #
# #Set the basic filters. Default to TRUE if NA
# #srat@meta.data$PASS_MT = !(srat@meta.data$percent.mt >= maxMT)
# srat@meta.data$PASS_MT = TRUE
# srat@meta.data$PASS_nGenes = !(srat@meta.data$nFeature_RNA <= minGenes)
# srat@meta.data$PASS_nCounts = !(srat@meta.data$nCount_RNA <= minUMIs)
# srat@meta.data$PASS_doublet = !(!is.na(srat@meta.data$scrubCall) & (srat@meta.data$scrubCall | srat@meta.data$scrubScore > scrubScoreMax))
# table(srat@meta.data$PASS_doublet)
# #Work out which ones don't pass because of clustering
# # this gives Percentage of cells in each clusters that did not pass at least one of these criteria
# tt = aggregate(!PASS_MT | !PASS_nGenes | !PASS_nCounts | !PASS_doublet ~ seurat_clusters,FUN=mean,data=srat@meta.data)
# srat@meta.data$fracBadClust = tt[match(srat@meta.data$seurat_clusters,tt[,1]),2]
# # Removing clusters containing too many bad cells
# srat@meta.data$PASS_cluster = !(srat@meta.data$fracBadClust > maxBadFrac)
# #Work out the final filter
# srat@meta.data$PASS = with(srat@meta.data,PASS_MT & PASS_nGenes & PASS_nCounts & PASS_doublet & PASS_cluster)
# #Make a reason for fail variable
# nBasicFail = with(srat@meta.data,4 - (PASS_MT + PASS_nGenes + PASS_nCounts + PASS_doublet))
# tmp = rep('Pass',length(nBasicFail))
# tmp[!srat@meta.data$PASS_MT] = 'highMT'
# tmp[!srat@meta.data$PASS_nGenes] = '#genesLow'
# tmp[!srat@meta.data$PASS_nCounts] = '#countsLow'
# tmp[!srat@meta.data$PASS_doublet] = 'doublet'
# tmp[nBasicFail>1] = 'multiple'
# #Don't want "multiple" to encompass cluster otherwise there'll be heaps where cluster + bad obscures the true reason
# tmp[tmp=='Pass' & !srat@meta.data$PASS_cluster] = 'badCluster'
# srat@meta.data$reasonForFail = tmp
# 
# preSoup_mdat = srat@meta.data
# preSoup_mdat$cellID = rownames(preSoup_mdat)
# 
# #### =================================================== #
# #### Run SoupX to remove ambiant reads/mRNA
# if(!skipSoup){
#   for(channel_directory in dataDirs){
#     srat = Read10X(channel_directory)
#     srat = CreateSeuratObject(srat)
#     srat$cellID = rownames(srat@meta.data)
#     m = preSoup_mdat[preSoup_mdat$cellID %in% rownames(srat@meta.data),]
#     m = m[match(srat$cellID,m$cellID),!colnames(m) %in% colnames(srat@meta.data)]
#     srat@meta.data = cbind(srat@meta.data,m)
#     preSoupSrat = subset(srat,subset = PASS)
#     soupedSrat = cleanWithSoupX(cleanSrat = preSoupSrat,channel_directory,scPath = scPath,is10X = is10X,numPCs_soupX = numPCs_soupX,plotDir = plotDir, cleanCountDir=cleanCountDir,rho_max_limit=rho_max_limit)
#   }
# }
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
# 
# 
# 
# ##---------   Re import all data together  ----------##
# if(verbose)
#   message('Loading data.')
# if(!is.null(dataDirs)){
#   srat = Read10X(dataDirs)
#   if(class(srat) == 'list' & 'Gene Expression' %in% names(srat)){
#     srat = srat[['Gene Expression']]
#   }
# }else{
#   if(!is.null(mtx)){
#     srat = mtx
#   }else{
#     stop('No data provided. Please check!')
#   }
# }
# 
# #if(ncol(srat) > 50000){
# #  srat = CreateSeuratObject(srat,min.cells = 3,min.features = minGenes)
# #}else{
# #  srat = CreateSeuratObject(srat)
# #}
# srat = CreateSeuratObject(srat)
# 
# 
# #### =================================================== #
# # Load and store the gene objects
# if(verbose)
#   message("Loading gene meta-data.")
# srat@misc$geneMap = lapply(dataDirs,function(e) {
#   if(file.exists(file.path(e,'features.tsv.gz'))){
#     x = read.table(file.path(e,'features.tsv.gz'),sep='\t',header=FALSE)
#   }else{
#     x = read.table(file.path(e,'genes.tsv'),sep='\t',header=FALSE)
#   }
#   #Do the Seurat gene name conversion.  Assume column 2.
#   x[,2] = make.unique(gsub('_','-',x[,2]))
#   return(x)
# })
# 
# #Check if they're all the same.
# srat@misc$geneMapCommon=TRUE
# for(i in seq_along(dataDirs)){
#   if(i==1){
#     base = srat@misc$geneMap[[i]]
#   }else{
#     tst = srat@misc$geneMap[[i]]
#     #Check they have the same dimensions
#     if(!all(dim(base)==dim(tst))){
#       srat@misc$geneMapCommon=FALSE
#       next
#     }
#     #Check they're all the same
#     for(j in seq(ncol(base)))
#       if(!all(base[,j]==tst[,j]))
#         srat@misc$geneMapCommon=FALSE
#   }
# }
# #If they're all the same, store the same one
# if(srat@misc$geneMapCommon){
#   srat@misc$geneMapCommon = base
#   rownames(srat@misc$geneMapCommon) = base[,2]
# }else{
#   srat@misc$geneMapCommon = NULL
# }
# 
# 
# 
# #### =================================================== #
# srat = NormalizeData(srat,verbose=FALSE)
# # Get cell cycle scores
# data('cc.genes.updated.2019',package='Seurat')
# sPhaseGenes = cc.genes.updated.2019$s.genes
# g2mPhaseGenes = cc.genes.updated.2019$g2m.genes
# #The seed parameter magic is needed because the Seurat authors are dicks
# srat = CellCycleScoring(srat, s.features=sPhaseGenes,g2m.features=g2mPhaseGenes,seed=sample(1e9,1))
# 
# # Get MT content
# if(!is.null(mtPattern)){
#   srat[['percent.mt']] = PercentageFeatureSet(srat, pattern = mtPattern)
# }
# 
# # Get content of any other gene sets
# if(is.null(names(geneSets))){
#   if(length(geneSets)==0){
#     message('No extra gene sets provided.')
#   }else if (length(geneSets) > 0){
#     warning("No names given for gene sets.  Setting to uninformative auto-generated values.")
#     names(geneSets) = paste0('geneSet',seq_along(geneSets))
#   }
# 
# }
# if(verbose)
#   message("Generating extra meta-data.")
# for(nom in names(geneSets)){
#   tmp = geneSets[[nom]]
#   tmp = tmp[tmp %in% rownames(srat@assays$RNA@data)]
#   if(length(tmp)==0)
#     next
#   srat@meta.data[,nom] = colSums(srat@assays$RNA@counts[tmp,,drop=FALSE])/srat@meta.data$nCount_RNA
# }
# 
# 
# #### =================================================== #
# # Call doublets with Scrublet
# if(verbose)
#   message("Loading or running scrublet.")
# srat@meta.data$scrubScore=NA
# srat@meta.data$scrubCall=NA
# srat@meta.data$scrubSim1=NA
# srat@meta.data$scrubSim2=NA
# 
# # Determine where to save scrublet output (with 10X data or locally somewhere)
# if(is.null(cleanCountDir)){
#   scrubDirs = dataDirs
# }else{
#   scrubDirs = file.path(cleanCountDir,names(dataDirs))
# }
# 
# if(skipScrub){
#   scrubRun = rep(FALSE,length(scrubDirs))
# }else{
#   if(!is.null(scrubPath)){
#     scrubRun = !file.exists(file.path(scrubDirs,scrubPath))
#   }else{
#     scrubRun = rep(TRUE,length(scrubDirs))
#   }
# 
#   for(i in seq_along(dataDirs)){
#     if(scrubRun[i]){
#       if(verbose)
#         message(sprintf('Running for channel %s',names(dataDirs)[i]))
# 
#       if(!file.exists(file.path(scrubDirs[i],ifelse(grepl('\\.tsv$',scrubPath),dirname(scrubPath),scrubPath)))){
#         if(verbose)
#           message(sprintf('Creating Clean Count dir for channel %s',names(dataDirs)[i]))
#         dir.create(file.path(scrubDirs[i],ifelse(grepl('\\.tsv$',scrubPath),dirname(scrubPath),scrubPath)),recursive = T)
#       }
# 
#       scrubScores = runScrublet(dataDirs[i],nPCs=numPCs)
#       write_delim(scrubScores,file.path(scrubDirs[i],scrubPath),delim = '\t')
#     }else{
#       scrubScores = read.table(file.path(scrubDirs[i],scrubPath),sep='\t',header=TRUE)
#     }
# 
#     #Save the relevant information in meta.data
#     m = match(paste0(names(dataDirs)[i],'_',scrubScores$barcode),rownames(srat@meta.data))
#     o = !is.na(m)
#     if(sum(o) == 0){
#       m = match(scrubScores$barcode,rownames(srat@meta.data))
#       o = !is.na(m)
#     }
#     if(sum(o) == 0){
#       warning(sprintf('Skipping Scrublet for channel %s as scrublet cell barcodes do not match with seurat object. Please check!',names(dataDirs)[i]))
#     }else{
#       srat@meta.data$scrubScore[m[o]] = scrubScores$score[o]
#       srat@meta.data$scrubCall[m[o]] = scrubScores$call[o]
#       srat@meta.data$scrubSim1[m[o]] = scrubScores$simScore1[o]
#       srat@meta.data$scrubSim2[m[o]] = scrubScores$simScore2[o]
#     }
#   }
# }
# 
# 
# #### =================================================== #
# #Do basic processing and clustering
# if(verbose)
#   message("Basic processing to complete QC.")
# # srat = FindVariableFeatures(srat,verbose=FALSE)
# # srat = ScaleData(srat,verbose=FALSE)
# # srat = RunPCA(srat,npcs=numPCs,approx=FALSE,verbose=FALSE)
# # srat = FindNeighbors(srat,dims=seq(numPCs),verbose=FALSE)
# # srat = FindClusters(srat,res=clusteringRes,verbose=FALSE)
# # srat = RunUMAP(srat,dims=seq(numPCs),verbose=FALSE)
# #
# # #### =================================================== #
# # #Set the basic filters. Default to TRUE if NA
# # #srat@meta.data$PASS_MT = !(srat@meta.data$percent.mt >= maxMT)
# # srat@meta.data$PASS_MT = TRUE
# # srat@meta.data$PASS_nGenes = !(srat@meta.data$nFeature_RNA <= minGenes)
# # srat@meta.data$PASS_nCounts = !(srat@meta.data$nCount_RNA <= minUMIs)
# # srat@meta.data$PASS_doublet = !(!is.na(srat@meta.data$scrubCall) & (srat@meta.data$scrubCall | srat@meta.data$scrubScore > scrubScoreMax))
# # table(srat@meta.data$PASS_doublet)
# # #Work out which ones don't pass because of clustering
# # # this gives Percentage of cells in each clusters that did not pass at least one of these criteria
# # tt = aggregate(!PASS_MT | !PASS_nGenes | !PASS_nCounts | !PASS_doublet ~ seurat_clusters,FUN=mean,data=srat@meta.data)
# # srat@meta.data$fracBadClust = tt[match(srat@meta.data$seurat_clusters,tt[,1]),2]
# # # Removing clusters containing too many bad cells
# # srat@meta.data$PASS_cluster = !(srat@meta.data$fracBadClust > maxBadFrac)
# # #Work out the final filter
# # srat@meta.data$PASS = with(srat@meta.data,PASS_MT & PASS_nGenes & PASS_nCounts & PASS_doublet & PASS_cluster)
# # #Make a reason for fail variable
# # nBasicFail = with(srat@meta.data,4 - (PASS_MT + PASS_nGenes + PASS_nCounts + PASS_doublet))
# # tmp = rep('Pass',length(nBasicFail))
# # tmp[!srat@meta.data$PASS_MT] = 'highMT'
# # tmp[!srat@meta.data$PASS_nGenes] = '#genesLow'
# # tmp[!srat@meta.data$PASS_nCounts] = '#countsLow'
# # tmp[!srat@meta.data$PASS_doublet] = 'doublet'
# # tmp[nBasicFail>1] = 'multiple'
# # #Don't want "multiple" to encompass cluster otherwise there'll be heaps where cluster + bad obscures the true reason
# # tmp[tmp=='Pass' & !srat@meta.data$PASS_cluster] = 'badCluster'
# # srat@meta.data$reasonForFail = tmp
# 
# preSoup_mdat = preSoup_mdat[match(srat$cellID,preSoup_mdat$cellID),]
# srat@meta.data = cbind(srat@meta.data,preSoup_mdat[,!colnames(preSoup_mdat) %in% colnames(srat@meta.data)])
# 
# #### =================================================== #
# #### Run SoupX to remove ambiant reads/mRNA
# if(!skipSoup){
#   preSoupSrat = subset(srat,subset = PASS)
#   soupedSrat = cleanWithSoupX(cleanSrat = preSoupSrat,dataDirs,scPath = scPath,is10X = is10X,numPCs_soupX = numPCs_soupX,plotDir = plotDir, cleanCountDir=cleanCountDir,rho_max_limit=rho_max_limit)
# 
#   # Add metadata to srat
#   srat@meta.data$PASS_soupX = F
#   m=match(rownames(soupedSrat@meta.data),rownames(srat@meta.data))
#   if(sum(is.na(m))>0){
#     # try renaming these cells
#     #soupedSrat$cellID = paste0(gsub('_\\d+$','',rownames(soupedSrat@meta.data)),'_',soupedSrat$orig.ident)
#     #srat$cellID = paste0(rownames(srat@meta.data),'_',srat$orig.ident)
#     m=match(gsub(rownames(soupedSrat@meta.data)),rownames(srat@meta.data))
#     m=match(soupedSrat$cellID,srat$cellID)
#     print(sum(is.na(m)))
#   }
#   srat@meta.data$PASS_soupX[m] = !(soupedSrat@meta.data$nCount_RNA <= minUMIs)
#   srat@meta.data$preSoupX_nCount=srat@meta.data$nCount_RNA
#   srat@meta.data$postSoupX_nCount = NA
#   srat@meta.data$postSoupX_nCount[m] = soupedSrat@meta.data$nCount_RNA
# 
#   m = which((!srat@meta.data$PASS_soupX) & (srat@meta.data$PASS))
#   srat@meta.data$PASS[m] = F
#   srat@meta.data$reasonForFail[m] = 'lowcount_postSoupX'
# 
# }else{
#   soupedSrat = srat
# }
# 
# # Add Real MT PASS filter
# srat@meta.data$PASS_withMT = srat@meta.data$PASS
# 
# srat@meta.data$PASS_MT = !(srat@meta.data$percent.mt >= maxMT)
# m = which((!srat@meta.data$PASS_MT) & (srat@meta.data$PASS))
# srat@meta.data$PASS[m] = F
# srat@meta.data$reasonForFail[m] = 'highMT'
# 
# #Move PASS column to the end
# w = which(colnames(srat@meta.data)=='PASS')
# srat@meta.data = srat@meta.data[,c(1:(w-1),(w+1):ncol(srat@meta.data),w)]
# 
# 
# # Add metadata if provided
# srat[['cellID']] = rownames(srat@meta.data)
# 
# if(all(unique(srat@meta.data$orig.ident)=='HCA')){
#   srat$orig.ident = paste0(srat$orig.ident,'_',sapply(strsplit(srat$cellID,split = '_'),'[',2))
# }
# 
# if(!is.null(metadata)){
#   if(!is.null(matchBy) & length(matchBy) == 2){
#     m = match(srat[[matchBy[1]]][,1],metadata[[matchBy[2]]])
#     o = !is.na(m)
#   }else{
#     # Default is to match metadata by orig.ident
#     m = match(srat@meta.data$orig.ident,metadata$orig.ident)
#     o = !is.na(m)
#   }
# 
#   if(sum(o)==0){
#     warning('No matching found between seurat object and metadata provided. Please check!')
#   }else{
#     srat@meta.data = merge(srat@meta.data,metadata,by.x=matchBy[1],by.y=matchBy[2],all.x=T)
#     rownames(srat@meta.data) = srat@meta.data$cellID
#   }
# }
# 
# 
# #### =================================================== #
# #### Filter
# if(verbose)
#   message("Filtering and finalizing.")
# goodCells = paste0(rownames(srat@meta.data[srat@meta.data$PASS_withMT,]),'_',srat$orig.ident[srat@meta.data$PASS_withMT])
# w = which(srat@meta.data$PASS_withMT)
# #Make the new version, store the old one in it
# srat@meta.data$preSoupX_nCount = srat@meta.data$nCount_RNA
# sratOld = srat
# #Restart without them
# mDat = srat@meta.data[w,]
# #Drop old clustering
# mDat = mDat[,!colnames(mDat) %in% c('seurat_clusters',grep('^RNA_snn_res',colnames(mDat),value=TRUE))]
# #And drop genes that we don't care about
# genesToKeep = rownames(srat@assays$RNA@counts)
# genesToKeep = genesToKeep[!genesToKeep %in% excludeGenes]
# srat = CreateSeuratObject(soupedSrat@assays$RNA@counts[genesToKeep,colnames(soupedSrat@assays$RNA@counts) %in%
#                                                          rownames(soupedSrat@meta.data)[paste0(gsub('_\\d+$','',rownames(soupedSrat@meta.data)),'_',soupedSrat$orig.ident) %in% goodCells]])
# #Merge in metadata
# m=match(rownames(srat@meta.data),rownames(mDat))
# srat@meta.data = cbind(srat@meta.data,mDat[m,!(colnames(mDat) %in% colnames(srat@meta.data)),drop=FALSE])
# 
# srat = NormalizeData(srat,verbose=FALSE)
# srat = FindVariableFeatures(srat,verbose=FALSE)
# srat@misc$preQC = sratOld
# srat@misc$geneMap = sratOld@misc$geneMap
# srat@misc$geneMapCommon = sratOld@misc$geneMapCommon
# 
# 
# #### =================================================== #
# #Make some QC plots before dropping things
# if(doPlot){
#   if(verbose)
#     message("Making diagnostic plots.")
# 
#   #Plot distributions of basic things
#   if(!is.null(plotDir)){
#     pdf(file.path(paste0(plotDir,'QCplots.pdf')),width = 9.5, height = 8, paper = 'a4r')
#   }
# 
#   df = sratOld@meta.data
#   df.out = df %>% group_by(orig.ident) %>% summarise(nCells_prefilter=n(),medMT=median(percent.mt),muMT=mean(percent.mt),nPASS = sum(PASS),nDoublet = sum(PASS_doublet==FALSE),medCnt = median(nCount_RNA),medGenes=median(nFeature_RNA),nBadClust = sum(PASS_cluster==FALSE),
#                                                      lowCnt = sum(PASS_nCounts==F),highMT = sum(PASS_MT==F),lowGenes = sum(PASS_nGenes==F))
# 
#   ## QC Part 1: How many cells prefiltered? ##
#   counts_preprocess = df %>% group_by(orig.ident) %>% summarise(nCells_prefilter = n())
#   p = ggplot(data=counts_preprocess, aes(x=orig.ident, y=nCells_prefilter)) +
#     geom_bar(stat="identity") +
#     #scale_fill_manual(values = c("#85a040","#8774ca","#ca7040","#4bae8d","#ca5688"))+
#     theme_bw() +
#     theme(text=element_text(size=11), axis.text.x=element_text(angle=45,hjust = 1)) +
#     ggtitle('Prefiltered Counts') + xlab('') + ylab('# cells')
#   plot(p)
#   #ggsave(file="/home/jovyan/Aneuploidy/Dec21/adrenal/prefiltered_counts.pdf", dpi=300, height=5, width=7)
# 
#   #QC Part 2: Mitochondrial Content ##
#   p=ggplot(data=df, aes(x=percent.mt)) +
#     geom_histogram(bins=500,color='black',fill='white') +
#     facet_grid(orig.ident~., scales="free") +
#     #scale_fill_manual(values = c("#85a040","#8774ca","#ca7040","#4bae8d","#ca5688"))+
#     theme_bw() + theme(text=element_text(size=10)) + xlab('')+
#     geom_vline(xintercept=maxMT, linetype="dashed", color="red")+
#     ggtitle('MT content')
#   plot(p)
#   #ggsave(file="/home/jovyan/Aneuploidy/Dec21/adrenal/MT_Distribution.pdf", dpi=300, height=10, width=15)
# 
#   #QC Part 3: nFeature and nUMI Counts ##
#   p=ggplot(data=df, aes(x=log10(nCount_RNA), y=nFeature_RNA, color=percent.mt)) +
#     geom_point(stat="identity",size=0.2) +
#     geom_hline(yintercept = minGenes, linetype="dashed", color="black") +
#     geom_vline(xintercept = log10(minUMIs), linetype="dashed", color="black") +
#     facet_grid(.~orig.ident) +
#     theme_bw() +
#     scale_y_log10()+
#     theme(text=element_text(size=15)) +
#     #geom_hline(yintercept=9000, linetype="dashed", color="red") +
#     scale_color_gradient(high="#cc2b5e", low="grey") +
#     ggtitle("UMI Count vs Features")
#   plot(p)
#   #ggsave(file="/home/jovyan/Aneuploidy/Dec21/adrenal/nFeature_vs_nUMI.pdf", dpi=300, height=10, width=15)
# 
#   p=ggplot(df,aes(x = scrubScore))+
#     geom_histogram(bins=300,color='black',fill='white') +
#     facet_grid(orig.ident~., scales="free") +
#     #scale_fill_manual(values = c("#85a040","#8774ca","#ca7040","#4bae8d","#ca5688"))+
#     theme_bw() + theme(text=element_text(size=10)) + xlab('')+
#     geom_vline(xintercept=scrubScoreMax, linetype="dashed", color="red")+
#     ggtitle('Scrub Score')
#   plot(p)
# 
#   #QC Part 4: How many cells are being removed from each filter? ##
#   nCount_removed <- rownames(df[df$nCount_RNA <= minUMIs,]) #1
#   nFeature_removed <- rownames(df[df$nFeature_RNA <= minGenes,]) #2
#   MT_Content_removed <- rownames(df[df$percent.mt > maxMT,]) #3
#   #badCluster <- rownames(df[df$PASS_cluster==F,]) #4
#   #doublets <- rownames(df[df$PASS_doublet==F,]) #4
# 
#   x <- list(nCount = nCount_removed, nFeature = nFeature_removed, MT_Content = MT_Content_removed)#,
#   #badClusters = badCluster, doublets = doublets)
# 
#   p=ggVennDiagram(x, label_size=4, color="black") + scale_fill_gradient(low = "white", high = "white")
#   plot(p)
# 
#   # barplot of #cells removed by each fileter per channel
#   filter_by_channel = as.data.frame(table(sratOld$orig.ident,sratOld$reasonForFail))
#   colnames(filter_by_channel) = c('orig.ident','reasonForFail','Freq')
# 
#   p = ggplot(data=filter_by_channel, aes(x=orig.ident, y=Freq,fill=reasonForFail)) +
#     geom_bar(stat="identity",position = 'stack') +
#     scale_fill_manual(values = c(brewer.pal(7,'Paired'),'grey'))+
#     theme_bw() +
#     theme(text=element_text(size=11), axis.text.x=element_text(angle=45,hjust = 1)) +
#     ggtitle('Reason for fail') + xlab('') + ylab('# cells')
#   plot(p)
#   # Normalized against total number of cells version
#   filter_by_channel = filter_by_channel %>% group_by(orig.ident) %>% mutate(total_orig_cells = sum(Freq))
#   p=ggplot(data=filter_by_channel, aes(x=orig.ident, y=Freq/total_orig_cells,fill=reasonForFail)) +
#     geom_bar(stat="identity",position = 'stack') +
#     scale_fill_manual(values = c(brewer.pal(7,'Paired'),'grey'))+
#     theme_bw() +
#     theme(text=element_text(size=11), axis.text.x=element_text(angle=45,hjust = 1)) +
#     ggtitle('Reason for fail',subtitle = 'Normalized against total # cells') + xlab('') + ylab('# cells')
#   plot(p)
# 
# 
#   ## QC Part 5: Plot UMI count before (post Scrublet) and after SoupX ##
#   df2 = sratOld@meta.data
#   p=ggplot(data=df2, aes(x=log10(nCount_RNA), y=log10(preSoupX_nCount))) +
#     geom_point(stat="identity",size=0.01) +
#     geom_hline(yintercept = log10(minUMIs), linetype="dashed", color="red",linewidth=0.1) +
#     geom_vline(xintercept = log10(minUMIs), linetype="dashed", color="red",linewidth=0.1) +
#     facet_grid(.~orig.ident) +
#     theme_bw() +
#     theme(text=element_text(size=15)) +
#     #geom_hline(yintercept=9000, linetype="dashed", color="red") +
#     scale_color_gradient(high="#cc2b5e", low="grey") +
#     ggtitle("UMI count pre/post SoupX")
#   plot(p)
# 
#   ## QC Part 6: How many cells pre vs post-filtered, Scrublet and SoupX? ##
#   counts_postprocess = srat@meta.data %>% group_by(orig.ident) %>% summarise(nCells_postfilter = n())
#   counts_postprocess = merge(counts_postprocess,counts_preprocess,by='orig.ident',all.y=T)
#   counts_postprocess = pivot_longer(counts_postprocess,cols = c(2:3),names_to = 'type',values_to = 'nCells')
#   counts_postprocess$type = factor(counts_postprocess$type,levels = c('nCells_prefilter','nCells_postfilter'))
#   p = ggplot(data=counts_postprocess, aes(x=orig.ident, y=nCells,fill=type)) +
#     geom_bar(stat="identity",position = 'dodge') +
#     scale_fill_manual(values = c("darkgrey","black"))+
#     theme_bw() +
#     theme(text=element_text(size=11), axis.text.x=element_text(angle=45,hjust = 1)) +
#     ggtitle('Postfiltered Cell Counts') + xlab('') + ylab('# cells')
#   plot(p)
# 
#   #Standard QC things
#   plot(DimPlot(sratOld,group.by='seurat_clusters',label=TRUE,repel=TRUE)+guides(colour="none"))
# 
#   #A hack for another of Seurat's stupid decisions
#   if('labels' %in% names(as.list(args(DimPlot)))){
#     plot(DimPlot(sratOld,group.by=c('Phase','orig.ident'),labels=c('Phase','Source')))
#   }else{
#     gg = DimPlot(sratOld,group.by=c('Phase','orig.ident'),combine=FALSE)#,cols = c(brewer.pal(12,'Paired'),brewer.pal(12,'Dark2'),brewer.pal(12,'Set1')))
#     labs = c('Phase','Source')
#     for(i in seq_along(labs)){
#       gg[[i]] = gg[[i]]  + ggtitle(labs[i])
#       plot(gg[[i]])
#     }
#   }
#   plot(FeaturePlot(sratOld,c('percent.mt','nFeature_RNA','nCount_RNA','fracBadClust')))
#   #Which cells pass based on what
#   if('labels' %in% names(as.list(args(DimPlot)))){
#     plot(DimPlot(sratOld,group.by=c('PASS_MT','PASS_nGenes','PASS_nCounts','PASS_doublet','PASS_cluster','PASS_soupX','PASS'),labels=c('Pass MT Frac?','Pass #Genes?','Pass #Counts?','Pass Doublet?','Pass cluster','Pass soupX?','PASS'),legend='none'))
#   }else{
#     gg = DimPlot(sratOld,group.by=c('PASS_MT','PASS_nGenes','PASS_nCounts','PASS_doublet','PASS_cluster','PASS'),combine=FALSE)
#     labs =c('Pass MT Frac?','Pass #Genes?','Pass #Counts?','Pass Doublet?','Pass cluster?','PASS?')
#     for(i in seq_along(labs)){
#       gg[[i]] = gg[[i]] +guides(colour="none") + ggtitle(labs[i])
#     }
#     plot(patchwork::wrap_plots(gg))
#   }
#   plot(DimPlot(sratOld,group.by=c('reasonForFail'),cols = c(brewer.pal(9,'Paired')[-c(1,6)],'lightgrey')))
# 
#   if(!is.null(plotDir)){
#     dev.off()
#   }
# }
# 
# 
# if(!is.null(outPath)){
#   if(keepMTCells){
#     # Save the results with all highMT cells
#     message(sprintf('Saving cleanSrat with hightMT cells to %s',paste0(outPath,'_clean_withMTCells.RDS')))
#     saveRDS(srat,paste0(outPath,'_clean_withMTCells.RDS'))
#   }
# 
#   # Save the results without highMT cells
#   message(sprintf('Saving cleanSrat to %s',paste0(outPath,'_clean_noMTCells.RDS')))
#   print(dim(srat))
#   print(table(srat$PASS))
#   srat = subset(srat, subset = PASS)
#   print(dim(srat))
#   print(table(srat$PASS))
#   # srat = NormalizeData(srat,verbose=FALSE)
#   # srat = FindVariableFeatures(srat,verbose=FALSE)
#   srat = standard_clustering(srat)
#   saveRDS(srat,paste0(outPath,'_clean_noMTCells.RDS'))
# }
# 
# 
# 
# 
