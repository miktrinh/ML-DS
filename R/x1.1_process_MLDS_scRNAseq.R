##----   Processing leukaemia scRNAseq datasets    -----##

outDir = "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/cleanSrat_scRNAseq/Results/1_scProcessing_leukaemia"

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
doPlot=F
verbose = T
skipIfExists=T
keepMTCells=T



##----------------------##
##     TAM / ML-DS    ####
##----------------------##

sub_outDir = file.path(outDir,'cleanSrat')
plotDir = file.path(sub_outDir,'cleanSrat_')
outPath = file.path(sub_outDir,'cleanSrat')
metadata = NULL
matchBy = NULL
cleanCountDir = file.path(sub_outDir,'cleanCount')



cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/cleanSrat_clean_withMTCells.RDS'),paste0(sub_outDir,'/cleanSrat_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)
}else{
  if(!dir.exists(sub_outDir)){
    message(sprintf('Creating output directory'))
    dir.create(sub_outDir,recursive = T)
  }

  setwd(sub_outDir)

  # List location of cellranger outputs
  dataDirs = c(list.dirs('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/cleanSrat_GEX',full.names = T,recursive = T),
               list.dirs('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/cleanSrat_scRNAseq/cleanSrat_scRNAseq_Data/',recursive = T,full.names = T))
  dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
  dataDirs = dataDirs[!grepl('GRCh38-1_2_0',dataDirs)]
  
  
  names(dataDirs) = basename(dirname(dataDirs))
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_ALeuk_','',names(dataDirs)))
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',names(dataDirs)))
  names(dataDirs) = gsub('_','.',gsub('_GRCh38-2020-A.*$','',gsub('^.*_MY_','MY_',names(dataDirs))))
  
  ## Remove L062 - LPD samples
  dataDirs = dataDirs[!names(dataDirs) %in% c("Leuk13234200", "Leuk13234201", "Leuk13234202", "NB14406184",   "NB14406185")]
  length(names(dataDirs))


  dataDirs=dataDirs[file.exists(dataDirs)]

  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))
  
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose)

  cleanSrat = QC.output[[1]]

  # df.out = QC.output[[2]]
  # 
  # write.csv(df.out,paste0(sub_outDir,'/cleanSrat_240422_qc_summary.csv'))
  # cleanSrat = standard_clustering(cleanSrat)
  # saveRDS(cleanSrat,file.path(sub_outDir,'cleanSrat_clean_noMTcells_2404.RDS'))
  
}















# ##----------------##
# ##     B-ALL    ####
# ##----------------##
# numPCs = 50
# sub_outDir = file.path(outDir,'BALL')
# plotDir = file.path(sub_outDir,'BALL_')
# outPath = file.path(sub_outDir,'BALL')
# 
# cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/BALL_clean_withMTCells.RDS'),paste0(sub_outDir,'/BALL_clean_noMTCells.RDS'))
# if(file.exists(cleanSrat_fp) & skipIfExists){
#   cleanSrat = readRDS(cleanSrat_fp)
# }else{
#   if(!dir.exists(sub_outDir)){
#     message(sprintf('Creating output directory'))
#     dir.create(sub_outDir,recursive = T)
#   }
# 
#   setwd(sub_outDir)
# 
#   # List location of cellranger outputs
#   sampleSheet = read_excel('~/lustre_mt22/Data/B_ALL/Manifest_2023_09_13.xlsx',sheet = 'B-ALL (Study 6918_7153_3030)')
#   sampleSheet = sampleSheet[sampleSheet$Experiment == 'GEX' &
#                               grepl('BALL',sampleSheet$Sample_ID),]
#   sampleSheet = sampleSheet[!sampleSheet$Cellgen_ID %in% c(sampleSheet$Cellgen_ID[sampleSheet$Patient_ID %in% c('L007','L071','L024','L031','L096','L057','L045','L073','L033','L059','L043','L072','L023','L087','L085','L032','L004','L006','L025','L048','L003') & sampleSheet$Timepoint == 'D0' | sampleSheet$Subtype == 'Other'],'SB_Leuk13652352'),]
# 
#   dataDirs = list.dirs('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/GEX/',recursive = T,full.names = T)
#   dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
#   names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',dataDirs))
# 
#   dataDirs=dataDirs[file.exists(dataDirs)]
#   dataDirs = dataDirs[names(dataDirs) %in% gsub('SB_','',sampleSheet$Cellgen_ID)]
# 
#   # Find brain Lymphoma samples
#   brainLPD_dirs = list.dirs('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/cleanSrat_scRNAseq/cleanSrat_scRNAseq_Data/',recursive = T,full.names = T)
#   brainLPD_dirs = brainLPD_dirs[grepl('filtered_feature_bc_matrix$',brainLPD_dirs) & grepl('SB_NB14406184|SB_NB14406185',brainLPD_dirs)]
#   names(brainLPD_dirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',brainLPD_dirs))
# 
#   dataDirs = c(dataDirs,brainLPD_dirs)
#   dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
#   print(n_distinct(dataDirs))
# 
# 
#   metadata = NULL
#   matchBy = NULL
#   cleanCountDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/cleanSrat_scRNAseq/Results/1_scProcessing_cleanSrat/BALL/cleanCount'
#   # Run basicQC
#   plotDir = file.path(sub_outDir,'BALL_')
#   outPath = file.path(sub_outDir,'BALL')
# 
#   message('\nPerforming scRNAseq QC...')
#   QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
#                       clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
#                       skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
#                       metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
#                       doPlot=doPlot,plotDir=plotDir,verbose=verbose)
# 
#   cleanSrat = QC.output[[1]]
# 
#   cleanSrat = standard_clustering(cleanSrat)
#   saveRDS(cleanSrat,'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/cleanSrat_scRNAseq/Results/1_scProcessing_cleanSrat/BALL/BALL_clean_noMTcells_2404.RDS')
# }








# ##---------------------##
# ##     MDS + pAML    ####
# ##---------------------##
# numPCs = 75
# sub_outDir = file.path(outDir,'pAML')
# plotDir = file.path(sub_outDir,'pAML_')
# outPath = file.path(sub_outDir,'pAML')
# 
# cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/pAML_clean_withMTCells.RDS'),paste0(sub_outDir,'/pAML_clean_noMTCells.RDS'))
# if(file.exists(cleanSrat_fp) & skipIfExists){
#   cleanSrat = readRDS(cleanSrat_fp)
# }else{
#   if(!dir.exists(sub_outDir)){
#     message(sprintf('Creating output directory'))
#     dir.create(sub_outDir,recursive = T)
#   }
#   
#   setwd(sub_outDir)
#   
#   # List location of cellranger outputs
#   sampleSheet = read_excel('~/lustre_mt22/Data/B_ALL/Manifest_2023_09_13.xlsx',sheet = 'B-ALL (Study 6918_7153_3030)')
#   sampleSheet = sampleSheet[sampleSheet$Experiment == 'GEX' &
#                               grepl('AML|MDS',sampleSheet$Sample_ID),]
#   
#   dataDirs = list.dirs('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/GEX/',recursive = T,full.names = T)
#   dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
#   names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_|^.*_ALeuk_','',dataDirs))
#   
#   dataDirs=dataDirs[file.exists(dataDirs)]
#   dataDirs = dataDirs[names(dataDirs) %in% gsub('SB_','',sampleSheet$Cellgen_ID)]
#   
#   dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
#   print(n_distinct(dataDirs))
#   
#   
#   
#   metadata = NULL
#   matchBy = NULL
#   cleanCountDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/cleanSrat_scRNAseq/Results/1_scProcessing_cleanSrat/pAML/cleanCount'
#   # Run basicQC
#   plotDir = file.path(sub_outDir,'pAML_')
#   outPath = file.path(sub_outDir,'pAML')
#   
#   message('\nPerforming scRNAseq QC...')
#   QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
#                       clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
#                       skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
#                       metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
#                       doPlot=doPlot,plotDir=plotDir,verbose=verbose)
#   
#   cleanSrat = QC.output[[1]]
#   
#   cleanSrat = standard_clustering(cleanSrat)
#   saveRDS(cleanSrat,'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/cleanSrat_scRNAseq/Results/1_scProcessing_cleanSrat/pAML/pAML.MDS_clean_noMTcells_2404.RDS')
# }
