##----   Processing leukaemia scRNAseq datasets    -----##

outDir = "Results/02_MLDS_scQC"

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}



##---------------##
#### Libraries ####
##---------------##
library(tidyverse)
#source("R/utils/misc.R")
source("R/utils/sc_basicQC.R")



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
metadata = NULL
matchBy = NULL
cleanCountDir = file.path(sub_outDir,'cleanCount')



cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/MLDS_clean_withMTCells.RDS'),paste0(sub_outDir,'/MLDS_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)
}else{
  if(!dir.exists(sub_outDir)){
    message(sprintf('Creating output directory'))
    dir.create(sub_outDir,recursive = T)
  }
  
  #setwd(sub_outDir)
  
  # List location of cellranger outputs
  dataDirs = c(list.dirs('Data/MLDS_scRNAseq/',full.names = T,recursive = T))
  dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
  dataDirs = dataDirs[!grepl('GRCh38-1_2_0',dataDirs)]
  
  
  names(dataDirs) = basename(dirname(dataDirs))
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_ALeuk_','',names(dataDirs)))
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',names(dataDirs)))
  names(dataDirs) = gsub('_','.',gsub('_GRCh38-2020-A.*$','',gsub('^.*_MY_','MY_',names(dataDirs))))
  
  ## Remove some low quality / not relevant samples 
  dataDirs = dataDirs[!names(dataDirs) %in% c("Leuk13234200", "Leuk13234201", "Leuk13234202", # L041_D
                                              "NB14406184", "NB14406185")] # L062 - LPD samples
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
  print(n_distinct(QC.output[[1]]$orig.ident))
}