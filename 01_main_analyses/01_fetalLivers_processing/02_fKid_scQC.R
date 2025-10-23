##----   Processing fetal Livers scRNAseq datasets    -----##

outDir = "Results/02_fKid_scQC"

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
minUMIs=500
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
keepMTCells=F
excludeGenes = c()
tissue = 'kidney'



##----------------------##
##    fetal Kidneys   ####
##----------------------##



sub_outDir = file.path(outDir,'fKid')
plotDir = file.path(sub_outDir,'fKid_')
outPath = file.path(sub_outDir,'fKid')
metadata = NULL
matchBy = NULL


cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/fKid_clean_withMTCells.RDS'),paste0(sub_outDir,'/fKid_clean_noMTCells.RDS'))

if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)
}else{
  if(!dir.exists(sub_outDir)){
    message(sprintf('Creating output directory'))
    dir.create(sub_outDir,recursive = T)
  }

  #setwd(sub_outDir)

  # List location of cellranger outputs
  dataDirs = c(list.dirs('Data/fKidneys_scRNAseq/',full.names = T,recursive = T))
  dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
  #dataDirs = dataDirs[grepl(genomeVersion,dataDirs)]

  names(dataDirs) = gsub('_GRCh.*$','',gsub('.*_MY','MY',dataDirs))
  names(dataDirs) = gsub('_','.',names(dataDirs))

  length(names(dataDirs))
  dataDirs=dataDirs[file.exists(dataDirs)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))

  cleanCountDir = dataDirs
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose)
  print(n_distinct(QC.output[[1]]$orig.ident))
}

