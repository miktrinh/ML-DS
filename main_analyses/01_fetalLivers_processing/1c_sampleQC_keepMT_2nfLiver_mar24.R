# This script was extracted from 1_sampleQC_keepMT.R, orignially written in Aug22

# Dec2023 QC 10X 2n fetal liver data
# We decided to do scRNAseq on a new 2n fetal liver data for validation. Need to process this

setwd('~/lustre_mt22/Aneuploidy/')

#############
# Libraries #
#############
# Load libraries
library(tidyverse)
library(Seurat)
library(readxl)
library(SoupX)
library(RColorBrewer)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/sc_basicQC.R")


#############
#   Params  #
#############

maxMT = 30
minGenes = 300
minUMIs=500

maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup = F
scrubScoreMax = 0.5
scrubPath='../cleanCounts/scrubletScores.tsv'
scPath="../cleanCounts/strainedCounts"
doPlot=T
verbose = T
skipIfExists=T
keepMTCells=T
excludeGenes = c()
tissue = 'liver'

outDir = file.path('~/lustre_mt22/Aneuploidy/Results/1_sc_processing',tissue,'jun24') 
if(!dir.exists(outDir)){
  message(sprintf('Creating output directory for %s',tissue))
  dir.create(outDir,recursive = T)
}

setwd(outDir)

dataDirs = list.dirs(c('~/lustre_mt22/Aneuploidy/Data/Hsb38/',
                       '~/lustre_mt22/Aneuploidy/Data/Hsb39/',
                       '~/lustre_mt22/Aneuploidy/Data/Hsb40/'),recursive = T,full.names = T)
dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
names(dataDirs) = gsub('_GRCh.*$','',gsub('.*_MY','MY',dataDirs))
names(dataDirs) = gsub('_','.',names(dataDirs))
dataDirs=dataDirs[file.exists(dataDirs)]
print(n_distinct(dataDirs))

# Run basicQC
# plotDir = file.path(outDir,paste0(tissue,'_kidney_'))
# outPath = file.path(outDir,paste0(tissue,'_kidney'))
# plotDir = file.path(outDir,paste0(tissue,'_adrenal_'))
# outPath = file.path(outDir,paste0(tissue,'_adrenal'))
plotDir = file.path(outDir,paste0(tissue,'_'))
outPath = file.path(outDir,paste0(tissue,'_'))

metadata = NULL
matchBy = NULL

QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                    clusteringRes=clusteringRes,
                    skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                    metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                    doPlot=doPlot,plotDir=plotDir,verbose=verbose)

cleanSrat = QC.output[[1]]

df.out = QC.output[[2]]
qc.summary=rbind(qc.summary,df.out)

write.csv(df.out,paste0('~/lustre_mt22/Aneuploidy/Results/1_sc_processing/',tissue,'/jun24/',tissue,'_qc_summary.csv'))









