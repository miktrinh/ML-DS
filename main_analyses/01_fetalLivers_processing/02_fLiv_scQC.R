##----   Processing fetal Livers scRNAseq datasets    -----##

outDir = "Results/02_fLiv_scQC"

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
tissue = 'liver'



##----------------------##
##     fetal Livers   ####
##----------------------##




# for(genomeVersion in c('GRCh38-2020-A','GRCh38-1_2_0')){
#     sub_outDir = file.path(outDir,'fLiv',genomeVersion)
#     plotDir = file.path(sub_outDir,'fLiv_')
#     outPath = file.path(sub_outDir,'fLiv')
#     metadata = NULL
#     matchBy = NULL

    
#     cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/fLiv_clean_withMTCells.RDS'),paste0(sub_outDir,'/fLiv_clean_noMTCells.RDS'))
    
#     if(file.exists(cleanSrat_fp) & skipIfExists){
#       cleanSrat = readRDS(cleanSrat_fp)
#     }else{
#       if(!dir.exists(sub_outDir)){
#         message(sprintf('Creating output directory'))
#         dir.create(sub_outDir,recursive = T)
#       }
      
#       #setwd(sub_outDir)
      
#       # List location of cellranger outputs
#       dataDirs = c(list.dirs('Data/fLivers_scRNAseq/',full.names = T,recursive = T))
#       dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
#       dataDirs = dataDirs[grepl(genomeVersion,dataDirs)]
      
#       names(dataDirs) = gsub('_GRCh.*$','',gsub('.*_MY','MY',dataDirs))
#       names(dataDirs) = gsub('_','.',names(dataDirs))
      
#       length(names(dataDirs))
#       dataDirs=dataDirs[file.exists(dataDirs)]
#       dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
#       print(n_distinct(dataDirs))
      
#       cleanCountDir = dataDirs
#       # Run basicQC
#       message('\nPerforming scRNAseq QC...')
#       QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
#                           clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
#                           skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
#                           metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
#                           doPlot=doPlot,plotDir=plotDir,verbose=verbose)
#       print(n_distinct(QC.output[[1]]$orig.ident))
#     }
# }


## Merge the two datasets
cleanSrat_fps = c()
if(keepMTCells){
    for(genomeVersion in c('GRCh38-2020-A','GRCh38-1_2_0')){
          cleanSrat_fps = c(cleanSrat_fps,
                            file.path(outDir,'fLiv',genomeVersion,'fLiv_clean_withMTCells.RDS'))
    }
}else{
  for(genomeVersion in c('GRCh38-2020-A','GRCh38-1_2_0')){
          cleanSrat_fps = c(cleanSrat_fps,
                            file.path(outDir,'fLiv',genomeVersion,'fLiv_clean_noMTCells.RDS'))
    }
}

fLiver_1 = readRDS(cleanSrat_fps[1])
fLiver_1@misc = fLiver_1@misc[names(fLiver_1@misc) != 'preQC']
fLiver_2 = readRDS(cleanSrat_fps[2])
fLiver_2@misc = fLiver_2@misc[names(fLiver_2@misc) != 'preQC']

fLiver = merge_seurat_objects(fLiver_1, fLiver_2, keepAllGenes=F, genomeVersions = c('v38','v38'))
fLiver = standard_clustering(fLiver)

out_fp = ifelse(keepMTCells,
                file.path(outDir,'fLiv','fLiver_clean_withMTCells.RDS'),
                file.path(outDir,'fLiv','fLiver_clean_noMTCells.RDS'))
saveRDS(fLiver,out_fp)

