# Processing of MLDS single cell RNAseq data

# Location of the samples
#/lustre/scratch119/casm/team274sb/bl10/B-ALL/Data/GEX/cellranger700_count_45842_SB_Leuk13104278_GRCh38-2020-A
#/lustre/scratch119/casm/team274sb/bl10/B-ALL/Data/GEX/cellranger700_count_45842_SB_Leuk13104279_GRCh38-2020-A

# sudo apt-get install bcftools
# sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
# sudo chmod +x /usr/local/bin/alleleCounter

setwd("/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/")

#############
# Libraries #
#############
# Load libraries
library(tidyverse)
library(alleleIntegrator)
library(readxl)
library(Seurat)
library(SoupX)
library(RColorBrewer)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")



##############################################################################
# Genotyping to check if the samples are from the same or different patients #
##############################################################################

#---------------------------------------------------------------------------------------------------------------------------------


# Install relevant packages
#install.packages('/nfs/users/nfs_m/my4/alleleIntegrator_0.7.3.tar.gz',repos = NULL,type='source')
#BiocManager::install("VariantAnnotation")
#BiocManager::install('SNPRelate')

#sudo apt-get install bcftools
#sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
#sudo chmod +x /usr/local/bin/alleleCounter

# We don't have DNA data so just look at RNA only

# refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
# liftChain = '~/lustre_mt22/hg19ToHg38_noChr.over.chain'
# nParallel=24
# 
# 
# outDir = file.path('MLDS_scRNAseq/genotypeCheck')
# if(!dir.exists(outDir)){
#   dir.create(outDir,recursive = T)
# }


#### Get list of RNA and DNA BAM files ####

#----- RNA BAMs
# bams10X = list.files('/lustre/scratch119/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX',full.names = T,recursive = T,pattern = 'possorted_genome_bam.bam$')
# 
# names(bams10X) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',bams10X))
# bams10X = bams10X[file.exists(bams10X)]
# 
# #----- DNA BAMs 
# dnaBAMs = list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3030',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)
# names(dnaBAMs) = gsub('\\..*$','',basename(dnaBAMs))
# if(length(dnaBAMs) > 0){
#   dnaBAMs = dnaBAMs[file.exists(dnaBAMs)]
# }
# 
# # Check that each dnaBAM file has a unique name
# if(length(unique(names(dnaBAMs))) != length(dnaBAMs)){
#   stop(sprintf('Duplicated bams10X names detected: %s',names(dnaBAMs)[duplicated(names(dnaBAMs))]))
# }
# 
# 
# #############################
# # Check genotype consistency
# #Are all the BAMs you're going to use from the same individual?  Check before you start
# genoCheck = matchBAMs(BAMs = BAMs,
#                       refGenomes = rep(c(refGenome,refGenome10X),c(length(dnaBAMs),length(bams10X))),
#                       outputs = file.path(outDir,paste0(c(names(dnaBAMs),names(bams10X)),'_genotypeCheck.tsv')),
#                       liftOvers=rep(c(NA,liftChain),c(length(bams10X))),
#                       is10X=rep(c(FALSE,TRUE),c(length(dnaBAMs),length(bams10X))),
#                       nParallel=nParallel,nMaxSim=10,nChunks=4,skipIfExists=skipIfExists)
# 
# 
# 
# #If anything is less than 0.8 and you should be concerned...
# message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))












##################################################################################################################################


###############################
# Preprocessing scRNAseq data #
###############################



#---------------------------------------------------------------------------------------------------------------------------------


source("~/lustre_mt22/generalScripts/utils/sc_basicQC.R")

message('\n--------  Hello! --------\n')

##############
##  Params  ##
##############
maxMT = 30
minGenes = 300
minUMIs = 500  
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

# 
# projMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'MLDS_GOSH')
# projMani = projMani[!is.na(projMani$Disease),]
# projMani = projMani[!is.na(projMani$dataLocation),]
# projMani = projMani[grepl('B-ALL',projMani$Disease),]
# projMani = projMani[projMani$assay == "5' V2 Dual Index",]


### 1. Import remapped data
###    Run SoupX
###    Subset to keep only cells present in the original publications
###    Add cell labels (as published)

outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/sept23/ek12_infantALL') 
plotDir = file.path(outDir,'infantALL_')
outPath = file.path(outDir,'infantALL')

cleanSrat_fp = ifelse(keepMTCells,paste0(outDir,'/infantALL_clean_withMTCells.RDS'),paste0(outDir,'/infantALL_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)  
}else{
  if(!dir.exists(outDir)){
    message(sprintf('Creating output directory'))
    dir.create(outDir,recursive = T)
  }
  
  setwd(outDir)
  
  dataDirs = list.dirs('/lustre/scratch125/casm/team274sb/mt22/Data/ek12_infantALL_paper/',recursive = T,full.names = T)
  dataDirs = dataDirs[grepl('filtered_gene_bc_matrices$|filtered_feature_bc_matrix',dataDirs)]
  
  names(dataDirs) = ifelse(grepl('_SB_',dataDirs),gsub('^.*_SB_','',gsub('_GRCh38-1_2_0.*$','',dataDirs)),
                           gsub('^.*_','',gsub('_GRCh38-1_2_0.*$','',dataDirs)))
  
  dataDirs=dataDirs[file.exists(dataDirs)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))
  
  metadata = NULL
  matchBy = NULL
  cleanCountDir = NULL
  # Run basicQC
  plotDir = file.path(outDir,'infantALL_')
  outPath = file.path(outDir,'infantALL')
  
  
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose)
  
  cleanSrat = QC.output[[1]]
  
  df.out = QC.output[[2]]
  
  write.csv(df.out,paste0(outDir,'/infantALL_0923_qc_summary.csv'))
}








##################################################################################################################################


#############################
##       Annotation      ####
#############################
old_infant = readRDS('~/lustre_mt22/Aneuploidy/ek12_infantALL_scRNAseq/ek12_8iALL_1iAMKL.RDS')
old_infant$percent.mt = PercentageFeatureSet(old_infant,pattern = '^MT-')
old_infant_mdat = old_infant@meta.data
old_infant_mdat$cellID_bc = gsub('/.*$','',old_infant_mdat$cellID)

infantALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/sept23/ek12_infantALL/infantALL_clean_withMTCells.RDS')
infantALL_preQC = infantALL@misc$preQC
infantALL_preQC$cellID_bc = gsub('.*_','',infantALL_preQC$cellID)
infantALL_preQC$cellID_bc[infantALL_preQC$orig.ident %in% c('4602STDY7920961','4602STDY7920960')] = old_infant_mdat$cellID[match(infantALL_preQC$cellID_bc[infantALL_preQC$orig.ident %in% c('4602STDY7920961','4602STDY7920960')],old_infant_mdat$cellID_bc[grepl('P1_iALL|P2_iALL',old_infant_mdat$patient_cancer)])]
infantALL_preQC$cellID_bc[is.na(infantALL_preQC$cellID_bc)] = gsub('.*_','',infantALL_preQC$cellID[is.na(infantALL_preQC$cellID_bc)])

infantALL_preQC$donorID = gsub('.*/','',infantALL_preQC$cellID_bc)
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('4602STDY7920965')] = 'P9_iAML'
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('4602STDY7920966')] = 'P9_iAML_TP1'
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('BON10188628')] = 'P5_iALL'
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('BON10188629')] = 'P6_iALL'
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('BON10188630')] = 'P8_iALL_ETV6'
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('BON10188631')] = 'P7_iALL_NUTM1'
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('NB8791864')] = 'P3_iALL' # only NB8791864 was used in the paper
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('NB8791865')] = 'P3_iALL_maybe' # only NB8791864 was used in the paper

infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('NB8791866')] = 'P4_iALL_run1'
infantALL_preQC$donorID[infantALL_preQC$orig.ident %in% c('NB8791867')] = 'P4_iALL_run2'
#infantALL_preQC$donorID[!infantALL_preQC$orig.ident %in% c('4602STDY7920961','4602STDY7920960')] = as.character(infantALL$donorID[match(infantALL_preQC$orig.ident[!infantALL_preQC$orig.ident %in% c('4602STDY7920961','4602STDY7920960')],infantALL$orig.ident)])

infantALL_preQC$cellID2 = infantALL_preQC$cellID_bc
infantALL_preQC$cellID2[!infantALL_preQC$orig.ident %in% c('4602STDY7920961','4602STDY7920960')] = paste0(infantALL_preQC$cellID_bc[!infantALL_preQC$orig.ident %in% c('4602STDY7920961','4602STDY7920960')],'/',
                                                                                                          infantALL_preQC$donorID[!infantALL_preQC$orig.ident %in% c('4602STDY7920961','4602STDY7920960')])

# Check that all published cells are included
table(old_infant_mdat$cellID %in% infantALL_preQC$cellID2)
table(duplicated(infantALL_preQC$cellID2))


# Remove low QC cells (but included those published cells)
infantALL = subset(infantALL_preQC,subset = cellID %in% infantALL_preQC$cellID[infantALL_preQC$cellID2 %in% old_infant_mdat$cellID | infantALL_preQC$PASS == T])
infantALL = standard_clustering(infantALL)
infantALL$tissue = 'pBM'

## Annotation
infantALL$patient_cancer = as.character(old_infant_mdat$patient_cancer[match(infantALL$cellID2,old_infant_mdat$cellID)])
infantALL$patient_cancer[is.na(infantALL$patient_cancer) & infantALL$donorID == 'P9_iAML_TP1'] = 'P9_iAML_TP1'
infantALL$patient_cancer[is.na(infantALL$patient_cancer) ] = 'NA'
DimPlot(infantALL,group.by = 'patient_cancer',label = T,repel = T,label.box = T,cols = col25) + NoLegend()
DimPlot(infantALL,group.by = 'seurat_clusters',label = T,repel = T,label.box = T) + NoLegend()
infantALL$donorID[grepl('-',infantALL$donorID) & infantALL$seurat_clusters %in% c(6) & infantALL$orig.ident %in% c('4602STDY7920961','4602STDY7920960')] = 'P1_iALLM'
infantALL$donorID[grepl('-',infantALL$donorID) & infantALL$seurat_clusters %in% c(13) & infantALL$orig.ident %in% c('4602STDY7920961','4602STDY7920960')] = 'P2_iALL'
infantALL$donorID[grepl('-',infantALL$donorID)] = 'unknown'

DimPlot(infantALL,cells.highlight = infantALL$cellID[infantALL$seurat_clusters %in% c(9) & infantALL$annot == 'NA'])
infantALL$annot = gsub('.*/','',infantALL$patient_cancer)
infantALL$annot[infantALL$seurat_clusters %in% c(26,30)] = 'EE'
infantALL$annot[infantALL$seurat_clusters %in% c(18)] = 'MEMP_MEP'
infantALL$annot[infantALL$seurat_clusters %in% c(9) & infantALL$annot == 'P9_iAML_TP1'] = 'MK'

FeaturePlot(infantALL,c('GATA1','PLEK','PF4','HBB','KLF1','ALAS2'))



outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/ek12_infantALL/'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

saveRDS(infantALL,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/ek12_infantALL/ek12_infantALL_clean_annotated.RDS')
