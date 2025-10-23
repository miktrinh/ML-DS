## Processing of "Other Leukaemia" single cell RNAseq data
## This includes the following samples:
# 1. B-ALL: 10 cases (8 diploid + 2 T21)
# 2. LPD_T21: 1 case (L062)
# 3. AML: 4 cases (all diploid)
# 4. AEL: 1 case (L069)
# 5. MDS: 1 case (L067)
# 6. AMKL(infant): 1 case (P9)


# sudo apt-get install bcftools
# sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
# sudo chmod +x /usr/local/bin/alleleCounter

setwd("/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia")

##------------------##
##    Libraries   ####
##------------------##
library(tidyverse)
#library(alleleIntegrator)
library(readxl)
library(Seurat)
#library(SoupX)
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
outDir = "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia"

source("~/lustre_mt22/generalScripts/utils/sc_basicQC.R")

message('\n--------  Hello! --------\n')

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
skipIfExists=F
keepMTCells=T


##----------------------------------------##
##   Import list of samples to process  ####
##----------------------------------------##
projMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'GOSH_others_included')
projMani = projMani[!is.na(projMani$assay) & projMani$assay == "5' V2 Dual Index" & projMani$`Point in treatment` != 'D29',]
projMani = projMani[projMani$donorID %in% c('L027','L062','L007','L001','L002','L080','L014','L016','L050','L068',
                                            'L051','L058','L044','L100',
                                            'L067','L069','L010'),]
projMani$Disease = ifelse(grepl('B-ALL',projMani$Disease),'BALL',
                          ifelse(projMani$Disease == 'Brain Lymphoma','LPD',
                                 ifelse(projMani$donorID == 'L067','MDS',
                                        ifelse(projMani$donorID == 'L069','AEL','AML'))))
#table(projMani$Disease,projMani$donorID)
#table(projMani$`Point in treatment`,projMani$donorID)
#table(projMani$sangerSampleID)
#length(projMani$sangerSampleID)
projMani$sampleID = gsub('.*SB_|.*ALeuk_','',projMani$sangerSampleID)

### 1. Import data

sub_outDir = file.path(outDir,'otherLeuk')
plotDir = file.path(sub_outDir,'otherLeuk_')
outPath = file.path(sub_outDir,'otherLeuk')
metadata = NULL
matchBy = NULL
cleanCountDir = file.path(sub_outDir,'cleanCount')



cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/otherLeuk_clean_withMTCells.RDS'),paste0(sub_outDir,'/otherLeuk_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)
}else{

# cleanSrat_fp = ifelse(keepMTCells,paste0(outDir,'/BALL_clean_withMTCells.RDS'),paste0(outDir,'/BALL_clean_noMTCells.RDS'))
# if(file.exists(cleanSrat_fp) & skipIfExists){
#   cleanSrat = readRDS(cleanSrat_fp)  
# }else{
  if(!dir.exists(sub_outDir)){
    message(sprintf('Creating output directory'))
    dir.create(sub_outDir,recursive = T)
  }
  
  setwd(sub_outDir)
  
  

  # List location of cellranger outputs
  dataDirs = list.dirs(c('/lustre/scratch126/casm/team274sb/project_folders/GOSH_Leuk/sc_raw_data/',
                         '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data'),recursive = T,full.names = T)
  dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_|^.*_ALeuk_','',dataDirs))
  ## Only keep samples of interest
  table(names(dataDirs) %in% projMani$sampleID)
  table(projMani$sampleID %in% names(dataDirs))
  projMani$sampleID[!projMani$sampleID %in% names(dataDirs)]
  dataDirs = dataDirs[names(dataDirs) %in% projMani$sampleID]
  dataDirs=dataDirs[file.exists(dataDirs)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  dataDirs = dataDirs[names(dataDirs) != 'Leuk13652354']
  
  ## Add Ellie's AMKL case from the infantALL paper
  amkl = c('/lustre/scratch126/casm/team274sb/project_folders/GOSH_Leuk/sc_raw_data/cellranger700_count_29622_4602STDY7920965_GRCh38-2020-A/filtered_feature_bc_matrix',
           '/lustre/scratch126/casm/team274sb/project_folders/GOSH_Leuk/sc_raw_data/cellranger700_count_29622_4602STDY7920966_GRCh38-2020-A/filtered_feature_bc_matrix')
  names(amkl) = c('4602STDY7920965','4602STDY7920966')
  dataDirs = c(dataDirs,amkl)
  
  print(n_distinct(dataDirs))
  
  
  
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose)
  print(n_distinct(QC.output[[1]]$orig.ident))
  #cleanSrat = QC.output[[1]]
  
  #df.out = QC.output[[2]]
  
  #write.csv(df.out,paste0(outDir,'/BALL_0923_qc_summary.csv'))
}








##################################################################################################################################


#############################
##       Annotation      ####
#############################
cleanSrat_fp = ifelse(keepMTCells,paste0(sub_outDir,'/otherLeuk_clean_withMTCells.RDS'),paste0(sub_outDir,'/otherLeuk_clean_noMTCells.RDS'))
cleanSrat = readRDS(paste0(sub_outDir,'/otherLeuk_clean_noMTCells.RDS'))
n_distinct(cleanSrat$orig.ident)

cleanSrat@misc = cleanSrat@misc[names(cleanSrat@misc) != 'preQC']
#cleanSrat = standard_clustering(cleanSrat)
DimPlot(cleanSrat,cols = col25,label = T,repel = T,label.box = T) + NoLegend() + NoAxes()

mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS')

##-------- Add metadata  --------##
colnames(mlds@meta.data)[!colnames(mlds@meta.data) %in% colnames(cleanSrat@meta.data)]
cleanSrat$donorID = projMani$donorID[match(cleanSrat$orig.ident,projMani$sampleID)]
cleanSrat$donorID[cleanSrat$orig.ident %in% c('4602STDY7920965','4602STDY7920966')] = 'P9'

cleanSrat$age_yrs = projMani$`age (years)`[match(cleanSrat$orig.ident,projMani$sampleID)]
cleanSrat$age_yrs[cleanSrat$donorID == 'P9'] = 'infant'

avgExpr = AverageExpression(cleanSrat,group.by = 'donorID',features = c('XIST','RPS4Y1'))
avgExpr = as.data.frame(t(avgExpr$RNA))
avgExpr$sex = ifelse(avgExpr$XIST > avgExpr$RPS4Y1,'F','M')
cleanSrat$sex = avgExpr$sex[match(cleanSrat$donorID,rownames(avgExpr))]

cleanSrat$tissue = projMani$Tissue[match(cleanSrat$orig.ident,projMani$sampleID)]
cleanSrat$tissue[cleanSrat$donorID == 'P9'] = 'BM'
cleanSrat$tissue[cleanSrat$tissue == 'Bone Marrow'] = 'BM'

cleanSrat$timePoint = projMani$`Point in treatment`[match(cleanSrat$orig.ident,projMani$sampleID)]
cleanSrat$timePoint[cleanSrat$orig.ident %in% c('4602STDY7920965')] = 'Diagnostic'
cleanSrat$timePoint[cleanSrat$orig.ident %in% c('4602STDY7920966')] = 'TP1'
cleanSrat$timePoint[cleanSrat$timePoint == 'D0'] = 'Diagnostic'

cleanSrat$blastPerc = projMani$`Blast Percentage`[match(cleanSrat$orig.ident,projMani$sampleID)]
cleanSrat$blastPerc[cleanSrat$donorID == 'P9'] = 'unknown'

cleanSrat$clinicalOutcome = projMani$`Clinical outcome`[match(cleanSrat$orig.ident,projMani$sampleID)]
cleanSrat$clinicalOutcome[cleanSrat$donorID == 'P9'] = 'Unknown'

cleanSrat$assay = 'GEX5p'
cleanSrat$assay[cleanSrat$donorID == 'P9'] = 'unknown'

cleanSrat$dataset = 'GOSH'
cleanSrat$dataset[cleanSrat$donorID == 'P9'] = 'unknown'

cleanSrat$Genotype = ifelse(cleanSrat$donorID %in% c('L027','L062'),'T21','diploid')
cleanSrat$disease = ifelse(cleanSrat$donorID %in% c('L062') & cleanSrat$orig.ident %in% c('NB14406184','NB14406185'),'LPD',
                           ifelse(cleanSrat$donorID %in% c('L067'),'MDS',
                                  ifelse(cleanSrat$donorID %in% c('L044','L051','L058','L069','L100','L010'),'pAML',
                                         ifelse(cleanSrat$donorID %in% c('P9'),'AMKL','pBALL'))))


## Add annotation
##---- Holly's annotation
mdat_pipeline = read.csv('/lustre/scratch126/casm/team274sb/project_folders/GOSH_Leuk/sample_metadata/manual_annotation_B-ALL_20240709.csv')
mdat_pipeline$cellID = gsub('::','_',gsub('ALeuk_|SB_','',mdat_pipeline$CellID))
table(cleanSrat$donorID,cleanSrat$cellID %in% mdat_pipeline$cellID)
cleanSrat$annot_pipeline = mdat_pipeline$manual_ann_VK[match(cleanSrat$cellID,mdat_pipeline$cellID)]
cleanSrat$annot_pipeline[is.na(cleanSrat$annot_pipeline)] = 'NA'
##---- Previous annotation
mdat_previousversion = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/combinedLeuk_mdat_2405.csv')
table(cleanSrat$cellID %in% mdat_previousversion$cellID)
table(mdat_previousversion$cellID %in% cleanSrat$cellID,mdat_previousversion$disease)
table(cleanSrat$orig.ident[!cleanSrat$cellID %in% mdat_previousversion$cellID])
cleanSrat$annot_apr24 = mdat_previousversion$finalAnn_broad[match(cleanSrat$cellID,mdat_previousversion$cellID)]
cleanSrat$annot_apr24[is.na(cleanSrat$annot_apr24)] = 'NA'

dd = as.data.frame(table(cleanSrat$annot_apr24,cleanSrat$annot_pipeline))
dd = dd[dd$Freq > 0,]

cleanSrat$annot_aug24 = '?'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'B-ALL' & cleanSrat$annot_apr24 == 'Tumour'] = 'Tumour'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'NA' & cleanSrat$annot_apr24 == 'Tumour'] = 'Tumour'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'B-ALL' & cleanSrat$annot_apr24 == 'NA'] = 'Tumour'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'B cells' & cleanSrat$annot_apr24 == 'naive.B'] = 'naive.B'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'T cells' & cleanSrat$annot_apr24 %in% c('T_CD4','T_CD8','T_gd','T_MAIT','NK','NK_T')] = cleanSrat$annot_apr24[cleanSrat$annot_pipeline == 'T cells' & cleanSrat$annot_apr24 %in% c('T_CD4','T_CD8','T_gd','T_MAIT','NK','NK_T')]
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'Myeloid' & cleanSrat$annot_apr24 %in% c('Mono_CD14')] = 'Mono_CD14'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'Myeloid' & cleanSrat$annot_apr24 %in% c('Mono_CD16')] = 'Mono_CD16'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'Myeloid' & cleanSrat$annot_apr24 %in% c('Neutrophil')] = 'Neutrophil'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'Myeloid' & cleanSrat$annot_apr24 %in% c('Myelocytes')] = 'Myelocytes'
cleanSrat$annot_aug24[cleanSrat$annot_pipeline == 'Erythroid' & !cleanSrat$annot_apr24 %in% c('NA')] = cleanSrat$annot_apr24[cleanSrat$annot_pipeline == 'Erythroid' & !cleanSrat$annot_apr24 %in% c('NA')]


DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$annot_apr24 == 'Tumour' & cleanSrat$annot_pipeline == 'NA'])
DimPlot(cleanSrat,group.by = 'donorID',cols = c(col25,pal34H),label.box = T,label = T,repel = T,label.size = 3) + NoAxes() + NoLegend()
DimPlot(cleanSrat,group.by = 'seurat_clusters',label.box = T,label = T,repel = T,label.size = 3) + NoAxes() + NoLegend()

## Annotation
library(SoupX)
qm = quickMarkers(cleanSrat@assays$RNA@counts,cleanSrat$seurat_clusters)
FeaturePlot(cleanSrat,'percent.mt')
cleanSrat = FindClusters(cleanSrat, resolution = 0.5)
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$seurat_clusters %in% c(45)])
DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$donorID =='P9' & cleanSrat$timePoint == 'Diagnostic'])
table(cleanSrat$annot_apr24[cleanSrat$seurat_clusters %in% c(43)])
table(cleanSrat$seurat_clusters[cleanSrat$annot_aug24 == '?'])


cleanSrat$annot_aug24[cleanSrat$seurat_clusters == 36] = 'Neuron_?'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(22,0,5,27,1,38,39)] = 'Tumour'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(44)] = 'Fibroblast_?'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(37)] = 'Neuron2_?'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(42)] = 'Glial_cells?'

cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(18,24,20,3,
                                                       25,33,40,2,24,12,4,
                                                       31,6,9,20,26,10,
                                                       3,21,14)] = 'Tumour'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(15) & cleanSrat$annot_apr24 == 'pDC'] = 'pDC'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(15) & cleanSrat$annot_aug24 == '?'] = 'naive.B'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(15) & cleanSrat$annot_aug24 == '?'] = 'naive.B'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(41)] = 'Plasma.cell'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(7)] = 'T_cells'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(7)] = 'T_cells'

cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(8) & cleanSrat$annot_aug24 == '?'] = 'doublets'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(11) & cleanSrat$annot_aug24 == '?'] = 'Tumour'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(13) & cleanSrat$annot_aug24 == '?' & !cleanSrat$annot_apr24 %in% c('18','25','Normal','doublets','NA','B.cell','Plasma.cell')] = cleanSrat$annot_apr24[cleanSrat$seurat_clusters %in% c(13) & cleanSrat$annot_aug24 == '?' & !cleanSrat$annot_apr24 %in% c('18','25','Normal','doublets','NA','B.cell','Plasma.cell')]
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(13) & cleanSrat$annot_apr24 %in% c('18','25','Normal','doublets','NA','B.cell','Plasma.cell')] = 'T_cells'

cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(19) & cleanSrat$annot_aug24 == '?' & cleanSrat$annot_apr24 %in% c('EE','HSC_MPP','ME','LE','MEMP_MEP','MEP','MK')] = cleanSrat$annot_apr24[cleanSrat$seurat_clusters %in% c(19) & cleanSrat$annot_aug24 == '?' & cleanSrat$annot_apr24 %in% c('EE','HSC_MPP','ME','LE','MEMP_MEP','MEP','MK')]
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(29) & cleanSrat$annot_aug24 == '?' & cleanSrat$annot_apr24 %in% c('NA','naive.B')] = 'T_cells'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(29) & cleanSrat$annot_aug24 == '?' & !cleanSrat$annot_apr24 %in% c('NA','naive.B')] = cleanSrat$annot_apr24[cleanSrat$seurat_clusters %in% c(29) & cleanSrat$annot_aug24 == '?' & !cleanSrat$annot_apr24 %in% c('NA','naive.B')]
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(30)] = 'T_cells'

cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(34) & cleanSrat$annot_aug24 == '?' & cleanSrat$annot_apr24 %in% c('activated_neutrophil','Neutrophil','neutrophils')] = 'Neutrophil'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(34) & cleanSrat$annot_aug24 == '?' & cleanSrat$annot_apr24 %in% c('25','8','doublets','NK_T','T_CD4','T_CD8','T/NK')] = 'doublets'

cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(43) & cleanSrat$annot_aug24 == '?'] = 'Tumour_?'



##----------------------------------##
##    Subclustering some cells    ####
##----------------------------------##
s = subset(cleanSrat,subset = cellID %in% unique(c(cleanSrat$cellID[cleanSrat$annot_aug24 %in% c('?','MEP','EE','ME','LE','HSC_MPP','MK','MEMP_MEP','naive.B','Plasma.cell')],
                                            cleanSrat$cellID[cleanSrat$seurat_clusters %in% c(16,23,17,19,32,45)],
                                            cleanSrat$cellID[cleanSrat$donorID %in% c('L051','L016','P9')])))

s = standard_clustering(s,clusteringRes = 0.3)
s = FindClusters(s, resolution = 1)
DimPlot(s,group.by = 'annot_aug24',cols = c(col25,pal34H),label.box = T,label = T,repel = T,label.size = 3) + NoAxes() + NoLegend()
DimPlot(s,cells.highlight = s$cellID[s$seurat_clusters== 23 & s$annot_aug24 == '?' & s$donorID != 'P9'])

qm = quickMarkers(s@assays$RNA@counts,s$seurat_clusters)

table(s$seurat_clusters[s$annot_aug24 == '?'])
table(s$donorID[s$seurat_clusters %in% c(23) & s$annot_aug24 == '?'])

s$annot_aug24[s$annot_aug24 == 'Tumour_?'] = 'Tumour'
s$annot_aug24[s$seurat_clusters == 11] = 'LE'
s$annot_aug24[s$seurat_clusters == 9 & s$annot_aug24 == '?'] = 'EE'
s$annot_aug24[s$seurat_clusters == 18 & s$annot_aug24 == 'Tumour'] = 'ME'
s$annot_aug24[s$seurat_clusters == 23 & s$annot_aug24 == '?' & s$timePoint != 'Diagnostic'] = 'MK'
s$annot_aug24[s$seurat_clusters %in% c(0,3,5,7,8,10,12,14,15,16,24,28,34,37,21,22,25,29) & s$annot_aug24 == '?'] = 'unknown'
s$annot_aug24[s$seurat_clusters %in% c(4,10,17,27,31) & s$annot_aug24 == '?'] = 'Tumour'
s$annot_aug24[s$seurat_clusters %in% c(20,32) & s$annot_aug24 == '?'] = 'T_cells'
s$annot_aug24[s$seurat_clusters %in% c(9) & s$annot_aug24 == '?' & s$donorID == 'P9'] = 'Tumour'
s$annot_aug24[s$seurat_clusters %in% c(9) & s$annot_aug24 == '?' & s$donorID != 'P9'] = 'unknown'
s$annot_aug24[s$seurat_clusters %in% c(23) & s$annot_aug24 == '?' & s$donorID != 'P9'] = 'MK'

table(s$annot_aug24)

## Add to cleanSrat
cleanSrat$annot_aug24[cleanSrat$cellID %in% s$cellID] = s$annot_aug24[match(cleanSrat$cellID[cleanSrat$cellID %in% s$cellID],s$cellID)]


## Add annotation for L010
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_withUnknown_2408_mdat.csv')
table(cleanSrat$cellID %in% mdat$cellID)
cleanSrat$annot_aug24[cleanSrat$cellID %in% mdat$cellID] = mdat$annot_aug24[match(cleanSrat$cellID[cleanSrat$cellID %in% mdat$cellID],mdat$cellID)]

DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$annot_aug24 == 'unknown'])
table(cleanSrat$seurat_clusters[cleanSrat$annot_aug24 == 'unknown' & cleanSrat$donorID == 'L010'])
table(cleanSrat$donorID[cleanSrat$seurat_clusters == 23],cleanSrat$annot_aug24[cleanSrat$seurat_clusters == 23])
table(cleanSrat$annot_aug24[cleanSrat$seurat_clusters == 29])
FeaturePlot(cleanSrat,'CD3D')

cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(10,13,3) & cleanSrat$donorID == 'L010' & cleanSrat$annot_aug24=='?'] = 'T_cells'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(11,17,30,31,32,34,57) & cleanSrat$donorID == 'L010' & cleanSrat$annot_aug24=='?'] = 'unknown'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(15) & cleanSrat$donorID == 'L010' & cleanSrat$annot_aug24=='?'] = 'Tumour'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(20,41) & cleanSrat$donorID == 'L010' & cleanSrat$annot_aug24=='?'] = 'naive.B'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(29) & cleanSrat$donorID == 'L010' & cleanSrat$annot_aug24=='?'] = 'Mono_CD14'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters %in% c(38) & cleanSrat$donorID == 'L010' & cleanSrat$annot_aug24=='?'] = 'EE'
cleanSrat$annot_aug24[cleanSrat$seurat_clusters == 17] = 'T_cells'

library(SoupX)
qm = quickMarkers(cleanSrat@assays$RNA@counts[,cleanSrat$cellID[cleanSrat$seurat_clusters %in% c(2,23,18)]],
                  cleanSrat$seurat_clusters[cleanSrat$seurat_clusters %in% c(2,23,18)])

table(cleanSrat$donorID[cleanSrat$annot_aug24 == '?'])
##--------------------------##
##    Save final object   ####
##--------------------------##

DimPlot(cleanSrat,group.by = 'annot_aug24',cols = c(col25,pal34H),label.box = T,label = T,repel = T,label.size = 3) + NoAxes() + NoLegend()
DimPlot(cleanSrat,label.box = T,label = T,repel = T,label.size = 3) + NoAxes() + NoLegend()
for(d in unique(cleanSrat$donorID)){
  p = DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$donorID == d & cleanSrat$annot_aug24 == 'Tumour']) + ggtitle(d)  
  print(p)
}

DimPlot(cleanSrat,cells.highlight = cleanSrat$cellID[cleanSrat$orig.ident == 'Leuk13104275'])
## Define broad lineages
cleanSrat$annot_aug24[cleanSrat$annot_aug24 == 'MEMP_MEP']='MEP'
cleanSrat$annot_aug24[cleanSrat$annot_aug24 == 'Tumour_?']='Tumour'
cleanSrat$broadLineage = '?'
table(cleanSrat$annot_aug24,cleanSrat$broadLineage)
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('naive.B','Plasma.cell','pre.B.cell','pro.B.cell')] = 'B lineage'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('DC1','DC2','pDC')] = 'Dendritic cells'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('EE','ME','LE')] = 'Erythroblasts'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('CMP_GMP','HSC_MPP','MEP')] = 'HSC & prog.'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('MK')] = 'Megakaryocytes'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('Mono_CD14','Mono_CD16','Macrophage')] = 'Monocyte/Macrophage'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('Myelocytes','activated_neutrophil','Neutrophil')] = 'Neutrophil'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('NK','NK_T','T_CD4','T_CD8','T_gd','T_reg','T_cells','T_MAIT','T/NK')] = 'T/NK lineage'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('Tumour')] = 'Tumour'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('doublets','NA')] = 'lowQual'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('unknown','?')] = 'unknown'
cleanSrat$broadLineage[cleanSrat$annot_aug24 %in% c('Fibroblast_?','Glial_cells?','Neuron_?','Neuron2_?')] = 'others'


cleanSrat$finalAnn_broad = as.character(cleanSrat$annot_aug24)
cleanSrat$annot = as.character(cleanSrat$annot_aug24)
cleanSrat$finalAnn = as.character(cleanSrat$annot_aug24)

mdat = cbind(cleanSrat@meta.data,cleanSrat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_withUnknown_2408_mdat.csv')

##------    Remove doublets / unknown   --------##
cleanSrat = subset(cleanSrat,subset = annot_aug24 %in% unique(cleanSrat$annot_aug24[!cleanSrat$annot_aug24 %in% c('?','doublets','unknown')]))
cleanSrat = standard_clustering(cleanSrat)
mdat = cbind(cleanSrat@meta.data,cleanSrat@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2408_mdat.csv')
saveRDS(cleanSrat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2408.RDS')





##--------------------------------------##
##    Subclustering each cancer type  ####
##--------------------------------------##

##---------- pAML
s = subset(srat,subset = disease == 'pAML')
s = standard_clustering(s)
DimPlot(s,group.by = 'seurat_clusters',cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()
qm = quickMarkers(s@assays$RNA@counts[,s$cellID[s$donorID == 'L010']],s$annot[s$donorID == 'L010'])

##---------- pBALL
s = subset(srat,subset = disease == 'pBALL')
s = standard_clustering(s)
DimPlot(s,group.by = 'annot_aug24',cols = c(col25,pal34H),label = T,repel = T,label.box = T) + NoLegend()
qm = quickMarkers(s@assays$RNA@counts[,s$cellID[s$donorID == 'L010']],s$annot[s$donorID == 'L010'])






## Investigate L010 - weird 2 clusters of blasts

l10 = subset(srat,subset = donorID == 'L010')
l10 = standard_clustering(l10)
DimPlot(l10,group.by = 'annot')
DimPlot(l10,cells.highlight = mdat$cellID[mdat$broadLineage == 'B lineage'])
FeaturePlot(l10,c('CD34','CD14','CD3D','MLLT10','KMT2A','CD5'))
FeaturePlot(mlds,c('CD3D','MLLT10','KMT2A'))
qm = quickMarkers(l10@assays$RNA@counts,l10$annot)
FeaturePlot(l10,)


avgExpr = AverageExpression(l10,features = genes$geneSym,group.by = 'annot')
avgExpr = avgExpr$RNA
genes = geneMap[geneMap$chr == 'chr10',]
genes$coord = start(gns[match(genes$ensID,gns$gene_id)])
genes$arm = ifelse(genes$coord >= 126000000,'q','p')
table(genes$arm)
genes = cbind(genes,avgExpr[match(genes$geneSym,rownames(avgExpr)),])
ggplot(genes,aes(coord/1000,1,col=arm))+
  geom_point()+
  xlab('')+
  theme(axis.text = element_blank())
ggplot(genes,aes(log10(Tumour),log10(naive.B),col=arm))+
  geom_point(size=0.1)+
  xlab('')+
  theme(axis.text = element_blank())

genes_long = pivot_longer(genes,cols = colnames(genes)[c(6:10)],names_to = 'annot',values_to = 'expr')
ggplot(genes_long,aes(annot,log2(expr)))+
  geom_boxplot()


DimPlot(l10,group.by = 'Phase')
FeaturePlot(l10,c('ZFPM2','CD5'))


DotPlot(l10,group.by = 'annot',features = unique(c("FCN1", "CD14", "FCGR3A", "FCER1A", "CLEC10A",
                                                   "HBA1", "HBB", "GYPA",
                                                   "GNLY", "NKG7", "KLRG1",
                                                   "MS4A1",
                                                   "ITGB1",
                                                   "COL4A4",
                                                   "PRDM1",
                                                   "IRF4",
                                                   "PAX5",
                                                   "BCL11A",
                                                   "BLK",
                                                   "IGHD",
                                                   "IGHM",
                                                   "MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN",
                                                   "CD4", "IL7R", "TRBC2", "CD8A", "CD8B", "GZMK", "LEF1", "CCR7", "TCF7",
                                                   'TRDC', 'TRAC', 'TRBC2', 'CD8A', 'CD8B',
                                                   "GZMB", "IL3RA", "COBLL1", "TCF4", "LILRA4", "CLEC4C",
                                                   "HLA-DRA", "CD38", "CD44", "CD19", "CD79A", "CD34", "PTPRC", "CSPG4", "CD22", "CD24", "MME", "DNTT")
                                                 
))+
  RotatedAxis()



#---------------------------------------------------------------------------------------------------------------------------------
# outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL'
# if(!dir.exists(outDir)){
#   dir.create(outDir)
# }
# 
# 
# ## 1. Import relevant objects
# sampleSheet = read_excel('~/lustre_mt22/Data/B_ALL/Manifest_2023_09_13.xlsx',sheet = 'B-ALL (Study 6918_7153_3030)')
# sampleSheet = sampleSheet[sampleSheet$Experiment == 'GEX' & 
#                             grepl('BALL',sampleSheet$Sample_ID),]
# sampleSheet = sampleSheet[!sampleSheet$Cellgen_ID %in% sampleSheet$Cellgen_ID[sampleSheet$Patient_ID %in% c('L007','L071','L024','L031','L096','L057','L045','L073','L033','L059','L043','L072','L023','L087','L085','L032','L004','L006','L025','L048','L003') & sampleSheet$Timepoint == 'D0' | sampleSheet$Subtype == 'Other'],]
# 
# 
# 
# ## 2. Import object
# tgt.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/sept23/BALL/BALL_clean_noMTCells.RDS')
# tgt.srat = standard_clustering(tgt.srat)
# DimPlot(tgt.srat,label.box = T,label = T,repel = T)
# DimPlot(tgt.srat,group.by = 'orig.ident',label.box = T,label = T,repel = T)
# 
# FeaturePlot(tgt.srat,c('GATA1','KLF1'))
# 
# ## 2. Add metadata
# m = match(tgt.srat$orig.ident,gsub('SB_','',sampleSheet$Cellgen_ID))
# table(is.na(m))
# tgt.srat$assay = 'GEX5p'
# tgt.srat$donorID = sampleSheet$Patient_ID[m]
# tgt.srat$age_yrs = 'paed'
# tgt.srat$sex = sampleSheet$Sex[m]
# tgt.srat$tissue = sampleSheet$Tissue[m]
# tgt.srat$mutation = sampleSheet$Subtype[m]
# tgt.srat$timePoint = sampleSheet$Timepoint[m]
# tgt.srat$blastPerc = sampleSheet$`Flow_blast_
# percent`[m]
# tgt.srat$clinicalOutcome = sampleSheet$Response[m]
# 
# DimPlot(tgt.srat,group.by = 'seurat_clusters',label.box = T,label = T,repel = T)
# DimPlot(tgt.srat,group.by = 'timePoint',label.box = T,label = F,repel = T,cols = col25)





## 3. Broad annotation
keyMarkers = unique(c('HLA-DRA','CLEC9A','CD34','CD38','HLF','SPINK2','MLLT3','PRSS57', # HSC_MPP
                      'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1',#MEMP
                      'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
                      #'KIT',
                      'CTSG',	'PRTN3', # CMP
                      'AZU1','MPO','FLT3','PTPRC', # GMP
                      'ZBTB16','LTB', 'CD52',# ILC precursor
                      'IL7R','CD3D','GZMA','CD4','CD8A', #Early.lymphoid_T.lymphocyte
                      'TRDV2','TRGV9', # gamma delta T-cell
                      'SLC4A10','TRAV1-2', #MAIT t-cell
                      'PRF1', # effector T cells
                      'FOXP3',	'CDH1', # regulatory T cells
                      'NKG7','KLRD1', #NK
                      'IGLL1','CD99', # Pre-pro
                      'DNTT','CD79B','VPREB1','EBF1','CD19','RAG1',# pro-b-cells
                      'MME','CD79A',# pre-b-cells
                      'TCL1A','MME','RAG1','MS4A1',  # B-cells
                      'MAFB','ITGB2','SELL','ITGAM','CD14','CCR2','FCGR3A','CX3CR1',# Monocytes
                      'CD68','MSR1',#Macrophages
                      'IRF8',	'CLEC10A', #DC.precursor
                      'ITGAX','CD1C','CLEC9A','THBD',#DC1
                      'CLEC10A', # DC2 
                      'IL3RA', 'CLEC4C',#pDC
                      'CD27',#plasma cell?
                      'CSF2RB','HDC','SERPINB1','TPSAB1','KIT', # Mast.cell
                      'PF4','ITGA2B', #Megakaryocyte
                      'PPBP','TUBB1', # Platelets
                      'GATA1','KLF1', # early Erythroid
                      'ALAS2', # mid.erythroid
                      'HBA1','BPGM', # late.erythroid
                      'CSF3R','FPR1','FCGR3B','NAMPT','MNDA', # NEUTROPHILS
                      'DEFA3','DEFA4','CAMP', 'LCN2',
                      'ESAM','PECAM1', 'KDR', 'PTPRB', 'PLVAP', # Endothelium
                      'DCN','SERPINF1','COL3A1','BGN','COL1A1','VIM','ECM1', # endosteal fibroblast / fibroblast
                      'APOA1','SCD','ALB','TTR' # Hepatocyte
))
# DotPlot(tgt.srat,
#         features = keyMarkers)+RotatedAxis()
# 
# DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$seurat_clusters %in% c(17,19,20,29,36,7)])
# VlnPlot(tgt.srat,group.by = 'cluster_ann',features =c('percent.mt','nFeature_RNA','nCount_RNA'),ncol=1)
# 
# DimPlot(tgt.srat,group.by = 'donorID',label = T,label.box = T,repel = T,cols = col25) 
# 
# tgt.srat$cluster_ann = as.character(tgt.srat$seurat_clusters)
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(7,17,41,47)] = 'Erythroblasts'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(26)] = 'MEMP_MEP'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(44)] = 'MK'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(19)] = 'neutrophils'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(0:6,8,10,11,13,15,18,20,22,23,27:37,39,40,42,43,46)] = 'Tumour'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(9,12,14,16,21)] = 'T.cell'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(14)] = 'T_CD8'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(25)] = 'NK'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(24)] = 'CMP_GMP'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(32)] = 'MonoMac'
# tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(45)] = 'B.cell'
# 
# write.csv(tgt.srat@meta.data,'MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_clean_annotated_0923_mdat.csv')
# 
# 
# library(SoupX)
# qm = quickMarkers(tgt.srat@assays$RNA@counts,tgt.srat$cluster_ann)
# write.csv(qm,'MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_clean_annotated_0923_quickMarkers.csv')
# 
# 
# 
# 
# ##------------------------------------------##
# ##    Sub-clustering non-tumour cells     ####
# ##------------------------------------------##
# 
# normCells = subset(tgt.srat,subset = cellID %in% tgt.srat$cellID[tgt.srat$cluster_ann != 'Tumour'])
# normCells = standard_clustering(normCells)
# DimPlot(normCells,group.by = 'cluster_ann',label = T,repel = T,label.box = T,cols=c(col25,pal34H)) + NoLegend()
# FeaturePlot(normCells,'scrubScore')
# 
# qm.2 = quickMarkers(normCells@assays$RNA@counts,normCells$seurat_clusters)
# 
# DotPlot(normCells,group.by = 'seurat_clusters',#idents = unique(s$ann_tmp[grepl('B.cell',s$ann_tmp)]),
#         features = keyMarkers)+RotatedAxis()
# 
# table(normCells$cluster_ann[normCells$seurat_clusters  %in% c(19,16,13,9,7,41,0)])
# DimPlot(normCells,cells.highlight = memp)
# DimPlot(normCells,cells.highlight = normCells$cellID[normCells$cluster_ann %in% c('neutrophils')])
# 
# 
# normCells$cluster_ann[normCells$seurat_clusters %in% c(32,18)] = 'LE'
# normCells$cluster_ann[normCells$seurat_clusters %in% c(3,10,5)] = 'ME'
# normCells$cluster_ann[normCells$seurat_clusters %in% c(2,22)] = 'EE'
# memp = WhichCells(normCells,cells = normCells$cellID[normCells$seurat_clusters %in% c(2,22)],expression = (PRSS57 > 0.5 & HBB < 1))
# normCells$cluster_ann[normCells$seurat_clusters %in% c(2,22) & normCells$cellID %in% memp] = 'MEMP_MEP'
# normCells$cluster_ann[normCells$seurat_clusters %in% c(2,22)] = 'EE'
# normCells$cluster_ann[normCells$seurat_clusters %in% c(6)] = 'neutrophil'
# normCells$cluster_ann[normCells$seurat_clusters %in% c(23)] = 'myelocyte'
# 
# gd_Tcells = WhichCells(normCells,cells = normCells$cellID[normCells$seurat_clusters == 11],expression = (TRDV2 > 0.5 | TRGV9 > 0.5))
# mait_Tcells = WhichCells(normCells,cells = normCells$cellID[normCells$seurat_clusters == 11],expression = (SLC4A10 > 0.5))
# normCells$cluster_ann[normCells$cellID %in% c(gd_Tcells)] = 'T_gd'
# normCells$cluster_ann[normCells$cellID %in% c(mait_Tcells)] = 'T_MAIT'
# normCells$cluster_ann[normCells$seurat_clusters %in% c(15,24)] = 'T_CD8'
# normCells$cluster_ann[normCells$seurat_clusters %in% c(19,16,13,9,7,41,0) & normCells$cluster_ann == 'T_CD8'] = 'T.cell'
# tcell = WhichCells(normCells,cells = normCells$cellID[normCells$cluster_ann %in% c('T_CD8') & !normCells$seurat_clusters %in% c(15,24)],expression = (CD8A <  1 | CD8B < 1))
# normCells$cluster_ann[normCells$cellID %in% c(tcell)] = 'T.cell'
# 
# 
# normCells$cluster_ann[normCells$seurat_clusters %in% c(14)] = 'Mono_CD14'
# normCells$cluster_ann[normCells$seurat_clusters %in% c(28)] = 'Mono_CD16'
# 
# normCells$cluster_ann[normCells$seurat_clusters %in% c(25) & !normCells$cluster_ann %in% c('Erythroblasts')] = 'HSC_MPP'
# #ilc = WhichCells(normCells,cells = normCells$cellID[normCells$cluster_ann %in% c('HSC_MPP')],expression = (ZBTB16 > 0.5 & LTB >  0.5 & CD52 > 0.5))
# normCells$cluster_ann[normCells$seurat_clusters %in% c(26)] = 'MEMP_MEP'
# 
# normCells$cluster_ann[normCells$seurat_clusters %in% c(20)] = 'Mono_weird'
# 
# normCells$cluster_ann[normCells$seurat_clusters %in% c(35,36,30,31,37,34,33,29)] = 'doublets'
# normCells$cluster_ann[normCells$cluster_ann %in% c('Erythroblasts','MonoMac','neutrophils')] = 'doublets'
# 
# 
# ## Add this annotation back to tgt.srat
# tgt.srat$cluster_ann[tgt.srat$cellID %in% normCells$cellID] = normCells$cluster_ann[match(tgt.srat$cellID[tgt.srat$cellID %in% normCells$cellID],normCells$cellID)]
# normCells$finalAnn = normCells$cluster_ann
# 
# saveRDS(normCells,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_normCellsOnly_clean_annotated_0923.RDS')
# write.csv(normCells@meta.data,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_normCellsOnly_clean_annotated_0923_mdat.csv')
# 
# ##------------------------------##
# ##    Previous annotation     ####
# ##------------------------------##
# goshALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/may23/BALL/DSLeuk_clean_noMTCells_annotated_may23.RDS')
# table(goshALL$finalAnn_broad[goshALL$cellID %in% normCells$cellID])
# mdat = goshALL@meta.data[goshALL@meta.data$cellID %in% normCells$cellID,]
# mdat$newLab = normCells$cluster_ann[match(mdat$cellID,normCells$cellID)]
# View(table(mdat$finalAnn_broad,mdat$newLab))
# 
# # Looks alright - not much changed...
# 
# ## Import tgt.srat mdat
# mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_clean_annotated_0923_mdat.csv')
# table(tgt.srat$cellID %in% mdat$cellID)
# tgt.srat$cluster_ann = mdat$cluster_ann[match(tgt.srat$cellID,mdat$cellID)] 
# 
# DimPlot(tgt.srat,group.by = 'cluster_ann',label = T,repel = T,label.box = T,cols = col25) + NoLegend()
# 
# saveRDS(tgt.srat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_clean_annotated_0923.RDS')
# 
# 
# 
# ##---------------------------------------------##
# ##    May 2023 session - merged with MLDS    ####
# ##---------------------------------------------##
# 
# ### Import previous annotation
# ds.prev_ann = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/may23/MLDS_DSLeuk_clean_annotated_may23.RDS')
# m = match(tgt.srat$cellID,ds.prev_ann$cellID)
# table(is.na(m))
# tgt.srat$may23_ann = 'NA'
# tgt.srat$may23_ann[!is.na(m)] = ds.prev_ann$finalAnn_broad[m[!is.na(m)]]
# View(table(tgt.srat$may23_ann,tgt.srat$cluster_ann))
# 
# DimPlot(tgt.srat,cells.highlight = tgt.srat$cellID[tgt.srat$may23_ann == 'naive_B_cell'])
# DimPlot(ds.prev_ann,group.by = 'finalAnn_broad')
# 
# 
# mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/may23/MLDS_clean_LRwCTannotated_may23.RDS')
# 
# mlds.ball.merged = merge_seurat_objects(mlds,tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# mlds.ball.merged = standard_clustering(mlds.ball.merged)
# 
# DimPlot(mlds.ball.merged, group.by = 'cluster_ann_merged',label = T,label.box = T,cols = c(col25[-6],col22),repel = T) + NoLegend()
# DimPlot(mlds.ball.merged, group.by = 'seurat_clusters',label = T,label.box = T,repel = T) + NoLegend()
# 
# DotPlot(mlds.ball.merged,#idents = c(18,40,57),
#         #group.by = 'cluster_ann',scale = F,
#         group.by = 'cluster_ann_merged', scale=F,
#         features = keyMarkers)+RotatedAxis()
# 13,25,26,
# 15,29,30,32,34,38,47,53
# DimPlot(mlds.ball.merged,cells.highlight = mlds.ball.merged$cellID[mlds.ball.merged$seurat_clusters %in% c(10) & mlds.ball.merged$cluster_ann_merged %in% c('Tumour')])
# table(mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(10)])
# mlds.ball.merged$cluster_ann_merged = as.character(mlds.ball.merged$seurat_clusters)
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(0,24)] = 'ME_LE'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(0,24) & mlds.ball.merged$cluster_ann %in% c('pro_B_cells','pre_B_cells','MEMP_MEP','Neutrophil','Myelocyte')] = 'doublets'
# 
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(18,40,57) & mlds.ball.merged$cluster_ann %in% c('Ery','ME')] = 'ME'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(18,40,57) & mlds.ball.merged$cluster_ann %in% c('LE')] = 'LE'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(18,40,57) & mlds.ball.merged$cluster_ann %in% c('EE','CD4_Tcell')] = 'EE'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(18,40,57) & mlds.ball.merged$cluster_ann %in% c('Neutrophil','HSC_MPP','MEMP_MEP','Tumour')] = 'doublets'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(41)] = 'Myelocyte'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(41) & mlds.ball.merged$cluster_ann %in% c('Neutrophil','Myelocyte','34')] = 'Myelocyte'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(41) & mlds.ball.merged$cluster_ann %in% c('CMP_GMP','CD4_Tcell','prog')] = 'Myelocyte'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(41) & mlds.ball.merged$cluster_ann %in% c('48')] = 'Neutrophil'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(41) & mlds.ball.merged$cluster_ann %in% c('pro_B_cells','pre_B_cells')] = 'doublets'
# 
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(17,19,2,21,22,3,
#                                                                             33)] = 'T_cell'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(17,19,2,21,22,3,33) & mlds.ball.merged$cluster_ann %in% c('CD4_Tcell','CD8_Tcell','gd_Tcell')] = mlds.ball.merged$cluster_ann[mlds.ball.merged$seurat_clusters %in% c(17,19,2,21,22,3,33) & mlds.ball.merged$cluster_ann %in% c('CD4_Tcell','CD8_Tcell','gd_Tcell')]
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(17,19,2,21,22,3,33) & mlds.ball.merged$cluster_ann %in% c('Tcell','40')] = 'CD4_Tcell'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(17,19,2,21,22,3,33) & mlds.ball.merged$cluster_ann %in% c('34')] = 'CD8_Tcell'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(17,19,2,21,22,3,33) & mlds.ball.merged$cluster_ann %in% c('NK','NK_T')] = 'NK_T'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(17,19,2,21,22,3,33) & mlds.ball.merged$cluster_ann %in% c('unknown')] = 'clonal_NK_T?'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(17,19,2,21,22,3,33) & mlds.ball.merged$cluster_ann %in% c('pro_B_cells','pre_B_cells','B_cells','DC1','MK','Myelocyte','Tumour','Tumour?')] = 'doublets'
# 
# 
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(32,11,47) & mlds.ball.merged$cluster_ann %in% c('EE','LE','Ery','ME','MK','CD4_Tcell','HSC_MPP','pDC','pro_B_cells','pre_B_cells','Tcell','Tumour?','weird_MEMP')] = 'doublets'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(32,11,47) & mlds.ball.merged$cluster_ann %in% c('CD14_Mono','CMP_GMP','DC2','Mono.Mac','prog','Tumour')] = 'CD14_Mono'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(32,11,47) & mlds.ball.merged$cluster_ann %in% c('CD16_Mono')] = 'CD16_Mono'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(32,11,47) & mlds.ball.merged$cluster_ann %in% c('Neutrophil')] = 'Neutrophil'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(32,11,47) & mlds.ball.merged$cluster_ann %in% c('Myelocyte')] = 'Myelocyte'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(32,11,47) & mlds.ball.merged$cluster_ann %in% c('Myelocyte')] = 'Myelocyte'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$cluster_ann_merged %in% c(11)] = 'CD14_Mono'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('B_cells','pro_B_cells','pre_B_cells','CD14_Mono','DC1')] = 'doublets'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('CD4_Tcell','Tcell')] = 'CD4_Tcell'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('CD8_Tcell')] = 'NK_T'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('CMP_GMP')] = 'CMP_GMP'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('LE')] = 'LE'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('ME')] = 'ME'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('EE')] = 'EE_MEP'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('Ery')] = 'EE_MEP'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('MEMP_MEP','HSC_MPP')] = 'MEMP_MEP?'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('MK')] = 'MK?'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('Myelocyte','Neutrophil')] = 'Myelocyte'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('gd_Tcell','NK','NK_T','unknown')] = 'NK'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('weird_MEMP')] = 'weird_MEMP'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('Tumour?','Tumour')] = 'Tumour'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('40')] = 'Endo'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13,14,15,16,30,31,45,47,5,51,53,54,6,7) & mlds.ball.merged$cluster_ann %in% c('42')] = 'Plasma_cell'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(54)] = 'Plasma_cell'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$cluster_ann_merged %in% c('14')] = 'weird_MEMP'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$cluster_ann_merged %in% c('16')] = 'Tumour'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$cluster_ann_merged %in% c('31')] = 'doublets'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(28) & mlds.ball.merged$cluster_ann %in% c('HSC_MPP','CMP_GMP','CD14_Mono','pDC','Myelocyte','48','30','6','23','24')] = 'CMP_GMP'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(28) & mlds.ball.merged$cluster_ann %in% c('Neutrophil','pro_B_cells','pre_B_cells','Ery','43','DC2')] = 'doublets'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(28) & mlds.ball.merged$cluster_ann %in% c('CD4_Tcell')] = 'CD4_Tcell'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13) & mlds.ball.merged$cluster_ann %in% c('Tumour','Tumour?')] = 'Tumour'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(13) & mlds.ball.merged$cluster_ann %in% c('EE')] = 'doublets'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(7) & mlds.ball.merged$cluster_ann_merged %in% c('7')] = 'Tumour'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('0','1','11','13','14','15','16','18',2,23,25,26,28,29,3,30,32,35,37,39,4,42,43,49,5,8,9,33,34,6,22,20,40,44,46)] = 'Tumour'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c(7,19,24,27,41,47,48,29,20,36,41,43,'CD4_Tcell','CD8_Tcell','Ery','weird_MEMP','Tumour','Tumour?','pDC','Plasma_cell','NK','MK','MEMP_MEP')] = 'doublets'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('CD14_Mono')] = 'CD14_Mono'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('CD16_Mono')] = 'CD16_Mono'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('CMP_GMP')] = 'DC_precursor'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('pre_B_cells')] = 'pre_B_cells'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('pro_B_cells')] = 'pro_B_cells'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('B_cells')] = 'B_cells'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('DC1')] = 'DC1'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('DC2')] = 'DC2'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(1,10,12,20,23,25,26,27,35,36,37,39,43,44,46,48,49,50,52,55,56,58,9,8) & mlds.ball.merged$cluster_ann %in% c('HSC_MPP')] = 'HSC_MPP'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(10) & mlds.ball.merged$cluster_ann_merged %in% c('10','Tumour')] = 'B_cells'
# 
# 13,15,29,30,32,34,38,4,42,45,47,5,53,
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(50) & mlds.ball.merged$cluster_ann_merged %in% c('50')] = 'HSC_MPP'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(5) & mlds.ball.merged$cluster_ann_merged %in% c('5')] = 'NK'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(45) & mlds.ball.merged$cluster_ann_merged %in% c('45')] = 'MK'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(42) & mlds.ball.merged$cluster_ann %in% c('24',47,'CD14_Mono','CD16_Mono','30','DC1')] = 'CD16_Mono'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(42) & mlds.ball.merged$cluster_ann %in% c('weird_MEMP','pro_B_cells','Neutrophil','pDC')] = 'doublets'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(42) & mlds.ball.merged$cluster_ann_merged %in% c('42')] = 'CD16_Mono'
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(4) & mlds.ball.merged$cluster_ann %in% c('CD16_Mono')] = 'CD16_Mono'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(4) & mlds.ball.merged$cluster_ann %in% c('CD14_Mono','CMP_GMP','23','24','27')] = 'CD14_Mono'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(4) & mlds.ball.merged$cluster_ann %in% c('DC1')] = 'DC1'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(4) & mlds.ball.merged$cluster_ann %in% c('DC2')] = 'DC2'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(4) & mlds.ball.merged$cluster_ann %in% c('CD4_Tcell','34','B_cells')] = 'doublets'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(4) & mlds.ball.merged$cluster_ann_merged %in% c('4')] = 'CD14_Mono'
# 
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(38,34) & mlds.ball.merged$cluster_ann_merged %in% c('4')] = 'CD14_Mono'
# 
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(15) & mlds.ball.merged$cluster_ann %in% c('19','22','23','29','NK_T','7','6','18','17','14')] = 'CD4_Tcell'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(15) & mlds.ball.merged$cluster_ann %in% c('11','20','4','34','B_cells','CD14_Mono','pro_B_cells','pre_B_cells','NK','NK_T','Tumour?')] = 'doublets'
# 
# 
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(29,34,38,35) & mlds.ball.merged$cluster_ann %in% c('weird_MEMP')] = 'MEMP_MEP'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(29,34,38,35) & mlds.ball.merged$cluster_ann %in% c('pro_B_cells')] = 'pro_B_cells'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(29,34,38,35) & mlds.ball.merged$cluster_ann %in% c('NK')] = 'NK'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(29,34,38,35) & mlds.ball.merged$cluster_ann %in% c('pDC','6')] = 'pDC'
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(29,34,38,35) & mlds.ball.merged$cluster_ann %in% c('CD14_Mono','CD16_Mono','DC2','DC1','HSC_MPP','ILC','MEMP_MEP')] = mlds.ball.merged$cluster_ann[mlds.ball.merged$seurat_clusters %in% c(29,34,38,35) & mlds.ball.merged$cluster_ann %in% c('CD14_Mono','CD16_Mono','DC2','DC1','HSC_MPP','ILC','MEMP_MEP')]
# mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(29,34,38,35) & mlds.ball.merged$cluster_ann %in% c('Tumour?','Ery','EE','B_cells','46','36','44','34','33','30','27','19','13','0')] = 'doublets'
# 
# 
# DimPlot(mlds.ball.merged,cells.highlight = mlds.ball.merged$cellID[mlds.ball.merged$seurat_clusters %in% c(29,34,38,35) & 
#                                                                      mlds.ball.merged$cluster_ann %in% c('0')])
# 
# View(table(mlds.ball.merged$cluster_ann_merged[mlds.ball.merged$seurat_clusters %in% c(10)],
#       mlds.ball.merged$cluster_ann[mlds.ball.merged$seurat_clusters %in% c(10)]) 
# )
# 
# table(mlds.ball.merged$cluster_ann[mlds.ball.merged$seurat_clusters %in% c(38)])

keyMarkers_2 = unique(c('HLA-DRA','CLEC9A','CD34','CD38','HLF','SPINK2','MLLT3','PRSS57', # HSC_MPP
                        'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1',#MEMP
                        'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
                        #'KIT',
                        'CTSG',	'PRTN3', # CMP
                        'AZU1','MPO','FLT3','PTPRC', # GMP
                        'ZBTB16','LTB', 'CD52',# ILC precursor
                        'IL7R','CD3D','GZMA','CD4','CD8A', #Early.lymphoid_T.lymphocyte
                        'TRDV2','TRGV9', # gamma delta T-cell
                        'SLC4A10','TRAV1-2', #MAIT t-cell
                        'PRF1', # effector T cells
                        'FOXP3',	'CDH1', # regulatory T cells
                        'NKG7','KLRD1', #NK
                        'IGLL1','CD99', # Pre-pro
                        'DNTT','CD79B','VPREB1','EBF1','CD19','RAG1',# pro-b-cells
                        'MME','CD79A',# pre-b-cells
                        'TCL1A','MME','RAG1','MS4A1',  # B-cells
                        'ITGB2','SELL','ITGAM','CD14','CCR2','FCGR3A','CX3CR1',# Monocytes
                        'CD68','MSR1',#Macrophages
                        'IRF8',	'CLEC10A', #DC.precursor
                        'ITGAX','CD1C','CLEC9A','THBD',#DC1
                        'CLEC10A', # DC2 
                        'IL3RA', 'CLEC4C',#pDC
                        'CD27',#plasma cell?
                        'CSF2RB','HDC','SERPINB1','TPSAB1','KIT', # Mast.cell
                        'PF4','ITGA2B', #Megakaryocyte
                        'PPBP','TUBB1', # Platelets
                        'GATA1','KLF1', # early Erythroid
                        'ALAS2', # mid.erythroid
                        'HBA1','BPGM', # late.erythroid
                        'CSF3R','FPR1','FCGR3B','NAMPT','MNDA', # NEUTROPHILS
                        'DEFA3','DEFA4','CAMP', 'LCN2',
                        'ESAM','PECAM1', 'KDR', 'PTPRB', 'PLVAP', # Endothelium
                        'DCN','SERPINF1','COL3A1','BGN','COL1A1','VIM','ECM1', # endosteal fibroblast / fibroblast
                        'APOA1','SCD','ALB','TTR' # Hepatocyte
))
# DotPlot(mlds,#idents = c(23),
#         group.by = 'finalAnn_broad',scale = F,
#         #group.by = 'cluster_ann_merged', scale=F,
#         features = keyMarkers_2)+RotatedAxis()
# 
# 
# 
# write.csv(mlds.ball.merged@meta.data,'MLDS_scRNAseq/1_annotation_scRNA_MLDS/may23/MLDS_DSLeuk_clean_annotated_may23_tmp.csv')
# 
# sampleSheet = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'MLDS_GOSH')
# sampleSheet = sampleSheet[!is.na(sampleSheet$Disease),]
# sampleSheet = sampleSheet[!is.na(sampleSheet$dataLocation),]
# sampleSheet = sampleSheet[sampleSheet$assay == "5' V2 Dual Index",]
# 
# m = match(mlds.ball.merged$donorID,sampleSheet$donorID)
# table(is.na(m))
# mlds.ball.merged$disease = sampleSheet$Disease[m]
# 
# saveRDS(mlds.ball.merged,file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/may23/MLDS_DSLeuk_clean_noMTCells_annotated_may23.RDS'))
# 
# 
# 
# m = match(tgt.srat$cellID,mlds.ball.merged$cellID)
# table(is.na(m))
# tgt.srat$cluster_ann = mlds.ball.merged$cluster_ann_merged[m]
# tgt.srat$finalAnn_broad = tgt.srat$cluster_ann
# tgt.srat$finalAnn_broad[tgt.srat$cluster_ann %in% c('T_cell','CD8_Tcell','CD4_Tcell','gd_Tcell')] = 'T_cell'
# tgt.srat$finalAnn_broad[tgt.srat$cluster_ann %in% c('LE','ME','ME_LE','EE')] = 'Ery'
# #tgt.srat$finalAnn_broad[tgt.srat$cluster_ann %in% c('30','31','?')] = '?'
# DimPlot(tgt.srat,group.by = 'finalAnn_broad',cols = c(col25[-6],col22),label = T,label.box = T,repel = T) + NoLegend()
# DimPlot(tgt.srat,group.by = 'seurat_clusters',label = T,label.box = T,repel = T) + NoLegend()
# DimPlot(mlds.ball.merged,cells.highlight = tgt.srat$cellID[tgt.srat$seurat_clusters == 23 & tgt.srat$finalAnn_broad == 'Tumour'])
# DimPlot(mlds.ball.merged,cells.highlight = mlds$cellID[mlds$finalAnn_broad == 'unknown'])
# saveRDS(tgt.srat,file.path(outDir,'DSLeuk_clean_noMTCells_annotated_may23.RDS'))
# 
# ## 4. Merge with MLDS
# mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/feb23/MLDS_clean_LRwCTannotated_feb23.RDS')
# big.srat = merge_seurat_objects(mlds,tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# big.srat = standard_clustering(big.srat)
# 
# big.srat$cluster_ann[match(tgt.srat$cellID,big.srat$cellID)] = tgt.srat$cluster_ann
# big.srat$finalAnn_broad[match(tgt.srat$cellID,big.srat$cellID)] = tgt.srat$finalAnn_broad
# big.srat$finalAnn_broad_level2 = as.character(big.srat$finalAnn_broad)
# big.srat$finalAnn_broad_level2[grepl('Tumour',big.srat$finalAnn_broad) & big.srat$donorID %in% c('L027','L062')] = 'Tumour_B.ALL'
# big.srat$finalAnn_broad_level2[grepl('Tumour',big.srat$finalAnn_broad) & !big.srat$donorID %in% c('L027','L062')] = 'Tumour_MLDS'
# big.srat$finalAnn_broad_level2[grepl('yelocyte',big.srat$finalAnn_broad)] = 'Myelocyte'
# big.srat$finalAnn_broad_level2[grepl('lowCnt|\\?|doublet',big.srat$finalAnn_broad)] = '?'
# DimPlot(big.srat,group.by = 'finalAnn_broad_level2',label = T,label.box = T,label.size = 4,repel = T,cols=sample(c(col25[-6],col22),n_distinct(big.srat$finalAnn_broad))) + NoLegend()
# DimPlot(big.srat,group.by = 'Phase',pt.size = 0.8,label = F,label.box = T,label.size = 4,repel = T,cols=col25[-6])
# DimPlot(big.srat,label = T,label.box = T,repel = T) + NoLegend()
# DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$finalAnn_broad=='MEMP_MEP'])
# 
# saveRDS(big.srat,file.path(outDir,'MLDS_DSLeuk_clean_annotated_may23.RDS'))
# 
# ### Subclustering MK lineage ####
# mkLin = subset(big.srat,subset = finalAnn_broad %in% unique(big.srat$finalAnn_broad[grepl('Tumour|MK|MEMP',big.srat$finalAnn_broad)]))
# mkLin = standard_clustering(mkLin)
# mkLin$finalAnn_broad_level2 = big.srat$finalAnn_broad_level2[match(mkLin$cellID,big.srat$cellID)]
# mkLin$ann = paste0(mkLin$finalAnn_broad,':',mkLin$donorID)
# DimPlot(mkLin,label = T,label.box = T,repel = T) + NoLegend()
# DimPlot(mkLin,group.by = 'finalAnn_broad_level2',label = T,label.box = T,label.size = 3,repel = T,cols=c(col25[-6],col22)) + NoLegend()
# DimPlot(mkLin,group.by = 'donorID',label = T,label.box = T,repel = T,label.size = 3,cols=c(col25[-6],col22))
# DimPlot(mkLin,group.by = 'Phase',label = F,label.box = T,repel = T,label.size = 3,cols=c(col25[-6],col22))
# DimPlot(mkLin,cells.highlight = mkLin$cellID[mkLin$donorID=='L019'])
# 
# m.15.19 = FindMarkers(mkLin,ident.1 = 15,ident.2 = 19)
# m.15.19$geneSym = rownames(m.15.19)
# m.15.19 = m.15.19[m.15.19$p_val_adj<0.05 & !grepl('^RPL|^RPS|MALAT1',m.15.19$geneSym),]
# 
# FeaturePlot(mkLin,c('HSP90B1','HSP90AB1','LY6E','PPIB','NRGN','TUBB1','PPBP','MYL9'),ncol = 4)
# FeaturePlot(mkLin,c('HSP90B1'))
# DotPlot(mkLin,group.by = 'finalAnn_broad_level2',features = c('FAM162A', 'SLC25A37', 'DDIT4', 'BNIP3', 'KRT18', 'P4HB', 'HSP90B1', 'GYPB', 'MT1X', 'HBA1'))
# DotPlot(mkLin,group.by = 'finalAnn_broad_level2',features = c('TUBA4A', 'TUBB1', 'PPBP', 'PTCRA', 'C19orf33', 'SH3BP5', 'MYL9', 'CCL5', 'ACRBP', 'HIST1H2AC', 'MYLK', 'GFI1B', 'NFE2', 'GP1BA'))
# 
# 
# 
# 
# 
# 
# 
# ##################################################################################################################################
# 
# ##---------------------------------##
# ##  Preprocessing Brain Lymphoma ####
# ##---------------------------------##
# 
# outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/nov23/DS_brainLymph') 
# plotDir = file.path(outDir,'dsBrain_')
# outPath = file.path(outDir,'dsBrain')
# 
# cleanSrat_fp = ifelse(keepMTCells,paste0(outDir,'/dsBrain_clean_withMTCells.RDS'),paste0(outDir,'/dsBrain_clean_noMTCells.RDS'))
# if(file.exists(cleanSrat_fp) & skipIfExists){
#   cleanSrat = readRDS(cleanSrat_fp)  
# }else{
#   if(!dir.exists(outDir)){
#     message(sprintf('Creating output directory'))
#     dir.create(outDir,recursive = T)
#   }
#   
#   setwd(outDir)
#   
#   # Find brain Lymphoma samples 
#   dataDirs = list.dirs('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/',recursive = T,full.names = T)
#   dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs) & grepl('SB_NB14406184|SB_NB14406185',dataDirs)]
#   names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',dataDirs))
#   
#   
#   dataDirs=dataDirs[file.exists(dataDirs)]
#   
#   dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
#   print(n_distinct(dataDirs))
#   
#   metadata = NULL
#   matchBy = NULL
#   cleanCountDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/nov23/cleanCount'
#   # Run basicQC
#   plotDir = file.path(outDir,'dsBrain_')
#   outPath = file.path(outDir,'dsBrain')
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
#   df.out = QC.output[[2]]
#   
#   write.csv(df.out,paste0(outDir,'/dsBrain_qc_summary.csv'))
# }
# 
# 
# 
# 
# 
# 
# 
# ##--------------------------------------------##
# ##       Feb 2024: Annotation - dsBrain     ####
# ##--------------------------------------------##
# outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24'
# if(!dir.exists(outDir)){
#   dir.create(outDir)
# }
# 
# 
# ## 1. Import object
# tgt.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/nov23/DS_brainLymph/dsBrain_clean_noMTCells.RDS')
# 
# ## 2. Add metadata
# tgt.srat$assay = 'GEX5p'
# tgt.srat$donorID = 'L062'
# tgt.srat$age_yrs = '?'
# tgt.srat$sex = 'M'
# tgt.srat$tissue = 'Brain'
# tgt.srat$mutation = '?'
# tgt.srat$timePoint = '?'
# tgt.srat$blastPerc = '?'
# tgt.srat$clinicalOutcome = '?'
# tgt.srat$dataset = 'GOSH'
# tgt.srat$disease = 'Lymphoma'
# tgt.srat$Genotype = 'T21'
# 
# DimPlot(tgt.srat,group.by = 'seurat_clusters',label.box = T,label = T,repel = T)
# DimPlot(tgt.srat,group.by = 'Phase',label.box = T,label = F,repel = T,cols = col25)
# 
# 
# library(SoupX)
# qm = quickMarkers(tgt.srat@assays$RNA@counts,tgt.srat$seurat_clusters)
# 
# ## 3. Import DS-BALL
# dsBALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_clean_annotated_0923.RDS')
# dsBALL$Genotype = ifelse(dsBALL$donorID %in% c('L027','L062'),'T21','diploid')
# dsBALL$dataset = 'GOSH'
# dsBALL$disease = 'BALL'
# 
# ## 4. Combine DS-BALL and BrainLymphoma
# dsBALL = merge_seurat_objects(dsBALL,tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# dsBALL = standard_clustering(dsBALL)
# 
# apply(dsBALL@meta.data,2,function(x){sum(is.na(x))})
# 
# DimPlot(dsBALL,group.by = 'cluster_ann',label = T,repel = T,label.box = T)+NoLegend()
# 
# 
# 
# ## 5. Combine DS-BALL and MLDS
# mlds@meta.data = mlds@meta.data[,!colnames(mlds@meta.data) %in% c('group','RNA_snn_res.2','RNA_snn_res.10')]
# colnames(mlds@meta.data)[!colnames(mlds@meta.data) %in% colnames(dsBALL@meta.data)]
# 
# big.srat = merge_seurat_objects(mlds,dsBALL,keepAllGenes = F,genomeVersions = c('v38','v38'))
# big.srat = standard_clustering(big.srat)
# 
# DimPlot(big.srat,group.by = 'disease',label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
# DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$disease == 'Lymphoma'])
# 
# 
# ## 6. Combine DS-BALL and MLDS and MDS
# big.srat = merge_seurat_objects(big.srat,mds,keepAllGenes = F,genomeVersions = c('v38','v38'))
# big.srat = standard_clustering(big.srat)
# 
# DimPlot(big.srat,group.by = 'annot_feb24',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
# DimPlot(dsBALL,group.by = 'cluster_ann',label.size = 3,label = T,repel = T,label.box = T,cols = c(col25,pal34H))+NoLegend()
# DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$disease == 'Lymphoma'])
# 
# big.srat$annot_feb24 = as.character(big.srat$annot_jan24)
# 
# 
# ##--------     Annotate    ---------
# annotDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24'
# if(!dir.exists(annotDir)){
#   dir.create(annotDir,recursive = T)
# }
# 
# # saveRDS(mlds,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_clean_annotated_noUnknowns_jan24_tmp.RDS')
# # mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_clean_annotated_noUnknowns_jan24_tmp.RDS')
# 
# 
# 
# 
# 
# 
# 
# DimPlot(mlds,group.by = 'annot_feb24',cols = c(col25,pal34H),label = T,repel = T,label.box = T)+NoLegend()
# DimPlot(mlds,cells.highlight = mlds$cellID[mlds$donorID == 'L156'])
# 
# 
# ## Annotate each cluster with NA cells
# for(clust in unique(as.character(big.srat$seurat_clusters[is.na(big.srat$annot_feb24) | big.srat$annot_feb24 == 'NA']))){
#   nNewCells = length(big.srat$cellID[is.na(big.srat$annot_feb24) & as.character(big.srat$seurat_clusters) == clust])
#   #nNewCells = length(big.srat$cellID[big.srat$annot_feb24 == 'NA' & as.character(big.srat$seurat_clusters) == clust])
#   nTot = length(big.srat$cellID[as.character(big.srat$seurat_clusters) == clust]) 
#   p = DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$seurat_clusters == clust]) + ggtitle(clust) + NoLegend() 
#   p1 = DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$seurat_clusters == clust & is.na(big.srat$annot_feb24)]) + ggtitle(clust,subtitle = 'new cells')+ NoLegend()
#   #p1 = DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$seurat_clusters == clust & big.srat$annot_feb24 == 'NA']) + ggtitle(clust,subtitle = 'new cells')+ NoLegend()
#   print(p+p1)
#   
#   current_annot = unique(big.srat$annot_feb24[!is.na(big.srat$annot_feb24) & as.character(big.srat$seurat_clusters) == clust])
#   #current_annot = unique(big.srat$annot_feb24[big.srat$annot_feb24 != 'NA' & as.character(big.srat$seurat_clusters) == clust])
#   print(sprintf('Cluster %s: %d new cells out of %d total cell count, representing a fraction of %f',clust,nNewCells,nTot,nNewCells/nTot))
#   if(length(current_annot) == 1){
#     assign = readline(sprintf('Cluster %s: current annot is %s. Assigned cluster annot to this now? ', clust, current_annot))
#   }else{
#     print(sprintf('Cluster %s current annotation:',clust))
#     print(table(big.srat$annot_feb24[!is.na(big.srat$annot_feb24) & as.character(big.srat$seurat_clusters) == clust]))
#     #print(table(big.srat$annot_feb24[big.srat$annot_feb24 != 'NA' & as.character(big.srat$seurat_clusters) == clust]))
#     assign = readline(sprintf('Cluster %s: Assigned cluster annot now? ', clust))
#   }
#   
#   if(assign == 'n'){
#     big.srat$annot_feb24[is.na(big.srat$annot_feb24) & big.srat$seurat_clusters == clust] = 'NA'
#     #big.srat$annot_feb24[big.srat$annot_feb24 =='NA' & big.srat$seurat_clusters == clust] = 'NA'
#   }else if(assign == 'y'){
#     big.srat$annot_feb24[is.na(big.srat$annot_feb24) & big.srat$seurat_clusters == clust] = current_annot
#     #big.srat$annot_feb24[big.srat$annot_feb24 =='NA' & big.srat$seurat_clusters == clust] = current_annot
#   }else{
#     big.srat$annot_feb24[is.na(big.srat$annot_feb24) & big.srat$seurat_clusters == clust] = assign
#     #big.srat$annot_feb24[big.srat$annot_feb24 =='NA' & big.srat$seurat_clusters == clust] = assign
#   }
# }
# 
# big.srat$annot_feb24[big.srat$annot_feb24 == 'Cacner'] = 'Cancer'
# big.srat$annot_feb24[big.srat$annot_feb24 == ''] = 'NA'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'EE' ] = 'EE'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'ME' ] = 'ME'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NK_T' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'T_CD8' ] = 'T_CD8'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NK_T' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'T_gd' ] = 'T_gd'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'Mono_CD14' ] = 'Mono_CD14'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'MK' ] = 'MK'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'MEMP_MEP' ] = 'MEP'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'Cancer' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'doublets' ] = 'doublets'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NK_T' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'T_MAIT' ] = 'T_MAIT'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'HSC_MPP' ] = 'HSC_MPP'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'Mono_CD16' ] = 'Mono_CD16'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NK_T' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'NK' ] = 'NK'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'myelocyte' ] = 'Myelocytes'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NK' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'T_gd' ] = 'T_gd'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NK' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'LE' ] = 'LE'
# big.srat$annot_feb24[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann == 'T_CD8' ] = 'T_CD8'
# 
# View(as.data.frame(table(big.srat$annot_feb24[!is.na(big.srat$cluster_ann)],big.srat$cluster_ann[!is.na(big.srat$cluster_ann)])))
# DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$annot_feb24 == 'NA' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann =='T_CD8'])
# DimPlot(dsBALL,cells.highlight = big.srat$cellID[big.srat$annot_feb24 == 'naive.B' & !is.na(big.srat$cluster_ann) & big.srat$cluster_ann =='Tumour'])
# 
# 
# big.srat@meta.data = big.srat@meta.data[,colnames(big.srat@meta.data) != 'tmp']
# umap = big.srat@reductions$umap@cell.embeddings %>% as.data.frame()
# umap$cellID = big.srat$cellID
# mdat = merge(big.srat@meta.data,umap,by='cellID',all=T)
# write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_combined_feb24_mdat.csv')
# 
# mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_combined_feb24_mdat.csv')
# dsBALL$annot_feb24 = mdat$annot_feb24[match(dsBALL$cellID,mdat$cellID)]
# saveRDS(dsBALL,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/BALL_clean_annotated_feb24.RDS')
# 
# 
# 
# 
# 
# big.srat$tmp = paste0(big.srat$annot_feb24,':',big.srat$cluster_ann)
# Idents(big.srat) = big.srat$tmp
# DotPlot(big.srat,idents = c('NK:T_gd'),
#         features = keyMarkers_2)+RotatedAxis()



##
