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

refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
liftChain = '~/lustre_mt22/hg19ToHg38_noChr.over.chain'
nParallel=24


outDir = file.path('MLDS_scRNAseq/genotypeCheck')
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}


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
scrubPath='./scrubletScores.tsv'
scPath="./strainedCounts"
doPlot=F
verbose = T
skipIfExists=T
keepMTCells=T


# projMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'MLDS_GOSH')
# projMani = projMani[!is.na(projMani$Disease),]
# projMani = projMani[!is.na(projMani$dataLocation),]
# projMani = projMani[grepl('B-ALL',projMani$Disease),]
# projMani = projMani[projMani$assay == "5' V2 Dual Index",]

sampleSheet = read_excel('~/lustre_mt22/Data/B_ALL/Manifest_2023_09_13.xlsx',sheet = 'B-ALL (Study 6918_7153_3030)')
sampleSheet = sampleSheet[sampleSheet$Experiment == 'GEX' & 
                            grepl('AML',sampleSheet$Sample_ID),]
#sampleSheet = sampleSheet[!sampleSheet$Cellgen_ID %in% sampleSheet$Cellgen_ID[sampleSheet$Patient_ID %in% c('L007','L071','L024','L031','L096','L057','L045','L073','L033','L059','L043','L072','L023','L087','L085','L032','L004','L006','L025','L048','L003') & sampleSheet$Timepoint == 'D0' | sampleSheet$Subtype == 'Other'],]
#samplesToKeep = unique(c(projMani$sangerSampleID,sampleSheet$Cellgen_ID))

### 1. Import remapped data
###    Run SoupX
###    Subset to keep only cells present in the original publications
###    Add cell labels (as published)

outDir = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/sept23/AML') 
plotDir = file.path(outDir,'AML_')
outPath = file.path(outDir,'AML')

cleanSrat_fp = ifelse(keepMTCells,paste0(outDir,'/AML_clean_withMTCells.RDS'),paste0(outDir,'/AML_clean_noMTCells.RDS'))
if(file.exists(cleanSrat_fp) & skipIfExists){
  cleanSrat = readRDS(cleanSrat_fp)  
}else{
  if(!dir.exists(outDir)){
    message(sprintf('Creating output directory'))
    dir.create(outDir,recursive = T)
  }
  
  setwd(outDir)
  
  dataDirs = list.dirs('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/GEX/',recursive = T,full.names = T)
  dataDirs = dataDirs[grepl('filtered_feature_bc_matrix$',dataDirs)]
  
  names(dataDirs) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_','',dataDirs))
  
  dataDirs=dataDirs[file.exists(dataDirs)]
  dataDirs = dataDirs[names(dataDirs) %in% gsub('SB_','',sampleSheet$Cellgen_ID)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))
  
  metadata = NULL
  matchBy = NULL
  cleanCountDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/sept23/cleanCount'
  # Run basicQC
  plotDir = file.path(outDir,'AML_')
  outPath = file.path(outDir,'AML')
  
  message('\nPerforming scRNAseq QC...')
  QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,
                      clusteringRes=clusteringRes,cleanCountDir=cleanCountDir,
                      skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                      metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                      doPlot=doPlot,plotDir=plotDir,verbose=verbose)
  
  cleanSrat = QC.output[[1]]
  
  df.out = QC.output[[2]]
  
  write.csv(df.out,paste0(outDir,'/AML_0923_qc_summary.csv'))
}








##################################################################################################################################


#############################
##       Annotation      ####
#############################


#---------------------------------------------------------------------------------------------------------------------------------
outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML'
if(!dir.exists(outDir)){
  dir.create(outDir)
}


## 1. Import relevant objects
sampleSheet = read_excel('~/lustre_mt22/Data/B_ALL/Manifest_2023_09_13.xlsx',sheet = 'B-ALL (Study 6918_7153_3030)')
sampleSheet = sampleSheet[sampleSheet$Experiment == 'GEX' & 
                            grepl('AML',sampleSheet$Sample_ID),]


## 2. Import object
tgt.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/0_preprocessing_scRNA_MLDS/sept23/AML/AML_clean_noMTCells.RDS')
tgt.srat = standard_clustering(tgt.srat)
DimPlot(tgt.srat,label.box = T,label = T,repel = T)
DimPlot(tgt.srat,group.by = 'orig.ident',label.box = T,label = T,repel = T)

FeaturePlot(tgt.srat,c('GATA1','KLF1'))

## 2. Add metadata
m = match(tgt.srat$orig.ident,gsub('SB_','',sampleSheet$Cellgen_ID))
table(is.na(m))
tgt.srat$assay = 'GEX5p'
tgt.srat$donorID = sampleSheet$Patient_ID[m]
tgt.srat$age_yrs = 'paed'
tgt.srat$sex = sampleSheet$Sex[m]
tgt.srat$tissue = sampleSheet$Tissue[m]
tgt.srat$mutation = sampleSheet$Subtype[m]
tgt.srat$timePoint = sampleSheet$Timepoint[m]
tgt.srat$blastPerc = sampleSheet$`Flow_blast_
percent`[m]
tgt.srat$clinicalOutcome = sampleSheet$Response[m]

DimPlot(tgt.srat,group.by = 'seurat_clusters',label.box = T,label = T,repel = T)
DimPlot(tgt.srat,group.by = 'Phase',label.box = T,label = F,repel = T,cols = col25)





## 3. Broad annotation

DotPlot(tgt.srat,group.by = 'cluster_ann',
        features = unique(c('ANPEP','FUT4','CD33','CD64','CD36','VIM','HLA-DRA','CLEC9A','CD34','CD38','HLF','SPINK2','MLLT3','PRSS57', # HSC_MPP
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
        )))+RotatedAxis()

DimPlot(tgt.srat,cells.highlight = nk)
VlnPlot(tgt.srat,group.by = 'cluster_ann',features =c('percent.mt','nFeature_RNA','nCount_RNA'),ncol=1)

DimPlot(tgt.srat,group.by = 'cluster_ann',label = T,label.box = T,repel = T,cols = c(col25,pal34H)) 

tgt.srat$cluster_ann = as.character(tgt.srat$seurat_clusters)
tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(4,13,15,16,5,3,25)] = 'Tumour'
tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(30)] = 'MK'
tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(19,24)] = 'neutrophils'

tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(1,6,7,10,20,21)] = 'T/NK'
t_cd8 = WhichCells(tgt.srat,cells = tgt.srat$cellID[tgt.srat$cluster_ann == 'T/NK'],expression = (CD8A > 1))
nk = WhichCells(tgt.srat,cells = tgt.srat$cellID[tgt.srat$cluster_ann == 'T/NK'],expression = (CD3D > 0.5 & KLRD1 > 1))
tgt.srat$cluster_ann[tgt.srat$cellID %in% t_cd8] = 'T_CD8'
tgt.srat$cluster_ann[tgt.srat$cellID %in% nk] = 'NK'
tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(9,27)] = 'B.cell'
tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(36)] = 'doublets'

tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(24)] = 'CMP_GMP'
tgt.srat$cluster_ann[tgt.srat$seurat_clusters %in% c(32)] = 'MonoMac'

tgt.srat$cluster_ann = mlds.aml.merged$finalAnn[match(tgt.srat$cellID,mlds.aml.merged$cellID)]

write.csv(tgt.srat@meta.data,'MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML/AML_clean_noMTcells_annotated_0923_mdat.csv')


library(SoupX)
qm = quickMarkers(tgt.srat@assays$RNA@counts,tgt.srat$cluster_ann)
write.csv(qm,'MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_clean_annotated_0923_quickMarkers.csv')




##----------------------------##
##    Combine with MLDS     ####
##----------------------------##

mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/may23/MLDS_clean_LRwCTannotated_may23.RDS')

mlds.aml.merged = merge_seurat_objects(mlds,tgt.srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
mlds.aml.merged = standard_clustering(mlds.aml.merged)
table(is.na(mlds.aml.merged$cluster_ann),mlds.aml.merged$dataset)
mlds.aml.merged$cluster_ann[is.na(mlds.aml.merged$cluster_ann)] =  mlds.aml.merged$finalAnn_broad[is.na(mlds.aml.merged$cluster_ann)]

DimPlot(mlds.aml.merged, group.by = 'finalAnn', label = T,label.box = T,cols = c(col25[-6],col22,pal34H),repel = T) + NoLegend()
DimPlot(mlds.aml.merged, group.by = 'seurat_clusters',label = T,label.box = T,repel = T) + NoLegend()



DotPlot(mlds.aml.merged,idents = c(8),
        group.by = 'finalAnn',scale = T,
        #group.by = 'cluster_ann_merged', scale=F,
        features = unique(c('HLA-DRA','CLEC9A','CD34','CD38','HLF','SPINK2','MLLT3','PRSS57', # HSC_MPP
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
                            'DEFA3','DEFA4','CAMP', 'LCN2', #myelocytes
                            'ESAM','PECAM1', 'KDR', 'PTPRB', 'PLVAP', # Endothelium
                            'DCN','SERPINF1','COL3A1','BGN','COL1A1','VIM','ECM1', # endosteal fibroblast / fibroblast
                            'APOA1','SCD','ALB','TTR' # Hepatocyte
        )))+RotatedAxis()

qm = quickMarkers(tgt.srat@assays$RNA@counts,tgt.srat$seurat_clusters)
DimPlot(tgt.srat,cells.highlight = mlds.aml.merged$cellID[mlds.aml.merged$seurat_clusters %in% c(48) & mlds.aml.merged$finalAnn %in% c('31')])

DimPlot(mlds.aml.merged,cells.highlight = tgt.srat$cellID[tgt.srat$seurat_clusters %in% c('40','38')])

mlds.aml.merged$finalAnn = mlds.aml.merged$cluster_ann

table(mlds.aml.merged$cluster_ann[mlds.aml.merged$seurat_clusters %in% c(8)])
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(19,16,34) & mlds.aml.merged$cluster_ann == 'l076_bm'] = 'B.cell' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(1,32) & mlds.aml.merged$cluster_ann == 'l076_bm'] = 'LE' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(26) & mlds.aml.merged$cluster_ann %in% c('l076_bm','37')] = 'ME' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(15) & mlds.aml.merged$cluster_ann %in% c('l076_bm','37')] = 'EE' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(29,47) & mlds.aml.merged$cluster_ann %in% c('18')] = 'Tumour' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(29,47) & mlds.aml.merged$cluster_ann %in% c('T_CD8','T/NK')] = 'doublets' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(24) & mlds.aml.merged$cluster_ann %in% c('25','3','5')] = 'Tumour' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(2,23,2,13,12,10,0) & mlds.aml.merged$cluster_ann %in% c('l076_bm')] = 'Tumour' 
t_cd8 = WhichCells(mlds.aml.merged,cells = mlds.aml.merged$cellID[mlds.aml.merged$seurat_clusters %in% c(2,23,2,13,12,10,0) & mlds.aml.merged$cluster_ann %in% c('l076_bm','0','18','36','31')],expression = (CD3D > 0.5 & CD8A > 0.5))
t_cell = WhichCells(mlds.aml.merged,cells = mlds.aml.merged$cellID[mlds.aml.merged$seurat_clusters %in% c(2,23,2,13,12,10,0) & mlds.aml.merged$cluster_ann %in% c('l076_bm','0','18','36','31')],expression = (CD3D > 0.5 & CD8A < 0.5))
mlds.aml.merged$finalAnn[mlds.aml.merged$cellID %in% t_cd8] = 'T_CD8'
mlds.aml.merged$finalAnn[mlds.aml.merged$cellID %in% t_cell] = 'T.cell'

nk = WhichCells(mlds.aml.merged,cells = mlds.aml.merged$cellID[mlds.aml.merged$seurat_clusters %in% c(2,23,2,13,12,10,0,4,20,33) & mlds.aml.merged$cluster_ann %in% c('l076_bm','0','18','36','31','12','36')],expression = (CD3D < 0.5 & KLRD1 > 0.5))
mlds.aml.merged$finalAnn_broad[mlds.aml.merged$cellID %in% nk] = 'NK'

mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(11,17,45,6) & mlds.aml.merged$cluster_ann %in% c('l076_bm','5','22','8')] = 'Mono_CD14' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(11,17,45,6) & mlds.aml.merged$cluster_ann %in% c('28')] = 'Mono_CD16' 

mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(46,44) & mlds.aml.merged$cluster_ann %in% c('18')] = 'Tumour' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(30) & mlds.aml.merged$cluster_ann %in% c('11','12','14','8') & mlds.aml.merged$donorID == 'L058'] = 'Tumour' 

mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(14) & mlds.aml.merged$cluster_ann %in% c('l076_bm')] = 'pro.B.cell' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(14) & mlds.aml.merged$cluster_ann %in% c('26','31')] = 'doublets' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(8) & mlds.aml.merged$cluster_ann %in% c('29','8')] = 'doublets' 
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(8) & !mlds.aml.merged$cluster_ann %in% c('29','8')] = 'Tumour' 

t_cd8 = WhichCells(mlds.aml.merged,cells = mlds.aml.merged$cellID[mlds.aml.merged$finalAnn %in% c('l076_bm','T/NK')],expression = (CD3D > 0.5 & CD8A > 0.5))
t_cell = WhichCells(mlds.aml.merged,cells = mlds.aml.merged$cellID[mlds.aml.merged$finalAnn %in% c('l076_bm','T/NK')],expression = (CD3D > 0.5 & CD8A < 0.5))
nk = WhichCells(mlds.aml.merged,cells = mlds.aml.merged$cellID[mlds.aml.merged$finalAnn %in% c('l076_bm','T/NK')],expression = (CD3D < 0.5 & KLRD1 > 0.5))
mlds.aml.merged$finalAnn[mlds.aml.merged$cellID %in% t_cd8] = 'T_CD8'
mlds.aml.merged$finalAnn[mlds.aml.merged$cellID %in% t_cell] = 'T_CD4'
mlds.aml.merged$finalAnn[mlds.aml.merged$cellID %in% nk] = 'NK'

mlds.aml.merged$finalAnn[mlds.aml.merged$finalAnn %in% c('T/NK','36','l076_bm') & mlds.aml.merged$seurat_clusters %in% c(10,23)] = 'T_CD8'
mlds.aml.merged$finalAnn[mlds.aml.merged$finalAnn %in% c('T/NK','36','l076_bm') & mlds.aml.merged$seurat_clusters %in% c(2,12,0,13)] = 'T_CD4'
mlds.aml.merged$finalAnn[mlds.aml.merged$seurat_clusters %in% c(48) & mlds.aml.merged$finalAnn == '31'] = 'Plasma.cell'
mlds.aml.merged$finalAnn[mlds.aml.merged$cellID %in% tgt.srat$cellID[tgt.srat$seurat_clusters %in% c(40,38,31)]] = 'Plasma.cell'


saveRDS(mlds.aml.merged,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML/AML_mergedMLDS_clean_noMTcells_annotated_0923.RDS')


tgt.srat$cluster_ann = mlds.aml.merged$finalAnn[match(tgt.srat$cellID,mlds.aml.merged$cellID)]
tgt.srat$finalAnn_broad = tgt.srat$cluster_ann

write.csv(tgt.srat@meta.data,'MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML/AML_clean_noMTcells_annotated_0923_mdat.csv')
saveRDS(tgt.srat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML/AML_clean_noMTcells_annotated_0923.RDS')



