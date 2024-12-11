##-------- Prepare Reference foetal objects (fLiver, fAdr, fKidney)  --------##
## 1. Reprocess remapped foetal references scRNA data for kidney, liver, adrenal
## [230324]: Not completed - need to revise properly...



setwd('~/lustre_mt22/Aneuploidy/Results/')

#------------------------#
##      Libraries     ####
#------------------------#
library(tidyverse)
library(readxl)
library(Seurat)
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/misc.R")
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/sc_utils.R")
source("/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/sc_basicQC.R")



#---------------------------------------------------#
##    0. Checking the published annotation       ####
#---------------------------------------------------#

##----  Kidney  ------##
kidREF = readRDS('/lustre/scratch125/casm/team274sb/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/preProcess/scData/finalisedSeurat_fKid.RDS')
kidREF@meta.data$tissue = 'kidney'
kidREF@meta.data$percent.mt = 100*kidREF@meta.data$mtGenes


##----  Adrenal  ------##
adrREF = readRDS('/lustre/scratch125/casm/team274sb/mt22/abnormal_karyotypes/adrenal/adrREF.rds')
adrREF@meta.data$annot = adrREF@meta.data$cell_type
adrREF@meta.data$tissue = 'adrenal'
adrREF@meta.data$percent.mt = 100*adrREF@meta.data$mtGenes
adrREF = subset(adrREF, subset = annot %in% unique(adrREF@meta.data$annot[adrREF@meta.data$annot != 'Other']))

DimPlot(adrREF,group.by = 'annot', label = T)
DimPlot(adrREF,group.by = 'annot', label = T,label.box = T,repel = T,cols = col25[-6])
### Subclustering of adrenal Cortex
adr_cortex = subset(adrREF, subset = annot =='Cortex')
adr_cortex = standard_clustering(adr_cortex)
DimPlot(adr_cortex,group.by = 'sample_name', label = T,label.box = T,repel = T)
FeaturePlot(adr_cortex,c('CYP11B2', 'DACH1', 'ANO4')) # ZG (outermost layer)
FeaturePlot(adr_cortex,c('CYP17A1', 'CYP11B1', 'HSD3B2')) # ZF        
FeaturePlot(adr_cortex,c('CYB5A', 'SULT2A1', 'GSTA1')) # ZR     
library(SoupX)
markers = quickMarkers(adr_cortex@assays$RNA@counts,clusters = adr_cortex@meta.data$seurat_clusters)



##----  Liver  ------##
## Old non-SoupXed fetal Liver
# livREFold = readRDS('/lustre/scratch125/casm/team274sb/mt22/REF_datasets/Muz_fLiver_REF.rds') # 10X indexed GRCh38 1.2.0 reference
# livREF@meta.data$annot = livREF@meta.data$cell.labels
# livREF@meta.data$tissue = 'liver'
# livREF@meta.data$percent.mt = 100*livREF@meta.data$percent.mito
livREF = readRDS('~/lustre_mt22/REF_datasets/Muz_fLiver_SoupXed.RDS')

## To obtain a more fine-grained annotation,especially for the HSC_MPP population, we will use some information from "cell.labels_prePub" column
## "cell.labels" = published labels
a = as.data.frame(table(livREF@meta.data$cell.labels_prePub,livREF@meta.data$cell.labels))
a = a[a$Freq >0,]
a$Var1 = as.character(a$Var1)
a$Var2 = as.character(a$Var2)
View(a[a$Var1!=a$Var2,])


DimPlot(livREF,group.by = 'cell.labels',label = T)+NoLegend()
DimPlot(livREF,label = T,cells.highlight = rownames(livREF@meta.data[livREF@meta.data$cell.labels == 'Monocyte precursor',]))+NoLegend()
DimPlot(livREF,label = T,cells.highlight = rownames(livREF@meta.data[livREF@meta.data$cell.labels_prePub == 'promonocyte' & livREF@meta.data$cell.labels == 'Neutrophil-myeloid progenitor',]))+NoLegend()
DimPlot(livREF,label = F,cells.highlight = rownames(livREF@meta.data[livREF@meta.data$cell.labels_prePub == 'Early Erythroid' & livREF@meta.data$cell.labels == 'MEMP',]))+NoLegend()
DimPlot(livREF,cells.highlight = rownames(livREF@meta.data[livREF@meta.data$cell.labels.groupped == 'DC precursor',]))+NoLegend()
DimPlot(livREF,cells.highlight = rownames(livREF@meta.data[livREF@meta.data$cell.labels_prePub == 'macrophage',]))+NoLegend()
FeaturePlot(livREF,c('VCAM1','CD7','CD14','CD34'))

livREF@meta.data$cell.labels.groupped = as.character(livREF@meta.data$cell.labels_prePub)
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub == 'sinusoidal EC') & (livREF@meta.data$cell.labels == 'Endothelial cell')] = 'Endothelial cell'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub == 'MOP')] = 'GMP'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub == 'promyelocyte') & (livREF@meta.data$cell.labels == 'Neutrophil-myeloid progenitor')] = 'GMP'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub == 'MEMP') & (livREF@meta.data$cell.labels == 'HSC_MPP')] = 'MPP'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub %in% c('pro-B cell','pre-B cell'))] = 'pre_pro_B_cell'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub %in% c('ELP','LMPP'))] = 'LMPP_ELP'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub %in% c('myeloid DC progenitor','promyelocyte'))] = 'GMP'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub %in% c('macrophage'))] = 'MPP'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub %in% c('pDC progenitor'))] = 'MPP'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub %in% c('erythroid-MPP hybrid'))] = 'MPP'
livREF@meta.data$cell.labels.groupped[(livREF@meta.data$cell.labels_prePub %in% c('eo/baso/mast precursor'))] = 'Mast cell'


a = as.data.frame(table(livREF@meta.data$cell.labels.groupped,livREF@meta.data$cell.labels))
a = a[a$Freq >0,]

DimPlot(livREF,group.by = 'cell.labels.groupped',label = T,repel = T)+NoLegend()

saveRDS(livREF,'~/lustre_mt22/REF_datasets/Muz_fLiver_SoupXed.RDS')
write.csv(livREF@meta.data,'~/lustre_mt22/Aneuploidy/Results/0_reprocessing_fetalREF/Muz_fLiver_SoupXed_sratObj_metadata.csv')


##---- assess proportion of B lineage in fLiver -----##
a = livREF@meta.data %>% group_by(donorID,cell.labels,gestationalAge) %>% summarise(n=n()) 
a$lin = ifelse(a$cell.labels %in% c('B cell','Pre pro B cell','pre-B cell','pro-B cell'),'B',
               ifelse(a$cell.labels %in% c('Fibroblast','Endothelial cell','Hepatocyte'),'stromal','others'))
a = a[a$lin != 'stromal',] %>% group_by(donorID) %>% mutate(totalCell = sum(n))

a = a %>% group_by(donorID,lin,totalCell,gestationalAge) %>% summarise(nCell_lin = sum(n)) %>% mutate(frac = 100*nCell_lin/totalCell)
View(a)
ggplot(a[a$lin == 'B',],aes(donorID,frac))+
  geom_col(aes(fill=lin))+
  facet_wrap(vars(gestationalAge),scales = 'free_x')





#-----------------------------------------------------#
##    1. Reprocess remapped foetal REF dataset     ####
#-----------------------------------------------------#

##----- Set parameters ------##
maxMT = 30
minGenes = 300

maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup = F
scrubScoreMax = 0.5
scrubPath='../cleanCounts_2208/scrubletScores.tsv'
scPath="../cleanCounts_2208/strainedCounts"
doPlot=T
verbose = T
skipIfExists=T
keepMTCells=T


##----- ReQC ------##
###    Import remapped data
###    Run SoupX
###    Subset to keep only cells present in the original publications
###    Add cell labels (as published)

for(tissue in c('adrenal','kidney','liver')){
  outDir = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF',tissue,'oct22') 
  plotDir = file.path(outDir,paste0(tissue,'_'))
  outPath = file.path(outDir,paste0(tissue))
  
  cleanSrat_fp = ifelse(keepMTCells,paste0(outDir,'/',tissue,'_clean_withMTCells.RDS'),paste0(outDir,'/',tissue,'_clean_noMTCells.RDS'))
  if(file.exists(cleanSrat_fp) & skipIfExists){
    cleanSrat = readRDS(cleanSrat_fp)  
  }else{
    if(!dir.exists(outDir)){
      message(sprintf('Creating output directory for %s',tissue))
      dir.create(outDir,recursive = T)
    }
    
    if(tissue == 'liver'){
      minUMIs = 500  
    }else{
      minUMIs = 1000
    }
    
    setwd(outDir)
    
    dataDirs = list.files(paste0('~/lustre119_mt22/fetalREF_scRNAseq_cr3.0.2_v38-1.2.0/',tissue),full.names = T)
    dataDirs = paste0(dataDirs,'/filtered_feature_bc_matrix')
    names(dataDirs) = basename(dirname(dataDirs))
    names(dataDirs) = gsub('_','.',names(dataDirs))
    dataDirs=dataDirs[file.exists(dataDirs)]
    dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
    print(n_distinct(dataDirs))
    # for(i in 1:length(dataDirs)){
    #   if(!dir.exists(sprintf('%s/filtered_feature_bc_matrix',dataDirs[i]))){
    #     system(sprintf('mkdir %s/filtered_feature_bc_matrix',dataDirs[i]))  
    #     system(sprintf('mv %s/* %s/filtered_feature_bc_matrix/',dataDirs[i],dataDirs[i]))  
    #   }
    #   
    #   if(dir.exists(sprintf('%s/filtered_feature_bc_matrix',dataDirs[i]))){
    #     
    #     system(sprintf('mv %s/filtered_feature_bc_matrix/* %s/',dataDirs[i],dataDirs[i]))  
    #     system(sprintf('rm %s/filtered_feature_bc_matrix -r',dataDirs[i]))
    #   }
    # }
    
    
    #metadata = mani[,c("donorID","PDID","gestationalAge","Sex","Genotype","Tissue","chanelID","assay")]
    #metadata$orig.ident = metadata$chanelID
    metadata = NULL
    matchBy = NULL
    cleanCountDir = NULL
    # Run basicQC
    plotDir = file.path(outDir,paste0(tissue,'_'))
    outPath = file.path(outDir,paste0(tissue))
    
    QC.output = basicQC(dataDirs = dataDirs,maxMT = maxMT, minGenes=minGenes,minUMIs=minUMIs,maxBadFrac=maxBadFrac,numPCs=numPCs,clusteringRes=clusteringRes,
                        skipScrub=skipScrub,skipSoup=skipSoup,scrubScoreMax=scrubScoreMax,scrubPath=scrubPath,
                        metadata=metadata,matchBy=matchBy,scPath=scPath,outPath=outPath,skipIfExists=skipIfExists,
                        doPlot=doPlot,plotDir=plotDir,verbose=verbose)
    
    cleanSrat = QC.output[[1]]
    #cleanSrat = sratOld
    #srat.outPath = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results',tissue,paste0(tissue,'_clean.RDS'))
    #message(sprintf('[Tissue: %s] cleanSrat exists? %s',tissue,file.exists(srat.outPath)))
    #saveRDS(cleanSrat,srat.outPath)
    
    df.out = QC.output[[2]]
    qc.summary=rbind(qc.summary,df.out)
    
    write.csv(df.out,paste0('~/lustre_mt22/Aneuploidy/Results/0_reprocessing_fetalREF/',tissue,'/oct22/',tissue,'_qc_summary.csv'))
  }
  
  
  #### Add metadata to new object #####
  ## Key metadata includes: donorID, Genotype, Gestational_age, Sex, assay, termination
  ## Import normal samples manifest
  message('\n\nImporting normal samples manifest')
  normSampleMani = read_excel("/lustre/scratch125/casm/team274sb/mt22/projectManifest.xlsx",sheet = "Normal samples channelID")
  normSampleMani = normSampleMani[!is.na(normSampleMani$sampleID),]
  
  normSampleMani$donorID[normSampleMani$Tissue == 'Kidney'] = paste0(normSampleMani$donorID[normSampleMani$Tissue == 'Kidney'],'_kidney')
  normSampleMani$donorID[normSampleMani$Tissue == 'Adrenal'] = paste0(normSampleMani$donorID[normSampleMani$Tissue == 'Adrenal'],'_adrenal')
  normSampleMani$donorID[normSampleMani$Tissue == 'Liver'] = paste0(normSampleMani$donorID[normSampleMani$Tissue == 'Liver'],'_liver')
  
  normSampleMani$channelID[normSampleMani$Tissue == 'Kidney'] = paste0(normSampleMani$channelID[normSampleMani$Tissue == 'Kidney'],'_kidney')
  normSampleMani$channelID[normSampleMani$Tissue == 'Adrenal'] = paste0(normSampleMani$channelID[normSampleMani$Tissue == 'Adrenal'],'_adrenal')
  normSampleMani$channelID[normSampleMani$Tissue == 'Liver'] = paste0(normSampleMani$channelID[normSampleMani$Tissue == 'Liver'],'_liver')
  
  
  normSampleMani2 = read_excel("/lustre/scratch125/casm/team274sb/mt22/projectManifest.xlsx",sheet = "existingNormalSamples")
  normSampleMani2 = normSampleMani2[!is.na(normSampleMani2$donorID),]
  normSampleMani2$donorID[normSampleMani2$tissue == 'Kidney'] = paste0(normSampleMani2$donorID[normSampleMani2$tissue == 'Kidney'],'_kidney')
  normSampleMani2$donorID[normSampleMani2$tissue == 'Adrenal'] = paste0(normSampleMani2$donorID[normSampleMani2$tissue == 'Adrenal'],'_adrenal')
  normSampleMani2$donorID[normSampleMani2$tissue == 'Liver'] = paste0(normSampleMani2$donorID[normSampleMani2$tissue == 'Liver'],'_liver')
  
  m = match(normSampleMani$donorID,normSampleMani2$donorID)
  sum(is.na(m))
  normSampleMani = merge(normSampleMani,normSampleMani2[,1:5],by = 'donorID')
  normSampleMani$termination[is.na(normSampleMani$termination)] = '??'
  
  cleanSrat$channelID = as.character(cleanSrat$orig.ident)
  cleanSrat$channelID = paste0(cleanSrat$channelID,'_',tissue)
  m = match(gsub('_','.',cleanSrat$channelID),gsub('_','.',normSampleMani$channelID))  
  sum(is.na(m))
  if(length(m) > 0 & sum(is.na(m)) == 0){
    cleanSrat@meta.data$donorID = normSampleMani$donorID[m]
    cleanSrat@meta.data$Sex = normSampleMani$sex[m]
    cleanSrat@meta.data$assay = normSampleMani$technology[m]
    cleanSrat@meta.data$termination = normSampleMani$termination[m]
    cleanSrat@meta.data$gestationalAge = normSampleMani$age[m]
  }
  
  cleanSrat@meta.data$Genotype = 'diploid'
  
  
  
  #### Add published annotation to object ####
  
  published_ann_fp = paste0('~/lustre_mt22/Aneuploidy/Results/0_reprocessing_fetalREF/',tissue,'/oct22/',tissue,'_published_metadata.csv')
  if(!file.exists(published_ann_fp)){
    ### Import current annotation
    ref_srat_fp = ifelse(tissue == 'kidney','/lustre/scratch125/casm/team274sb/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/preProcess/scData/finalisedSeurat_fKid.RDS',
                         ifelse(tissue == 'liver','~/lustre_mt22/REF_datasets/Muz_fLiver_SoupXed.RDS',
                                '/lustre/scratch125/casm/team274sb/mt22/abnormal_karyotypes/adrenal/adrREF.rds'))
    ref.srat = readRDS(ref_srat_fp)
    ref.srat@meta.data$Genotype = 'diploid'
    if(tissue == 'adrenal'){
      ref.srat@meta.data$annot = ref.srat@meta.data$cell_type  
      ref.srat@meta.data$tissue = 'adrenal'
      ref.srat@meta.data$percent.mt = 100*ref.srat@meta.data$mtGenes
      
    }else if(tissue == 'kidney'){
      #ref.srat = readRDS('/lustre/scratch125/casm/team274sb/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/preProcess/scData/finalisedSeurat_fKid.RDS')
      ref.srat@meta.data$tissue = 'kidney'
      ref.srat@meta.data$percent.mt = 100*ref.srat@meta.data$mtGenes
      ref.srat@meta.data$Genotype = 'diploid'
    }#else if(tissue == 'liver'){
    
    #}
    
    
    mdat = ref.srat@meta.data
    
    write.csv(mdat,published_ann_fp)
  }else{
    mdat = read.csv(published_ann_fp)
    mdat = column_to_rownames(mdat,var = 'X')
  }
  
  # Match up channel IDs
  mdat$cellID = rownames(mdat)
  #cleanSrat = cleanSrat@misc$preQC
  
  if(tissue == 'adrenal'){
    mdat$cell_bc = gsub('^.*_','',rownames(mdat))  
    mdat$channelID = as.character(mdat$orig.ident)
    mdat$channelID[as.character(mdat$orig.ident) == 'babyAdrenal1'] = '?WSSS8012017'
    mdat$channelID[as.character(mdat$orig.ident) == 'babyAdrenal2'] = '?WSSS8011223'
    mdat$channelID[as.character(mdat$orig.ident) == 'cellranger302'] = gsub('^.*_WSSS','WSSS',mdat$cellID[as.character(mdat$orig.ident) == 'cellranger302'])
    mdat$channelID[as.character(mdat$orig.ident) == 'cellranger302'] = gsub('_GRCh38.*$','',mdat$channelID[as.character(mdat$orig.ident) == 'cellranger302'])
    mdat$channelID[as.character(mdat$orig.ident) == 'cellranger302'] = gsub('_','.',mdat$channelID[as.character(mdat$orig.ident) == 'cellranger302'])
    mdat$cellID = paste0(mdat$channelID,'_',mdat$cell_bc)
    
    cleanSrat$channelID = as.character(cleanSrat$orig.ident)
  }else if (tissue == 'kidney'){
    fetalID_toKeep = c('F16','F17','F35','F38','F41','F45')
    mdat = mdat[mdat$orig.ident %in% fetalID_toKeep,]
    #DimPlot(ref.srat,group.by = 'annot',label=T,repel=T)
    #DimPlot(ref.srat,cells.highlight = rownames(ref.srat@meta.data[!ref.srat@meta.data$orig.ident %in% fetalID_toKeep,]))
    
    cleanSrat$channelID = as.character(cleanSrat$orig.ident)
    cleanSrat$channelID[cleanSrat$orig.ident == '4834STDY7002876'] = 'F16_Kid_N_neg_1_1'
    cleanSrat$channelID[cleanSrat$orig.ident == '4834STDY7002875'] = 'F16_Kid_N_pos_1_1'
    cleanSrat$channelID[cleanSrat$orig.ident == '4834STDY7002881'] = 'F17_Kid_N_ldf_1_1'
    cleanSrat$channelID[cleanSrat$orig.ident == '4834STDY7002886'] = 'F17_Kid_N_neg_1_1'
    cleanSrat$channelID[cleanSrat$orig.ident == '4834STDY7002885'] = 'F17_Kid_N_pos_1_1'
    cleanSrat$channelID[cleanSrat$orig.ident == 'FCAImmP7462243'] = 'F35_KI_45N'
    cleanSrat$channelID[cleanSrat$orig.ident == 'FCAImmP7462242'] = 'F35_KI_45P'
    cleanSrat$channelID[cleanSrat$orig.ident == 'FCAImmP7528293'] = 'F38_KI_45N'
    cleanSrat$channelID[cleanSrat$orig.ident == 'FCAImmP7528292'] = 'F38_KI_45P'
    cleanSrat$channelID[cleanSrat$orig.ident == 'FCAImmP7555850'] = 'F41_KI_45N'
    cleanSrat$channelID[cleanSrat$orig.ident == 'FCAImmP7555849'] = 'F41_KI_45P'
    cleanSrat$channelID[cleanSrat$orig.ident == 'FCAImmP7579215'] = 'F45_KI_45N'
    cleanSrat$channelID[cleanSrat$orig.ident == 'FCAImmP7579214'] = 'F45_KI_45P'
    
    
  }else if(tissue == 'liver'){
    cleanSrat$channelID = as.character(cleanSrat$orig.ident)
    mdat$annot = as.character(mdat$annot)
  }
  
  
  
  cleanSrat@meta.data$cellID = rownames(cleanSrat@meta.data)
  
  cleanSrat@meta.data$cellID = gsub('-\\d*$','',cleanSrat@meta.data$cellID)
  cleanSrat@meta.data$cellID = paste0(cleanSrat@meta.data$channelID,'_',sapply(strsplit(cleanSrat@meta.data$cellID,split='_'),'[',2))
  #cleanSrat@meta.data$cellID_bc = sapply(strsplit(cleanSrat@meta.data$cellID,split='_'),'[',2)
  
  table(mdat$cellID %in% cleanSrat$cellID)
  table(mdat$channel[!mdat$cellID %in% cleanSrat$cellID])
  table(mdat$annot[!mdat$cellID %in% cleanSrat$cellID])
  table(cleanSrat$reasonForFail[gsub('-\\d.*$','',rownames(cleanSrat@meta.data)) %in% mdat$cellID[!mdat$cellID %in% cleanSrat$cellID]])
  message(sprintf('%d cells in current object are missing from new object',sum(!mdat$cellID %in% cleanSrat$cellID)))
  
  
  cleanSrat$published_ann = NA
  m = match(cleanSrat$cellID,mdat$cellID)
  sum(is.na(m))
  cleanSrat$published_ann[!is.na(m)] = mdat$annot[m[!is.na(m)]]
  library(Seurat)
  DimPlot(cleanSrat,group.by = 'published_ann')
  
  ## Look at the reason why these published cells were removed
  ## If it's bad clusters --> keep
  ## if it's doublets/low counts/low genes/highMT or multiple --> remove
  cleanSrat@meta.data$PASS_withMT = ifelse(cleanSrat@meta.data$scrubScore > 0.2,F,cleanSrat@meta.data$PASS_withMT)
  table(cleanSrat@meta.data$reasonForFail[cleanSrat@meta.data$PASS_withMT == F & !is.na(cleanSrat@meta.data$published_ann)])
  cellID_tokeep = c(rownames(cleanSrat@meta.data[!is.na(cleanSrat@meta.data$published_ann) & cleanSrat@meta.data$reasonForFail !='doublet',]),
                    rownames(cleanSrat@meta.data[is.na(cleanSrat@meta.data$published_ann) & cleanSrat@meta.data$PASS_withMT,]))
  cleanSrat$rownames = rownames(cleanSrat@meta.data)
  cleanSrat.sub = subset(cleanSrat,rownames %in% cellID_tokeep)
  cleanSrat.sub = subset(cleanSrat.sub,subset = scrubScore <= 0.2)
  if(tissue=='kidney'){
    ## remove cells from 5GEX channel
    cellID_tokeep = rownames(cleanSrat.sub@meta.data[!is.na(cleanSrat.sub@meta.data$channelID),])
    cleanSrat.sub = subset(cleanSrat.sub,rownames %in% cellID_tokeep)
  }
  cleanSrat.sub = standard_clustering(cleanSrat.sub)
  DimPlot(cleanSrat.sub,group.by = 'published_ann',label = T)
  DimPlot(cleanSrat.sub,group.by = 'Phase',label = T)
  DimPlot(cleanSrat.sub,group.by = 'channelID',label = T)
  DimPlot(cleanSrat.sub,label = T)
  FeaturePlot(cleanSrat.sub,c('PTPRC','HBB'))
  
  #cleanSrat.sub = subset(cleanSrat,rownames %in% srat$cellID)
  saveRDS(cleanSrat.sub,paste0(outPath,'_clean_withMTCells_annotated_tmp.RDS')) 
  
}



