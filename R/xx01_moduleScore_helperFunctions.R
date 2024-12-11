##--- [Helper functions] Module scoring in bulk / sc / sn RNAseq dataset ---##


##----------------##
##   Libraries  ####
##----------------##
library(tidyverse)
library(Seurat)
library(UCell)
source('~/lustre_mt22/generalScripts/utils/pseudobulk.R')

scRNA_moduleScore_UCell = function(geneList,outDir=NULL,module_type,datasets = c('fLiver','fLiver_Muz','MLDS','MDS','BALL','fBM_2n','fBM_T21','infantALL','pAML','otherLeuk'),ncores=1){
  known_datasets = c('fLiver','fLiver_Muz','MLDS','MDS','BALL','fBM_2n','fBM_T21','infantALL','pAML','otherLeuk')
  
  ## check that datasets is NOT null
  if(is.null(datasets)){
    stop('Please provide some datasets for module scoring')
  }
  
  
  if(!all(datasets %in% known_datasets)){
    warning(sprintf('The following datasets are not recognised - will ignore them for now: %s',paste(datasets[!datasets %in% known_datasets],collapse = ', ')))
    datasets = datasets[datasets %in% known_datasets]
  }
  
  if(is.null(outDir)){
    stop('Please provide output directory')
  }
  
  if(!dir.exists(outDir)){
    dir.create(outDir,recursive = T)
  }
  
  
  ##-------------------------------##
  ##    Score module in fLiver   ####
  ##-------------------------------##
  if('fLiver' %in% datasets){
    print('Scoring fLiver...')
    # fLiver = filter_sratObj(tissue = 'liver',ageMatch=T,remove_cyclingCells=F)
    # fLiver = fLiver[[1]]

    akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'
    akLiv_mdat_fp = gsub('_0824.RDS','_0824_mdat.csv',akLiv_srat_fp)
    
    fLiver = readRDS(akLiv_srat_fp)
    fLiver$annot = fLiver$annot_aug24
    fLiver$finalAnn_broad = fLiver$annot_aug24
    
    
    fLiver$tissue = 'fLiver'
    fLiver$dataset = 'fLiver'
    set.seed(12234)
    # cells_toKeep_fLiver = c(fLiver$cellID[fLiver$finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','EE','earlyMK','MK','Mast.cell')],
    #                         #sample(fLiver$cellID[fLiver$finalAnn_broad == 'B.cell' & fLiver$Genotype == 'T21'],1000),
    #                         fLiver$cellID[fLiver$finalAnn_broad == 'B.cell' & fLiver$Genotype == 'T21'],
    #                         sample(fLiver$cellID[fLiver$finalAnn_broad == 'B.cell' & fLiver$Genotype == 'diploid'],1000),
    #                         sample(fLiver$cellID[fLiver$finalAnn_broad == 'Endo' & fLiver$Genotype == 'T21'],1000),
    #                         sample(fLiver$cellID[fLiver$finalAnn_broad == 'Endo' & fLiver$Genotype == 'diploid'],1000),
    #                         sample(fLiver$cellID[fLiver$finalAnn_broad == 'NK_T'& fLiver$Genotype == 'T21'],1000),
    #                         sample(fLiver$cellID[fLiver$finalAnn_broad == 'NK_T'& fLiver$Genotype == 'diploid'],1000))
    # table(fLiver$Genotype[fLiver$cellID %in% cells_toKeep_fLiver])
    # fLiver = subset(fLiver,subset = cellID %in% cells_toKeep_fLiver)
    
    fLiver <- UCell::AddModuleScore_UCell(fLiver, features = geneList,ncores = ncores)
    
    write.csv(fLiver@meta.data,file.path(outDir,paste0('fLiver_',module_type,'_moduleScore.csv')))
    
    rm(fLiver)
    print('Scoring fLiver... - COMPLETED!')
    
  }
  
  
  if('fLiver_Muz' %in% datasets){
    print('Scoring published sc fLiver atlas ...')
    #fLiver = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/0_reprocessing_fetalREF/liver/oct22/liver_clean_withMTCells.RDS')
    fLiver = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver_2n/jan24/liver_liverREFmerged_clean_processed_annotated_noUnknowns_0124.RDS')
    ## Only keep published cells
    fLiver = subset(fLiver,subset = cellID %in% fLiver$cellID[fLiver$published_ann_2 != 'NA'])
    fLiver$tissue = 'fLiver'
    fLiver$dataset = 'fLiverMuz'
    
    set.seed(12234)
    
    fLiver <- UCell::AddModuleScore_UCell(fLiver, features = geneList,ncores = ncores)
    
    write.csv(fLiver@meta.data,file.path(outDir,paste0('fLiverMuz_',module_type,'_moduleScore.csv')))
    
    rm(fLiver)
    print('Scoring fLiver... - COMPLETED!')
    
  }
  
  ##-----------------------------##
  ##    Score module in MLDS   ####
  ##-----------------------------##
  if('MLDS' %in% datasets){
    cat('\n')
    print('Scoring MLDS...')
    # mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_231127.RDS')
    # # reannotated L038 tumour cells
    # mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_mdat_231221')
    # mlds$annot_dec23 = mdat$annot_dec23[match(mlds$cellID,mdat$cellID)]
    # mlds$finalAnn_broad = mlds$annot_dec23
    
    #mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_clean_annotated_noUnknowns_jan24.RDS')
    #mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
    mlds_srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS'
    mlds_mdat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv'
    
    mlds = readRDS(mlds_srat_fp)
    mlds$annot = mlds$annot_aug24
    mlds$finalAnn_broad = as.character(mlds$annot_aug24)
    
    set.seed(12234)
    cells_toKeep = c(mlds$cellID[mlds$finalAnn_broad %in% c('EE','HSC_MPP','MEP','MK','MK_WT','Tum_MK?')],
                     mlds$cellID[grepl('Tumour',mlds$finalAnn_broad)],
                     sample(mlds$cellID[mlds$finalAnn_broad == 'naive.B'],1000),
                     sample(mlds$cellID[mlds$finalAnn_broad == 'Mono_CD14'],1000),
                     sample(mlds$cellID[mlds$finalAnn_broad == 'NK'],1000))
    mlds = subset(mlds,subset = cellID %in% cells_toKeep)
    
    
    mlds <- UCell::AddModuleScore_UCell(mlds, features = geneList,ncores = ncores)
    
    write.csv(mlds@meta.data,file.path(outDir,paste0('MLDS_',module_type,'_moduleScore.csv')))
    
    rm(mlds)
    print('Scoring MLDS - COMPLETED!')
  }
  
  ##------------------------------------------##
  ##    Score module in "otherLeukaemia"    ####
  ##------------------------------------------##
  if('otherLeuk' %in% datasets){
    cat('\n')
    print('Scoring Other Leukaemia...')
    ## other leukaemia
    otherLeuk = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2408.RDS')
    
    otherLeuk <- UCell::AddModuleScore_UCell(otherLeuk, features = geneList,ncores = ncores)
    
    write.csv(otherLeuk@meta.data,file.path(outDir,paste0('otherLeuk_',module_type,'_moduleScore.csv')))
    
    rm(otherLeuk)
    print('Scoring Other Leukaemia... - COMPLETED!')  
  }
  
  
  ##------------------------------##
  ##    Score module in MDS    ####
  ##-----------------------------##
  if('MDS' %in% datasets){
    cat('\n')
    print('Scoring MDS...')
    ## MDS
    mds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MDS/MDS_clean_annotated_tmp.RDS')
    mds$tissue = 'BM'
    mds$dataset = 'MDS'
    
    mds <- UCell::AddModuleScore_UCell(mds, features = geneList,ncores = ncores)
    
    write.csv(mds@meta.data,file.path(outDir,paste0('MDS_',module_type,'_moduleScore.csv')))
    
    rm(mds)
    print('Scoring MDS... - COMPLETED!')  
  }
  
  
  
  
  ##-----------------------------##
  ##    Score module in bALL   ####
  ##-----------------------------##
  if('BALL' %in% datasets){
    cat('\n')
    print('Scoring gosh B-ALL...')
    ## Import goshBALL
    #bALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/BALL/BALL_clean_annotated_0923.RDS')
    bALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/BALL_clean_annotated_jan24.RDS')
    bALL$finalAnn_broad = bALL$annot_feb24
    bALL$dataset = 'goshBALL'
    #bALL$Genotype = ifelse(grepl('risomy',bALL$mutation),'T21','diploid')
    bALL = subset(bALL,subset = cellID %in% bALL$cellID[bALL$finalAnn_broad %in% c('HSC_MPP','MEMP_MEP','MEP','MK','EE','NK','Mono_CD14','Tumour','Cancer')])
    
    bALL <- UCell::AddModuleScore_UCell(bALL, features = geneList,ncores = ncores)
    
    write.csv(bALL@meta.data,file.path(outDir,paste0('bALL_',module_type,'_moduleScore.csv')))
    
    rm(bALL)
    
    print('Scoring gosh B-ALL... - COMPLETED!')  
  }
  
  
  
  
  ##--------------------------------##
  ##    Score module in fBM-2n    ####
  ##--------------------------------##
  if('fBM_2n' %in% datasets){
    cat('\n')
    print('Scoring fBM 2n...')
    ## Import fBM object
    fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_diploid_Laura_sratObj_HARM_0923.RDS'
    fBM = readRDS(fBM_fp)
    fBM$tissue = 'fBM'
    fBM$dataset = 'fBM'
    cells_toKeep = c(fBM$cellID[fBM$finalAnn_broad %in% c('pre B progenitor','pro B progenitor','pre pro B progenitor','immature B cell','naive B cell','ELP','LMPP','HSC/MPP',
                                                          'endothelium',
                                                          'CD56 bright NK','mature NK','NK progenitor','NK T cell','transitional NK cell',
                                                          'early erythroid','earlyMK','early MK','MK','MEMP','MEP','mast cell')])
    
    fBM = subset(fBM, subset = cellID %in% cells_toKeep)
    
    fBM <- UCell::AddModuleScore_UCell(fBM, features = geneList,ncores = ncores)
    
    write.csv(fBM@meta.data,file.path(outDir,paste0('fBM_',module_type,'_moduleScore.csv')))
    
    rm(fBM)
    print('Scoring fBM 2n... - COMPLETED')  
  }
  
  
  
  
  ##---------------------------------##
  ##    Score module in fBM-T21    ####
  ##---------------------------------##
  if('fBM_T21' %in% datasets){
    cat('\n')
    print('Scoring fBM T21...')
    # Import fetal bone marrow T21
    t21_fBM_fp = '~/lustre_mt22/REF_datasets/fBM_Jardine21/fBM_T21_Laura_sratObj_HARM_0923.RDS'
    t21_fBM = readRDS(t21_fBM_fp)
    t21_fBM$tissue = 'fBM'
    t21_fBM$dataset = 't21_fBM'
    cells_toKeep = c(t21_fBM$cellID[t21_fBM$finalAnn %in% c('early B cell','pre B cell','mature B cell','pre B progenitor','pro B progenitor','pre pro B progenitor','immature B cell','naive B cell','ELP','LMPP','HSC/MPP',
                                                            'endothelium',
                                                            'CD56 bright NK','mature NK','NK progenitor','NK T cell','transitional NK cell',
                                                            'early erythroid','earlyMK','early MK','MK','MEMP','MEP','mast cell')])
    t21_fBM = subset(t21_fBM, subset = cellID %in% cells_toKeep)
    
    t21_fBM <- UCell::AddModuleScore_UCell(t21_fBM, features = geneList,ncores = ncores)
    
    write.csv(t21_fBM@meta.data,file.path(outDir,paste0('fBM.T21_',module_type,'_moduleScore.csv')))
    
    rm(t21_fBM)
    print('Scoring fBM T21... - COMPLETED')  
  }
  
  
  
  
  ##-------------------------------------##
  ##    Score module in Ellie's ALL    ####
  ##-------------------------------------##
  if('infantALL' %in% datasets){
    cat('\n')
    print('Scoring B-aLL... ')
    ## Add AMKL + B-ALL
    ## Ellie's infant ALL paper - 8 diploid ALL + 1 AMKL
    infantALL = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/ek12_infantALL/ek12_infantALL_clean_annotated.RDS')
    infantALL$finalAnn_broad = infantALL$annot
    infantALL$tissue = 'BM'
    infantALL$dataset = 'infantALL'
    #infantALL.cancer = subset(infantALL,subset = cellID %in% infantALL$cellID[infantALL$finalAnn_broad == 'Cancer'])
    
    infantALL <- UCell::AddModuleScore_UCell(infantALL, features = geneList,ncores = ncores)
    
    write.csv(infantALL@meta.data,file.path(outDir,paste0('infantALL_',module_type,'_moduleScore.csv')))
    
    rm(infantALL)
    print('Scoring b-ALL ... - COMPLETED')  
  }
  
  
  
  ##-----------------------------##
  ##    Score module in AML    ####
  ##-----------------------------##
  if('pAML' %in% datasets){
    cat('\n')
    print('Scoring pAML...')
    ## Add AML
    pAML = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/AML/AML_clean_noMTcells_annotated_0923.RDS')
    #pAML$tissue = 'BM'
    pAML$dataset = 'pAML'
    
    pAML <- UCell::AddModuleScore_UCell(pAML, features = geneList,ncores = ncores)
    
    write.csv(pAML@meta.data,file.path(outDir,paste0('pAML_',module_type,'_moduleScore.csv')))
    
    rm(pAML)
    
    print('Scoring pAML... - COMPLETED')  
  }
  
  
}



import_UCell_result = function(outDir,module_type){
  ## Import UCELL score across different dataset
  data = data.frame()
  for(f in list.files(outDir,full.names = T,pattern = paste0(module_type,'_*moduleScore.csv'))){
    if(!grepl('csv$',f)){next}
    d = read.csv(f)
    print(unique(d$dataset))
    d$dataset_ID = gsub('_.*$','',basename(f))
    if(unique(d$dataset) %in% c('fLiver','fLiverMuz')){
      d$sex = d$Sex
      d$disease = 'fLiver'
    }else if(unique(d$dataset) == 'MLDS'){
      d$finalAnn_broad = d$annot
      d$disease = d$disease
    }else if(unique(d$dataset) == 'infantALL'){
      d$donorID = d$patient_cancer
      d$sex = 'NA'
      d$Genotype = 'diploid'
      d$disease = 'infantALL'
    }else if(unique(d$dataset) == 'MDS'){
      d$donorID = 'L067'
      d$sex = 'NA'
      d$finalAnn_broad = d$finalAnn
      d$Genotype = 'diploid'
      d$disease = 'MDS'
    }else if(unique(d$dataset) == 'pAML'){
      d$Genotype = 'diploid'
      d$disease = 'pAML'
    }else if(unique(d$dataset) == 'otherLeuk'){
      mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_withUnknown_2408_mdat.csv',row.names = 1)
      d = cbind(d[,c('orig.ident','cellID','dataset_ID',colnames(d)[grepl('UCell',colnames(d))])],mdat[match(d$cellID,mdat$cellID),!colnames(mdat) %in% colnames(d)])
      d$finalAnn_broad = d$annot
      #d$finalAnn_broad = mdat$annot[match(d$cellID,mdat$cellID)]
    }else if(unique(d$dataset %in% c('fBM','t21_fBM'))){
      d$disease = 'fBM'
    }
    d = d[,c('cellID',"orig.ident","donorID","sex",'tissue',"finalAnn_broad",'Genotype','dataset','dataset_ID','disease',colnames(d)[grepl('UCell',colnames(d))])]
    if(nrow(data) == 0){
      data = d
    }else{
      data = rbind(data,d)
    }
  }
  
  return(data)
}





import_HenningsBulkSamples = function(gns=NULL,rm_henning_healthy=F,oxfordBulk=T,filter_lowExpr_gene=T,
                                      tpm_fp = c('~/lustre_mt22/Down_Leukemia/Down_Leukemia_Klusmann_StarSalmon_tpm.csv')){
  library(edgeR)
  if(is.null(gns)){
    # Create geneMap 
    library(GenomicFeatures)
    #Define genomic coordinates
    gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf.gz'
    txdb = makeTxDbFromGFF(gtf)
    gns = genes(txdb)
  }
  
  
  ##------  Create Gene Map -------##
  rawCnt = read_excel('~/lustre_mt22/Down_Leukemia/Down_Leukemia_Klusmann_Dec22.xlsx',sheet = 'StarSalmon_Counts')
  bulk_geneMap = rawCnt[,colnames(rawCnt) %in% c('gene_id','gene_name','EffectiveLength')]
  m = match(bulk_geneMap$gene_id,gns$gene_id)
  bulk_geneMap$chr = 'unknown'
  bulk_geneMap$chr[!is.na(m)] = as.character(seqnames(gns[m[!is.na(m)]]))
  
  
  
  ##------  Import Hennings bulk samples metadata -------##
  message('Importing bulk samples metadata....')
  ## Import bulk sample metadata
  mdat = read_excel('~/lustre_mt22/Down_Leukemia/Down_Leukemia_Klusmann_Dec22.xlsx',sheet = 'Annotation')
  mdat = mdat[,c('Sample','Source','Subgroup')]
  mdat$Sample = gsub('-','.',mdat$Sample)
  
  # Remove healthy control from Hennings data
  if(rm_henning_healthy){
    mdat = mdat[!grepl('healthy',mdat$Sample),] 
  }
  
  # remove samples with problems - from Konstantin email
  toExclude = c('Patient_TMD_1020','Patient_TMD_1037','Patient_TMD_1303','Patient_TMD_BCKR',
                'Patient_TMD_KHLR','Patient_MLL_1569','Patient_MLL_1569_2','Patient_MLL_615',
                'Patient_MLL_DEE7','Patient_MLL_DEEF3','Patient_MLL_DEHS1','Patient_MLL_DEMK11','Patient_MLL_DES10',
                'Patient_MLL_642')
  mdat = mdat[!mdat$Sample %in% toExclude,]
  
  
  ## Import survival data ##
  #survival_mdat = read_excel('~/lustre_mt22/mt22/Down_Leukemia_Klusmann_Apr23.xlsx',sheet = 'QC_survival')
  survival_mdat = read_excel('~/lustre_mt22/Hennings_bulkRNAseq_metadata.xlsx')
  ## Only keep samples with event information
  #survival_mdat = survival_mdat[survival_mdat$Event != '-',]
  survival_mdat$Sample = gsub('-','.',survival_mdat$Sample)
  table(survival_mdat$Sample %in% mdat$Sample)
  
  # # Remove samples with high duplication rate
  # mdat = mdat[mdat$Sample %in% survival_mdat$Sample[as.numeric(survival_mdat$Duplication_rate) < 0.65 | survival_mdat$Duplication_rate == '-'],]
  survival_mdat = survival_mdat[survival_mdat$Sample %in% mdat$Sample,]
  # Add survival data to mdat
  mdat = merge(mdat,survival_mdat[,!grepl('\\.\\.',colnames(survival_mdat))],by=c('Sample','Source','Subgroup'),all=T)
  
  
  ##------------------------------------##
  ##   Compute bulkRNAdata CPM counts ####
  ##------------------------------------##
  rawCnt = column_to_rownames(rawCnt,var ='gene_id')
  rawCnt = rawCnt[,colnames(rawCnt) %in% mdat$Sample]
  
  ## Filter out lowly expressed genes
  bulk_dge = DGEList(counts = rawCnt, genes = rownames(rawCnt))
  cpmCnt = edgeR::cpm(bulk_dge)
  
  
  ##-----------------------------------##
  ##   Import bulkRNAdata TPM counts ####
  ##-----------------------------------##
  message('Importing TPM count....')
  
  if(filter_lowExpr_gene & !oxfordBulk){
    
    prop_expressed = rowMeans(edgeR::cpm(bulk_dge) > 1)
    keep = prop_expressed > 0.5
    op = par(no.readonly = TRUE)
    par(mfrow = c(1, 2))
    hist(edgeR::cpm(bulk_dge, log = TRUE), main = 'Unfiltered', xlab = 'logCPM')
    abline(v = log(1), lty = 2, col = 2)
    hist(edgeR::cpm(bulk_dge[keep, ], log = TRUE), main = 'Filtered', xlab = 'logCPM')
    abline(v = log(1), lty = 2, col = 2)
    par(op)
    
    ## Subset the count matrix
    bulk_dge = bulk_dge[keep, , keep.lib.sizes = FALSE]
    rawCnt = rawCnt[keep, ]
    
    
    ##---- Calculate TPM ----##
    
    geneLen = as.vector(bulk_geneMap$EffectiveLength[match(rownames(rawCnt),bulk_geneMap$gene_id)])

    rpk = apply(rawCnt[,!colnames(rawCnt) %in% c('Length','gene_name','Geneid')], 2, function(x) x/(as.vector(geneLen)/1000))
    #normalize by the sample size using rpk values
    tpm = apply(rpk, 2, function(x) x / (sum(as.numeric(x)) / 1e6))  %>% as.data.frame()
    colNames = colnames(tpm)
    tpm$ensID = rownames(rawCnt)
    tpm$geneLength = bulk_geneMap$EffectiveLength[match(tpm$ensID,bulk_geneMap$gene_id)]
    tpm_count = cbind(tpm[,c('ensID','geneLength',colNames)])
  }else{
    ## Import Hennings TPM count
    
    tpm_count = data.frame()
    
    for(file in tpm_fp){
      tpm = read.csv(file,row.names = 1)
      #tpm$countMethod = sapply(strsplit(basename(file),split='_'),'[',4)
      colnames(tpm) = gsub('-','.',colnames(tpm))
      colnames(tpm) = gsub('\\.markdup\\.sorted\\.bam','',colnames(tpm))
      rownames(tpm) = tpm$ensID
      if(ncol(tpm_count) == 0){
        tpm_count = tpm
      }else{
        m = match(colnames(tpm_count),colnames(tpm))
        print(sum(is.na(m)))
        tpm_count = rbind(tpm_count,tpm[,m])
      }
      
    }
    tpm_count = tpm_count[,!grepl('^X\\.',colnames(tpm_count))]  
  }
  
  # Check that all samples have got metadata
  #if(!all(colnames(tpm_count)[!colnames(tpm_count) %in% c('ensID','geneLength','countMethod','geneName')] %in% mdat$Sample)){
  if(!all(mdat$Sample %in% colnames(tpm_count))){
    stop(sprintf('Not all samples are found in metadata sheet. Please check!'))
  }
  # subset TMP count to include only samples in metadata
  tpm_count = tpm_count[,c('ensID',colnames(tpm_count)[colnames(tpm_count) %in% mdat$Sample])]
  
  
  
  
  ##--- Import oxford bulk dataset -----##
  if(oxfordBulk){
    ## 1. Meta data
    mdat_oxfordBulk = read.delim('~/lustre_mt22/Down_Leukemia/Oxford_flowSorted_Blineage_bulkRNA/GSE122982_series_matrix.txt',skip = 28,sep = '\t')
    mdat_oxfordBulk = as.data.frame(t(mdat_oxfordBulk[c(10,11),-1]))
    colnames(mdat_oxfordBulk) = c('celltype','uniqueID')
    mdat_oxfordBulk$uniqueID = gsub('unique identifier: ','',mdat_oxfordBulk$uniqueID)
    mdat_oxfordBulk$Source = gsub('population: ','',mdat_oxfordBulk$celltype)
    mdat_oxfordBulk$Subgroup = mdat_oxfordBulk$Source
    mdat_oxfordBulk$Sample = paste0(mdat_oxfordBulk$Source,'_',mdat_oxfordBulk$uniqueID)
    ## Merge metadata
    mdat_oxfordBulk[,colnames(mdat)[!colnames(mdat) %in% colnames(mdat_oxfordBulk)]] = 'NA'
    mdat = rbind(mdat,mdat_oxfordBulk[,colnames(mdat)])
    
    
    ## 2. TPM count
    
    ## Import Oxford bulk
    oxford_bulk = read.csv('~/lustre_mt22/Down_Leukemia/Oxford_flowSorted_Blineage_bulkRNA/GSE122982_FBM_TPM.csv')
    colnames(oxford_bulk)[1] = 'ensID'
    colnames(oxford_bulk) = c('ensID',mdat_oxfordBulk$Sample[match(colnames(oxford_bulk)[-1],mdat_oxfordBulk$uniqueID)])
    ## Merge TPM count table
    tpm_count = merge(tpm_count,oxford_bulk,by='ensID')
    tpm_count = column_to_rownames(tpm_count,'ensID')
  }
  
  return(list('tpm_count' = tpm_count,'raw_count' = rawCnt,'cpm_count'=cpmCnt,
              'mdat' = mdat,'bulk_geneMap'=bulk_geneMap))
}


bulkRNA_moduleScore_singScore = function(moduleList,module_type,gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T){
  library(edgeR)
  library(singscore)
  
  ##----- 1. Import bulk counts -----##
  bulkSamples = import_HenningsBulkSamples(gns=gns,rm_henning_healthy=rm_henning_healthy,oxfordBulk=oxfordBulk,filter_lowExpr_gene=filter_lowExpr_gene,
                                           tpm_fp = c('~/lustre_mt22/Down_Leukemia/Down_Leukemia_Klusmann_StarSalmon_tpm.csv'))
  tpm_count = bulkSamples[['tpm_count']]
  bulk_mdat = bulkSamples[['mdat']]
  bulk_geneMap = bulkSamples[['bulk_geneMap']]
  
  
  ##----- 3. Prepare mdat -----##
  bulk_mdat$os = bulk_mdat$OS
  bulk_mdat$os_group = '-'
  bulk_mdat$os_group[!is.na(bulk_mdat$os) & bulk_mdat$os != '-'] = ifelse(as.numeric(bulk_mdat$os[bulk_mdat$os != '-'])<2,'low',
                                                                          ifelse(as.numeric(bulk_mdat$os[bulk_mdat$os != '-'])<5,'mid',
                                                                                 ifelse(as.numeric(bulk_mdat$os[bulk_mdat$os != '-'])<10,'mid2','high')))
  bulk_mdat$category = ifelse(bulk_mdat$Subgroup %in% c('CMP','GMP','MEP','HSC','MPP','Bcell','ELP','LMPP','PreProB','ProB'),'Normal',
                              ifelse(bulk_mdat$Subgroup %in% c('TMD','MLDS'),'TAM / MLDS','Other leukaemia'))
  
  bulk_mdat$category = factor(bulk_mdat$category,c('Normal','TAM / MLDS','Other leukaemia'))
  bulk_mdat$sampleGroup = as.character(bulk_mdat$Subgroup)
  bulk_mdat$sampleGroup[grepl('^t',bulk_mdat$sampleGroup)] = 'MLL_rearrangement'
  bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD'] = 'TAM'
  bulk_mdat$sampleGroup = factor(bulk_mdat$sampleGroup,c('TAM','MLDS','AMKL','MLL','MLL_rearrangement',
                                                         'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))
  
  
  
  ##----- 4.  Score the module -----##
  mtx = tpm_count[,colnames(tpm_count) %in% bulk_mdat$Sample[!grepl('^PDX',bulk_mdat$Sample)]]
  ## apply the rankGenes method
  bulk_ranked = rankGenes(mtx)
  
  
  ## apply the scoring function
  allScore = data.frame()
  for(i in 1:length(moduleList)){
    if(names(moduleList)[i] == 'all'){
      if(module_type == 'GATA1s_topGenes'){
        moduleScores = simpleScore(bulk_ranked,
                                   upSet = moduleList[['TAM.MLDS.up']],
                                   downSet = moduleList[['TAM.MLDS.down']])  
      }else if(module_type == 'MLDS_topGenes'){
        moduleScores = simpleScore(bulk_ranked,
                                   upSet = moduleList[['MLDS_up']],
                                   downSet = moduleList[['MLDS_down']])  
      }else{
        print('Unrecognized module type, please check!')
      }
      
      # moduleScores = simpleScore(bulk_ranked,
      #                          upSet = moduleList[['all_up']],
      #                          downSet = moduleList[['all_down']])
      # moduleScores = simpleScore(bulk_ranked,
      #                            upSet = moduleList[['L075_down']],
      #                            downSet = moduleList[['L075_up']])
      moduleScores = moduleScores[,c('TotalScore', 'TotalDispersion')]
    }else{
      moduleScores = simpleScore(bulk_ranked,upSet = moduleList[[i]])
    }
    
    # create a dataframe with the data required: scores and sample group
    scoredf = merge(bulk_mdat,moduleScores,by.x='Sample',by.y=0)
    scoredf$moduleType = paste0(module_type,'_',names(moduleList)[i])
    
    ## Add to allScore
    allScore = rbind(allScore,scoredf)
  }
  
  print(table(allScore$moduleType))
  
  return(allScore)
}

