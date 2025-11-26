## conda activate alleleIntegrator
## L075 - unmatched WGS analysis ##


donorID = 'L075'
outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L075'
setwd(outDir)

##-------------------------------------------##
##    Define CN regions - TP1 sample       ####
##-------------------------------------------##
cnSegs = GRanges(c('chr5','chr17','chr21'),
                 IRanges(c(rep(1,3)),c(1e9,27.5e6,1e9)))

# Define sample metadata
samples_manifest = data.frame(PDID=c('PD62336a'),
                              projectid=c(3030),
                              haveCNV=c(F),
                              timepoint=c('Diagnostic'),
                              donorID=c('L075'))

# Pipeline for UNMATCHED subs analysis
unmatched_wrapper = function(outDir,cnSegs=NULL,samples_manifest,donorID=NULL,sampleID_name='PDID',insilico_ref='PDv38is_wgs',
                             run_cgpVaf_step1=T,
                             cgpVaf_script = '~/lustre_mt22/Aneuploidy/scripts/activeScritps/WGS_analysis/x10_MLDS_cgpVaf.sh',
                             gender = 'female',genotype = 'T21',
                             intv = 10,
                             doPlot=T){
  
  if(!dir.exists(outDir)){
    dir.create(outDir,recursive = T)
  }
  
  if(sum(samples_manifest$haveCNV) > 0 & is.null(cnSegs)){
    stop('Please provide CNV info')
  }
  
  library(VariantAnnotation)
  library(tidyverse)
  library(GenomicRanges)
  source('~/lustre_mt22/generalScripts/caveman_readvcf.R')
  source('~/lustre_mt22/generalScripts/binom_mix_model_Tim.R')
  source('~/lustre_mt22/generalScripts/utils/wgs_analysis_helperFunctions.R')
  source('~/lustre_mt22/generalScripts/utils/misc.R')
  
  
  #-----------------------------------##
  #    Import Caveman output        ####
  #-----------------------------------##
  var4cgpVaf = import_multiple_caveman(samples_manifest,sampleID_name='PDID')
  # if(is.null(donorID)){
  #   if('donorID' %in% colnames(samples_manifest)){
  #     donorID = samples_manifest$donorID  
  #   }else{
  #     donorID = samples_manifest$PDID
  #   }
  #   
  # }
  
  if(run_cgpVaf_step1){
    ##--------------------------------------------------##
    ##      cgpVAF of the somatic variants            ####
    ##--------------------------------------------------##
    # write bed file of the position of interest
    # requried format:
    #   chr pos ref_allele alt_allele
    #   filename ends with .bed, delim='\t'
    
    bed = data.frame(snv_ID = unique(var4cgpVaf$snv_ID))
    bed = merge(bed,var4cgpVaf[,c('Chr','Pos','Ref','Alt','snv_ID')],by='snv_ID',all.x=T)
    bed = bed[!duplicated(bed),]
    nrow(bed) == n_distinct(var4cgpVaf$snv_ID)
    bed = bed[order(bed$Chr),]
    write_delim(bed[,colnames(bed) != 'snv_ID'],file.path(outDir,paste0(donorID,'_unmatched_cavemanVar.bed')),col_names = F,delim = '\t')
    
    
    ##--- prepare for cgpVAF
    cgpVAF_dir = file.path(outDir,c('cgpVaf_input','cgpVaf_output'))
    for(d in cgpVAF_dir){
      if(!dir.exists(d)){
        dir.create(d,recursive = T)
      }
    }
    
    # Copy relevant BAM files to cgpVAF input directory
    files_to_copy = c(file.path('/nfs/cancer_ref01/nst_links/live',samples_manifest$projectid,samples_manifest$PDID,paste0(samples_manifest$PDID,'.sample.dupmarked.bam')),
                      file.path('/nfs/cancer_ref01/nst_links/live',samples_manifest$projectid,samples_manifest$PDID,paste0(samples_manifest$PDID,'.sample.dupmarked.bam.bai')),
                      '/nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam',
                      '/nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam.bai')
    
    
    for(f in files_to_copy){
      system(sprintf('cp -s %s  %s',f,file.path(outDir,'cgpVaf_input')))
    }
    
    ##--- Run cgpVAF
    cmd = sprintf('bash %s %s %s %s %s %s',cgpVaf_script, donorID, paste(unique(samples_manifest$PDID,collapse='_')),
                  cgpVAF_dir[grepl('input',cgpVAF_dir)],cgpVAF_dir[grepl('output',cgpVAF_dir)],
                  file.path(outDir,paste0(donorID,'_unmatched_cavemanVar.bed')))
    print(cmd)
    #bash ~/lustre_mt22/Aneuploidy/scripts/activeScritps/WGS_analysis/x10_MLDS_cgpVaf.sh L038 PD60301a_PD61846a ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/cgpVaf_input ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/cgpVaf_output ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/L038_unmatched_cavemanVar.bed
    return()
  }
  
  
  
  ##---- Import cgpVAF output with read filters applied
  output_wFilter_fp = list.files(file.path(outDir,'cgpVaf_output'),pattern = '_snp_vaf.tsv$',full.names = T)
  output_noFilter_fp = list.files(file.path(outDir,'cgpVaf_output_nofilter'),pattern = '_snp_vaf.tsv$',full.names = T)
  cgpVaf_QC_og = import_cgpVAF_output(output_wFilter_file = output_wFilter_fp,
                                      caveman_samples=samples_manifest,sampleID_name='PDID',
                                      output_noFilter_file = output_noFilter_fp)
  ## Add annotation info from caveman
  columns_toAdd = c('Gene','Impact','Type','AAchange','Snp')
  cgpVaf_QC_og = cbind(cgpVaf_QC_og[,!colnames(cgpVaf_QC_og) %in% columns_toAdd],
                       var4cgpVaf[match(cgpVaf_QC_og$snv_ID,var4cgpVaf$snv_ID),columns_toAdd])
  
  

  ##----------------------------------------------------------------##
  ##      cgpVAF depth filter - consider sex chromosome           ####
  ##----------------------------------------------------------------##
  
  
  # Plot distribution of depth by chromosome
  if(all(grepl('chr',cgpVaf_QC_og$Chrom))){
    cgpVaf_QC_og$Chrom = factor(cgpVaf_QC_og$Chrom,levels =  paste0('chr',c(1:22,'X','Y')))
  }else{
    cgpVaf_QC_og$Chrom = factor(cgpVaf_QC_og$Chrom,levels =  c(1:22,'X','Y'))
  }
  
  if(doPlot){
    pdf(file.path(outDir,paste0(donorID,'_unmatched_analysis.pdf')),width = 10,height = 10)
  }
  cgpVaf_QC_og$sample = paste0(cgpVaf_QC_og$PDID,' - ',cgpVaf_QC_og$timePoint)
  
  p = ggplot(cgpVaf_QC_og,aes(Chrom,wFilter_sample_DEP))+
    geom_boxplot(outlier.colour = 'white')+
    facet_wrap(vars(sample),ncol=1)+
    theme_classic(base_size = 11) + ggtitle('cgpVAF DEPTH distribution from both samples',
                                            subtitle = paste0(n_distinct(cgpVaf_QC_og$snv_ID),' variants')) + xlab('')+
    theme(panel.border = element_rect(fill=F,colour = 'black',size=1),axis.line = element_blank(),
          strip.background = element_blank(),axis.text = element_text(color='black'),
          strip.text = element_text(size = 13))
  
  
  print(p)


  # Calculate the average depth of each variant across samples of the same individual
  # exclude the copy number regions in TP1 sample
  cgpVaf_QC_og_gr = df2granges(data=cgpVaf_QC_og, chr='Chrom', start='Pos', end='Pos',id='snv_ID',sampleID='PDID')
  cgpVaf_QC_og_TP1_CNsegs = subsetByOverlaps(cgpVaf_QC_og_gr[cgpVaf_QC_og_gr$PDID == 'PD61846a'],cnSegs)
  toExclude = c(which(cgpVaf_QC_og$Chrom == 'chr21'),
                which(cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID & cgpVaf_QC_og$sample == 'PD61846a'))
  data_forDepthFilter = cgpVaf_QC_og[-toExclude,]
  
  data_forDepthFilter = data_forDepthFilter %>% group_by(snv_ID,Chrom) %>% summarise(meanDEP = mean(wFilter_sample_DEP))
  colnames(data_forDepthFilter)[3] = 'sample_DEP'
  
  cgpVaf_var_avgDepth = depthFilter(data = data_forDepthFilter,gender,genotype,useChr=NULL,upperCutOff=F)
  
  table(cgpVaf_QC_og$Chrom[cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID])
  cgpVaf_QC_og$depthFilter = cgpVaf_var_avgDepth$depthFilter[match(cgpVaf_QC_og$snv_ID,cgpVaf_var_avgDepth$snv_ID)]
  
  ## Remove variants with no coverage in one of the samples
  noCovVar = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$wFilter_sample_DEP == 0]
  cgpVaf_QC_og$depthFilter[cgpVaf_QC_og$snv_ID %in% noCovVar] = 'tooLow_DEP'
  cgpVaf_QC_og$depthFilter[is.na(cgpVaf_QC_og$depthFilter)] = 'NA'
  
  table(cgpVaf_QC_og$depthFilter,cgpVaf_QC_og$sample)



  ##---------------------------------------------------------------##
  ##      cgpVAF fraction of high-quality read filter            ####
  ##---------------------------------------------------------------##
  
  cgpVaf_QC_og$highQualRatio = cgpVaf_QC_og$wFilter_sample_DEP / cgpVaf_QC_og$noFilter_sample_DEP
  cgpVaf_QC_og$highQualRatio_filter = ifelse(cgpVaf_QC_og$highQualRatio >= 0.75,'PASS','low')
  
  ## Remove variants with low high-quality read fraction in at least one of the samples
  table(cgpVaf_QC_og$highQualRatio_filter[cgpVaf_QC_og$depthFilter != 'tooLow_DEP'],cgpVaf_QC_og$sample[cgpVaf_QC_og$depthFilter != 'tooLow_DEP'])
  
  lowHiQualVar = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$highQualRatio_filter == 'low']
  cgpVaf_QC_og$highQualRatio_filter[cgpVaf_QC_og$snv_ID %in% lowHiQualVar] = 'low'
  
  
  cgpVaf_QC_og$reasonForFail = ifelse(cgpVaf_QC_og$depthFilter == 'tooLow_DEP' | cgpVaf_QC_og$wFilter_sample_DEP == 0,'low_depth',
                                      ifelse(cgpVaf_QC_og$highQualRatio_filter == 'low','low_highQualRatio','PASS'))
  table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$sample)
  
  write.csv(cgpVaf_QC_og,file.path(outDir,'cgpVaf_output',paste0(donorID,'_cgpVafQCog_depth.highQualRatio_filtered.csv')))






  ##----------------------------------##
  ##      Pindel filter             ####
  ##----------------------------------##
  # Including FAILED Pindel variants too
  
  indels = import_multiple_pindel(samples_manifest,sampleID_name='PDID',intv=intv)

  
  for(s in samples_manifest$PDID[samples_manifest$haveCNV==T]){
    # Remove copy number regions in TP1 sample only
    cnSegs_indels = subsetByOverlaps(indels[indels$PDID == s],cnSegs[cnSegs$PDID == s])
    indels = indels[!(indels$PDID == s & indels$varID %in% cnSegs_indels$varID)]  
  }
  

  # find mutations within 10bp of indels
  #cgpVaf_QC_og = read.csv(file.path(outDir,'cgpVaf_output/L038_cgpVafQCog_depth.highQualRatio_filtered.csv'))
  cgpVaf_QC = cgpVaf_QC_og[cgpVaf_QC_og$reasonForFail == 'PASS',]
  
  ## Convert cgpVaf_QC into Granges
  cgpVaf_QC_gr = GRanges(cgpVaf_QC$Chrom,IRanges(cgpVaf_QC$Pos,cgpVaf_QC$Pos),
                         snv_ID = cgpVaf_QC$snv_ID,sample=cgpVaf_QC$PDID,timepoint=cgpVaf_QC$timePoint)
  mcols(cgpVaf_QC_gr) = cbind(mcols(cgpVaf_QC_gr),cgpVaf_QC[match(paste0(cgpVaf_QC_gr$snv_ID,':',cgpVaf_QC_gr$sample,':',cgpVaf_QC_gr$timepoint),
                                                                  paste0(cgpVaf_QC_gr$snv_ID,':',cgpVaf_QC_gr$sample,':',cgpVaf_QC_gr$timepoint)),])
  length(cgpVaf_QC_gr) == nrow(cgpVaf_QC)
  cgpVaf_QC = cgpVaf_QC_gr
  
  # Identify variants overlapping with indels
  cgpVaf_QC_overlapped_wIndels = subsetByOverlaps(cgpVaf_QC,indels)
  
  cgpVaf_QC$overlapped_with_indels = ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC_overlapped_wIndels$snv_ID,T,F)
  table(cgpVaf_QC$overlapped_with_indels,cgpVaf_QC$reasonForFail,cgpVaf_QC$timepoint)
  
  cgpVaf_QC = as.data.frame(mcols(cgpVaf_QC))

  cgpVaf_QC_og$pindelFilter = '?'
  for(s in unique(cgpVaf_QC_og$PDID)){
    cgpVaf_QC_og$pindelFilter[cgpVaf_QC_og$PDID == s & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == s & cgpVaf_QC$overlapped_with_indels == F]]  = 'PASS'
    cgpVaf_QC_og$pindelFilter[cgpVaf_QC_og$PDID == s & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == s & cgpVaf_QC$overlapped_with_indels == T]]  = 'overlapped_wINDELs'
  }
  
  # cgpVaf_QC_og$pindelFilter = ifelse(cgpVaf_QC_og$PDID == 'PD60301a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$overlapped_with_indels == F],'PASS',
  #                                    ifelse(cgpVaf_QC_og$PDID == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$overlapped_with_indels == F],'PASS',
  #                                           ifelse(cgpVaf_QC_og$PDID == 'PD60301a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$overlapped_with_indels == T],'overlapped_wINDELs',
  #                                                  ifelse(cgpVaf_QC_og$PDID == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$overlapped_with_indels == T],'overlapped_wINDELs','removed'))))


  table(cgpVaf_QC_og$pindelFilter,cgpVaf_QC_og$sample)
  table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$pindelFilter,cgpVaf_QC_og$sample)
  
  write.csv(cgpVaf_QC,file.path(outDir,'cgpVaf_output',paste0(donorID,'_cgpVaf_depth.highQualRatio.Pindel_filtered.csv')))




  # Extract only PASS variants for germline filter
  cgpVaf_QC=cgpVaf_QC[cgpVaf_QC$overlapped_with_indels == F,]
  table(cgpVaf_QC$sample)


  ##------------------------------------------------------##
  ##      Binomial test for germline filter             ####
  ##------------------------------------------------------##
  
  # Run binomial filter
  # data is a data.frame where: rownames = snv_ID (chr1:pos_REF/ALT), required columns are refCnt, altCnt, DP (total depth)
  
  cgpVaf_QC_og$refCnt = cgpVaf_QC_og$wFilter_sample_DEP - cgpVaf_QC_og$wFilter_sample_MTR
  cgpVaf_QC_og$altCnt = cgpVaf_QC_og$wFilter_sample_MTR
  cgpVaf_QC_og$DP = cgpVaf_QC_og$wFilter_sample_DEP
  
  # Exclude SNPs from CN regions in TP1 sample
  cgpVaf_QC_og$CN_tot2major = '2:1'
  if(genotype == 'T21'){
    cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$Chrom == 'chr21'] = '3:2'
  }
  if(donorID == 'L038'){
    cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$PDID == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID[cgpVaf_QC_og_TP1_CNsegs$Chrom == 'chr5']] = '1:0'
    cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$PDID == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID[cgpVaf_QC_og_TP1_CNsegs$Chrom == 'chr17']] = '1:0'
  }
  
  
  ## Perform binomial test filter
  cgpVaf_QC_og.sub = snvBinomTest(gender=gender,data=cgpVaf_QC_og[cgpVaf_QC_og$pindelFilter == 'PASS',],
                                  genotype=genotype,pCut=0.05)
  table(cgpVaf_QC_og.sub$isGermline,cgpVaf_QC_og.sub$sample)
  table(cgpVaf_QC_og.sub$isGermline,cgpVaf_QC_og.sub$Chrom)
  
  cgpVaf_QC_og$isGermline = cgpVaf_QC_og.sub$isGermline[match(paste0(cgpVaf_QC_og$snv_ID,':',cgpVaf_QC_og$sample),
                                                              paste0(cgpVaf_QC_og.sub$snv_ID,':',cgpVaf_QC_og.sub$sample))]
  
  
  cgpVaf_QC_og$germlineFilter = ifelse(cgpVaf_QC_og$PDID == 'PD60301a' &
                                         cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$PDID == 'PD60301a' & cgpVaf_QC_og.sub$isGermline == T],'germline',
                                       ifelse(cgpVaf_QC_og$PDID == 'PD61846a' &
                                                cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$PDID == 'PD61846a' & cgpVaf_QC_og.sub$isGermline == T],'germline',
                                              ifelse(cgpVaf_QC_og$PDID == 'PD60301a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$PDID == 'PD60301a' & cgpVaf_QC_og.sub$isGermline == F],'likely_somatic',
                                                     ifelse(cgpVaf_QC_og$PDID == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$PDID == 'PD61846a' & cgpVaf_QC_og.sub$isGermline == F],'likely_somatic','others'))))
  
  
  table(cgpVaf_QC_og$germlineFilter,cgpVaf_QC_og$pindelFilter)


  ## Summarise reason for fails
  cgpVaf_QC_og$reasonForFail = ifelse(cgpVaf_QC_og$CN_tot2major != '2:1' & cgpVaf_QC_og$germlineFilter != 'germline','CN_region',
                                      ifelse(cgpVaf_QC_og$germlineFilter == 'germline','germline',
                                             ifelse(cgpVaf_QC_og$germlineFilter == 'likely_somatic','likely_somatic',
                                                    ifelse((!is.na(cgpVaf_QC_og$depthFilter) & cgpVaf_QC_og$depthFilter == 'tooLow_DEP') | cgpVaf_QC_og$wFilter_sample_DEP == 0,'low_depth',
                                                           ifelse(cgpVaf_QC_og$highQualRatio_filter == 'low','low_highQualRatio',
                                                                  ifelse(cgpVaf_QC_og$pindelFilter == 'overlapped_wINDELs','overlapped_wINDELs','others'))))))
  
  
  table(is.na(cgpVaf_QC_og$reasonForFail))
  table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$sample)
  table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$germlineFilter)


  p = ggplot(cgpVaf_QC_og,aes(sample,fill = reasonForFail))+
    geom_bar(position='fill')+
    scale_fill_manual(values = c(col25[1],grey(0.8),col25[-1]))+
    theme_classic(base_size = 13)+
    theme(panel.border = element_rect(fill=F,colour = 'black',size=1),axis.line = element_blank(),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+xlab('')
  
  print(p)
  
  write.csv(cgpVaf_QC_og,file.path(outDir,paste0(donorID,'_cgpVaf_allFiltersApplied.csv')))
  #cgpVaf_QC_og = read.csv(file.path(outDir,'L038_cgpVaf_allFiltersApplied.csv'))

  cgpVaf_QC = cgpVaf_QC_og[cgpVaf_QC_og$germlineFilter == 'likely_somatic',]
  cgpVaf_QC = cgpVaf_QC_og

  ##----------------------------------------------------##
  ##      Plot - number of shared variants            ####
  ##----------------------------------------------------##
  table(cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD60301a'] %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD61846a'])
  
  varList = list('D_germline' = cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline == T & cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$wFilter_sample_VAF > 0],
                 'TP1_germline' = cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline == T & cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$wFilter_sample_VAF > 0],
                 'D_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$reasonForFail == 'likely_somatic' & !is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline != T & cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$wFilter_sample_VAF > 0],
                 'TP1_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$reasonForFail == 'likely_somatic' & !is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline != T & cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$wFilter_sample_VAF > 0])
  table(varList[[1]] %in% varList[[3]])
  
  library(UpSetR)
  upset(fromList(varList),nsets = 20,text.scale = 2)

  ##-------------------------------------------------------------------##
  ##      Plot - distribution of VAF / across chromosomes            ####
  ##-------------------------------------------------------------------##
  varToKeep = cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline == F & !cgpVaf_QC$snv_ID %in% unique(cgpVaf_QC$snv_ID[cgpVaf_QC$isGermline==T])]
  cgpVaf_QC$snv_type = ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$isGermline == F & cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$isGermline == T]],'likely_germline_D',
                              ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$isGermline == F & cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$isGermline == T]],'likely_germline_TP1',
                                     ifelse(cgpVaf_QC$CN_tot2major !='2:1' & cgpVaf_QC$reasonForFail != 'germline','CN_region',cgpVaf_QC$reasonForFail)))
  
  
  cgpVaf_QC$snv_type[cgpVaf_QC$snv_type %in% c('low_depth','low_highQualRatio','overlapped_wINDELs')] = 'removed'
  table(cgpVaf_QC$snv_type,cgpVaf_QC$reasonForFail)
  table(cgpVaf_QC$snv_type,cgpVaf_QC$CN_tot2major,cgpVaf_QC$sample)
  
  p = ggplot(cgpVaf_QC[cgpVaf_QC$snv_type != 'removed',],aes(wFilter_sample_VAF))+
    geom_density()+
    facet_grid(snv_type ~ timePoint)+
    geom_vline(xintercept = 0.5,lty=2,lwd=0.7,col='grey')+
    geom_vline(xintercept = c(0.13),lty=2,lwd=0.7,col=c('red'))+
    #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
    theme_classic(base_size = 12)+
    theme(axis.line = element_blank(),panel.border = element_rect(fill=F,colour = 'black',size=1),
          strip.background = element_blank(),axis.text = element_text(color='black'))
  print(p)

  table(cgpVaf_QC$snv_type,cgpVaf_QC$sample)


  ## Before shearwater-filter
  cgpVaf_QC$snv_type = factor(cgpVaf_QC$snv_type,c('germline',"likely_germline_D", "likely_germline_TP1",'CN_region','likely_somatic','removed'))
  cgpVaf_QC = cgpVaf_QC[order(cgpVaf_QC$snv_type),]
  cgpVaf_QC$Chrom = factor(cgpVaf_QC$Chrom,paste0('chr',c(1:22,'X','Y')))
  
  p = ggplot(cgpVaf_QC,aes(Pos,wFilter_sample_VAF))+
    geom_point(size=0.1,aes(col=snv_type,size=wFilter_sample_DEP))+
    facet_grid(sample ~ Chrom,scales = 'free_x')+
    scale_size_manual(values = c(0.01,3))+
    scale_color_manual(values = c(grey(0.8),col25))+
    geom_hline(yintercept = 0.5,lty=2,lwd=0.4,col=grey(0.4))+
    geom_hline(yintercept = 0.13,lty=2,lwd=0.4,col='black')+
    #geom_vline(xintercept = c(0.13),lty=2,lwd=0.7,col=c('red'))+
    #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
    theme_classic(base_size = 15)+
    theme(panel.border = element_rect(fill=F),
          axis.line = element_blank(),axis.text.x = element_blank(),axis.ticks.length.x = unit(0,'cm')) +
    xlab('Genomic position') + ylab('VAF')
  
  
  
  
  
  
  
  ##--------------------------------------------------------------------------------------------------------##
  ##      Tim's shearwater-like filter: cgpVAF of the somatic variants in normal blood samples            ####
  ##--------------------------------------------------------------------------------------------------------##
  
  ##--- Import list of somatic variants for shearwater
  # outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038'
  # cgpVaf_QC_og = read.csv(file.path(outDir,'L038_cgpVaf_allFiltersApplied.csv'))
  cgpVaf_QC = cgpVaf_QC_og
  ##------ TO BE SORTED OUT
  cgpVaf_QC$snv_type = ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$isGermline == F & cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$isGermline == T]],'likely_germline_D',
                              ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$isGermline == F & cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$isGermline == T]],'likely_germline_TP1',
                                     ifelse(cgpVaf_QC$CN_tot2major !='2:1' & cgpVaf_QC$reasonForFail != 'germline','CN_region',cgpVaf_QC$reasonForFail)))
  
  
  cgpVaf_QC$snv_type[cgpVaf_QC$snv_type %in% c('low_depth','low_highQualRatio','overlapped_wINDELs')] = 'removed'
  
  
  somaticVar = cgpVaf_QC_og[cgpVaf_QC_og$germlineFilter == 'likely_somatic',]
  
  
  outDir_sw = file.path(outDir,'shearwater')
  
  ##--- Upset plot
  varList = list('D_somatic' = somaticVar$snv_ID[somaticVar$PDID == 'PD60301a' & somaticVar$wFilter_sample_VAF > 0],
                 'TP1_somatic' = somaticVar$snv_ID[somaticVar$PDID == 'PD61846a' & somaticVar$wFilter_sample_VAF > 0])
  library(UpSetR)
  upset(fromList(varList),nsets = 20,text.scale = 2)
  
  ##--------------------------------------------------------------------------##
  ##      cgpVAF of the somatic variants in normal blood samples            ####
  ##--------------------------------------------------------------------------##
  ## STEP1: cgpVAF of variants of interst in panel of normal samples ##
  print('cgpVaf on normal panel')
  cgpVaf_onNormPanel(outDir_sw,somaticVar=somaticVar,
                     normPanel_fp = "/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt")
  
  
  ## STEP2: alleleCount at each variant position across panel of normal samples ##
  allelecount_dir = file.path(outDir_sw,'2_allelecount')
  columns = c('sample','Chrom','Pos')
  
  
  ##---- Import cgpVAF output with read filters applied
  
  cgpVaf_withReadFilters_output = import_cgpVAF_output(output_wFilter_file = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/cgpVaf_output/PDv38is_wgs_PD60301a_snp_vaf.tsv',
                                                       caveman_samples=caveman_samples,sampleID_name='timePoint',columns='all')
  
  bed = cgpVaf_QC_og[cgpVaf_withReadFilters_output$snv_ID %in% somaticVar$snv_ID,!grepl('PDv38is',colnames(cgpVaf_withReadFilters_output))]
  bed$Count_A = bed[,grepl('FAZ',colnames(bed))] + bed[,grepl('RAZ',colnames(bed))]
  bed$Count_C = bed[,grepl('FCZ',colnames(bed))] + bed[,grepl('RCZ',colnames(bed))]
  bed$Count_T = bed[,grepl('FTZ',colnames(bed))] + bed[,grepl('RTZ',colnames(bed))]
  bed$Count_G = bed[,grepl('FGZ',colnames(bed))] + bed[,grepl('RGZ',colnames(bed))]
  
  bed = split(bed,bed$PDID)
  bed = lapply(bed,function(i){
    sample = unique(i$PDID)
    print(sample)
    rownames(i) = paste0(gsub('chr','',i$Chrom),'_',i$Pos)
    write.table(i[,c('Chrom','Pos','Count_A','Count_C','Count_G','Count_T')], file = file.path(allelecount_dir,paste0(sample, "_allelecounts.txt")),
                sep = "\t", row.names = F, quote = F, col.names = T)
    
    return(i)
  })
  
  
  allele_count = bed[[1]][,c('Chrom','Pos','Count_A','Count_C','Count_G','Count_T')]
  
  # Run the same thing again but this time for the normal panel.
  alleleCount_onNormPanel(outDir_sw,allele_count)
  
  
  ## STEP3: Run Tim's shearwater-like filter to remove false mutations ##
  
  ## Input for shearwater filter
  # Vector of all sample names
  samples_ID=c('PD60301a','PD61846a')
  
  mutations = runShearWaterLike(samples_ID=c('PD60301a','PD61846a'),outDir=outDir_sw,patient='L038')
  
  
  # Tim uses <0.001 as the mutations that passed the shearwater filter and Anna uses <0.005
  table(mutations$PD60301a < 0.001)
  table(mutations$PD60301a < 0.005)
  mutations$PD60301a_shearwaterPASS = ifelse(mutations$PD60301a < 0.005,T,F)
  mutations$PD61846a_shearwaterPASS = ifelse(mutations$PD61846a < 0.005,T,F)
  mutations$snv_ID = paste0(mutations$Chr,':',mutations$Pos,'_',mutations$Ref,'/',mutations$Alt)
  table(somaticVar$snv_ID %in% mutations$snv_ID)
  
  
  ##  STEP4: Add shearwater-like filter results to main object ##
  
  cgpVaf_QC$snv_type = as.character(cgpVaf_QC$snv_type)
  cgpVaf_QC$snv_type[cgpVaf_QC$PDID == 'PD60301a' &
                       cgpVaf_QC$snv_type != 'removed' &
                       cgpVaf_QC$snv_ID %in% mutations$snv_ID[mutations$PD60301a_shearwaterPASS == F]] = 'shearwater_Failed'
  cgpVaf_QC$snv_type[cgpVaf_QC$PDID == 'PD61846a' &
                       cgpVaf_QC$snv_type != 'removed' &
                       cgpVaf_QC$snv_ID %in% mutations$snv_ID[mutations$PD61846a_shearwaterPASS == F]] = 'shearwater_Failed'
  
  
  cgpVaf_QC$snv_type = factor(cgpVaf_QC$snv_type,c('germline','likely_somatic',"likely_germline_D", "likely_germline_TP1",'CN_region','removed','shearwater_Failed'))
  cgpVaf_QC = cgpVaf_QC[order(cgpVaf_QC$snv_type),]
  cgpVaf_QC$Chrom = factor(cgpVaf_QC$Chrom,paste0('chr',c(1:22,'X','Y')))
  ggplot(cgpVaf_QC,aes(Pos,wFilter_sample_VAF))+
    geom_point(size=0.1,aes(col=snv_type,size=wFilter_sample_DEP))+
    facet_grid(sample ~ Chrom,scales = 'free_x')+
    scale_size_manual(values = c(0.01,3))+
    scale_color_manual(values = c(grey(0.8),'black',col25[1:3],grey(0.9),'purple'))+
    geom_hline(yintercept = 0.5,lty=2,lwd=0.4,col=grey(0.4))+
    geom_hline(yintercept = 0.13,lty=2,lwd=0.4,col='black')+
    #geom_vline(xintercept = c(0.13),lty=2,lwd=0.7,col=c('red'))+
    #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
    theme_classic(base_size = 12)+
    theme(panel.border = element_rect(fill=F,color='black'),
          axis.line = element_blank(),axis.text.x = element_blank(),
          axis.text = element_text(color='black'),
          axis.ticks.length.x = unit(0,'cm')) +
    xlab('Genomic position') + ylab('VAF')
  
  
  ##--- Upset plot
  table(cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD60301a'] %in% cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD61846a'])
  
  varList = list('D_germline' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type %in% c('germline','likely_germline_D') & cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$wFilter_sample_VAF >0],
                 'TP1_germline' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type %in% c('germline','likely_germline_TP1') & cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$wFilter_sample_VAF >0],
                 'D_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'likely_somatic' & cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$wFilter_sample_VAF > 0],
                 'TP1_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'likely_somatic' & cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$wFilter_sample_VAF > 0])
  table(varList[[1]] %in% varList[[3]])
  
  library(UpSetR)
  upset(fromList(varList),nsets = 20,text.scale = 2)
  
  
  
  
  ##------ Extract somatic variants after shearwater filter
  somaticVar = cgpVaf_QC[cgpVaf_QC$snv_type == 'likely_somatic',]
  somaticVar$somaticVar_type = ifelse(somaticVar$snv_ID %in% intersect(somaticVar$snv_ID[somaticVar$PDID != 'PD60301a' & somaticVar$wFilter_sample_VAF > 0],
                                                                       somaticVar$snv_ID[somaticVar$PDID == 'PD60301a' & somaticVar$wFilter_sample_VAF > 0]),'shared',
                                      ifelse(somaticVar$PDID == 'PD60301a' & somaticVar$wFilter_sample_VAF > 0,'unique_D',
                                             ifelse(somaticVar$PDID == 'PD61846a' & somaticVar$wFilter_sample_VAF > 0,'unique_TP1','others')))
  
  table(somaticVar$somaticVar_type,somaticVar$sample)
  
  # What are the reason for the unique variants to be unique
  cgpVaf_QC_tp1 = cgpVaf_QC[cgpVaf_QC$PDID != 'PD60301a',]
  cgpVaf_QC_tp1$snv_type = as.character(cgpVaf_QC_tp1$snv_type)
  
  cgpVaf_QC_d = cgpVaf_QC[cgpVaf_QC$PDID == 'PD60301a',]
  somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_D'] = paste0(somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_D'],':',
                                                                                cgpVaf_QC_tp1$snv_type[match(somaticVar$snv_ID[somaticVar$somaticVar_type == 'unique_D'],cgpVaf_QC_tp1$snv_ID)])
  somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_TP1'] = paste0(somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_TP1'],':',
                                                                                  cgpVaf_QC_d$snv_type[match(somaticVar$snv_ID[somaticVar$somaticVar_type == 'unique_TP1'],cgpVaf_QC_d$snv_ID)])
  
  
  table(somaticVar$somaticVar_type,somaticVar$sample)
  
  
  p = ggplot(somaticVar,aes(wFilter_sample_VAF))+
    geom_density()+
    facet_grid(somaticVar_type ~ timePoint)+
    geom_vline(xintercept = 0.5,lty=2,lwd=0.7,col='grey')+
    #geom_vline(xintercept = c(expected_vaf_d),lty=2,lwd=0.7,col=c('red'))+
    #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
    theme_classic(base_size = 14)+theme(axis.line = element_line())
  print(p)
  
  
  ## pivot somaticVar wider to plot vaf of D vs TP1
  df = cgpVaf_QC_og[cgpVaf_QC_og$snv_ID %in% somaticVar$snv_ID,]
  df$somaticVar_type = somaticVar[somaticVar$somaticVar_type != 'others',]$somaticVar_type[match(df$snv_ID,somaticVar$snv_ID[somaticVar$somaticVar_type != 'others'])]
  
  df = pivot_wider(df,id_cols = c('snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect','Impact','AAchange','somaticVar_type'),
                   names_from = 'PDID',values_from = 'wFilter_sample_VAF')
  
  cgpVaf_QC_og$reasonForFail[cgpVaf_QC_og$reasonForFail == 'PASS' & cgpVaf_QC_og$pindelFilter == 'overlapped_wINDELs'] = 'overlapped_wINDELs'
  cgpVaf_QC_og$reasonForFail[cgpVaf_QC_og$reasonForFail == 'PASS' & cgpVaf_QC_og$germlineFilter == 'germline'] = 'germline'
  df$snv_removed_PD60301a = cgpVaf_QC_og[cgpVaf_QC_og$PDID == 'PD60301a',]$reasonForFail[match(df$snv_ID,cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$PDID == 'PD60301a'])]
  df$snv_removed_PD61846a = cgpVaf_QC_og[cgpVaf_QC_og$PDID == 'PD61846a',]$reasonForFail[match(df$snv_ID,cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$PDID == 'PD61846a'])]
  
  
  p = ggplot(df,aes(PD60301a,PD61846a,col=somaticVar_type))+
    geom_point(size=0.7)+
    scale_color_manual(values = c(grey(0.8),grey(0.2),col25))+
    theme_classic(base_size = 14) + ggtitle('cgpVAF output of somatic variants')
  print(p)
  
  write.csv(df, '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/L038_D.vs.TP1_finalSomaticVariants.csv')
  
  
  
  
  ##--------------------------------------------------##
  ##  Generate jbrowse images of these mutations    ####
  ##--------------------------------------------------##
  # df_old = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_D.vs.TP1_finalSomaticVariants.csv')
  # colnames(df_old) = paste0(colnames(df_old),'_old')
  # df_old = cbind(df_old,df[match(df_old$snv_ID_old,df$snv_ID),])
  # table(df_old$somaticVar_type_old,df_old$somaticVar_type)
  df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/L038_D.vs.TP1_finalSomaticVariants.csv')
  ## 1. Define JBROWSE URL
  # L038 samples
  jbrowse_url_l038 = '# JBROWSE https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F3030&loc=1%3A1..247467938&tracks=PD60301a_bwa%2CPD61846a_bwa'
  # normal panel samples from 3010 project
  normal_samples=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt", header=T)
  normalBlood_sampleIDs = normal_samples$Sample_ID
  jbrowse_url_normalblood = paste0('# JBROWSE https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F3010&loc=1%3A1..247467938&tracks=',
                                   paste0(normalBlood_sampleIDs,collapse = '_bwa%2C'))
  
  
  ## 2. Prepare bed file containing location of these mutations, with +/ 50 bases around the mutation position
  jbrowse_outDir = file.path(outDir,'jbrowse')
  jbrowse_script = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/jbrowse_rasterise.sh'
  
  df$start = df$Pos - 50
  df$end = df$Pos + 50
  df$Chrom = gsub('^chr','',df$Chrom)
  cmd_list = c()
  for(mutCat in unique(df$somaticVar_type)){
    var = df[df$somaticVar_type == mutCat,]
    d = file.path(jbrowse_outDir,mutCat)
    if(!dir.exists(d)){
      dir.create(d,recursive = T)
      for(d.sub in file.path(d,c(donorID,'normalPanel_3010'))){
        dir.create(d.sub,recursive = T)
      }
    }
    
    for(sampleGroup in c(donorID,'normalPanel_3010')){
      ## Write bed file containing mutation of this category
      bed = var[,c('Chrom','start','end')]
      bedFile = file.path(d,sampleGroup,'L038_mutation.bed')
      if(sampleGroup == donorID){
        cat(jbrowse_url_l038,file = bedFile,sep = '\n')
      }else{
        cat(jbrowse_url_normalblood,file = bedFile,sep = '\n')
      }
      write.table(bed,file = bedFile,sep = '\t',quote = F,row.names = F,col.names = F,append = T)
      
      
      ## call jbrowse bash script
      cmd_list = c(cmd_list,sprintf('bash %s %s %s',jbrowse_script,bedFile,file.path(d,sampleGroup)))
    }
    
  }
  
  write.table(cmd_list,file = file.path(jbrowse_outDir,'jbrowse_command.sh'),quote = F,col.names = F,row.names = F)
  
  
  ##------------------------------##
  ##    Plot binom mix model    ####
  ##------------------------------##
  # df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/L038_D.vs.TP1_finalSomaticVariants.csv')
  # df_long = pivot_longer(df,cols = c('PD60301a', 'PD61846a'),names_to = 'sampleID',values_to = 'VAF')
  # ggplot(df_long,aes(VAF))+
  #   geom_density()+
  #   facet_grid(.~sampleID)+
  #   theme_classic(base_size = 14)+
  #   theme(panel.border = element_rect(fill=F),axis.line = element_blank())+xlab('Variant allele frequency')
  # 
  # ## subset for mutations that has been called
  # cgpVaf_QC_og = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/L038_cgpVaf_allFiltersApplied.csv')
  # cgpVaf_QC = cgpVaf_QC[cgpVaf_QC$snv_ID %in% df$snv_ID,]
  # res_D = binom_mix(x=cgpVaf_QC$wFilter_sample_MTR[cgpVaf_QC$PDID == 'PD60301a'& cgpVaf_QC$noFilter_sample_DEP > 0],
  #                   size = cgpVaf_QC$noFilter_sample_DEP[cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$noFilter_sample_DEP > 0],mode = 'Full')
  # 
  # data_D = data.frame(NV = cgpVaf_QC$wFilter_sample_MTR[cgpVaf_QC$PDID == 'PD60301a'& cgpVaf_QC$noFilter_sample_DEP > 0],
  #                     NR = cgpVaf_QC$noFilter_sample_DEP[cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$noFilter_sample_DEP > 0],
  #                     snvID = cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD60301a' & cgpVaf_QC$noFilter_sample_DEP > 0],
  #                     cluster = res_D[['Which_cluster']]
  # )
  # data_D$cat = df$somaticVar_type[match(data_D$snvID,df$snv_ID)]
  # table(data_D$cat,data_D$cluster)
  # plot_binomMixModel(data_D,res_D,sampleID=NULL,mode='Full')
  # 
  # res_TP1 = binom_mix(x=cgpVaf_QC$wFilter_sample_MTR[cgpVaf_QC$PDID == 'PD61846a'& cgpVaf_QC$noFilter_sample_DEP > 0],
  #                     size = cgpVaf_QC$noFilter_sample_DEP[cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$noFilter_sample_DEP > 0],mode = 'Full',nrange = 1:2)
  # data_TP1 = data.frame(NV = cgpVaf_QC$wFilter_sample_MTR[cgpVaf_QC$PDID == 'PD61846a'& cgpVaf_QC$noFilter_sample_DEP > 0],
  #                       NR = cgpVaf_QC$noFilter_sample_DEP[cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$noFilter_sample_DEP > 0],
  #                       snvID = cgpVaf_QC$snv_ID[cgpVaf_QC$PDID == 'PD61846a' & cgpVaf_QC$noFilter_sample_DEP > 0],
  #                       cluster = res_TP1[['Which_cluster']])
  # data_TP1$cat = df$somaticVar_type[match(data_TP1$snvID,df$snv_ID)]
  # table(data_TP1$cat,data_TP1$cluster)
  # plot_binomMixModel(data_TP1,res_TP1,sampleID=NULL,mode='Full')
  
}


