## L076 - process caveman output between Diagnostic and TP1 samples ##

outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

library(VariantAnnotation)
library(tidyverse)
source('~/lustre_mt22/generalScripts/caveman_readvcf.R')
source('~/lustre_mt22/generalScripts/binom_mix_model_Tim.R')
source('~/lustre_mt22/generalScripts/utils/wgs_analysis_helperFunctions.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')

##-------------------------------------------##
##    Define CN regions - TP1 sample       ####
##-------------------------------------------##
cnSegs = GRanges(c('chr8','chr13','chr21'),
                 IRanges(c(rep(1,3)),c(2e9,27.5e6,2e9)))

# specify ChrX and chrY pairing regions
chrXY_par = GRanges(c('chrX','chrX','chrY','chrY'),
                    IRanges(c(10001,155701383,10001,56887903),c(2781479,156030895,2781479,57217415)),
                    region = c('PAR1','PAR2','PAR1','PAR1'))
  
## hg38 PAR regions as specified on ensembl website  
# chromosome:GRCh38:Y:1 - 10000 is unique to Y but is a string of 10000 Ns
# chromosome:GRCh38:Y:10001 - 2781479 is shared with X: 10001 - 2781479 (PAR1)
# chromosome:GRCh38:Y:2781480 - 56887902 is unique to Y
# chromosome:GRCh38:Y:56887903 - 57217415 is shared with X: 155701383 - 156030895 (PAR2)
# chromosome:GRCh38:Y:57217416 - 57227415 is unique to Y



#------------------------------------------------------##
#    Import Caveman output - Blood Diagnostic        ####
#------------------------------------------------------##
b_caveman = caveman2(sample='PD62331c',projectid = 3030)
b_caveman$tissue = 'Blood_D'
b_caveman = b_caveman[b_caveman$Flag == T,]
b_caveman$snv_ID = paste0(b_caveman$Chr,':',b_caveman$Pos,'_',b_caveman$Ref,'/',b_caveman$Alt)

##---------------------------------------------------##
##    Import Caveman output - BM Diagnostic        ####
##---------------------------------------------------##
bm_caveman = caveman2(sample='PD62331a',projectid = 3030)
bm_caveman$tissue = 'BM_D'
bm_caveman = bm_caveman[bm_caveman$Flag == T,]
bm_caveman$snv_ID = paste0(bm_caveman$Chr,':',bm_caveman$Pos,'_',bm_caveman$Ref,'/',bm_caveman$Alt)


var4cgpVaf = rbind(b_caveman,bm_caveman)
# 
# 
# ##--------------------------------------------------##
# ##      cgpVAF of the somatic variants            ####
# ##--------------------------------------------------##
# # write bed file of the position of interest
# # requried format:
# #   chr pos ref_allele alt_allele
# #   filename ends with .bed, delim='\t'
# 
# bed = data.frame(snv_ID = unique(var4cgpVaf$snv_ID))
# bed = merge(bed,var4cgpVaf[,c('Chr','Pos','Ref','Alt','snv_ID')],by='snv_ID',all.x=T)
# bed = bed[!duplicated(bed),]
# dim(bed)
# n_distinct(var4cgpVaf$snv_ID)
# #rownames(bed) = bed$snv_ID
# bed = bed[order(bed$Chr),]
# write_delim(bed[,colnames(bed) != 'snv_ID'],'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L076_unmatched_cavemanVar.bed',col_names = F,delim = '\t')
# 
# 
# ##--- prepare for cgpVAF
# # cd ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis
# # mkdir cgpVaf_input/L076
# # mkdir cgpVaf_output/L076
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD62331c/PD62331c.sample.dupmarked.bam cgpVaf_input/L076/
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD62331c/PD62331c.sample.dupmarked.bam.bai cgpVaf_input/L076
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD62331c/PD62331c.caveman_c.annot.vcf.gz cgpVaf_input/L076
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD62331a/PD62331a.sample.dupmarked.bam cgpVaf_input/L076
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD62331a/PD62331a.sample.dupmarked.bam.bai cgpVaf_input/L076
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD62331a/PD62331a.caveman_c.annot.vcf.gz cgpVaf_input/L076
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam cgpVaf_input/L076/
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam.bai cgpVaf_input/L076/
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam cgpVaf_input/L076/
# 
# ##--- Run cgpVAF
# # bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/x10_MLDS_cgpVaf.sh L076 PD62331c
# # bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/x10_MLDS_cgpVaf.sh L076 PD62331a
# 
# 
# 
# 
# 
##---- Import cgpVAF output with read filters applied
# Diagnostic
b_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076/PDv38is_wgs_PD62331c_snp_vaf.tsv',sep = '\t',skip = 36)
b_cgpVaf$tissue = 'Blood_D'
b_cgpVaf$snv_ID = paste0(b_cgpVaf$Chrom,':',b_cgpVaf$Pos,'_',b_cgpVaf$Ref,'/',b_cgpVaf$Alt)
b_cgpVaf$sample = 'PD62331c'
b_cgpVaf$sample_MTR = b_cgpVaf[,grepl("MTR",colnames(b_cgpVaf))&!colnames(b_cgpVaf)%in%"PDv38is_wgs_MTR"]
b_cgpVaf$sample_VAF = b_cgpVaf[,grepl("VAF",colnames(b_cgpVaf))&!colnames(b_cgpVaf)%in%"PDv38is_wgs_VAF"]
b_cgpVaf$sample_DEP = b_cgpVaf[,grepl("DEP",colnames(b_cgpVaf))&!colnames(b_cgpVaf)%in%"PDv38is_wgs_DEP"]

# TP1
bm_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076/PDv38is_wgs_PD62331a_snp_vaf.tsv',sep = '\t',skip = 36)
bm_cgpVaf$tissue = 'BM_D'
bm_cgpVaf$snv_ID = paste0(bm_cgpVaf$Chrom,':',bm_cgpVaf$Pos,'_',bm_cgpVaf$Ref,'/',bm_cgpVaf$Alt)
bm_cgpVaf$sample = 'PD62331a'
bm_cgpVaf$sample_MTR = bm_cgpVaf[,grepl("MTR",colnames(bm_cgpVaf))&!colnames(bm_cgpVaf)%in%"PDv38is_wgs_MTR"]
bm_cgpVaf$sample_VAF = bm_cgpVaf[,grepl("VAF",colnames(bm_cgpVaf))&!colnames(bm_cgpVaf)%in%"PDv38is_wgs_VAF"]
bm_cgpVaf$sample_DEP = bm_cgpVaf[,grepl("DEP",colnames(bm_cgpVaf))&!colnames(bm_cgpVaf)%in%"PDv38is_wgs_DEP"]



columns = c('sample','tissue','snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect','sample_MTR','sample_VAF','sample_DEP')
cgpVaf_withReadFilters_output = rbind(b_cgpVaf[,columns],bm_cgpVaf[,columns])
dim(cgpVaf_withReadFilters_output)



##---- Import cgpVAF output with NO read filters applied
# Diagnostic
b_cgpVaf_noFilter = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_nofilter/PDv38is_wgs_PD62331c_snp_vaf.tsv',sep = '\t',skip = 36)
b_cgpVaf_noFilter$tissue = 'Blood_D'
b_cgpVaf_noFilter$snv_ID = paste0(b_cgpVaf_noFilter$Chrom,':',b_cgpVaf_noFilter$Pos,'_',b_cgpVaf_noFilter$Ref,'/',b_cgpVaf_noFilter$Alt)
b_cgpVaf_noFilter$sample = 'PD62331c'
b_cgpVaf_noFilter$sample_MTR = b_cgpVaf_noFilter[,grepl("MTR",colnames(b_cgpVaf_noFilter))&!colnames(b_cgpVaf_noFilter)%in%"PDv38is_wgs_MTR"]
b_cgpVaf_noFilter$sample_VAF = b_cgpVaf_noFilter[,grepl("VAF",colnames(b_cgpVaf_noFilter))&!colnames(b_cgpVaf_noFilter)%in%"PDv38is_wgs_VAF"]
b_cgpVaf_noFilter$sample_DEP = b_cgpVaf_noFilter[,grepl("DEP",colnames(b_cgpVaf_noFilter))&!colnames(b_cgpVaf_noFilter)%in%"PDv38is_wgs_DEP"]

# TP1
bm_cgpVaf_noFilter = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_nofilter/PDv38is_wgs_PD62331a_snp_vaf.tsv',sep = '\t',skip = 36)
bm_cgpVaf_noFilter$tissue = 'BM_D'
bm_cgpVaf_noFilter$snv_ID = paste0(bm_cgpVaf_noFilter$Chrom,':',bm_cgpVaf_noFilter$Pos,'_',bm_cgpVaf_noFilter$Ref,'/',bm_cgpVaf_noFilter$Alt)
bm_cgpVaf_noFilter$sample = 'PD62331a'
bm_cgpVaf_noFilter$sample_MTR = bm_cgpVaf_noFilter[,grepl("MTR",colnames(bm_cgpVaf_noFilter))&!colnames(bm_cgpVaf_noFilter)%in%"PDv38is_wgs_MTR"]
bm_cgpVaf_noFilter$sample_VAF = bm_cgpVaf_noFilter[,grepl("VAF",colnames(bm_cgpVaf_noFilter))&!colnames(bm_cgpVaf_noFilter)%in%"PDv38is_wgs_VAF"]
bm_cgpVaf_noFilter$sample_DEP = bm_cgpVaf_noFilter[,grepl("DEP",colnames(bm_cgpVaf_noFilter))&!colnames(bm_cgpVaf_noFilter)%in%"PDv38is_wgs_DEP"]



columns = c('sample','tissue','snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect','sample_MTR','sample_VAF','sample_DEP')
cgpVaf_noReadFilters_output = rbind(b_cgpVaf_noFilter[,columns],bm_cgpVaf_noFilter[,columns])
dim(cgpVaf_noReadFilters_output)
colnames(cgpVaf_noReadFilters_output)[colnames(cgpVaf_noReadFilters_output) %in% c('sample_MTR','sample_VAF','sample_DEP')] = paste0(colnames(cgpVaf_noReadFilters_output)[colnames(cgpVaf_noReadFilters_output) %in% c('sample_MTR','sample_VAF','sample_DEP')],':noFilter')

cgpVaf_QC_og = merge(cgpVaf_withReadFilters_output,cgpVaf_noReadFilters_output,by=c('sample','tissue','snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect'),all=T)
#cgpVaf_QC = merge(cgpVaf_withReadFilters_output,cgpVaf_noReadFilters_output,by=c('sample','tissue','snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect'),all=T)



##----------------------------------------------------------------##
##      cgpVAF depth filter - consider sex chromosome           ####
##----------------------------------------------------------------##

gender = 'male'
genotype = 'T21'

# Plot distribution of depth by chromosome
cgpVaf_QC_og$Chrom = factor(cgpVaf_QC_og$Chrom,levels = paste0('chr',c(1:22,'X','Y')))
ggplot(cgpVaf_QC_og,aes(Chrom,sample_DEP))+
  geom_boxplot(outlier.colour = 'white')+
  facet_wrap(vars(sample),ncol=1)+
  theme_classic() + ggtitle('cgpVAF DEPTH distribution from both samples',
                            subtitle = paste0(n_distinct(cgpVaf_QC_og$snv_ID),' variants')) + xlab('')






# Calculate the average depth of each variant across samples of the same individual
# exclude the copy number regions in TP1 sample
cgpVaf_QC_og_gr = df2granges(data=cgpVaf_QC_og, chr='Chrom', start='Pos', end='Pos',id='snv_ID',sampleID='sample')
cgpVaf_QC_og_TP1_CNsegs = subsetByOverlaps(cgpVaf_QC_og_gr,cnSegs)
toExclude = c(which(cgpVaf_QC_og$Chrom == 'chr21'),
              which(cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID))
data_forDepthFilter = cgpVaf_QC_og[-toExclude,]

data_forDepthFilter = data_forDepthFilter %>% group_by(snv_ID,Chrom) %>% summarise(meanDEP = mean(sample_DEP))
colnames(data_forDepthFilter)[3] = 'sample_DEP'

cgpVaf_var_avgDepth = depthFilter(data = data_forDepthFilter,gender,genotype,useChr=NULL,upperCutOff=F)

table(cgpVaf_QC_og$Chrom[cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID])
cgpVaf_QC_og$depthFilter = cgpVaf_var_avgDepth$depthFilter[match(cgpVaf_QC_og$snv_ID,cgpVaf_var_avgDepth$snv_ID)]

## Remove variants with no coverage in one of the samples
noCovVar = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample_DEP == 0]
cgpVaf_QC_og$depthFilter[cgpVaf_QC_og$snv_ID %in% noCovVar] = 'tooLow_DEP'
cgpVaf_QC_og$depthFilter[is.na(cgpVaf_QC_og$depthFilter)] = 'NA'
table(cgpVaf_QC_og$depthFilter,cgpVaf_QC_og$sample)


##---------------------------------------------------------------##
##      cgpVAF fraction of high-quality read filter            ####
##---------------------------------------------------------------##

cgpVaf_QC_og$highQualRatio = cgpVaf_QC_og$sample_DEP / cgpVaf_QC_og$`sample_DEP:noFilter`
cgpVaf_QC_og$highQualRatio_filter = ifelse(cgpVaf_QC_og$highQualRatio >= 0.75,'PASS','low')

## Remove variants with low high-quality read fraction in at least one of the samples
lowHiQualVar = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$highQualRatio_filter == 'low']
cgpVaf_QC_og$highQualRatio_filter[cgpVaf_QC_og$snv_ID %in% lowHiQualVar] = 'low'

table(cgpVaf_QC_og$highQualRatio_filter[cgpVaf_QC_og$depthFilter != 'tooLow_DEP'],cgpVaf_QC_og$sample[cgpVaf_QC_og$depthFilter != 'tooLow_DEP'])

cgpVaf_QC_og$reasonForFail = ifelse(cgpVaf_QC_og$depthFilter == 'tooLow_DEP' | cgpVaf_QC_og$sample_DEP == 0,'low_depth',
                                    ifelse(cgpVaf_QC_og$highQualRatio_filter == 'low','low_highQualRatio','PASS'))
table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$sample)


write.csv(cgpVaf_QC_og,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_cgpVafQCog_depth.highQualRatio_filtered.csv')


df = pivot_wider(cgpVaf_QC_og,id_cols = c('snv_ID','reasonForFail'),names_from = sample,values_from = sample_VAF)
ggplot(df,aes(PD62331a,PD62331c,col = reasonForFail))+
  geom_point(size=0.01)+
  theme_classic()
table(df$reasonForFail[df$PD62331a == 0 | df$PD62331c == 0])



##----------------------------------##
##      Pindel filter             ####
##----------------------------------##
# Including FAILED Pindel variants too

intv = 10

# pindel filter
library(data.table)
indels = pindel(sample = 'PD62331c',projectid=3030,filter_lowQual=F)
dim(indels)
n_distinct(indels$varID)
indels$End=indels$Pos+nchar(indels$Ref)-1
indels$sampleID = 'PD62331c'

tp1_indels = pindel(sample = 'PD62331a',projectid=3030,filter_lowQual=F)
dim(tp1_indels)
n_distinct(tp1_indels$varID)
tp1_indels$sampleID = 'PD62331a'


## combine indels from Diagnostic and TP1 samples
indels = rbind(indels,tp1_indels)
indels_gr = GRanges(indels$Chr,IRanges(indels$Pos - intv,indels$End + intv),
                    varID = indels$varID,sampleID=indels$sampleID)
mcols(indels_gr) = cbind(mcols(indels_gr),indels[match(paste0(indels_gr$varID,':',indels_gr$sampleID),
                                                       paste0(indels$varID,':',indels$sampleID)),])
length(indels_gr) == nrow(indels)

# Remove copy number regions in both sample
cnSegs_indels = subsetByOverlaps(indels_gr,cnSegs)

indels_gr = indels_gr[!indels_gr$varID %in% cnSegs_indels$varID,]
indels = indels_gr

# find mutations within 10bp of indels
intv = 10
cgpVaf_QC_og = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_cgpVafQCog_depth.highQualRatio_filtered.csv')
cgpVaf_QC = cgpVaf_QC_og[cgpVaf_QC_og$reasonForFail == 'PASS',]


## Convert cgpVaf_QC into Granges
cgpVaf_QC_gr = GRanges(cgpVaf_QC$Chrom,IRanges(cgpVaf_QC$Pos,cgpVaf_QC$Pos),
                       snv_ID = cgpVaf_QC$snv_ID,sample=cgpVaf_QC$sample,tissue=cgpVaf_QC$tissue)
mcols(cgpVaf_QC_gr) = cbind(mcols(cgpVaf_QC_gr),cgpVaf_QC[match(paste0(cgpVaf_QC_gr$snv_ID,':',cgpVaf_QC_gr$sample,':',cgpVaf_QC_gr$timepoint),
                                                                paste0(cgpVaf_QC_gr$snv_ID,':',cgpVaf_QC_gr$sample,':',cgpVaf_QC_gr$timepoint)),])
length(cgpVaf_QC_gr) == nrow(cgpVaf_QC)
cgpVaf_QC = cgpVaf_QC_gr



# Identify variants overlapping with indels
cgpVaf_QC_overlapped_wIndels = subsetByOverlaps(cgpVaf_QC,indels)

cgpVaf_QC$overlapped_with_indels = ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC_overlapped_wIndels$snv_ID,T,F)
table(cgpVaf_QC$overlapped_with_indels,cgpVaf_QC$reasonForFail,cgpVaf_QC$sample)

cgpVaf_QC = as.data.frame(mcols(cgpVaf_QC))

# Taryn's way of filtering
# cgpVaf_QC$overlapped_with_indels = F
# 
# library(doParallel)
# 
# parallel = T
# #Setup backend to use many processors
# totalCores = detectCores()
# #Leave some core to avoid overload your computer
# if(parallel){
#   cluster <- makeCluster(totalCores[1]-10)
# }else{
#   cluster <- makeCluster(1)
# }
# 
# registerDoParallel(cluster)
# 
# 
# #Load foreach library
# start_time <- Sys.time()
# cor_df <- foreach(i = 1:nrow(cgpVaf_QC), .combine=rbind) %dopar% {
# #for (i in 1:nrow(cgpVaf_QC)){
#   print(i)
#   chr = cgpVaf_QC$Chrom[i]
#   pos = cgpVaf_QC$Pos[i]
#   #sample = cgpVaf_QC$sample[i]
#   
#   overlapped_indels = indels[#indels$sampleID == sample & 
#                                indels$Chr == chr & indels$Pos-intv < pos & indels$End+intv > pos,]
#   if(nrow(overlapped_indels) > 0){
#     cgpVaf_QC$overlapped_with_indels[i] = T
#   }
#   
#   return(cgpVaf_QC[i,])
# }
# 
#   
# #Stop cluster
# stopCluster(cluster)
# end_time <- Sys.time()
# print(end_time - start_time)



write.csv(cgpVaf_QC,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_cgpVaf_depth.highQualRatio.Pindel_filtered.csv')
cgpVaf_QC = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_cgpVaf_depth.highQualRatio.Pindel_filtered.csv')

# cgpVaf_QC_og$pindelFilter = ifelse(cgpVaf_QC_og$sample == 'PD62331c' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$overlapped_with_indels == F],'PASS',
#                                    ifelse(cgpVaf_QC_og$sample == 'PD62331a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$overlapped_with_indels == F],'PASS',
#                                           ifelse(cgpVaf_QC_og$sample == 'PD62331c' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$overlapped_with_indels == T],'overlapped_wINDELs',
#                                                  ifelse(cgpVaf_QC_og$sample == 'PD62331a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$overlapped_with_indels == T],'overlapped_wINDELs','removed'))))

cgpVaf_QC_og$pindelFilter = ifelse(cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$overlapped_with_indels == F],'PASS',
                                   ifelse(cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$overlapped_with_indels == T],'overlapped_wINDELs', 
                                          ifelse(cgpVaf_QC_og$reasonForFail != 'PASS','removed','others')))
                                                 

table(cgpVaf_QC_og$pindelFilter,cgpVaf_QC_og$sample)

# Extract only PASS variants for germline filter
cgpVaf_QC=cgpVaf_QC[cgpVaf_QC$overlapped_with_indels == F,]
table(cgpVaf_QC$sample)





##------------------------------------------------------##
##      Binomial test for germline filter             ####
##------------------------------------------------------##

# Run binomial filter
# data is a data.frame where: rownames = snv_ID (chr1:pos_REF/ALT), required columns are refCnt, altCnt, DP (total depth)
cgpVaf_QC_og$refCnt = NA
cgpVaf_QC_og$refCnt[cgpVaf_QC_og$sample == 'PD62331c'] =   b_cgpVaf$PD62331c_WTR[match(cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD62331c'],b_cgpVaf$snv_ID)]
cgpVaf_QC_og$refCnt[cgpVaf_QC_og$sample == 'PD62331a'] =   bm_cgpVaf$PD62331a_WTR[match(cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD62331a'],bm_cgpVaf$snv_ID)]

#colnames(cgpVaf_withReadFilters_output) = gsub('\\.sample_MTR$','',colnames(cgpVaf_withReadFilters_output))
cgpVaf_QC_og$altCnt = cgpVaf_withReadFilters_output$sample_MTR[match(paste0(cgpVaf_QC_og$snv_ID,'_',cgpVaf_QC_og$sample),
                                                                     paste0(cgpVaf_withReadFilters_output$snv_ID,'_',cgpVaf_withReadFilters_output$sample))]

cgpVaf_QC_og$DP = cgpVaf_withReadFilters_output$sample_DEP[match(paste0(cgpVaf_QC_og$snv_ID,'_',cgpVaf_QC_og$sample),
                                                                 paste0(cgpVaf_withReadFilters_output$snv_ID,'_',cgpVaf_withReadFilters_output$sample))]

cgpVaf_QC_og$sample_MTR = cgpVaf_withReadFilters_output$sample_MTR[match(paste0(cgpVaf_QC_og$snv_ID,'_',cgpVaf_QC_og$sample),
                                                                         paste0(cgpVaf_withReadFilters_output$snv_ID,'_',cgpVaf_withReadFilters_output$sample))]


# Exclude SNPs from CN regions in both samples
cgpVaf_QC_og$CN_tot2major = '2:1'
cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$Chrom == 'chr21'] = '3:2'
cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID[cgpVaf_QC_og_TP1_CNsegs$Chrom == 'chr8']] = '3:1'
cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID[cgpVaf_QC_og_TP1_CNsegs$Chrom == 'chr13']] = '1:1'


## Perform binomial test filter

cgpVaf_QC_og.sub = snvBinomTest(gender='male',data=cgpVaf_QC_og[cgpVaf_QC_og$pindelFilter == 'PASS',],
                                genotype='T21',pCut=0.05)
table(cgpVaf_QC_og.sub$isGermline,cgpVaf_QC_og.sub$sample)
table(cgpVaf_QC_og.sub$isGermline,cgpVaf_QC_og.sub$Chrom)
table(!is.na(cgpVaf_QC_og.sub$isGermline),cgpVaf_QC_og.sub$Chrom)

cgpVaf_QC_og$isGermline = cgpVaf_QC_og.sub$isGermline[match(paste0(cgpVaf_QC_og$snv_ID,':',cgpVaf_QC_og$sample),
                                                            paste0(cgpVaf_QC_og.sub$snv_ID,':',cgpVaf_QC_og.sub$sample))]
cgpVaf_QC_og$germlineFilter = ifelse(cgpVaf_QC_og$pindelFilter != 'PASS','removed',
                                     ifelse(cgpVaf_QC_og$sample == 'PD62331c' &
                                              cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$sample == 'PD62331c' & cgpVaf_QC_og.sub$isGermline == T],'germline',
                                       ifelse(cgpVaf_QC_og$sample == 'PD62331a' &
                                                cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$sample == 'PD62331a' & cgpVaf_QC_og.sub$isGermline == T],'germline',
                                              ifelse(cgpVaf_QC_og$sample == 'PD62331c' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$sample == 'PD62331c' & cgpVaf_QC_og.sub$isGermline == F],'likely_somatic',
                                                     ifelse(cgpVaf_QC_og$sample == 'PD62331a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$sample == 'PD62331a' & cgpVaf_QC_og.sub$isGermline == F],'likely_somatic','others')))))


table(cgpVaf_QC_og$germlineFilter,cgpVaf_QC_og$pindelFilter)


## Summarise reason for fails
cgpVaf_QC_og$reasonForFail = ifelse(cgpVaf_QC_og$CN_tot2major != '2:1' & cgpVaf_QC_og$germlineFilter == 'others','CN_region',
                                    ifelse(cgpVaf_QC_og$germlineFilter == 'germline','germline',
                                           ifelse(cgpVaf_QC_og$germlineFilter == 'likely_somatic','likely_somatic',
                                              ifelse((!is.na(cgpVaf_QC_og$depthFilter) & cgpVaf_QC_og$depthFilter == 'tooLow_DEP') | cgpVaf_QC_og$sample_DEP == 0,'low_depth',
                                                     ifelse(cgpVaf_QC_og$highQualRatio_filter == 'low','low_highQualRatio',
                                                            ifelse(cgpVaf_QC_og$pindelFilter == 'overlapped_wINDELs','overlapped_wINDELs','others'))))))
                                                           

table(is.na(cgpVaf_QC_og$reasonForFail))
table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$sample)
table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$germlineFilter)

## Remove SNPs within PAR regions on sex chromosomes
cgpVaf_QC_gr = GRanges(cgpVaf_QC$Chrom,IRanges(cgpVaf_QC$Pos,cgpVaf_QC$Pos),
                       snv_ID = cgpVaf_QC$snv_ID,sample=cgpVaf_QC$sample,tissue=cgpVaf_QC$tissue)
cgpVaf_QC_overlapped_wPAR = subsetByOverlaps(cgpVaf_QC_gr,chrXY_par)
cgpVaf_QC$reasonForFail[cgpVaf_QC$snv_ID %in% cgpVaf_QC_overlapped_wPAR$snv_ID & cgpVaf_QC$reasonForFail %in% c('likely_somatic')] = 'PAR'

cgpVaf_QC_og$reasonForFail[cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_overlapped_wPAR$snv_ID & cgpVaf_QC_og$reasonForFail %in% c('likely_somatic')] = 'PAR'

ggplot(cgpVaf_QC_og,aes(sample,fill = reasonForFail))+
  geom_bar(position='fill')+
  scale_fill_manual(values = c(col25[1],grey(0.8),col25[-1]))+
  theme_classic(base_size = 13)

write.csv(cgpVaf_QC_og,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_cgpVaf_filtersApplied.csv')
cgpVaf_QC_og = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_cgpVaf_filtersApplied.csv')

#cgpVaf_QC = cgpVaf_QC_og[cgpVaf_QC_og$germlineFilter == 'likely_somatic',]

cgpVaf_QC = cgpVaf_QC_og[!is.na(cgpVaf_QC_og$isGermline),]
cgpVaf_QC = cgpVaf_QC_og

##----------------------------------------------------##
##      Plot - number of shared variants            ####
##----------------------------------------------------##
table(cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331c'] %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331a'])

varList = list('B_germline' = cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline == T & cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$sample_VAF > 0],
               'BM_germline' = cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) &cgpVaf_QC$isGermline == T & cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$sample_VAF > 0],
               'B_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$reasonForFail == 'likely_somatic' & !is.na(cgpVaf_QC$isGermline) &cgpVaf_QC$isGermline != T & cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$sample_VAF > 0],
               'BM_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$reasonForFail == 'likely_somatic' & !is.na(cgpVaf_QC$isGermline) &cgpVaf_QC$isGermline != T & cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$sample_VAF > 0])
table(varList[[1]] %in% varList[[3]])
table(varList[[1]] %in% intersect(varList[[4]],intersect(varList[[2]],varList[[3]])))

library(UpSetR)
upset(fromList(varList),nsets = 20,text.scale = 2)

##-------------------------------------------------------------------##
##      Plot - distribution of VAF / across chromosomes            ####
##-------------------------------------------------------------------##
varToKeep = cgpVaf_QC$snv_ID[cgpVaf_QC$reasonForFail %in% c('likely_somatic') & !cgpVaf_QC$snv_ID %in% unique(cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline==T])]
cgpVaf_QC$snv_type = ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$sample_VAF > 0 & cgpVaf_QC$isGermline==F & cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$isGermline == T & cgpVaf_QC$sample_VAF >0]],'likely_germline_B',
                            ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$sample_VAF > 0 & cgpVaf_QC$isGermline==F & cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$isGermline == T & cgpVaf_QC$sample_VAF >0]],'likely_germline_BM',
                                   ifelse(cgpVaf_QC$CN_tot2major !='2:1' & cgpVaf_QC$reasonForFail != 'germline','CN_region',cgpVaf_QC$reasonForFail)))
cgpVaf_QC$snv_type[cgpVaf_QC$snv_type %in% c('low_depth','low_highQualRatio','overlapped_wINDELs')] = 'removed'
table(cgpVaf_QC$snv_type,cgpVaf_QC$reasonForFail)

cgpVaf_QC$snv_type = factor(cgpVaf_QC$snv_type,c('germline',"likely_germline_B", "likely_germline_BM",'CN_region','likely_somatic','PAR','removed'))


ggplot(cgpVaf_QC,aes(sample_VAF))+
  geom_density()+
  facet_grid(snv_type ~ tissue)+
  geom_vline(xintercept = 0.5,lty=2,lwd=0.7,col='grey')+
  geom_vline(xintercept = c(0.13),lty=2,lwd=0.7,col=c('red'))+
  #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
  theme_classic(base_size = 14)+theme(axis.line = element_line())

table(cgpVaf_QC$snv_type,cgpVaf_QC$sample)

## upset plot
varList = list('B_germline' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type %in% c('germline','likely_germline_B','likely_germline_BM') & cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$sample_VAF > 0],
               'BM_germline' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type %in% c('germline','likely_germline_B','likely_germline_BM') & cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$sample_VAF > 0],
               'B_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'likely_somatic' & cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$sample_VAF > 0],
               'BM_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'likely_somatic' & cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$sample_VAF > 0])

library(UpSetR)
upset(fromList(varList),nsets = 20,text.scale = 2)




cgpVaf_QC = cgpVaf_QC[order(cgpVaf_QC$snv_type),]
#cgpVaf_QC = cgpVaf_QC[cgpVaf_QC$snv_type != 'removed',]
cgpVaf_QC$Chrom = factor(cgpVaf_QC$Chrom,paste0('chr',c(1:22,'X','Y')))
ggplot(cgpVaf_QC,aes(Pos,sample_VAF))+
  geom_point(size=0.1,aes(col=snv_type,size=sample_DEP))+
  facet_grid(sample ~ Chrom,scales = 'free_x')+
  scale_size_manual(values = c(0.01,3))+
  scale_color_manual(values = c(grey(0.85),col25))+
  geom_hline(yintercept = 0.5,lty=2,lwd=0.4,col=grey(0.4))+
  geom_hline(yintercept = 0.13,lty=2,lwd=0.4,col='black')+
  #geom_vline(xintercept = c(0.13),lty=2,lwd=0.7,col=c('red'))+
  #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
  theme_classic(base_size = 12.5)+
  theme(panel.border = element_rect(fill=F),
        axis.line = element_blank(),axis.text.x = element_blank(),axis.ticks.length.x = unit(0,'cm')) +
  xlab('Genomic position') + ylab('VAF')







##--------------------------------------------------------------------------------------------------------##
##      Tim's shearwater-like filter: cgpVAF of the somatic variants in normal blood samples            ####
##--------------------------------------------------------------------------------------------------------##

outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L076'
somaticVar = cgpVaf_QC[cgpVaf_QC$snv_type == 'likely_somatic',]

## STEP1: cgpVAF of variants of interst in panel of normal samples ##

cgpVaf_onNormPanel(outDir,somaticVar)

## STEP2: alleleCount at each variant position across panel of normal samples ##
allelecount_dir = file.path(outDir,'2_allelecount')
columns = c('sample','Chrom','Pos')

##---- Import cgpVAF output with read filters applied
# Diagnostic
b_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076/PDv38is_wgs_PD62331c_snp_vaf.tsv',sep = '\t',skip = 36)
b_cgpVaf$tissue = 'Blood_D'
b_cgpVaf$snv_ID = paste0(b_cgpVaf$Chrom,':',b_cgpVaf$Pos,'_',b_cgpVaf$Ref,'/',b_cgpVaf$Alt)
b_cgpVaf$sample = 'PD62331c'

b_bed = b_cgpVaf[b_cgpVaf$snv_ID %in% somaticVar$snv_ID,c(columns,colnames(b_cgpVaf)[!grepl('PDv38is',colnames(b_cgpVaf)) & grepl('F*Z|R*Z',colnames(b_cgpVaf))])]
b_bed$Count_A = b_bed[,grepl('FAZ',colnames(b_bed))] + b_bed[,grepl('RAZ',colnames(b_bed))]
b_bed$Count_C = b_bed[,grepl('FCZ',colnames(b_bed))] + b_bed[,grepl('RCZ',colnames(b_bed))]
b_bed$Count_T = b_bed[,grepl('FTZ',colnames(b_bed))] + b_bed[,grepl('RTZ',colnames(b_bed))]
b_bed$Count_G = b_bed[,grepl('FGZ',colnames(b_bed))] + b_bed[,grepl('RGZ',colnames(b_bed))]
rownames(b_bed) = paste0(gsub('chr','',b_bed$Chrom),'_',b_bed$Pos)

write.table(b_bed[,c('Chrom','Pos','Count_A','Count_C','Count_G','Count_T')], file = file.path(allelecount_dir,paste0('PD62331c', "_allelecounts.txt")),
            sep = "\t", row.names = F, quote = F, col.names = T)

# TP1
bm_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076/PDv38is_wgs_PD62331a_snp_vaf.tsv',sep = '\t',skip = 36)
bm_cgpVaf$tissue = 'BM_D'
bm_cgpVaf$snv_ID = paste0(bm_cgpVaf$Chrom,':',bm_cgpVaf$Pos,'_',bm_cgpVaf$Ref,'/',bm_cgpVaf$Alt)
bm_cgpVaf$sample = 'PD62331a'


bm_bed = bm_cgpVaf[bm_cgpVaf$snv_ID %in% somaticVar$snv_ID,c(columns,colnames(bm_cgpVaf)[!grepl('PDv38is',colnames(bm_cgpVaf)) & grepl('F*Z|R*Z',colnames(bm_cgpVaf))])]
bm_bed$Count_A = bm_bed[,grepl('FAZ',colnames(bm_bed))] + bm_bed[,grepl('RAZ',colnames(bm_bed))]
bm_bed$Count_C = bm_bed[,grepl('FCZ',colnames(bm_bed))] + bm_bed[,grepl('RCZ',colnames(bm_bed))]
bm_bed$Count_T = bm_bed[,grepl('FTZ',colnames(bm_bed))] + bm_bed[,grepl('RTZ',colnames(bm_bed))]
bm_bed$Count_G = bm_bed[,grepl('FGZ',colnames(bm_bed))] + bm_bed[,grepl('RGZ',colnames(bm_bed))]
rownames(bm_bed) = paste0(gsub('chr','',bm_bed$Chrom),'_',bm_bed$Pos)
bm_bed = bm_bed[match(rownames(b_bed),rownames(bm_bed)),]
write.table(bm_bed[,c('Chrom','Pos','Count_A','Count_C','Count_G','Count_T')], file = file.path(allelecount_dir,paste0('PD62331a', "_allelecounts.txt")),
            sep = "\t", row.names = F, quote = F, col.names = T)



allele_count = b_bed[,c('Chrom','Pos','Count_A','Count_C','Count_G','Count_T')]


# Run the same thing again but this time for the normal panel.
alleleCount_onNormPanel(outDir,allele_count)
                                   

## STEP3: Run Tim's shearwater-like filter to remove false mutations ##

# Input for shearwater filter
samples_ID=c('PD62331a','PD62331c')

mutations = runShearWaterLike(samples_ID=c('PD62331a','PD62331c'),outDir=outDir,patient='L076')
  
# Tim uses <0.001 (FDR) as the mutations that passed the shearwater filter and I use <0.005 (FDR)
print('FDR = 0.001')
table(mutations[,5] < 0.001)
print('FDR = 0.005')
table(mutations[,5] < 0.005)
mutations$PD62331a_shearwaterPASS = ifelse(mutations$PD62331a < 0.001,T,F)
mutations$PD62331c_shearwaterPASS = ifelse(mutations$PD62331c < 0.001,T,F)
mutations$snv_ID = paste0(mutations$Chr,':',mutations$Pos,'_',mutations$Ref,'/',mutations$Alt)





##  STEP4: Add shearwater-like filter results to main object ##

cgpVaf_QC$snv_type = as.character(cgpVaf_QC$snv_type)
cgpVaf_QC$snv_type[cgpVaf_QC$sample == 'PD62331a' & 
                     cgpVaf_QC$snv_type != 'removed' &  
                     cgpVaf_QC$snv_ID %in% mutations$snv_ID[mutations$PD62331a_shearwaterPASS == F]] = 'shearwater_Failed'
cgpVaf_QC$snv_type[cgpVaf_QC$sample == 'PD62331c' & 
                     cgpVaf_QC$snv_type != 'removed' &  
                     cgpVaf_QC$snv_ID %in% mutations$snv_ID[mutations$PD62331c_shearwaterPASS == F]] = 'shearwater_Failed'
cgpVaf_QC$snv_type[cgpVaf_QC$reasonForFail == 'PAR'] = 'PAR'

cgpVaf_QC$snv_type = factor(cgpVaf_QC$snv_type,c('germline','likely_somatic',"likely_germline_B", "likely_germline_BM",'CN_region','shearwater_Failed','removed','PAR'))
cgpVaf_QC = cgpVaf_QC[order(cgpVaf_QC$snv_type),]
cgpVaf_QC$Chrom = factor(cgpVaf_QC$Chrom,paste0('chr',c(1:22,'X','Y')))
ggplot(cgpVaf_QC,aes(Pos,sample_VAF))+
  geom_point(size=0.1,aes(col=snv_type,size=sample_DEP))+
  facet_grid(sample ~ Chrom,scales = 'free_x')+
  scale_size_manual(values = c(0.01,3))+
  scale_color_manual(values = c(grey(0.8),'black',col25[1:3],'purple',grey(0.9),col25[5]))+
  geom_hline(yintercept = 0.5,lty=2,lwd=0.4,col=grey(0.4))+
  geom_hline(yintercept = 0.13,lty=2,lwd=0.4,col='black')+
  #geom_vline(xintercept = c(0.13),lty=2,lwd=0.7,col=c('red'))+
  #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
  theme_classic(base_size = 12)+
  theme(panel.border = element_rect(fill=F),
        axis.line = element_blank(),axis.text.x = element_blank(),axis.ticks.length.x = unit(0,'cm')) +
  xlab('Genomic position') + ylab('VAF')


##--- Upset plot
table(cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331a'] %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD62331c'])

varList = list('B_germline' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'germline' & cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$sample_VAF > 0],
               'BM_germline' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'germline' & cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$sample_VAF > 0],
               'B_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'likely_somatic' & cgpVaf_QC$sample == 'PD62331a' & cgpVaf_QC$sample_VAF > 0],
               'BM_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'likely_somatic' & cgpVaf_QC$sample == 'PD62331c' & cgpVaf_QC$sample_VAF > 0])
table(varList[[1]] %in% varList[[3]])

library(UpSetR)
upset(fromList(varList),nsets = 20,text.scale = 2)






##------ Extract somatic variants after shearwater filter
somaticVar = cgpVaf_QC[cgpVaf_QC$snv_type == 'likely_somatic',]
somaticVar$somaticVar_type = ifelse(somaticVar$snv_ID %in% intersect(somaticVar$snv_ID[somaticVar$sample != 'PD62331c' & somaticVar$sample_VAF > 0],
                                                                     somaticVar$snv_ID[somaticVar$sample == 'PD62331c' & somaticVar$sample_VAF > 0]),'shared',
                                    ifelse(somaticVar$sample == 'PD62331c' & somaticVar$sample_VAF > 0,'unique_B',
                                           ifelse(somaticVar$sample == 'PD62331a' & somaticVar$sample_VAF > 0,'unique_BM','others')))


# What are the reason for the unique variants to be unique
cgpVaf_QC$snv_type = as.character(cgpVaf_QC$snv_type)
cgpVaf_QC_b = cgpVaf_QC[cgpVaf_QC$sample == 'PD62331c',]
cgpVaf_QC_bm = cgpVaf_QC[cgpVaf_QC$sample == 'PD62331a',]

somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_B'] = paste0(somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_B'],':',
                                                                              cgpVaf_QC_bm$snv_type[match(somaticVar$snv_ID[somaticVar$somaticVar_type == 'unique_B'],cgpVaf_QC_bm$snv_ID)])
somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_BM'] = paste0(somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_BM'],':',
                                                                               cgpVaf_QC_b$snv_type[match(somaticVar$snv_ID[somaticVar$somaticVar_type == 'unique_BM'],cgpVaf_QC_b$snv_ID)])


table(somaticVar$somaticVar_type,somaticVar$sample)

# varList = list('B_somatic' = somaticVar$snv_ID[somaticVar$sample == 'PD62331c' & somaticVar$sample_VAF > 0],
#                'BM_somatic' = somaticVar$snv_ID[somaticVar$sample == 'PD62331a' & somaticVar$sample_VAF > 0 & somaticVar$somaticVar_type =='unique_TP1'],
#                'B_removed' = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD62331c' & cgpVaf_QC_og$germlineFilter == 'removed'],
#                'BM_removed' = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD62331a' & cgpVaf_QC_og$germlineFilter == 'removed'])
# 
# library(UpSetR)
# upset(fromList(varList),nsets = 20,text.scale = 2)
# 
# 
# ggplot(somaticVar,aes(sample_VAF))+
#   geom_density()+
#   facet_grid(somaticVar_type ~ tissue)+
#   geom_vline(xintercept = 0.5,lty=2,lwd=0.7,col='grey')+
#   #geom_vline(xintercept = c(expected_vaf_d),lty=2,lwd=0.7,col=c('red'))+
#   #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
#   theme_classic(base_size = 14)+theme(axis.line = element_line())

## Add variant annotation info from caveman
somaticVar$Gene = var4cgpVaf$Gene[match(somaticVar$snv_ID,var4cgpVaf$snv_ID)]
somaticVar$Impact = var4cgpVaf$Impact[match(somaticVar$snv_ID,var4cgpVaf$snv_ID)]
somaticVar$AAchange = var4cgpVaf$AAchange[match(somaticVar$snv_ID,var4cgpVaf$snv_ID)]

cgpVaf_withReadFilters_output$Gene = var4cgpVaf$Gene[match(cgpVaf_withReadFilters_output$snv_ID,var4cgpVaf$snv_ID)]
cgpVaf_withReadFilters_output$Impact = var4cgpVaf$Impact[match(cgpVaf_withReadFilters_output$snv_ID,var4cgpVaf$snv_ID)]
cgpVaf_withReadFilters_output$AAchange = var4cgpVaf$AAchange[match(cgpVaf_withReadFilters_output$snv_ID,var4cgpVaf$snv_ID)]

## pivot somaticVar wider to plot vaf of D vs TP1
df = cgpVaf_withReadFilters_output[cgpVaf_withReadFilters_output$snv_ID %in% somaticVar$snv_ID,]
#df = cgpVaf_withReadFilters_output[cgpVaf_withReadFilters_output$snv_ID %in% cgpVaf_QC$snv_ID,]
#df$somaticVar_type = cgpVaf_QC[cgpVaf_QC$somaticVar_type != 'others',]$somaticVar_type[match(df$snv_ID,cgpVaf_QC$snv_ID[cgpVaf_QC$somaticVar_type != 'others'])]
df$somaticVar_type = somaticVar[somaticVar$somaticVar_type != 'others',]$somaticVar_type[match(df$snv_ID,somaticVar$snv_ID[somaticVar$somaticVar_type != 'others'])]

df = pivot_wider(df,id_cols = c('snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect','Impact','AAchange','somaticVar_type'),
                 names_from = 'sample',values_from = 'sample_VAF')


# cgpVaf_QC_og$reasonForFail[cgpVaf_QC_og$reasonForFail == 'PASS' & cgpVaf_QC_og$pindelFilter == 'overlapped_wINDELs'] = 'overlapped_wINDELs'
# cgpVaf_QC_og$reasonForFail[cgpVaf_QC_og$reasonForFail == 'PASS' & cgpVaf_QC_og$germlineFilter == 'germline'] = 'germline'
# df$snv_removed_PD62331a = cgpVaf_QC_og[cgpVaf_QC_og$sample == 'PD62331a',]$reasonForFail[match(df$snv_ID,cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD62331a'])]
# df$snv_removed_PD62331c = cgpVaf_QC_og[cgpVaf_QC_og$sample == 'PD62331c',]$reasonForFail[match(df$snv_ID,cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD62331c'])]


#df$somaticVar_type[df$somaticVar_type=='unique_B/shearwater_Failed' & df$PD62331a == 0] = 'unique_B'
#df$somaticVar_type[df$somaticVar_type=='unique_BM/shearwater_Failed' & df$PD62331c == 0] = 'unique_BM'
#df$somaticVar_type = gsub('/',':',df$somaticVar_type)
#table(df$snv_removed,df$somaticVar_type)
#df$somaticVar_type_2 = ifelse(df$snv_removed == 'removed', 'removed',df$somaticVar_type)
ggplot(df,aes(PD62331c,PD62331a,col=somaticVar_type))+
  geom_point(size=0.7)+
  scale_color_manual(values = c(grey(0.8),grey(0.2),col25))+
  theme_classic(base_size = 14) + ggtitle('cgpVAF output of somatic variants')



# df$snv_removed = ifelse(df$snv_ID %in% cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$pindelFilter == 'removed' | cgpVaf_QC_og$germlineFilter == 'removed'],'removed','PASS')
# df$snv_removed = cgpVaf_QC_og$reasonForFail[match(df$snv_ID,cgpVaf_QC_og$snv_ID)]
# df$somaticVar_type_2 = df$snv_removed
# table(df$snv_removed,df$somaticVar_type)
# df$somaticVar_type_2 = ifelse(df$snv_removed == 'removed', 'removed',df$somaticVar_type)
# ggplot(df,aes(PD62331c,PD62331a,col=somaticVar_type_2))+
#   geom_point(size=0.7)+
#   facet_wrap(vars(Chrom))+
#   scale_color_manual(values = c(grey(0.3),grey(0.8),col25))+
#   theme_classic(base_size = 14) + ggtitle('cgpVAF output of somatic variants')
# 
# table(cgpVaf_QC_og[cgpVaf_QC_og$snv_ID %in% df$snv_ID[df$PD62331a == 0 & df$PD62331c > 0],]$sample,
#       cgpVaf_QC_og[cgpVaf_QC_og$snv_ID %in% df$snv_ID[df$PD62331a == 0 & df$PD62331c > 0],]$reasonForFail)


write.csv(df, '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_B.vs.BM_finalSomaticVariants.csv')
df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_B.vs.BM_finalSomaticVariants.csv')










##--------------------------------------------------##
##  Generate jbrowse images of these mutations    ####
##--------------------------------------------------##
df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076_B.vs.BM_finalSomaticVariants.csv')
## 1. Define JBROWSE URL
# L038 samples
jbrowse_url_l076 = '# JBROWSE https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F3030&loc=1%3A1..247467938&tracks=PD62331a_bwa%2CPD62331c_bwa'
# normal panel samples from 3010 project
normal_samples=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt", header=T) 
normalBlood_sampleIDs = normal_samples$Sample_ID
jbrowse_url_normalblood = paste0('# JBROWSE https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F3010&loc=1%3A1..247467938&tracks=',
                                 paste0(normalBlood_sampleIDs,collapse = '_bwa%2C'))


## 2. Prepare bed file containing location of these mutations, with +/ 50 bases around the mutation position
jbrowse_outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/jbrowse/L076/'
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
    for(d.sub in file.path(d,c('L076','normalPanel_3010'))){
      dir.create(d.sub,recursive = T)
    }
  }    
  
  for(sampleGroup in c('L076','normalPanel_3010')){
    ## Write bed file containing mutation of this category
    bed = var[,c('Chrom','start','end')]
    bedFile = file.path(d,sampleGroup,'L038_mutation.bed')
    if(sampleGroup == 'L076'){
      cat(jbrowse_url_l076,file = bedFile,sep = '\n') 
    }else{
      cat(jbrowse_url_normalblood,file = bedFile,sep = '\n') 
    }
    write.table(bed,file = bedFile,sep = '\t',quote = F,row.names = F,col.names = F,append = T)
    
    
    ## call jbrowse bash script
    cmd_list = c(cmd_list,sprintf('bash %s %s %s',jbrowse_script,bedFile,file.path(d,sampleGroup)))
  }
  
}

write.table(cmd_list,file = file.path(jbrowse_outDir,'jbrowse_command.sh'),quote = F,col.names = F,row.names = F)





##---------------------------------------------------##
## Finalising number of mutation in each sample    ####
##---------------------------------------------------##
## Go through jbrowse images and filter out "bad" variants
table(df$somaticVar_type)

shared_questionable = c('1_16588248-16588347','1_1702643-1702742','1_81390419-81390518',
                        '11_4344214-4344313','12_10421915-10422014','15_20187718-20187817',
                        '15_20201305-20201404','15_20206566-20206665','15_20235377-20235476','15_20255686-20255785',
                        '15_20263415-20263514','15_20345170-20345269','15_20414973-20415072','15_20440613-20440712','15_21257814-21257913',
                        '15_21317744-21317843','15_22612387-22612486','15_45071105-45071204','15_82478130-82478229',
                        '16_21529810-21529909','16_22583474-22583573','17_22137386-22137485','17_36159312-36159411','17_36211318-36211417',
                        '19_54903348-54903447','19_89255-89354','2_231832355-231832454','2_53233044-53233143',
                        '22_10705736-10705835','22_11285183-11285282','22_11328946-11329045','22_11328946-11329045','22_12205428-12205527','22_12389248-12389347','22_12390807-12390906',
                        '3_142513847-142513946','4_189635309-189635408','4_190112561-190112660','6_66670201-66670300','7_63579368-63579467','7_75105856-75105955',
                        '9_41036228-41036327','9_41043173-41043272','9_42388242-42388341',
                        'X_100550661-100550760')

# Chr15 is an acrocentric chromosome. Proper coverage only starts from ~25Mb
# NOTE: go through mutations on chr15,16,17 with someone - they look oddly high vaf and high depth? Normal panels sometimes also have higher coverage here, but not always... 
# is it a jbrowse thing where it randomly chooses to show more reads?
#   
# Check variants on Chr21 - should we just put this with the CN region and ignore them all?

  
## After manual inspection of jbrowse images, move bad variants to a seperate folder
for(mutType in c('shared','unique_B','unique_BM')){
  if(mutType == 'shared'){
    l076_files = list.files(file.path(jbrowse_outDir,mutType,'L076','JBROWSE'),pattern = '\\.png$',full.names = T)
    
    varToKeep = paste0(df$Chrom[df$somaticVar_type == 'shared'],'_',df$start[df$somaticVar_type == 'shared'] + 1,'-',df$end[df$somaticVar_type == 'shared'])
    filesToKeep = file.path(jbrowse_outDir,mutType,'L076','JBROWSE',paste0(varToKeep,'.png'))
    print(table(l076_files %in% filesToKeep))
    if(sum(l076_files %in% filesToKeep) < 10){
      stop('something is wrong....')
    }
    
    l076_files_toExcl = l076_files[!l076_files %in% filesToKeep]
    l076_files_toExcl_fn = gsub('\\.png$','',basename(l076_files_toExcl))
    cgpVaf_QC$tmp = paste0(gsub('chr','',cgpVaf_QC$Chrom),'_',cgpVaf_QC$Pos - 50 + 1,'-',cgpVaf_QC$Pos + 50)
    print(table(l076_files_toExcl_fn %in% cgpVaf_QC$tmp))
    if(sum(!l076_files_toExcl_fn %in% cgpVaf_QC$tmp) >0 | !all(l076_files_toExcl_fn %in% cgpVaf_QC$tmp)){
      stop('Something is wrong with the variants to exclulde')
    }
    
    t = cgpVaf_QC[cgpVaf_QC$tmp %in% l076_files_toExcl_fn,]
    for(type in unique(t$snv_type)){
      d = file.path(jbrowse_outDir,mutType,'L076','JBROWSE',type)
      dir.create(d)
      f = l076_files_toExcl[l076_files_toExcl_fn %in% t$tmp[t$snv_type == type]]
      system(sapply(f,function(s){sprintf('mv %s %s',s,d)}))
    }
    table(cgpVaf_QC$snv_type[cgpVaf_QC$tmp %in% l076_files_toExcl])
    
    
  }
}  


