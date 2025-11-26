## L038 - process caveman output between Diagnostic and TP1 samples ##
## for L038, we know that this case has: 
#       1. missense variant c.605G>A p.(Arg202Gln) in ETV6; 
#       2. missense variant c.2038G>A p.(Val680Met) in EZH2 
# we can use this information to determine tumour purity 
# vaf of EZH2 variant is 0.11 in Diagnostic sample --> tumour purity D is 0.11*100/0.5 = 22%
# vaf of EZH2 variant is 0.2040816 in TP1 sample --> tumour purity D is 0.2040816*100/0.5 = 40%
# Couldn't detect any ETV6 variants though...

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
cnSegs = GRanges(c('chr5','chr17','chr21'),
                 IRanges(c(rep(1,3)),c(1e9,27.5e6,1e9)))



##------------------------------------------------##
##    Import Caveman output - Diagnostic        ####
##------------------------------------------------##
d_caveman = caveman2(sample='PD60301a',projectid = 3030)
d_caveman$timepoint = 'Diagnostic'
d_caveman = d_caveman[d_caveman$Flag == T,]
d_caveman$snv_ID = paste0(d_caveman$Chr,':',d_caveman$Pos,'_',d_caveman$Ref,'/',d_caveman$Alt)



##-----------------------------------------------##
##      Import Caveman output - TP1            ####
##-----------------------------------------------##
tp1_caveman = caveman2(sample='PD61846a',projectid = 3030)
tp1_caveman$timepoint = 'TP1'
tp1_caveman = tp1_caveman[tp1_caveman$Flag == T,]
tp1_caveman$snv_ID = paste0(tp1_caveman$Chr,':',tp1_caveman$Pos,'_',tp1_caveman$Ref,'/',tp1_caveman$Alt)



var4cgpVaf = rbind(d_caveman,tp1_caveman)
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
# write_delim(bed[,colnames(bed) != 'snv_ID'],'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L038_unmatched_cavemanVar.bed',col_names = F,delim = '\t')
# write.csv(c('PD60301a','PD61846a'),'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L038_samples.csv',row.names = F,quote = F)
# 
# ##--- prepare for cgpVAF
# # cd ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis
# # mkdir cgpVaf_input/L038
# # mkdir cgpVaf_output/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD60301a/PD60301a.sample.dupmarked.bam cgpVaf_input/L038/  
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD60301a/PD60301a.sample.dupmarked.bam.bai cgpVaf_input/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD60301a/PD60301a.caveman_c.annot.vcf.gz cgpVaf_input/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD61846a/PD61846a.sample.dupmarked.bam cgpVaf_input/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD61846a/PD61846a.sample.dupmarked.bam.bai cgpVaf_input/L038  
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD61846a/PD61846a.caveman_c.annot.vcf.gz cgpVaf_input/L038
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam cgpVaf_input/L038/
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam.bai cgpVaf_input/L038/
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam cgpVaf_input/L038/
# 
# ##--- Run cgpVAF
# # bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/x10_MLDS_cgpVaf.sh L038 PD60301a  
# # bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/x10_MLDS_cgpVaf.sh L038 PD61846a
# 
# 
# 
# 
##---- Import cgpVAF output with read filters applied
# Diagnostic
d_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038/PDv38is_wgs_PD60301a_snp_vaf.tsv',sep = '\t',skip = 36)
d_cgpVaf$timepoint='Diagnostic'
d_cgpVaf$snv_ID = paste0(d_cgpVaf$Chrom,':',d_cgpVaf$Pos,'_',d_cgpVaf$Ref,'/',d_cgpVaf$Alt)
d_cgpVaf$sample = 'PD60301a'
d_cgpVaf$sample_MTR = d_cgpVaf[,grepl("MTR",colnames(d_cgpVaf))&!colnames(d_cgpVaf)%in%"PDv38is_wgs_MTR"]
d_cgpVaf$sample_VAF = d_cgpVaf[,grepl("VAF",colnames(d_cgpVaf))&!colnames(d_cgpVaf)%in%"PDv38is_wgs_VAF"]
d_cgpVaf$sample_DEP = d_cgpVaf[,grepl("DEP",colnames(d_cgpVaf))&!colnames(d_cgpVaf)%in%"PDv38is_wgs_DEP"]

# TP1
tp1_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038/PDv38is_wgs_PD61846a_snp_vaf.tsv',sep = '\t',skip = 36)
tp1_cgpVaf$timepoint='TP1'
tp1_cgpVaf$snv_ID = paste0(tp1_cgpVaf$Chrom,':',tp1_cgpVaf$Pos,'_',tp1_cgpVaf$Ref,'/',tp1_cgpVaf$Alt)
tp1_cgpVaf$sample = 'PD61846a'
tp1_cgpVaf$sample_MTR = tp1_cgpVaf[,grepl("MTR",colnames(tp1_cgpVaf))&!colnames(tp1_cgpVaf)%in%"PDv38is_wgs_MTR"]
tp1_cgpVaf$sample_VAF = tp1_cgpVaf[,grepl("VAF",colnames(tp1_cgpVaf))&!colnames(tp1_cgpVaf)%in%"PDv38is_wgs_VAF"]
tp1_cgpVaf$sample_DEP = tp1_cgpVaf[,grepl("DEP",colnames(tp1_cgpVaf))&!colnames(tp1_cgpVaf)%in%"PDv38is_wgs_DEP"]



columns = c('sample','timepoint','snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect','sample_MTR','sample_VAF','sample_DEP')
cgpVaf_withReadFilters_output = rbind(d_cgpVaf[,columns],tp1_cgpVaf[,columns])
dim(cgpVaf_withReadFilters_output)



##---- Import cgpVAF output with NO read filters applied
# Diagnostic
d_cgpVaf_noFilter = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_nofilter/PDv38is_wgs_PD60301a_snp_vaf.tsv',sep = '\t',skip = 36)
d_cgpVaf_noFilter$timepoint='Diagnostic'
d_cgpVaf_noFilter$snv_ID = paste0(d_cgpVaf_noFilter$Chrom,':',d_cgpVaf_noFilter$Pos,'_',d_cgpVaf_noFilter$Ref,'/',d_cgpVaf_noFilter$Alt)
d_cgpVaf_noFilter$sample = 'PD60301a'
d_cgpVaf_noFilter$sample_MTR = d_cgpVaf_noFilter[,grepl("MTR",colnames(d_cgpVaf_noFilter))&!colnames(d_cgpVaf_noFilter)%in%"PDv38is_wgs_MTR"]
d_cgpVaf_noFilter$sample_VAF = d_cgpVaf_noFilter[,grepl("VAF",colnames(d_cgpVaf_noFilter))&!colnames(d_cgpVaf_noFilter)%in%"PDv38is_wgs_VAF"]
d_cgpVaf_noFilter$sample_DEP = d_cgpVaf_noFilter[,grepl("DEP",colnames(d_cgpVaf_noFilter))&!colnames(d_cgpVaf_noFilter)%in%"PDv38is_wgs_DEP"]

# TP1
tp1_cgpVaf_noFilter = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_nofilter/PDv38is_wgs_PD61846a_snp_vaf.tsv',sep = '\t',skip = 36)
tp1_cgpVaf_noFilter$timepoint='TP1'
tp1_cgpVaf_noFilter$snv_ID = paste0(tp1_cgpVaf_noFilter$Chrom,':',tp1_cgpVaf_noFilter$Pos,'_',tp1_cgpVaf_noFilter$Ref,'/',tp1_cgpVaf_noFilter$Alt)
tp1_cgpVaf_noFilter$sample = 'PD61846a'
tp1_cgpVaf_noFilter$sample_MTR = tp1_cgpVaf_noFilter[,grepl("MTR",colnames(tp1_cgpVaf_noFilter))&!colnames(tp1_cgpVaf_noFilter)%in%"PDv38is_wgs_MTR"]
tp1_cgpVaf_noFilter$sample_VAF = tp1_cgpVaf_noFilter[,grepl("VAF",colnames(tp1_cgpVaf_noFilter))&!colnames(tp1_cgpVaf_noFilter)%in%"PDv38is_wgs_VAF"]
tp1_cgpVaf_noFilter$sample_DEP = tp1_cgpVaf_noFilter[,grepl("DEP",colnames(tp1_cgpVaf_noFilter))&!colnames(tp1_cgpVaf_noFilter)%in%"PDv38is_wgs_DEP"]



columns = c('sample','timepoint','snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect','sample_MTR','sample_VAF','sample_DEP')
cgpVaf_noReadFilters_output = rbind(d_cgpVaf_noFilter[,columns],tp1_cgpVaf_noFilter[,columns])
dim(cgpVaf_noReadFilters_output)
colnames(cgpVaf_noReadFilters_output)[colnames(cgpVaf_noReadFilters_output) %in% c('sample_MTR','sample_VAF','sample_DEP')] = paste0(colnames(cgpVaf_noReadFilters_output)[colnames(cgpVaf_noReadFilters_output) %in% c('sample_MTR','sample_VAF','sample_DEP')],':noFilter')

cgpVaf_QC_og = merge(cgpVaf_withReadFilters_output,cgpVaf_noReadFilters_output,by=c('sample','timepoint','snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect'),all=T)
#cgpVaf_QC = merge(cgpVaf_withReadFilters_output,cgpVaf_noReadFilters_output,by=c('sample','timepoint','snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect'),all=T)


##----------------------------------------------------------------##
##      cgpVAF depth filter - consider sex chromosome           ####
##----------------------------------------------------------------##
gender = 'female'
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
cgpVaf_QC_og_TP1_CNsegs = subsetByOverlaps(cgpVaf_QC_og_gr[cgpVaf_QC_og_gr$sample == 'PD61846a'],cnSegs)
toExclude = c(which(cgpVaf_QC_og$Chrom == 'chr21'),
              which(cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID & cgpVaf_QC_og$sample == 'PD61846a'))
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
table(cgpVaf_QC_og$highQualRatio_filter[cgpVaf_QC_og$depthFilter != 'tooLow_DEP'],cgpVaf_QC_og$sample[cgpVaf_QC_og$depthFilter != 'tooLow_DEP'])

lowHiQualVar = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$highQualRatio_filter == 'low']
cgpVaf_QC_og$highQualRatio_filter[cgpVaf_QC_og$snv_ID %in% lowHiQualVar] = 'low'


cgpVaf_QC_og$reasonForFail = ifelse(cgpVaf_QC_og$depthFilter == 'tooLow_DEP' | cgpVaf_QC_og$sample_DEP == 0,'low_depth',
                                    ifelse(cgpVaf_QC_og$highQualRatio_filter == 'low','low_highQualRatio','PASS'))
table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$sample)

write.csv(cgpVaf_QC_og,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_cgpVafQCog_depth.highQualRatio_filtered.csv')



##----------------------------------##
##      Pindel filter             ####
##----------------------------------##
# Including FAILED Pindel variants too

intv = 10

# pindel filter
library(data.table)
indels = pindel(sample = 'PD60301a',projectid=3030,filter_lowQual=F)
dim(indels)
n_distinct(indels$varID)
indels$End=indels$Pos+nchar(indels$Ref)-1
indels$sampleID = 'PD60301a'

tp1_indels = pindel(sample = 'PD61846a',projectid=3030,filter_lowQual=F)
dim(tp1_indels)
n_distinct(tp1_indels$varID)
tp1_indels$sampleID = 'PD61846a'
tp1_indels_gr = GRanges(tp1_indels$Chr,IRanges(tp1_indels$Pos - intv,tp1_indels$End + intv),
                        varID = tp1_indels$varID,sampleID=tp1_indels$sampleID)

mcols(tp1_indels_gr) = cbind(mcols(tp1_indels_gr),tp1_indels[match(paste0(tp1_indels_gr$varID,':',tp1_indels_gr$sampleID),
                                                                   paste0(tp1_indels$varID,':',tp1_indels$sampleID)),])
length(tp1_indels_gr) == nrow(tp1_indels)

# Remove copy number regions in TP1 sample
cnSegs_indels = subsetByOverlaps(tp1_indels_gr,cnSegs)

tp1_indels = tp1_indels[!tp1_indels$varID %in% cnSegs_indels$varID,]


## combine indels from Diagnostic and TP1 samples
indels = rbind(indels,tp1_indels)
indels_gr = GRanges(indels$Chr,IRanges(indels$Pos - intv,indels$End + intv),
                    varID = indels$varID,sampleID=indels$sampleID)
mcols(indels_gr) = cbind(mcols(indels_gr),indels[match(paste0(indels_gr$varID,':',indels_gr$sampleID),
                                                       paste0(indels$varID,':',indels$sampleID)),])
length(indels_gr) == nrow(indels)
indels = indels_gr

# find mutations within 10bp of indels
intv = 10
cgpVaf_QC_og = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_cgpVafQCog_depth.highQualRatio_filtered.csv')
cgpVaf_QC = cgpVaf_QC_og[cgpVaf_QC_og$reasonForFail == 'PASS',]

## Convert cgpVaf_QC into Granges
cgpVaf_QC_gr = GRanges(cgpVaf_QC$Chrom,IRanges(cgpVaf_QC$Pos,cgpVaf_QC$Pos),
                    snv_ID = cgpVaf_QC$snv_ID,sample=cgpVaf_QC$sample,timepoint=cgpVaf_QC$timepoint)
mcols(cgpVaf_QC_gr) = cbind(mcols(cgpVaf_QC_gr),cgpVaf_QC[match(paste0(cgpVaf_QC_gr$snv_ID,':',cgpVaf_QC_gr$sample,':',cgpVaf_QC_gr$timepoint),
                                                                paste0(cgpVaf_QC_gr$snv_ID,':',cgpVaf_QC_gr$sample,':',cgpVaf_QC_gr$timepoint)),])
length(cgpVaf_QC_gr) == nrow(cgpVaf_QC)
cgpVaf_QC = cgpVaf_QC_gr

# Identify variants overlapping with indels
cgpVaf_QC_overlapped_wIndels = subsetByOverlaps(cgpVaf_QC,indels)

cgpVaf_QC$overlapped_with_indels = ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC_overlapped_wIndels$snv_ID,T,F)
table(cgpVaf_QC$overlapped_with_indels,cgpVaf_QC$reasonForFail)

cgpVaf_QC = as.data.frame(mcols(cgpVaf_QC))

write.csv(cgpVaf_QC,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_cgpVaf_depth.highQualRatio.Pindel_filtered.csv')
cgpVaf_QC = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_cgpVaf_depth.highQualRatio.Pindel_filtered.csv')

cgpVaf_QC_og$pindelFilter = ifelse(cgpVaf_QC_og$sample == 'PD60301a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$overlapped_with_indels == F],'PASS',
                                   ifelse(cgpVaf_QC_og$sample == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$overlapped_with_indels == F],'PASS',
                                          ifelse(cgpVaf_QC_og$sample == 'PD60301a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$overlapped_with_indels == T],'overlapped_wINDELs',
                                                 ifelse(cgpVaf_QC_og$sample == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$overlapped_with_indels == T],'overlapped_wINDELs','removed'))))


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
cgpVaf_QC_og$refCnt[cgpVaf_QC_og$sample == 'PD60301a'] =   d_cgpVaf$PD60301a_WTR[match(cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD60301a'],d_cgpVaf$snv_ID)]
cgpVaf_QC_og$refCnt[cgpVaf_QC_og$sample == 'PD61846a'] =   tp1_cgpVaf$PD61846a_WTR[match(cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD61846a'],tp1_cgpVaf$snv_ID)]

#colnames(cgpVaf_withReadFilters_output) = gsub('\\.sample_MTR$','',colnames(cgpVaf_withReadFilters_output))
cgpVaf_QC_og$altCnt = cgpVaf_withReadFilters_output$sample_MTR[match(paste0(cgpVaf_QC_og$snv_ID,'_',cgpVaf_QC_og$sample),
                                                                  paste0(cgpVaf_withReadFilters_output$snv_ID,'_',cgpVaf_withReadFilters_output$sample))]

cgpVaf_QC_og$DP = cgpVaf_withReadFilters_output$sample_DEP[match(paste0(cgpVaf_QC_og$snv_ID,'_',cgpVaf_QC_og$sample),
                                                                  paste0(cgpVaf_withReadFilters_output$snv_ID,'_',cgpVaf_withReadFilters_output$sample))]

cgpVaf_QC_og$sample_MTR = cgpVaf_withReadFilters_output$sample_MTR[match(paste0(cgpVaf_QC_og$snv_ID,'_',cgpVaf_QC_og$sample),
                                                                      paste0(cgpVaf_withReadFilters_output$snv_ID,'_',cgpVaf_withReadFilters_output$sample))]

# Exclude SNPs from CN regions in TP1 sample
cgpVaf_QC_og$CN_tot2major = '2:1'
cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$Chrom == 'chr21'] = '3:2'
cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$sample == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID[cgpVaf_QC_og_TP1_CNsegs$Chrom == 'chr5']] = '1:0'
cgpVaf_QC_og$CN_tot2major[cgpVaf_QC_og$sample == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og_TP1_CNsegs$snv_ID[cgpVaf_QC_og_TP1_CNsegs$Chrom == 'chr17']] = '1:0'

## Perform binomial test filter
cgpVaf_QC_og.sub = snvBinomTest(gender='female',data=cgpVaf_QC_og[cgpVaf_QC_og$pindelFilter == 'PASS',],
                                genotype='T21',pCut=0.05)
table(cgpVaf_QC_og.sub$isGermline,cgpVaf_QC_og.sub$sample)
table(cgpVaf_QC_og.sub$isGermline,cgpVaf_QC_og.sub$Chrom)

cgpVaf_QC_og$isGermline = cgpVaf_QC_og.sub$isGermline[match(paste0(cgpVaf_QC_og$snv_ID,':',cgpVaf_QC_og$sample),
                                                            paste0(cgpVaf_QC_og.sub$snv_ID,':',cgpVaf_QC_og.sub$sample))]


cgpVaf_QC_og$germlineFilter = ifelse(cgpVaf_QC_og$sample == 'PD60301a' &
                                       cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$sample == 'PD60301a' & cgpVaf_QC_og.sub$isGermline == T],'germline',
                                     ifelse(cgpVaf_QC_og$sample == 'PD61846a' &
                                              cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$sample == 'PD61846a' & cgpVaf_QC_og.sub$isGermline == T],'germline',
                                            ifelse(cgpVaf_QC_og$sample == 'PD60301a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$sample == 'PD60301a' & cgpVaf_QC_og.sub$isGermline == F],'likely_somatic',
                                                   ifelse(cgpVaf_QC_og$sample == 'PD61846a' & cgpVaf_QC_og$snv_ID %in% cgpVaf_QC_og.sub$snv_ID[cgpVaf_QC_og.sub$sample == 'PD61846a' & cgpVaf_QC_og.sub$isGermline == F],'likely_somatic','others'))))


table(cgpVaf_QC_og$germlineFilter,cgpVaf_QC_og$pindelFilter)

## Summarise reason for fails
cgpVaf_QC_og$reasonForFail = ifelse(cgpVaf_QC_og$CN_tot2major != '2:1' & cgpVaf_QC_og$germlineFilter != 'germline','CN_region',
                                    ifelse(cgpVaf_QC_og$germlineFilter == 'germline','germline',
                                           ifelse(cgpVaf_QC_og$germlineFilter == 'likely_somatic','likely_somatic',
                                                  ifelse((!is.na(cgpVaf_QC_og$depthFilter) & cgpVaf_QC_og$depthFilter == 'tooLow_DEP') | cgpVaf_QC_og$sample_DEP == 0,'low_depth',
                                                         ifelse(cgpVaf_QC_og$highQualRatio_filter == 'low','low_highQualRatio',
                                                                ifelse(cgpVaf_QC_og$pindelFilter == 'overlapped_wINDELs','overlapped_wINDELs','others'))))))


table(is.na(cgpVaf_QC_og$reasonForFail))
table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$sample)
table(cgpVaf_QC_og$reasonForFail,cgpVaf_QC_og$germlineFilter)


ggplot(cgpVaf_QC_og,aes(sample,fill = reasonForFail))+
  geom_bar(position='fill')+
  scale_fill_manual(values = c(col25[1],grey(0.8),col25[-1]))+
  theme_classic(base_size = 13)


write.csv(cgpVaf_QC_og,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_cgpVaf_filtersApplied.csv')
cgpVaf_QC_og = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_cgpVaf_filtersApplied.csv')

cgpVaf_QC = cgpVaf_QC_og[cgpVaf_QC_og$germlineFilter == 'likely_somatic',]
cgpVaf_QC = cgpVaf_QC_og

##----------------------------------------------------##
##      Plot - number of shared variants            ####
##----------------------------------------------------##
table(cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD60301a'] %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD61846a'])

varList = list('D_germline' = cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline == T & cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$sample_VAF > 0],
               'TP1_germline' = cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline == T & cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$sample_VAF > 0],
               'D_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$reasonForFail == 'likely_somatic' & !is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline != T & cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$sample_VAF > 0],
               'TP1_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$reasonForFail == 'likely_somatic' & !is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline != T & cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$sample_VAF > 0])
table(varList[[1]] %in% varList[[3]])

library(UpSetR)
upset(fromList(varList),nsets = 20,text.scale = 2)

##-------------------------------------------------------------------##
##      Plot - distribution of VAF / across chromosomes            ####
##-------------------------------------------------------------------##
varToKeep = cgpVaf_QC$snv_ID[!is.na(cgpVaf_QC$isGermline) & cgpVaf_QC$isGermline == F & !cgpVaf_QC$snv_ID %in% unique(cgpVaf_QC$snv_ID[cgpVaf_QC$isGermline==T])]
cgpVaf_QC$snv_type = ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$isGermline == F & cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$isGermline == T]],'likely_germline_D',
                            ifelse(cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$isGermline == F & cgpVaf_QC$snv_ID %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$isGermline == T]],'likely_germline_TP1',
                                   ifelse(cgpVaf_QC$CN_tot2major !='2:1' & cgpVaf_QC$reasonForFail != 'germline','CN_region',cgpVaf_QC$reasonForFail)))


cgpVaf_QC$snv_type[cgpVaf_QC$snv_type %in% c('low_depth','low_highQualRatio','overlapped_wINDELs')] = 'removed'
table(cgpVaf_QC$snv_type,cgpVaf_QC$reasonForFail)
table(cgpVaf_QC$snv_type,cgpVaf_QC$CN_tot2major,cgpVaf_QC$sample)

ggplot(cgpVaf_QC[cgpVaf_QC$snv_type != 'removed',],aes(sample_VAF))+
  geom_density()+
  facet_grid(snv_type ~ timepoint)+
  geom_vline(xintercept = 0.5,lty=2,lwd=0.7,col='grey')+
  geom_vline(xintercept = c(0.13),lty=2,lwd=0.7,col=c('red'))+
  #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
  theme_classic(base_size = 12)+theme(axis.line = element_line())

table(cgpVaf_QC$snv_type,cgpVaf_QC$sample)


## Before shearwater-filter
cgpVaf_QC$snv_type = factor(cgpVaf_QC$snv_type,c('germline',"likely_germline_D", "likely_germline_TP1",'CN_region','likely_somatic','removed'))
cgpVaf_QC = cgpVaf_QC[order(cgpVaf_QC$snv_type),]
cgpVaf_QC$Chrom = factor(cgpVaf_QC$Chrom,paste0('chr',c(1:22,'X','Y')))

ggplot(cgpVaf_QC,aes(Pos,sample_VAF))+
  geom_point(size=0.1,aes(col=snv_type,size=sample_DEP))+
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

outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038'
somaticVar = cgpVaf_QC[cgpVaf_QC$snv_type == 'likely_somatic',]
##--- Upset plot
varList = list('D_somatic' = somaticVar$snv_ID[somaticVar$sample == 'PD60301a' & somaticVar$sample_VAF > 0],
               'TP1_somatic' = somaticVar$snv_ID[somaticVar$sample == 'PD61846a' & somaticVar$sample_VAF > 0])
library(UpSetR)
upset(fromList(varList),nsets = 20,text.scale = 2)


##--------------------------------------------------------------------------##
##      cgpVAF of the somatic variants in normal blood samples            ####
##--------------------------------------------------------------------------##
## STEP1: cgpVAF of variants of interst in panel of normal samples ##

cgpVaf_onNormPanel(outDir,somaticVar)


# 
# 
# ## Run cgpVaf on the normal samples
# ## 1. Create bed files
# # requried format:
# #   chr pos ref_allele alt_allele
# #   filename ends with .bed, delim='\t'
# somaticVar = cgpVaf_QC[cgpVaf_QC$snv_type == 'somatic',]
# 
# bed = data.frame(snv_ID = unique(somaticVar$snv_ID))
# bed = merge(bed,somaticVar[,c('Chrom','Pos','Ref','Alt','snv_ID')],by='snv_ID',all.x=T)
# bed = bed[!duplicated(bed),]
# dim(bed)
# n_distinct(somaticVar$snv_ID)
# bed = bed[order(bed$Chrom),]
# write_delim(bed[,colnames(bed) != 'snv_ID'],'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_potential_somaticVar_forShearwater.bed',col_names = F,delim = '\t')
# 
# 
# ## Import list of normal samples
# inDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_cgpVaf/cgpVaf_input/'
# 
# 
# normal_samples=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt", header=T) 
# for(i in 1:nrow(normal_samples)){
#   sampleID = normal_samples$Sample_ID[i]
#   projectID = normal_samples$Project_ID[i]
#   path = file.path('/nfs/cancer_ref01/nst_links/live',projectID,sampleID,paste0(sampleID,'.sample.dupmarked.bam'))
#   if(!file.exists(path)){
#     stop(sprintf('[%s] File does not exist: %s',sampleID,path))
#   }
#   
#   system(sprintf('cp -s %s %s',path,inDir))
#   system(sprintf('cp -s %s %s',paste0(path,'.bai'),inDir))
# 
# }
# 
# system('cp -s /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_input/L038/PDv38is_wgs.sample.dupmarked.bam /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_cgpVaf/cgpVaf_input/PDv38is_wgs.sample.dupmarked.bam')
# system('cp -s /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_input/L038/PDv38is_wgs.sample.dupmarked.bam.bai /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_cgpVaf/cgpVaf_input/PDv38is_wgs.sample.dupmarked.bam.bai')
# 
# script = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/x10_MLDS_cgpVaf_forShearwater.sh'
# cgpVaf_outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_cgpVaf/cgpVaf_output/'
# bedFile = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_potential_somaticVar_forShearwater.bed'
# sampleFile = "/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt"
# #sampleFile = "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/test_sample.txt"
# system(sprintf('bash %s %s %s %s %s',script,inDir,cgpVaf_outDir,bedFile,sampleFile))
# 
# 
# 
# 
# ## Create cgpVaf output folder
# cgpVaf_outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_cgpVaf_output'
# if(!dir.exists(cgpVaf_outDir)){
#   dir.create(cgpVaf_outDir,recursive = T)
# }
# 
# 
# 
# 
# 
# 
# 
# ## Run ShearWater on these potential somatic variants of interest
# sw_outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038'
# if(!dir.exists(sw_outDir)){
#   dir.create(sw_outDir,recursive = T)
#   subFolders = c('1_cgpVaf','1_cgpVaf/cgpVaf_output','1_cgpVaf/cgpVaf_input')
#   for(s in subFolders){
#     dir.create(file.path(sw_outDir,s),recursive = T)
#   }
# }






## STEP2: alleleCount at each variant position across panel of normal samples ##
allelecount_dir = file.path(outDir,'2_allelecount')
columns = c('sample','Chrom','Pos')


# allelecount_dir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/2_allelecount/'
# if(!dir.exists(allelecount_dir)){
#   dir.create(allelecount_dir,recursive = T)
# }



##---- Import cgpVAF output with read filters applied
# Diagnostic
d_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038/PDv38is_wgs_PD60301a_snp_vaf.tsv',sep = '\t',skip = 36)
d_cgpVaf$timepoint='Diagnostic'
d_cgpVaf$snv_ID = paste0(d_cgpVaf$Chrom,':',d_cgpVaf$Pos,'_',d_cgpVaf$Ref,'/',d_cgpVaf$Alt)
d_cgpVaf$sample = 'PD60301a'

d_bed = d_cgpVaf[d_cgpVaf$snv_ID %in% somaticVar$snv_ID,c(columns,colnames(d_cgpVaf)[!grepl('PDv38is',colnames(d_cgpVaf)) & grepl('F*Z|R*Z',colnames(d_cgpVaf))])]
d_bed$Count_A = d_bed[,grepl('FAZ',colnames(d_bed))] + d_bed[,grepl('RAZ',colnames(d_bed))]
d_bed$Count_C = d_bed[,grepl('FCZ',colnames(d_bed))] + d_bed[,grepl('RCZ',colnames(d_bed))]
d_bed$Count_T = d_bed[,grepl('FTZ',colnames(d_bed))] + d_bed[,grepl('RTZ',colnames(d_bed))]
d_bed$Count_G = d_bed[,grepl('FGZ',colnames(d_bed))] + d_bed[,grepl('RGZ',colnames(d_bed))]
rownames(d_bed) = paste0(gsub('chr','',d_bed$Chrom),'_',d_bed$Pos)

write.table(d_bed[,c('Chrom','Pos','Count_A','Count_C','Count_G','Count_T')], file = file.path(allelecount_dir,paste0('PD60301a', "_allelecounts.txt")),
            sep = "\t", row.names = F, quote = F, col.names = T)



# TP1
tp1_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038/PDv38is_wgs_PD61846a_snp_vaf.tsv',sep = '\t',skip = 36)
tp1_cgpVaf$timepoint='TP1'
tp1_cgpVaf$snv_ID = paste0(tp1_cgpVaf$Chrom,':',tp1_cgpVaf$Pos,'_',tp1_cgpVaf$Ref,'/',tp1_cgpVaf$Alt)
tp1_cgpVaf$sample = 'PD61846a'

tp1_bed = tp1_cgpVaf[tp1_cgpVaf$snv_ID %in% somaticVar$snv_ID,c(columns,colnames(tp1_cgpVaf)[!grepl('PDv38is',colnames(tp1_cgpVaf)) & grepl('F*Z|R*Z',colnames(tp1_cgpVaf))])]
tp1_bed$Count_A = tp1_bed[,grepl('FAZ',colnames(tp1_bed))] + tp1_bed[,grepl('RAZ',colnames(tp1_bed))]
tp1_bed$Count_C = tp1_bed[,grepl('FCZ',colnames(tp1_bed))] + tp1_bed[,grepl('RCZ',colnames(tp1_bed))]
tp1_bed$Count_T = tp1_bed[,grepl('FTZ',colnames(tp1_bed))] + tp1_bed[,grepl('RTZ',colnames(tp1_bed))]
tp1_bed$Count_G = tp1_bed[,grepl('FGZ',colnames(tp1_bed))] + tp1_bed[,grepl('RGZ',colnames(tp1_bed))]
rownames(tp1_bed) = paste0(gsub('chr','',tp1_bed$Chrom),'_',tp1_bed$Pos)
tp1_bed = tp1_bed[match(rownames(d_bed),rownames(tp1_bed)),]
write.table(tp1_bed[,c('Chrom','Pos','Count_A','Count_C','Count_G','Count_T')], file = file.path(allelecount_dir,paste0('PD61846a', "_allelecounts.txt")),
            sep = "\t", row.names = F, quote = F, col.names = T)


# alleleCount = rbind(d_bed[,c('sample','Chrom','Pos','Count_A','Count_C','Count_G','Count_T')],tp1_bed[,c('sample','Chrom','Pos','Count_A','Count_C','Count_G','Count_T')])
# for(s in unique(alleleCount$sample)){
#   d = alleleCount[alleleCount$sample == s,]
#   write.table(d[,colnames(d) != 'sample'], file = file.path(allelecount_dir,paste0(s, "_allelecounts.txt")),
#               sep = "\t", row.names = F, quote = F, col.names = T)
# }
#rownames(alleleCount)

allele_count = d_bed[,c('Chrom','Pos','Count_A','Count_C','Count_G','Count_T')]

# Run the same thing again but this time for the normal panel.
alleleCount_onNormPanel(outDir,allele_count)




# # This will be used for the shearwater analysis
# normal_samples=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt", header=T) 
# 
# normSamples_cgpVaf = read.table(paste0("~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_cgpVaf/cgpVaf_output/PDv38is_wgs_PD55504c_snp_vaf.tsv"), 
#                    header = TRUE, sep = '\t', stringsAsFactors = F, comment.char = "",skip = 64)
# colnames(normSamples_cgpVaf)[which(colnames(normSamples_cgpVaf)=="VariantID")]="ID"
# normSamples_cgpVaf=normSamples_cgpVaf[!duplicated(normSamples_cgpVaf[,c("Chrom", "Pos")]),]
# row.names(normSamples_cgpVaf)=paste0(sub('.*chr', '', normSamples_cgpVaf$Chrom), "_", normSamples_cgpVaf$Pos)
# #Reorder normSamples_cgpVaf to have the same order of rows as allele count
# table(row.names(allele_count) %in% row.names(normSamples_cgpVaf))
# normSamples_cgpVaf=normSamples_cgpVaf[match(row.names(allele_count), row.names(normSamples_cgpVaf)), ]
# 
# for (i in 1:nrow(normal_samples)) {
#   case=normal_samples$Case[i]
#   print(case)
#   num_samples=normal_samples$Sample_ID[startsWith(normal_samples$Sample_ID, case)]
#   
#   for (j in 1:length(num_samples)) {
#     sample=num_samples[j]
#     print(sample)
#     mut_temp=normSamples_cgpVaf[,grep(paste0(sample,"_"), colnames(normSamples_cgpVaf))]
#     #reset the allele counts to NA before starting a new sample. Not strictly necessary but doing it for my sanity/piece of mind
#     allele_count$Count_A=NA
#     allele_count$Count_C=NA
#     allele_count$Count_G=NA
#     allele_count$Count_T=NA
#     
#     allele_count$Count_A = mut_temp[,grepl('FAZ',colnames(mut_temp))] + mut_temp[,grepl('RAZ',colnames(mut_temp))]
#     allele_count$Count_C = mut_temp[,grepl('FCZ',colnames(mut_temp))] + mut_temp[,grepl('RCZ',colnames(mut_temp))]
#     allele_count$Count_G = mut_temp[,grepl('FGZ',colnames(mut_temp))] + mut_temp[,grepl('RGZ',colnames(mut_temp))]
#     allele_count$Count_T = mut_temp[,grepl('FTZ',colnames(mut_temp))] + mut_temp[,grepl('RTZ',colnames(mut_temp))]
#     
#     write.table(allele_count, file = file.path(allelecount_dir,paste0(sample, "_allelecounts.txt")),
#                 sep = "\t", row.names = F, quote = F, col.names = T)
#   }
# }





## STEP3: Run Tim's shearwater-like filter to remove false mutations ##

## Input for shearwater filter
# Vector of all sample names
samples_ID=c('PD60301a','PD61846a')

mutations = runShearWaterLike(samples_ID=c('PD60301a','PD61846a'),outDir=outDir,patient='L038')

# # Vector of normal panel (only blood from other project [3110])
# normal_samples=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt", header=T) 
# 
# # "Bed" file of all mutations to be considered (across all patients)
# # Format: Chr Ref Pos Alt
# variants = read.delim('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/1_potential_somaticVar_forShearwater.bed',header = F)
# #muts=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/cgpVAF/01_Input/all_mutations_all_pat.bed")
# variants=variants[!duplicated(variants[,1:2]),]
# coords=paste(variants$V1,variants$V2,sep="_")
# 
# # Read in data from AlleleCounter/cgpVAF
# all_counts = array(0,dim=c(length(samples_ID),length(coords),4),
#                    dimnames=list(samples_ID,coords,c("A","C","G","T")))
# print(length(samples_ID))
# for (k in 1:length(samples_ID)){
#   #Read in allele counts per sample
#   if(file.exists(paste0(allelecount_dir, samples_ID[k],"_allelecounts.txt"))){
#     print(samples_ID[k])
#     data=read.table(paste0(allelecount_dir, samples_ID[k],"_allelecounts.txt"),comment.char = '',header=T)
#     rownames(data)=paste(data$Chrom,data$Pos,sep="_")
#     all_counts[k,,]=as.matrix(data[coords,3:6])
#   }
# }
# 
# # Read in data from AlleleCounter/cgpVAF for all samples in the normal panel
# norm_all_counts = array(0,dim=c(nrow(normal_samples),length(coords),4),
#                         dimnames=list(normal_samples$Sample_ID,coords,c("A","C","G","T")))
# 
# for (k in 1:nrow(normal_samples)){
#   #Read in allele counts per sample
#   if(file.exists(paste0(allelecount_dir, normal_samples$Sample_ID[k],"_allelecounts.txt"))){
#     print(normal_samples$Sample_ID[k])
#     data=read.table(paste0(allelecount_dir, normal_samples$Sample_ID[k],"_allelecounts.txt"),comment.char = '',header=T)
#     rownames(data)=paste(data$Chrom,data$Pos,sep="_")
#     norm_all_counts[k,,]=as.matrix(data[coords,3:6])
#   }
# }
# 
# 
# 
# #Run shearwater filter
# setwd("~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/")
# source('~/lustre_mt22/generalScripts/shearwater/shearwater_flt_2019_WT_AW.R')
# #mutations=read.table(paste0("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/germline_flt/",pat_single,"_mut_germline_flt.txt"), header=T)
# #pval_mat=shearwater_probability(patient=pat_single)
# #patient = 'L038'
# colnames(variants) = c('Chr','Pos','Ref','Alt')
# pval_mat = shearwater_probability(patient, save=paste0(patient,"_shearwater_pval_mat.txt"),allVar=variants,
#                                 all_counts=all_counts,norm_all_counts=norm_all_counts,case_samples=samples_ID, rho=10^-3)
# #shearwater_probability(patient=patient,save=paste0(patient,"_shearwater_pval_mat.txt"))
# 
# # Separate files for correcting for multiple testing
# # One file with the variants where MTR is >1 (at least 2) in more than one sample
# #NOTE! I'm doing this because I have several thousands of private mutations in one sample which is skewing the
# # qvalue for shared mutations (ie inflating them) but you probably don't need to do that
# # NV=mutations[,c(paste0(samples_pat, "_MTR"))]
# # mut_several=pval_mat[which(rowSums(NV>1)>1),]
# # mut_private=pval_mat[which(rowSums(NV>1)<=1),]
# # qval_mat_several=apply(mut_several,2,function(x) p.adjust(x,method="BH",n = length(mut_several)))
# # qval_mat_private=apply(mut_private,2,function(x) p.adjust(x,method="BH",n = length(mut_private)))
# # qval_mat=rbind(qval_mat_private, qval_mat_several)
# 
# 
# qval_mat=apply(pval_mat,2,function(x) p.adjust(x,method="BH",n = length(pval_mat)))
# # Note that if MTR=DEP, the p-value will be so low that it is listed as NA(n). Change these values to 0
# qval_mat[which(qval_mat=="NaN")]=0
# row.names(qval_mat)=sub('.*chr', '', row.names(qval_mat))
# 
# # Add the qval info to the mutations object (make sure rows are in the same order first)
# mutations = variants
# row.names(mutations)=paste0(sub('.*chr', '', mutations$Chr), "_", mutations$Pos, "_", mutations$Ref,"_", mutations$Alt)
# table(row.names(mutations) %in% row.names(qval_mat))
# reorder_idx=match(row.names(mutations), row.names(qval_mat))
# qval_mat=qval_mat[reorder_idx, ]
# mutations=cbind(mutations, qval_mat)

# Tim uses <0.001 as the mutations that passed the shearwater filter and I use <0.005
table(mutations$PD60301a < 0.001)
table(mutations$PD60301a < 0.005)
mutations$PD60301a_shearwaterPASS = ifelse(mutations$PD60301a < 0.001,T,F)
mutations$PD61846a_shearwaterPASS = ifelse(mutations$PD61846a < 0.001,T,F)
mutations$snv_ID = paste0(mutations$Chr,':',mutations$Pos,'_',mutations$Ref,'/',mutations$Alt)
table(somaticVar$snv_ID %in% mutations$snv_ID)


##  STEP4: Add shearwater-like filter results to main object ##

cgpVaf_QC$snv_type = as.character(cgpVaf_QC$snv_type)
cgpVaf_QC$snv_type[cgpVaf_QC$sample == 'PD60301a' & 
                     cgpVaf_QC$snv_type != 'removed' &  
                     cgpVaf_QC$snv_ID %in% mutations$snv_ID[mutations$PD60301a_shearwaterPASS == F]] = 'shearwater_Failed'
cgpVaf_QC$snv_type[cgpVaf_QC$sample == 'PD61846a' & 
                     cgpVaf_QC$snv_type != 'removed' &  
                     cgpVaf_QC$snv_ID %in% mutations$snv_ID[mutations$PD61846a_shearwaterPASS == F]] = 'shearwater_Failed'


cgpVaf_QC$snv_type = factor(cgpVaf_QC$snv_type,c('germline','likely_somatic',"likely_germline_D", "likely_germline_TP1",'CN_region','removed','shearwater_Failed'))
cgpVaf_QC = cgpVaf_QC[order(cgpVaf_QC$snv_type),]
cgpVaf_QC$Chrom = factor(cgpVaf_QC$Chrom,paste0('chr',c(1:22,'X','Y')))
ggplot(cgpVaf_QC,aes(Pos,sample_VAF))+
  geom_point(size=0.1,aes(col=snv_type,size=sample_DEP))+
  facet_grid(sample ~ Chrom,scales = 'free_x')+
  scale_size_manual(values = c(0.01,3))+
  scale_color_manual(values = c(grey(0.8),'black',col25[1:3],grey(0.9),'purple'))+
  geom_hline(yintercept = 0.5,lty=2,lwd=0.4,col=grey(0.4))+
  geom_hline(yintercept = 0.13,lty=2,lwd=0.4,col='black')+
  #geom_vline(xintercept = c(0.13),lty=2,lwd=0.7,col=c('red'))+
  #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
  theme_classic(base_size = 12)+
  theme(panel.border = element_rect(fill=F),
        axis.line = element_blank(),axis.text.x = element_blank(),axis.ticks.length.x = unit(0,'cm')) +
  xlab('Genomic position') + ylab('VAF')


##--- Upset plot
table(cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD60301a'] %in% cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD61846a'])

varList = list('D_germline' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type %in% c('germline','likely_germline_D') & cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$sample_VAF >0],
               'TP1_germline' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type %in% c('germline','likely_germline_TP1') & cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$sample_VAF >0],
               'D_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'likely_somatic' & cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$sample_VAF > 0],
               'TP1_somatic' = cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type == 'likely_somatic' & cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$sample_VAF > 0])
table(varList[[1]] %in% varList[[3]])

library(UpSetR)
upset(fromList(varList),nsets = 20,text.scale = 2)




##------ Extract somatic variants after shearwater filter
somaticVar = cgpVaf_QC[cgpVaf_QC$snv_type == 'likely_somatic',]
somaticVar$somaticVar_type = ifelse(somaticVar$snv_ID %in% intersect(somaticVar$snv_ID[somaticVar$sample != 'PD60301a' & somaticVar$sample_VAF > 0],
                                                                     somaticVar$snv_ID[somaticVar$sample == 'PD60301a' & somaticVar$sample_VAF > 0]),'shared',
                                    ifelse(somaticVar$sample == 'PD60301a' & somaticVar$sample_VAF > 0,'unique_D',
                                           ifelse(somaticVar$sample == 'PD61846a' & somaticVar$sample_VAF > 0,'unique_TP1','others')))

table(somaticVar$somaticVar_type,somaticVar$sample)

# What are the reason for the unique variants to be unique
cgpVaf_QC_tp1 = cgpVaf_QC[cgpVaf_QC$sample != 'PD60301a',]
cgpVaf_QC_tp1$snv_type = as.character(cgpVaf_QC_tp1$snv_type)

cgpVaf_QC_d = cgpVaf_QC[cgpVaf_QC$sample == 'PD60301a',]
somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_D'] = paste0(somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_D'],':',
                                                                              cgpVaf_QC_tp1$snv_type[match(somaticVar$snv_ID[somaticVar$somaticVar_type == 'unique_D'],cgpVaf_QC_tp1$snv_ID)])
somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_TP1'] = paste0(somaticVar$somaticVar_type[somaticVar$somaticVar_type == 'unique_TP1'],':',
                                                                              cgpVaf_QC_d$snv_type[match(somaticVar$snv_ID[somaticVar$somaticVar_type == 'unique_TP1'],cgpVaf_QC_d$snv_ID)])


table(somaticVar$somaticVar_type,somaticVar$sample)

# removedVar = unique(c(cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$germlineFilter == 'removed'],
#                       cgpVaf_QC$snv_ID[cgpVaf_QC$snv_type != 'somatic']))
# varList = list('D_somatic' = somaticVar$snv_ID[somaticVar$sample == 'PD60301a' & somaticVar$sample_VAF > 0],
#                'TP1_somatic' = somaticVar$snv_ID[somaticVar$sample == 'PD61846a' & somaticVar$sample_VAF > 0])
# 
# 
# 
# varList = list('D_somatic' = somaticVar$snv_ID[somaticVar$sample == 'PD60301a' & somaticVar$sample_VAF > 0],
#                'TP1_somatic' = somaticVar$snv_ID[somaticVar$sample == 'PD61846a' & somaticVar$sample_VAF > 0],
#                'D_removed' = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD60301a' & cgpVaf_QC_og$germlineFilter == 'removed'],
#                'TP1_removed' = cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD61846a' & cgpVaf_QC_og$germlineFilter == 'removed'])
# 
# library(UpSetR)
# upset(fromList(varList),nsets = 20,text.scale = 2)


ggplot(somaticVar,aes(sample_VAF))+
  geom_density()+
  facet_grid(somaticVar_type ~ timepoint)+
  geom_vline(xintercept = 0.5,lty=2,lwd=0.7,col='grey')+
  #geom_vline(xintercept = c(expected_vaf_d),lty=2,lwd=0.7,col=c('red'))+
  #geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
  theme_classic(base_size = 14)+theme(axis.line = element_line())

## Add variant annotation info from caveman
somaticVar$Gene = var4cgpVaf$Gene[match(somaticVar$snv_ID,var4cgpVaf$snv_ID)]
somaticVar$Impact = var4cgpVaf$Impact[match(somaticVar$snv_ID,var4cgpVaf$snv_ID)]
somaticVar$AAchange = var4cgpVaf$AAchange[match(somaticVar$snv_ID,var4cgpVaf$snv_ID)]

cgpVaf_withReadFilters_output$Gene = var4cgpVaf$Gene[match(cgpVaf_withReadFilters_output$snv_ID,var4cgpVaf$snv_ID)]
cgpVaf_withReadFilters_output$Impact = var4cgpVaf$Impact[match(cgpVaf_withReadFilters_output$snv_ID,var4cgpVaf$snv_ID)]
cgpVaf_withReadFilters_output$AAchange = var4cgpVaf$AAchange[match(cgpVaf_withReadFilters_output$snv_ID,var4cgpVaf$snv_ID)]

## pivot somaticVar wider to plot vaf of D vs TP1
df = cgpVaf_withReadFilters_output[cgpVaf_withReadFilters_output$snv_ID %in% somaticVar$snv_ID,]
df$somaticVar_type = somaticVar[somaticVar$somaticVar_type != 'others',]$somaticVar_type[match(df$snv_ID,somaticVar$snv_ID[somaticVar$somaticVar_type != 'others'])]

df = pivot_wider(df,id_cols = c('snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect','Impact','AAchange','somaticVar_type'),
                 names_from = 'sample',values_from = 'sample_VAF')

cgpVaf_QC_og$reasonForFail[cgpVaf_QC_og$reasonForFail == 'PASS' & cgpVaf_QC_og$pindelFilter == 'overlapped_wINDELs'] = 'overlapped_wINDELs'
cgpVaf_QC_og$reasonForFail[cgpVaf_QC_og$reasonForFail == 'PASS' & cgpVaf_QC_og$germlineFilter == 'germline'] = 'germline'
df$snv_removed_PD60301a = cgpVaf_QC_og[cgpVaf_QC_og$sample == 'PD60301a',]$reasonForFail[match(df$snv_ID,cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD60301a'])]
df$snv_removed_PD61846a = cgpVaf_QC_og[cgpVaf_QC_og$sample == 'PD61846a',]$reasonForFail[match(df$snv_ID,cgpVaf_QC_og$snv_ID[cgpVaf_QC_og$sample == 'PD61846a'])]


# df$somaticVar_type[df$somaticVar_type=='unique_D/shearwater_Failed' & df$PD61846a == 0] = 'unique_D'
# df$somaticVar_type[df$somaticVar_type=='unique_TP1/shearwater_Failed' & df$PD60301a == 0] = 'unique_TP1'
# df$somaticVar_type = gsub('/',':',df$somaticVar_type)
#table(df$snv_removed,df$somaticVar_type)
#df$somaticVar_type_2 = ifelse(df$snv_removed == 'removed', 'removed',df$somaticVar_type)
ggplot(df,aes(PD60301a,PD61846a,col=somaticVar_type))+
  geom_point(size=0.7)+
  scale_color_manual(values = c(grey(0.8),grey(0.2),col25))+
  theme_classic(base_size = 14) + ggtitle('cgpVAF output of somatic variants')


write.csv(df, '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_D.vs.TP1_finalSomaticVariants.csv')









##--------------------------------------------------##
##  Generate jbrowse images of these mutations    ####
##--------------------------------------------------##
df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_D.vs.TP1_finalSomaticVariants.csv')

## 1. Define JBROWSE URL
# L038 samples
jbrowse_url_l038 = '# JBROWSE https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F3030&loc=1%3A1..247467938&tracks=PD60301a_bwa%2CPD61846a_bwa'
# normal panel samples from 3010 project
normal_samples=read.table("/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt", header=T) 
normalBlood_sampleIDs = normal_samples$Sample_ID
jbrowse_url_normalblood = paste0('# JBROWSE https://cgp-jbrowse.internal.sanger.ac.uk/Homo_sapiens/GRCh38/JBrowse/?data=auto%2F3010&loc=1%3A1..247467938&tracks=',
                                 paste0(normalBlood_sampleIDs,collapse = '_bwa%2C'))
  

## 2. Prepare bed file containing location of these mutations, with +/ 50 bases around the mutation position
jbrowse_outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/jbrowse/L038/'
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
    for(d.sub in file.path(d,c('L038','normalPanel_3010'))){
      dir.create(d.sub,recursive = T)
    }
  }    
  
  for(sampleGroup in c('L038','normalPanel_3010')){
    ## Write bed file containing mutation of this category
    bed = var[,c('Chrom','start','end')]
    bedFile = file.path(d,sampleGroup,'L038_mutation.bed')
    if(sampleGroup == 'L038'){
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
df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_D.vs.TP1_finalSomaticVariants.csv')
df_long = pivot_longer(df,cols = c('PD60301a', 'PD61846a'),names_to = 'sampleID',values_to = 'VAF')
ggplot(df_long,aes(VAF))+
  geom_density()+
  facet_grid(.~sampleID)+
  theme_classic(base_size = 14)+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank())+xlab('Variant allele frequency')

## subset for mutations that has been called
cgpVaf_QC_og = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038_cgpVaf_filtersApplied.csv')
cgpVaf_QC = cgpVaf_QC[cgpVaf_QC$snv_ID %in% df$snv_ID,]
res_D = binom_mix(x=cgpVaf_QC$sample_MTR[cgpVaf_QC$sample == 'PD60301a'& cgpVaf_QC$sample_DEP.noFilter > 0],size = cgpVaf_QC$sample_DEP.noFilter[cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$sample_DEP.noFilter > 0],mode = 'Full')
data_D = data.frame(NV = cgpVaf_QC$sample_MTR[cgpVaf_QC$sample == 'PD60301a'& cgpVaf_QC$sample_DEP.noFilter > 0],
                  NR = cgpVaf_QC$sample_DEP.noFilter[cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$sample_DEP.noFilter > 0],
                  snvID = cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD60301a' & cgpVaf_QC$sample_DEP.noFilter > 0],
                  cluster = res_D[['Which_cluster']]
                  )
data_D$cat = df$somaticVar_type[match(data_D$snvID,df$snv_ID)]
plot_binomMixModel(data_D,res_D,sampleID=NULL,mode='Full')

res_TP1 = binom_mix(x=cgpVaf_QC$sample_MTR[cgpVaf_QC$sample == 'PD61846a'& cgpVaf_QC$sample_DEP.noFilter > 0],
                    size = cgpVaf_QC$sample_DEP.noFilter[cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$sample_DEP.noFilter > 0],mode = 'Full',nrange = 1:2)
data_TP1 = data.frame(NV = cgpVaf_QC$sample_MTR[cgpVaf_QC$sample == 'PD61846a'& cgpVaf_QC$sample_DEP.noFilter > 0],
                  NR = cgpVaf_QC$sample_DEP.noFilter[cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$sample_DEP.noFilter > 0],
                  snvID = cgpVaf_QC$snv_ID[cgpVaf_QC$sample == 'PD61846a' & cgpVaf_QC$sample_DEP.noFilter > 0],
                  cluster = res_TP1[['Which_cluster']])
data_TP1$cat = df$somaticVar_type[match(data_TP1$snvID,df$snv_ID)]
table(data_TP1$cat,data_TP1$cluster)
plot_binomMixModel(data_TP1,res_TP1,sampleID=NULL,mode='Full')


# 
# 
# 
# 
# 
# 
# 
# View(df[df$somaticVar_type_2 == 'removed',])
# 
# View(somaticVar[somaticVar$somaticVar_type != 'shared',])
# View(d_cgpVaf[d_cgpVaf$snv_ID %in% somaticVar$snv_ID[somaticVar$somaticVar_type != 'shared'],])
# View(var4cgpVaf[var4cgpVaf$snv_ID %in% unique(somaticVar$snv_ID[somaticVar$somaticVar_type != 'shared']),])
# 


##---------------------------------------------------##
## Finalising number of mutation in each sample    ####
##---------------------------------------------------##
## Go through jbrowse images and filter out "bad" variants

## After manual inspection of jbrowse images, move bad variants to a seperate folder





















#
#
#
#
# # Binomial mix model to decompose clonal structure
#
# somaticVar$NV = somaticVar$refCnt + somaticVar$altCnt
# somaticVar$NR = somaticVar$altCnt
# Mode = 'Full'
# res = binom_mix(x = somaticVar$NV,size = somaticVar$NR,mode=Mode)
# plot_binomMixModel(d_caveman_somatic,res,sampleID = 'PD60301a')
#
#
#
#
#
#
#
# cgpVaf_withReadFilters_output = rbind(d_cgpVaf,tp1_cgpVaf)
# table(d_cgpVaf$snv_ID %in% somaticVar$snv_ID)
# d_cgpVaf = d_cgpVaf[d_cgpVaf$snv_ID %in% somaticVar$snv_ID,]
#
#
# table(tp1_cgpVaf$snv_ID %in% somaticVar$snv_ID)
# tp1_cgpVaf = tp1_cgpVaf[tp1_cgpVaf$snv_ID %in% somaticVar$snv_ID,]
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # Find EZH2 variant to determine tumour purity / expected vaf
# expected_vaf_d = d_caveman$VAF[!is.na(d_caveman$Gene) & d_caveman$Gene == 'EZH2']
#
# # Binomial test to filter variants around the expected vaf
# pval = sapply(seq(1:nrow(d_caveman_somatic)),function(s){
#     binom.test(x=d_caveman_somatic$altCnt[s],n=d_caveman_somatic$DP[s], p=expected_vaf_d,alt='greater')$p.value
# })
# qval = p.adjust(pval,method="BH")
# d_caveman_somatic$isEZH2clone = qval>pCut
#
#
#
# ##-----------------------------------------------##
# ##      Import Caveman output - TP1            ####
# ##-----------------------------------------------##
# tp1_caveman = caveman2(sample='PD61846a',projectid = 3030)
# tp1_caveman$timepoint = 'TP1'
# tp1_caveman = tp1_caveman[tp1_caveman$Flag == T,]
#
#
# # Run binomial filter
# # data is a data.frame where: rownames = snv_ID (chr1:pos_REF/ALT), required columns are refCnt, altCnt, DP (total depth)
# tp1_caveman$snv_ID = paste0(tp1_caveman$Chr,':',tp1_caveman$Pos,'_',tp1_caveman$Ref,'/',tp1_caveman$Alt)
# tp1_caveman$refCnt = tp1_caveman$NR - tp1_caveman$NV
# tp1_caveman$altCnt = tp1_caveman$NV
# tp1_caveman$DP = tp1_caveman$NR
#
# tp1_caveman = snvBinomTest(gender='female',data=tp1_caveman,genotype='diploid',pCut=0.05)
# table(tp1_caveman$isGermline)
#
#
# # Binomial mix model to decompose clonal structure
# tp1_caveman_somatic = tp1_caveman[tp1_caveman$isGermline == F,]
# res = binom_mix(tp1_caveman_somatic$NV,tp1_caveman_somatic$NR,mode=Mode)
# plot_binomMixModel(tp1_caveman_somatic,res,sampleID = 'PD61846a')
#
# # Find EZH2 variant to determine tumour purity / expected vaf
# expected_vaf_tp1 = tp1_caveman$VAF[!is.na(tp1_caveman$Gene) & tp1_caveman$Gene == 'EZH2']
#
# # Binomial test to filter variants around the expected vaf
# pval = sapply(seq(1:nrow(tp1_caveman_somatic)),function(s){
#   binom.test(x=tp1_caveman_somatic$altCnt[s],n=tp1_caveman_somatic$DP[s], p=expected_vaf_tp1,alt='greater')$p.value
# })
# qval = p.adjust(pval,method="BH")
# tp1_caveman_somatic$isEZH2clone = qval>pCut
#
#
#
#
# ##-------------------------------------------------##
# ##      Extract only Somatic variants            ####
# ##-------------------------------------------------##
# expected_vaf = 0.5*22/100
#
# somaticVar = rbind(d_caveman,tp1_caveman)
# varToKeep = c(tp1_caveman$snv_ID[tp1_caveman$isGermline == F & !tp1_caveman$snv_ID %in% d_caveman$snv_ID[d_caveman$isGermline == T]],
#               d_caveman$snv_ID[d_caveman$isGermline == F & !d_caveman$snv_ID %in% tp1_caveman$snv_ID[tp1_caveman$isGermline == T]])
# somaticVar$snv_type = ifelse(somaticVar$snv_ID %in% varToKeep,'somatic',
#                              ifelse(somaticVar$snv_ID %in% tp1_caveman$snv_ID[tp1_caveman$isGermline == F & tp1_caveman$snv_ID %in% d_caveman$snv_ID[d_caveman$isGermline == T]],'likely_germline_D',
#                                     ifelse(somaticVar$snv_ID %in% d_caveman$snv_ID[d_caveman$isGermline == F & d_caveman$snv_ID %in% tp1_caveman$snv_ID[tp1_caveman$isGermline == T]],'likely_germline_TP1','germline')))
#
# ggplot(somaticVar,aes(VAF))+
#   geom_density()+
#   facet_grid(snv_type ~ timepoint)+
#   geom_vline(xintercept = 0.5,lty=2,lwd=0.7,col='grey')+
#   geom_vline(xintercept = c(expected_vaf_d),lty=2,lwd=0.7,col=c('red'))+
#   geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
#   theme_classic(base_size = 14)+theme(axis.line = element_line())
#
# somaticVar = somaticVar[somaticVar$snv_ID %in% varToKeep,]
# somaticVar$shared_var = ifelse(somaticVar$snv_ID %in% intersect(d_caveman$snv_ID,tp1_caveman$snv_ID),'shared',
#                                ifelse(somaticVar$snv_ID %in% d_caveman$snv_ID & ! somaticVar$snv_ID %in% tp1_caveman$snv_ID,'unique_D','unique_TP1'))
# table(somaticVar$isGermline)
#
#
# ggplot(somaticVar,aes(VAF))+
#   geom_density()+
#   facet_grid(shared_var ~ timepoint)+
#   geom_vline(xintercept = 0.5,lty=2,lwd=0.7,col='grey')+
#   geom_vline(xintercept = c(expected_vaf_d),lty=2,lwd=0.7,col=c('red'))+
#   geom_vline(xintercept = c(expected_vaf_tp1),lty=2,lwd=0.7,col=c('blue'))+
#   theme_classic(base_size = 14)+theme(axis.line = element_line()) + ggtitle('Somatic Variants')
#
#
# #df = somaticVar[somaticVar$shared_var == 'shared',]
# df = pivot_wider(df,id_cols = 'snv_ID',names_from = 'timepoint',values_from = 'VAF')
# ggplot(df,aes(Diagnostic,TP1))+
#   geom_point(size=0.5)+theme_classic(base_size = 13)+geom_abline()+
#   ggtitle('Shared somatic variants')
#
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
# bed = data.frame(snv_ID = unique(somaticVar$snv_ID))
# bed = merge(bed,somaticVar[,c('Chr','Pos','Ref','Alt','snv_ID')],by='snv_ID',all.x=T)
# bed = bed[!duplicated(bed),]
# dim(bed)
# n_distinct(somaticVar$snv_ID)
# #rownames(bed) = bed$snv_ID
# bed = bed[order(bed$Chr),]
# write_delim(bed[,colnames(bed) != 'snv_ID'],'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L038_somaticVariants.bed',col_names = F,delim = '\t')
# write.csv(c('PD60301a','PD61846a'),'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L038_samples.csv',row.names = F,quote = F)
#
# ##--- prepare for cgpVAF
# # cd ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis
# # mkdir cgpVaf_input/L038
# # mkdir cgpVaf_output/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD60301a/PD60301a.sample.dupmarked.bam cgpVaf_input/L038/
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD60301a/PD60301a.sample.dupmarked.bam.bai cgpVaf_input/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD60301a/PD60301a.caveman_c.annot.vcf.gz cgpVaf_input/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD61846a/PD61846a.sample.dupmarked.bam cgpVaf_input/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD61846a/PD61846a.sample.dupmarked.bam.bai cgpVaf_input/L038
# # cp -s /lustre/scratch124/casm/team78pipelines/nst_links/live/3030/PD61846a/PD61846a.caveman_c.annot.vcf.gz cgpVaf_input/L038
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam cgpVaf_input/L038/
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam.bai cgpVaf_input/L038/
# # cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam cgpVaf_input/L038/
# ##--- Run cgpVAF
# # bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/x10_MLDS_cgpVaf.sh L038 PD60301a
# # bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/x10_MLDS_cgpVaf.sh L038 PD61846a
#
#
# ##---- Import cgpVAF output
# d_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038/PDv38is_wgs_PD60301a_snp_vaf.tsv',sep = '\t',skip = 36)
# d_cgpVaf$timepoint='Diagnostic'
# d_cgpVaf$snv_ID = paste0(d_cgpVaf$Chrom,':',d_cgpVaf$Pos,'_',d_cgpVaf$Ref,'/',d_cgpVaf$Alt)
# table(d_cgpVaf$snv_ID %in% somaticVar$snv_ID)
# d_cgpVaf = d_cgpVaf[d_cgpVaf$snv_ID %in% somaticVar$snv_ID,]
#
#
# tp1_cgpVaf = read.delim('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L038/PDv38is_wgs_PD61846a_snp_vaf.tsv',sep = '\t',skip = 36)
# tp1_cgpVaf$timepoint='TP1'
# tp1_cgpVaf$snv_ID = paste0(tp1_cgpVaf$Chrom,':',tp1_cgpVaf$Pos,'_',tp1_cgpVaf$Ref,'/',tp1_cgpVaf$Alt)
# table(tp1_cgpVaf$snv_ID %in% somaticVar$snv_ID)
# tp1_cgpVaf = tp1_cgpVaf[tp1_cgpVaf$snv_ID %in% somaticVar$snv_ID,]
#
# df = data.frame(snv_ID = unique(somaticVar$snv_ID))
# df$shared_var = somaticVar$shared_var[match(df$snv_ID,somaticVar$snv_ID)]
# df$vaf_D = d_cgpVaf$PD60301a_VAF[match(df$snv_ID,d_cgpVaf$snv_ID)]
# df$vaf_TP1 = tp1_cgpVaf$PD61846a_VAF[match(df$snv_ID,tp1_cgpVaf$snv_ID)]
#
# df = df[df$snv_ID %in% c(d_caveman$snv_ID[d_caveman$Flag == T],
#                          tp1_caveman$snv_ID[tp1_caveman$Flag == T]) &
#           !df$snv_ID %in% df$snv_ID[df$vaf_D == 0 & df$vaf_TP1 == 0],]
# ggplot(df,aes(vaf_D,vaf_TP1,col=shared_var))+
#   geom_point(size=0.2)+
#   scale_color_manual(values = col25)+
#   theme_classic() + ggtitle('cgpVAF output of somatic variants')
#
#
# ## Manual check on IGV for these mutations
# df = merge(df,somaticVar[,c('Chr','Pos','Ref','Alt','Gene','Impact','AAchange','ID','snv_ID')],by='snv_ID',all.x=T)
# df = df[!duplicated(df),]
# dim(df)
#
# df$caveman_vaf_D = d_caveman$VAF[match(df$snv_ID,d_caveman$snv_ID)]
# df$caveman_vaf_TP1 = tp1_caveman$VAF[match(df$snv_ID,tp1_caveman$snv_ID)]
#
# df$cgpVaf_DP_D = d_cgpVaf$PD60301a_DEP[match(df$snv_ID,d_cgpVaf$snv_ID)]
# df$cgpVaf_DP_TP1 = tp1_cgpVaf$PD61846a_DEP[match(df$snv_ID,tp1_cgpVaf$snv_ID)]
#
#
# View(df[df$snv_ID %in% c(tp1_caveman_somatic$snv_ID[!is.na(tp1_caveman_somatic$Impact) & tp1_caveman_somatic$Impact == 'missense'],
#                          d_caveman_somatic$snv_ID[!is.na(d_caveman_somatic$Impact) & d_caveman_somatic$Impact == 'missense']),c('snv_ID','shared_var','vaf_D','vaf_TP1','Gene','Impact','AAchange')])
#
# ## Filter for variants
#
# df.sub = df[df$cgpVaf_DP_D >= 10 | df$cgpVaf_DP_TP1 >= 10,]
# dim(df.sub)
# View(df.sub[df.sub$vaf_D == 0 | df.sub$vaf_TP1 == 0,])
