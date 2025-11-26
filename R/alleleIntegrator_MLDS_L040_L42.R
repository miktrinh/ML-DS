# Run alleleIntegrator on tumour samples

# sudo apt-get install bcftools
# sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
# sudo chmod +x /usr/local/bin/alleleCounter



outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/'
if(!dir.exists(outDir)){
  dir.create(outDir)
}
setwd(outDir)

#############
# Libraries #
#############
# Load libraries
library(alleleIntegrator)
library(GenomicFeatures)
library(tidyverse)

# library(RColorBrewer)
# library(circlize)
# library(ComplexHeatmap)
# library(Seurat)
# library(readxl)
# source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
# source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
# source("~/lustre_mt22/generalScripts/utils/misc.R")

#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
minSegLen=1e6
subCl.minSegLen=5e6
skipIfExists = T
normREF = T
mainDir = outDir

# Import chromInfo
#chromInfo_fp = '../chrom_abspos_kb.txt'
#chromInfo = read_delim(chromInfo_fp,delim = '\t',col_names = T)

refGenome = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa'
refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
liftChain = '~/lustre_mt22/hg19ToHg38_noChr.over.chain'
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
#txdb = makeTxDbFromGFF(gtf)
#gns = genes(txdb)
nParallel=48




## Run alleleIntegrator on MLDS samples


############
# PD60303a - L040 ####
############

#########################
# Sample specific params
# tumourDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD60303a/PD60303a.sample.dupmarked.bam'
# patientDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD61851a/PD61851a.sample.dupmarked.bam'
# bams10X = c('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234197_GRCh38-2020-A/possorted_genome_bam.bam',
#             '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234198_GRCh38-2020-A/possorted_genome_bam.bam',
#             '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234199_GRCh38-2020-A/possorted_genome_bam.bam',
#             '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415928_GRCh38-2020-A/possorted_genome_bam.bam',
#             '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415929_GRCh38-2020-A/possorted_genome_bam.bam',
#             '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415930_GRCh38-2020-A/possorted_genome_bam.bam',
#             '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47089_SB_Leuk13645524_GRCh38-2020-A/possorted_genome_bam.bam')
# bams10X = setNames(bams10X,gsub('.*SB_|_GRCh38-2020-A','',basename(dirname(bams10X))))
# PDID = 'PD60303a'
# 
# 
# 
# 
# 
# 
# source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
# ####------------------ Generate Battenberg CN summary file ----------------####
# # Battenberg .summary.csv file - only summarize Major Clone CNV, does not included CN states of minor clones
# btb.fp = '/nfs/cancer_ref01/nst_links/live/3030/PD60303a/PD60303a.battenberg.summary.csv'
# #----- Processing Battenberg data -------#
# segs = annotateBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID,tgtChrs=tgtChrs,removeBalancedSegs=T,longFormat = F,method = 'allelicRatio')
# 
# # Remove Chr21 and sex chr
# segs = segs[!segs$Chr %in% c(21,23),]
# segs$Chr = paste0('chr',segs$Chr)
# segs = GRanges(segs$Chr,IRanges(segs$Start,segs$Stop),
#               totCN=segs$tumTot, idx=segs$idx, chr=segs$Chr,
#               patNum=segs$patNum, matNum=segs$matNum,
#               tumFrac = segs$tumFrac,clonalType=segs$type)
# 
# message(sprintf('[PDID %s] final number of segments %d',PDID,nrow(mcols(segs))))
# if(nrow(mcols(segs)) == 0) {next()}
# 
# # #Define CN segments roughly
# # altChrs = c('chr6')
# # segs = GRanges(altChrs,IRanges(c(rep(1,1)),c(1e9,50e6,1e9)))
# #
# # segs$matNum = c(1,1,2)
# # segs$patNum = c(0,0,1)
# # segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
# # segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
# # names(segs) = altChrs
# 
# 
# 
# 
# if(skipIfExists & file.exists(file.path(outDir,paste0(PDID,'_phCnts.RDS')))){
#   phCnts = readRDS(file.path(outDir,paste0(PDID,'_phCnts.RDS')))
# 
# }else{
#   ######################
#   # Call and phase SNPs
# 
#   hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
# 
#   #Expectation is that we'll find ~ 3 million of them
#   message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
#   #Found 1,918,709 heterozygous SNPs
# 
#   ## Plot B-allele frequency of hSNPs on chr5 and chr7
#   # df = hSNPs[seqnames(hSNPs) %in% c('chr5','chr7','chr21','chr22')]
#   # df = as.data.frame(mcols(df))
#   # df$pos = gsub('.*:|_.*$','',rownames(df))
#   # df$chr = gsub(':.*$','',rownames(df))
#   #
#   # p = ggplot(df,aes(pos,BAF))+
#   #   geom_point(size=0.01,alpha=0.2)+
#   #   geom_hline(yintercept = 0.5)+
#   #   facet_wrap(vars(chr))+
#   #   theme_classic() + theme(axis.text.x = element_blank())
#   #
#   # print(p)
# 
# 
# 
#   #Use tumour DNA to phase them.
#   # As we are using tumour as patientDNA (i.e. no matched normal), we will switch off EM, and set minPhasable to be very low.
#   # This is because if nPhased/nSNPs < minPhasable, it will set phSNPs$passSanity to be False, as most likely there's no real CN changes
#   phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tBAM = tumourDNA,useEM = T,
#                            outPath=file.path(outDir,paste0(PDID,'_tumour_countAtHetSNPs.tsv')),
#                            FDR = 0.1,minPhasable=0.01,
#                            nParallel=nParallel)
#   #Liftover to GRCh38
#   #phSNPs38 = changeGenomeVersion(phSNPs,liftChain)
#   phSNPs38 = phSNPs
#   #Annotate SNPs using GTF
#   phSNPs38 = annotateSNPs(phSNPs38,gtf)
# 
#   ########################
#   # Integrate with 10X.
#   #If the majority of the high coverage SNPs don't look heterozygous, something has gone wrong...
#   phCnts = getAllelicExpression(loci=phSNPs38,refGenome = refGenome10X,bams = bams10X,
#                                 outputs=file.path(outDir,paste0(PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
#                                 nParallel=nParallel)
# 
#   # Save this object so that we can process it faster next time!
#   saveRDS(phCnts,file.path(outDir,paste0(PDID,'_phCnts.RDS')))
# }
# 
# 
# ##############
# # Add cell type annotation
# srat = readRDS(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_231127.RDS'))
# srat@meta.data$cellID = rownames(srat@meta.data)
# srat$annot = srat$finalAnn_broad
# 
# normREF = T
# PDID = 'L040'
# 
# passCellIDs = rownames(srat@meta.data)
# clusterIDs = setNames(srat@meta.data$annot,srat@meta.data$cellID)
# normIDs = setNames(srat@meta.data$cellID[srat$annot %in% c('Mono_CD14')],srat@meta.data$annot[srat$annot %in% c('Mono_CD14')])
# #clusterIDs = NULL
# 
# 
# #If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
# aggToClust=FALSE
# #dropUninformative = F
# if(aggToClust & !is.null(clusterIDs)){
#   gCnts = aggregateByClusters(phCnts,clusterIDs)
#   gCnts = filterCells(gCnts,passCellIDs=levels(clusterIDs),normIDs=normIDs)
# }else if (!aggToClust & !is.null(clusterIDs)){
#   # Not aggToClust but using clusterInfo, including normCells being Leukocytes
#   gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=normIDs)
# }else if (!aggToClust & is.null(clusterIDs)){
#   # No annotation info is passed
#   gCnts = filterCells(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs)
# }
# 
# 
# ##################
# #Calibrate model
# #Specify the error rate
# gCnts$errRate = c('Exonic'=0.01,'Intronic'=0.05,'Intergenic'=0.15)[gCnts$regionType]
# 
# #Detect allele specific expression
# gCnts = calcASE(gCnts,priorKappa = 60)
# 
# 
# #Get over-dispersion
# od = calcOverDispersion(gCnts)
# 
# 
# ############
# # Inference
# #gCnts@metadata$segs$tumFrac = gCnts@metadata$segs$matNum/(gCnts@metadata$segs$matNum + gCnts@metadata$segs$patNum)
# pp = abbSegProb(gCnts,od)
# 
# 
# #############
# # Validation
# if(normREF){
#   pdf(file.path(outDir,paste0(PDID,'_rawAI_output.pdf')))
# }else{
#   pdf(file.path(outDir,paste0(PDID,'_rawAI_output_noNormREF.pdf')))
# }
# names(gCnts)
# gCnts$regionID = gsub(':.*$','',gCnts$regionID)
# dat = plotRawData(gCnts,segs = gCnts@metadata$segs,returnData=TRUE)
# p = plotPosteriorHeatmap(pp,'nLL')
# p = plotPosteriorHeatmap(pp,'nLL',split='',km = 4)
# print(p)
# dev.off()





















############
# PD60304a - L042 ####
############


outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L042'
if(!dir.exists(outDir)){
  dir.create(outDir)
}
setwd(outDir)



#########################
# Sample specific params
tumourDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD60304a/PD60304a.sample.dupmarked.bam'
patientDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD61850a/PD61850a.sample.dupmarked.bam'
bams10X = c('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234203_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234204_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234205_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415934_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415935_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415936_GRCh38-2020-A/possorted_genome_bam.bam')
bams10X = setNames(bams10X,gsub('.*SB_|_GRCh38-2020-A','',basename(dirname(bams10X))))
PDID = 'PD60304a'






source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
####------------------ Generate Battenberg CN summary file ----------------####
# Battenberg .summary.csv file - only summarize Major Clone CNV, does not included CN states of minor clones
btb.fp = '/nfs/cancer_ref01/nst_links/live/3030/PD60304a/PD60304a.battenberg.summary.csv'
#----- Processing Battenberg data -------#
segs = annotateBTB(btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID,tgtChrs=tgtChrs,removeBalancedSegs=T,longFormat = F,method = 'allelicRatio')

# Remove Chr21 and sex chr
segs = segs[!segs$Chr %in% c(21,23),]
segs$Chr = paste0('chr',segs$Chr)
segs = GRanges(segs$Chr,IRanges(segs$Start,segs$Stop),
               totCN=segs$tumTot, idx=segs$idx, chr=segs$Chr,
               patNum=segs$patNum, matNum=segs$matNum,
               tumFrac = segs$tumFrac,clonalType=segs$type)

message(sprintf('[PDID %s] final number of segments %d',PDID,nrow(mcols(segs))))
if(nrow(mcols(segs)) == 0) {next()}

# #Define CN segments roughly
# altChrs = c('chr6')
# segs = GRanges(altChrs,IRanges(c(rep(1,1)),c(1e9,50e6,1e9)))
#
# segs$matNum = c(1,1,2)
# segs$patNum = c(0,0,1)
# segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
# segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
# names(segs) = altChrs




if(skipIfExists & file.exists(file.path(outDir,paste0(PDID,'_phCnts.RDS')))){
  phCnts = readRDS(file.path(outDir,paste0(PDID,'_phCnts.RDS')))

}else{
  ######################
  # Call and phase SNPs

  hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)

  #Expectation is that we'll find ~ 3 million of them
  message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
  #Found 2,180,314 heterozygous SNPs

  # # Plot B-allele frequency of hSNPs on chr5 and chr7
  # df = hSNPs[seqnames(hSNPs) %in% c('chr1','chr7','chr21')]
  # df = as.data.frame(mcols(df))
  # df$pos = gsub('.*:|_.*$','',rownames(df))
  # df$chr = gsub(':.*$','',rownames(df))
  #
  # p = ggplot(df,aes(pos,BAF))+
  #   geom_point(size=0.01,alpha=0.2)+
  #   geom_hline(yintercept = 0.5)+
  #   facet_wrap(vars(chr))+
  #   theme_classic() + theme(axis.text.x = element_blank())
  #
  # print(p)
  #
  #
  # df.sub = df[df$chr == 'chr1',]
  # df.sub$type = ifelse(df.sub$pos > 32375592 & df.sub$pos < 81247540,'segs','others')
  # p = ggplot(df.sub,aes(pos,BAF))+
  #   geom_point(size=0.01,alpha=0.2)+
  #   geom_hline(yintercept = 0.5)+
  #   facet_wrap(vars(type))+
  #   theme_classic() + theme(axis.text.x = element_blank())
  #
  # print(p)



  #Use tumour DNA to phase them.
  # As we are using tumour as patientDNA (i.e. no matched normal), we will switch off EM, and set minPhasable to be very low.
  # This is because if nPhased/nSNPs < minPhasable, it will set phSNPs$passSanity to be False, as most likely there's no real CN changes
  phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tBAM = tumourDNA,useEM = T,
                           outPath=file.path(outDir,paste0(PDID,'_tumour_countAtHetSNPs.tsv')),nParallel=nParallel)
  #Liftover to GRCh38
  #phSNPs38 = changeGenomeVersion(phSNPs,liftChain)
  phSNPs38 = phSNPs
  #Annotate SNPs using GTF
  phSNPs38 = annotateSNPs(phSNPs38,gtf)

  message(sprintf("Found %s phSNPs",prettyNum(length(phSNPs38),big.mark=',')))


  # Plot B-allele frequency of hSNPs - in Tumour samples
  # df = phSNPs38[seqnames(phSNPs38) %in% c('chr1','chr7','chr21','chr5')]
  # df = as.data.frame(mcols(df))
  # df$pos = as.numeric(gsub('.*:|_.*$','',rownames(df)))
  # df$chr = gsub(':.*$','',rownames(df))
  # df$BAF_tum = df$altCountTum/df$totCountTum
  # df$altIsMum[is.na(df$altIsMum)] = 'uninformative'
  # p = ggplot(df,aes(pos/1e6,BAF_tum,col=altIsMum))+
  #   geom_point(size=0.01,alpha=0.2)+
  #   scale_color_manual(values = c('black','red',grey(0.7)))+
  #   geom_hline(yintercept = 0.5)+
  #   facet_wrap(vars(chr),scales='free_x')+
  #   theme_classic() + #theme(axis.text.x = element_blank())+
  #   ggtitle(PDID,subtitle = 'L042') + xlab('Genomic position')
  # 
  # print(p)
  # 
  # 
  # df.sub = df[df$chr == 'chr1',]
  # df.sub$pos = as.numeric(df.sub$pos)
  # df.sub$type = ifelse(df.sub$pos > 32375592 & df.sub$pos < 81247540,'segs','others')
  # p = ggplot(df.sub,aes(pos/1e6,BAF_tum))+
  #   geom_point(size=0.01,alpha=0.2)+
  #   geom_hline(yintercept = 0.5)+
  #   facet_grid(vars(type),vars(informative),scales = 'free_x')+
  #   theme_classic() #+ theme(axis.text.x = element_blank())
  # 
  # print(p)


  ########################
  # Integrate with 10X.
  #If the majority of the high coverage SNPs don't look heterozygous, something has gone wrong...
  phCnts = getAllelicExpression(loci=phSNPs38,refGenome = refGenome10X,bams = bams10X,
                                outputs=file.path(outDir,paste0(PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                                nParallel=nParallel,nChunks = 24)
  
  # Save this object so that we can process it faster next time!
  saveRDS(phCnts,file.path(outDir,paste0(PDID,'_phCnts.RDS')))
}



##############
normREF = T
PDID = 'L042'

##------------------------------##
## Add single-cell annotation ####
##------------------------------##
#srat = readRDS(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_231127.RDS'))
srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
srat@meta.data$cellID = rownames(srat@meta.data)
srat$annot = srat$annot_mar24

# Subclustering L042
l042 = subset(srat,subset = donorID == PDID)
l042 = standard_clustering(l042)
DimPlot(l042,group.by = 'seurat_clusters',label = T,repel = T,label.box = T,cols=col25) + NoLegend()
DimPlot(l042,group.by = 'annot_mar24',label = T,repel = T,label.box = T,cols=c(col25,pal34H)) + NoLegend()

passCellIDs = rownames(l042@meta.data)
clusterIDs = setNames(l042@meta.data$cellID,paste0(l042@meta.data$annot,':',l042$timePoint,':',l042$seurat_clusters))

normIDs = setNames(l042@meta.data$cellID[l042$annot %in% c('Mono_CD14')],l042@meta.data$annot[l042$annot %in% c('Mono_CD14')])

#If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
aggToClust=F
#dropUninformative = F
if(aggToClust & !is.null(clusterIDs)){

  gCnts = aggregateByClusters(phCnts,clusters = clusterIDs,preserve = c("REF", "ALT", "geneID", "imprintingFDR",
                                                                        "imprinted", "errRate", "matASE", "matASE_postAlpha", "matASE_postBeta",
                                                                        "regionType", "altIsMum", "geneID", "geneName","informative"))
  gCnts = filterCells(gCnts,normIDs=unique(clusterIDs)[grepl('Mono|B',unique(clusterIDs))])
}else if (!aggToClust & !is.null(clusterIDs)){
  # Not aggToClust but using clusterInfo, including normCells being Leukocytes
  gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=normIDs)
}else if (!aggToClust & is.null(clusterIDs)){
  # No annotation info is passed
  gCnts = filterCells(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs)
  #gCnts = filterCells2(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs,verbose = 2,segs=segs)
}

##################
#Calibrate model
#Specify the error rate
gCnts$errRate = c('Exonic'=0.01,'Intronic'=0.05,'Intergenic'=0.15)[gCnts$regionType]

#Detect allele specific expression
gCnts = calcASE(gCnts,priorKappa = 50)


#Get over-dispersion
od = calcOverDispersion(gCnts)


############
# Inference
#gCnts@metadata$segs$tumFrac = gCnts@metadata$segs$matNum/(gCnts@metadata$segs$matNum + gCnts@metadata$segs$patNum)
pp = abbSegProb(gCnts,od)


#############
# Validation
if(normREF){
  pdf(file.path(outDir,paste0(PDID,'_rawAI_output.pdf')))
}else{
  pdf(file.path(outDir,paste0(PDID,'_rawAI_output_noNormREF.pdf')))
}
dat = plotRawData(gCnts,returnData=TRUE)
p = plotPosteriorHeatmap(pp,'nLL')
#p = plotPosteriorHeatmap(pp,'nLL',split='',km = 4)
print(p)
dev.off()

## Add single-cell call to l042 object
pp2 = pp[pp$cellID %in% l042@meta.data$cellID[l042@meta.data$donorID==PDID]]
m = match(pp2[seqnames(pp2) == 'genomeWide',]$cellID,l042@meta.data$cellID)
l042$AI_output = '?'
l042$AI_output[m] = ifelse(pp2[seqnames(pp2) == 'genomeWide',]$maxPostProb>0.8,pp2[seqnames(pp2) == 'genomeWide',]$mostLikelyState,'Uncalled')

View(table(l042$AI_output,l042$annot,l042$donorID))


DimPlot(l042,group.by = 'AI_output',label = T,repel = T,label.box = T,cols = col25) + NoLegend()
DimPlot(l042,cells.highlight = l042$cellID[l042$donorID == 'L042' & l042$annot == 'Tumour' & l042$AI_output =='abbFrac']) + NoLegend()


## Because of low cell number and small CN fragment, the method is not sensitive enough to call lots of cancer cells with high confidence 
#  --> should use this as validation methods for cancer cells, rather than denovo calling





##------ Generate BAF Copy Number plots  ----####
# mds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MDS/MDS_clean_annotated_tmp.RDS')
# mdat = mds@meta.data
# table(mdat$finalAnn)
mdat = l042@meta.data
normREF = T
# passCellIDs = mdat$cellID
# clusterIDs = setNames(mdat$cellID,mdat$finalAnn)


## Keep uninformative SNPs too
table(phCnts$informative)
phCnts$matCount[is.na(phCnts$matCount)] = ifelse(phCnts$altCountTum[is.na(phCnts$matCount)]>phCnts$refCountTum[is.na(phCnts$matCount)],phCnts$altCount[is.na(phCnts$matCount)],phCnts$refCount[is.na(phCnts$matCount)])
phCnts$patCount[is.na(phCnts$patCount)] = ifelse(phCnts$altCountTum[is.na(phCnts$patCount)]<phCnts$refCountTum[is.na(phCnts$patCount)],phCnts$refCount[is.na(phCnts$patCount)],phCnts$refCount[is.na(phCnts$patCount)])

# There is 1 cell with NA count for G and T???
table(is.na(phCnts$refCount))
phCnts[is.na(phCnts$refCount)]$refCount = 0

## Filter SNPs from imprinting regions, but keep uninformative SNPs
gCnts_og = gCnts
gCnts = filterCells(phCnts,dropUninformative = F,clusterIDs=clusterIDs,normIDs=normIDs,regionsToKeep = c('Exonic','Genic','Intronic'))
gCnts.sub = gCnts[gCnts$cellID %in% mdat$cellID]

gCnts.segs = gCnts.sub
gCnts.segs$regionID = names(gCnts.segs)
clusterIDs = setNames(l042@meta.data$cellID,ifelse(l042$annot == 'Tumour','Tumour','Normal'))
gCnts.segs$clusterID = names(clusterIDs)[match(gCnts.segs$cellID,clusterIDs)]
#gCnts.segs$clusterID = ifelse(grepl('Tumour',gCnts.segs$clusterID),gCnts.segs$clusterID,'Normal')

## Aggregate the allele count across all cells for each SNP 
clusterCnts = aggregateByLists(gCnts.segs, assays = c("altCount", "refCount"), gCnts.segs$clusterID, gCnts.segs$regionID)
clusterCnts.segs = clusterCnts
clusterCnts.segs$pos = as.numeric(gsub('.*:|_.*$','',clusterCnts.segs$regionID))
clusterCnts.segs$chr = gsub(':.*$','',clusterCnts.segs$regionID)
clusterCnts.segs$altFreq = clusterCnts.segs$altCount / (clusterCnts.segs$altCount + clusterCnts.segs$refCount)
clusterCnts.segs$totCount = clusterCnts.segs$altCount + clusterCnts.segs$refCount
clusterCnts.segs$cov = ifelse(clusterCnts.segs$totCount < 5,'<5',
                              ifelse(clusterCnts.segs$totCount < 10,'5-10',
                                     ifelse(clusterCnts.segs$totCount < 20,'10-20','>=20')))
clusterCnts.segs$cov = factor(clusterCnts.segs$cov,c('<5','5-10','10-20','>=20'))
clusterCnts.segs$timePoint = gsub('^.*:','',clusterCnts.segs$cellID)
clusterCnts.segs$celltype = clusterCnts.segs$cellID
clusterCnts.segs$altIsMum = gCnts$altIsMum[match(clusterCnts.segs$regionID,gCnts$regionID)]
clusterCnts.segs$altIsMum[is.na(clusterCnts.segs$altIsMum)] = 'uninformative'

chr_toPlot = paste0('chr',c(5,1,21,3))

chrom = chr_toPlot



plotDir='~/lustre_mt22/Aneuploidy/manuscriptDraft_0424/Plots/'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')

allSNP = T
figSupXX_alleleIntegrator_L067 = function(){
  if(file.exists(file.path(plotDir,paste0('SFigxx_L038Tum_rawBAF_',chrom,'_rawData.tsv')))){
    df = read.delim(file.path(plotDir,paste0('SFigxx_L038Tum_rawBAF_',chrom,'_rawData.tsv')),sep = '\t')
  }
  
  if(!allSNP){
    df_tum = clusterCnts.segs[clusterCnts.segs$totCount >= 5 & clusterCnts.segs$chr %in% chrom & clusterCnts.segs$cellID == 'Tumour',]
    df_norm = clusterCnts.segs[clusterCnts.segs$totCount >= 5 & clusterCnts.segs$chr %in% chrom & clusterCnts.segs$cellID != 'Tumour' & clusterCnts.segs$regionID %in% df_tum$regionID,]
    df = rbind(df_norm,df_tum)
  }else{
    df = clusterCnts.segs[clusterCnts.segs$totCount >= 10 & clusterCnts.segs$chr %in% chrom,]  
  }
  
  df$chr = factor(df$chr,c('chr1','chr5','chr3','chr21'))
  
  plotFun_rawBAF_L042Tum = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    p = ggplot(df,aes(pos/1e6,altFreq))+
      geom_point(aes(col=altIsMum),size=0.1,alpha=1)+
      #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
      facet_grid(cellID  ~ chr,scales = 'free_x')+
      #geom_hline(yintercept = 0.5,lty=1,lwd=0.5,col='black')+
      scale_color_manual(values = c('black','red',grey(0.85)))+
      #scale_color_manual(values = c(grey(0.7),col25[1:2]))+
      #geom_vline(xintercept = c(32375592/1e6,81247540/1e6),col='red')+
      #geom_vline(xintercept = 35e6/1e6,col='red')+
      geom_hline(yintercept = 0.5,lwd=0.3,col='black')+
      scale_y_continuous(breaks = c(0.0,0.5,1.0))+
      theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=F),
        axis.line = element_blank(),axis.text = element_text(size=8))+
      xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
    
    
    print(p)
    
  }
  
  if(!allSNP){
    saveFig(file.path(plotDir,'SFigxx_L042.MLDS_Tum_rawBAF_sameSNPs'),plotFun_rawBAF_L042Tum,rawData=df,width = 10,height = 3,res = 500)  
  }else{
    saveFig(file.path(plotDir,'SFigxx_L042.MLDS_Tum_rawBAF_allSNPs'),plotFun_rawBAF_L042Tum,rawData=df,width = 9,height = 2.7,res = 500)  
  }
  
  
  
  
  
  plotFun_majorAlleleFreq_perCell = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    p = ggplot(df,aes(pos/1e6,altFreq))+
      #geom_point(aes(col=cov),size=0.3,alpha=1)+
      #geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+
      geom_point(aes(col=cov),size=0.1,alpha=1)+
      #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
      facet_grid(cellID  ~ chr,scales = 'free_x')+
      geom_hline(yintercept = 0.5,lty=1,lwd=0.5,col='black')+
      #scale_color_manual(values = c('black','red',grey(0.5)))+
      scale_color_manual(values = c(grey(0.7),col25[1:2]))+
      #geom_vline(xintercept = 115,col='red')+
      #geom_vline(xintercept = 35e6/1e6,col='red')+
      scale_y_continuous(breaks = c(0.0,0.5,1.0))+
      theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=F),
        axis.line = element_blank(),axis.text = element_text(size=8))+
      xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
    
    
    print(p)
    
  }
}







