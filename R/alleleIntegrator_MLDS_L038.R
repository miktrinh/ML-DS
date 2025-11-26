# Run alleleIntegrator on tumour samples

# sudo apt-get install bcftools
# sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
# sudo chmod +x /usr/local/bin/alleleCounter



outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/'
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
library(Seurat)
library(tidyverse)
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/misc.R")

#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
#minSegLen=1e6
#subCl.minSegLen=5e6
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
nParallel=24


############
# PD60301a #
############
#########################
# Sample specific params
tumourDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD60301a/PD60301a.sample.dupmarked.bam'
patientDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD60301a/PD60301a.sample.dupmarked.bam'
bams10X = c('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234191_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234192_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234193_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415922_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415923_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415924_GRCh38-2020-A/possorted_genome_bam.bam')
bams10X = setNames(bams10X,gsub('.*SB_|_GRCh38-2020-A','',basename(dirname(bams10X))))


#Define CN segments roughly
altChrs = c('chr5','chr17','chr21')
segs = GRanges(altChrs,IRanges(c(rep(1,3)),c(1e9,50e6,1e9)))

segs$matNum = c(1,1,2)
segs$patNum = c(0,0,1)
segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
names(segs) = altChrs

PDID = 'PD60301a'


if(skipIfExists & file.exists(file.path(outDir,paste0(PDID,'_phCnts.RDS')))){
  phCnts = readRDS(file.path(outDir,paste0(PDID,'_phCnts.RDS')))

}else{
  ######################
  # Call and phase SNPs

  hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)

  #Expectation is that we'll find ~ 3 million of them
  message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
  #Found 1,861,469 heterozygous SNPs

  # ## Plot B-allele frequency of hSNPs on chr5 and chr7
  # df = hSNPs[seqnames(hSNPs) %in% c('chr5','chr7','chr21','chr22')]
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


  #pdf(file.path(outDir,paste0(PDID,'_AIplots.pdf')))
  #Use tumour DNA to phase them.
  # As we are using tumour as patientDNA (i.e. no matched normal), we will switch off EM, and set minPhasable to be very low.
  # This is because if nPhased/nSNPs < minPhasable, it will set phSNPs$passSanity to be False, as most likely there's no real CN changes
  phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tBAM = tumourDNA,
                           useEM = F,FDR=0.1,minPhasable = 1e-20,
                           outPath=file.path(outDir,paste0(PDID,'_tumour_countAtHetSNPs.tsv')),nParallel=nParallel)
  #pdf(file.path(outDir,paste0(PDID,'_AIplots.pdf')))

  #Liftover to GRCh38
  #phSNPs38 = changeGenomeVersion(phSNPs,liftChain)
  phSNPs38 = phSNPs
  #Annotate SNPs using GTF
  phSNPs38 = annotateSNPs(phSNPs38,gtf)

  ########################
  # Integrate with 10X.
  # If the majority of the high coverage SNPs don't look heterozygous, something has gone wrong...
  phCnts = getAllelicExpression(loci=phSNPs38,refGenome = refGenome10X,bams = bams10X,
                                outputs=file.path(outDir,paste0(PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                                nParallel=nParallel)
  # phCnts$matCount = ifelse(phCnts$altCnt < phCnts$refCnt,phCnts$refCount,phCnts$altCount)
  # phCnts$patCount = ifelse(phCnts$altCnt < phCnts$refCnt,phCnts$altCount,phCnts$refCount)
  #dev.off()
  
  # Save this object so that we can process it faster next time!
  saveRDS(phCnts,file.path(outDir,paste0(PDID,'_phCnts.RDS')))
}


##############
# Add cell type annotation
srat = readRDS(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_231127.RDS'))
srat@meta.data$cellID = rownames(srat@meta.data)
srat$annot = srat$finalAnn_broad

normREF = T
PDID = 'L038'

passCellIDs = rownames(srat@meta.data)
clusterIDs = setNames(srat@meta.data$cellID,srat@meta.data$annot)
normIDs = setNames(srat@meta.data$cellID[srat$annot %in% c('Mono_CD14')],srat@meta.data$annot[srat$annot %in% c('Mono_CD14')])
#clusterIDs = NULL


# gCnts.segs = phCnts[phCnts$cellID %in% srat$cellID[srat$donorID == 'L038'],]
# 
# clusterIDs = setNames(l038@meta.data$cellID,paste0(l038@meta.data$annot_dec23,':',l038$timePoint))
# 
# gCnts.segs$clusterID = names(clusterIDs)[match(gCnts.segs$cellID,clusterIDs)]
# clusterCnts = aggregateByLists(gCnts.segs, assays = c("altCount", "refCount"), gCnts.segs$clusterID, gCnts.segs$regionID)
# #clusterCnts.segs = clusterCnts[grepl('chr5|chr17|chr21|chr8',clusterCnts$regionID),]
# clusterCnts.segs = clusterCnts
# clusterCnts.segs$pos = as.numeric(gsub('.*:|_.*$','',clusterCnts.segs$regionID))
# clusterCnts.segs$chr = gsub(':.*$','',clusterCnts.segs$regionID)
# clusterCnts.segs$altFreq = clusterCnts.segs$altCount / (clusterCnts.segs$altCount + clusterCnts.segs$refCount)
# clusterCnts.segs$totCount = clusterCnts.segs$altCount + clusterCnts.segs$refCount
# clusterCnts.segs$cov = ifelse(clusterCnts.segs$totCount < 5,'<5',
#                               ifelse(clusterCnts.segs$totCount < 10,'<10',
#                                      ifelse(clusterCnts.segs$totCount < 20,'<20','>=20')))
# clusterCnts.segs$cov = factor(clusterCnts.segs$cov,c('<5','<10','<20','>=20'))
# clusterCnts.segs$timePoint = gsub('^.*:','',clusterCnts.segs$cellID)
# clusterCnts.segs$celltype = gsub(':.*$','',clusterCnts.segs$cellID)
# 
# chr_toPlot = paste0('chr',c(4,5,17,19))
# ct_toPlot = c('Tumour','LE','T_CD4','unsure_Tumour','Tumour_WT','unsureME?',unique(clusterCnts.segs$celltype[grepl('^Tumour',clusterCnts.segs$celltype)]))
# ct_toPlot = c(unique(clusterCnts.segs$celltype[grepl('Tumour|T_CD4',clusterCnts.segs$celltype)]))
# ggplot(clusterCnts.segs[clusterCnts.segs$totCount > 5 &
#                           clusterCnts.segs$celltype %in%  ct_toPlot &
#                           clusterCnts.segs$chr %in% chr_toPlot,],aes(pos,altFreq))+
#   geom_point(aes(col=cov),size=0.2,alpha=1)+
#   scale_color_manual(values = c(grey(0.7),col25[c(3,4)]))+
#   facet_grid(celltype + timePoint ~ chr,scales = 'free_x')+
#   theme_classic(base_size = 7)+theme(axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         panel.border = element_rect(fill=F),
#                         axis.line = element_blank())




#If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
aggToClust=FALSE
dropUninformative = F
if(aggToClust & !is.null(clusterIDs)){
  gCnts = aggregateByClusters(phCnts,clusterIDs)
  gCnts = filterCells(gCnts,passCellIDs=levels(clusterIDs),normIDs=normIDs)
}else if (!aggToClust & !is.null(clusterIDs)){
  # Not aggToClust but using clusterInfo, including normCells being Leukocytes
  gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=normIDs,dropUninformative = dropUninformative)
}else if (!aggToClust & is.null(clusterIDs)){
  # No annotation info is passed
  gCnts = filterCells(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs)
}


##################
#Calibrate model
#Specify the error rate
gCnts$errRate = c('Exonic'=0.01,'Intronic'=0.05,'Intergenic'=0.15)[gCnts$regionType]

#Detect allele specific expression
gCnts = calcASE(gCnts)


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
names(gCnts)
gCnts$regionID = gsub(':.*$','',gCnts$regionID)
dat = plotRawData(gCnts,segs = segs,returnData=TRUE)
p = plotPosteriorHeatmap(pp,'nLL')
p = plotPosteriorHeatmap(pp,'nLL',split='',km = 4)
print(p)
dev.off()























#--------------------------------      TP1 sample     --------------------------------####
############
# PD61846a #
############
#########################
# Sample specific params
tumourDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD61846a/PD61846a.sample.dupmarked.bam'
patientDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD61846a/PD61846a.sample.dupmarked.bam'
bams10X = c('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234191_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234192_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234193_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415922_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415923_GRCh38-2020-A/possorted_genome_bam.bam',
            '/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415924_GRCh38-2020-A/possorted_genome_bam.bam')
bams10X = setNames(bams10X,gsub('.*SB_|_GRCh38-2020-A','',basename(dirname(bams10X))))


#Define CN segments roughly
altChrs = c('chr5','chr17','chr21')
segs = GRanges(altChrs,IRanges(c(rep(1,3)),c(1e9,50e6,1e9)))

segs$matNum = c(1,1,2)
segs$patNum = c(0,0,1)
segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
names(segs) = altChrs

PDID = 'PD61846a'


if(skipIfExists & file.exists(file.path(outDir,paste0(PDID,'_phCnts.RDS')))){
  phCnts = readRDS(file.path(outDir,paste0(PDID,'_phCnts.RDS')))

}else{
  ######################
  # Call and phase SNPs

  hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel,
                      minDeviation = 0.25)

  #Expectation is that we'll find ~ 3 million of them
  message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
  #Found 1,784,108 heterozygous SNPs

  # ## Plot B-allele frequency of hSNPs on chr5 and chr17
  # df = hSNPs[seqnames(hSNPs) %in% c('chr5','chr17','chr21','chr22')]
  # df = as.data.frame(mcols(df))
  # df$pos = as.numeric(gsub('.*:|_.*$','',rownames(df)))
  # df$chr = gsub(':.*$','',rownames(df))
  #
  # p = ggplot(df,aes(pos,BAF))+
  #   geom_point(size=0.01,alpha=0.2)+
  #   geom_hline(yintercept = 0.5)+
  #   facet_wrap(vars(chr))+
  #   theme_classic() + theme(axis.text.x = element_blank())
  #
  # print(p)


  #pdf(file.path(outDir,paste0(PDID,'_AIplots.pdf')))
  #Use tumour DNA to phase them.
  # As we are using tumour as patientDNA (i.e. no matched normal), we will switch off EM, and set minPhasable to be very low.
  # This is because if nPhased/nSNPs < minPhasable, it will set phSNPs$passSanity to be False, as most likely there's no real CN changes
  phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tBAM = tumourDNA,
                           useEM = F,FDR=0.1,minPhasable = 1e-20,
                           outPath=file.path(outDir,paste0(PDID,'_tumour_countAtHetSNPs.tsv')),nParallel=nParallel)
  #Liftover to GRCh38
  #phSNPs38 = changeGenomeVersion(phSNPs,liftChain)
  phSNPs38 = phSNPs
  #Annotate SNPs using GTF
  phSNPs38 = annotateSNPs(phSNPs38,gtf)

  ########################
  # Integrate with 10X.
  #If the majority of the high coverage SNPs don't look heterozygous, something has gone wrong...
  phCnts = getAllelicExpression(loci=phSNPs38,refGenome = refGenome10X,bams = bams10X,
                                outputs=file.path(outDir,paste0(PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                                nParallel=nParallel)

  dev.off()
  # Save this object so that we can process it faster next time!
  saveRDS(phCnts,file.path(outDir,paste0(PDID,'_phCnts.RDS')))
}


# ##############
# # Add cell type annotation
# srat = readRDS(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_231127.RDS'))
# srat@meta.data$cellID = rownames(srat@meta.data)
# srat$annot = srat$finalAnn_broad
# 
# normREF = T
# PDID = 'L038'
# 
# passCellIDs = rownames(srat@meta.data[srat$donorID == 'L038',])
# 
# clusterIDs = setNames(srat@meta.data$cellID,paste0(srat@meta.data$annot,':',srat$timePoint))
# normIDs = setNames(srat@meta.data$cellID[srat$annot %in% c('Mono_CD14')],srat@meta.data$annot[srat$annot %in% c('Mono_CD14')])
# #clusterIDs = NULL
# #
# #
# # #If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
# # aggToClust=FALSE
# # dropUninformative = F
# # if(aggToClust & !is.null(clusterIDs)){
# #   gCnts = aggregateByClusters(phCnts,clusterIDs)
# #   gCnts = filterCells(gCnts,passCellIDs=levels(clusterIDs),normIDs=normIDs)
# # }else if (!aggToClust & !is.null(clusterIDs)){
# #   # Not aggToClust but using clusterInfo, including normCells being Leukocytes
# #   gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=normIDs,dropUninformative = dropUninformative)
# # }else if (!aggToClust & is.null(clusterIDs)){
# #   # No annotation info is passed
# #   gCnts = filterCells(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs)
# #   #gCnts = filterCells2(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs,verbose = 2,segs=segs)
# # }
# #
# # dat = plotRawData(gCnts,segs = segs,returnData=TRUE)
# #
# phCnts$clusterID = paste0(srat$annot[match(phCnts$cellID,srat$cellID)],':',
#                           srat$timePoint[match(phCnts$cellID,srat$cellID)],':',
#                           srat$donorID[match(phCnts$cellID,srat$cellID)])
# #max(start(phCnts[seqnames(phCnts) == 'chr1' & phCnts$clusterID == 'Tumour:TP1:L038']))
# 
# gCnts.segs = phCnts[phCnts$cellID %in% srat$cellID[srat$donorID == 'L038'],]
# 
# clusterIDs = setNames(l038@meta.data$cellID,paste0(l038@meta.data$categories,':',l038$timePoint))
# 
# gCnts.segs$clusterID = names(clusterIDs)[match(gCnts.segs$cellID,clusterIDs)]
# clusterCnts = aggregateByLists(gCnts.segs, assays = c("altCount", "refCount"), gCnts.segs$clusterID, gCnts.segs$regionID)
# #clusterCnts.segs = clusterCnts[grepl('chr5|chr17|chr21|chr8',clusterCnts$regionID),]
# clusterCnts.segs = clusterCnts
# clusterCnts.segs$pos = as.numeric(gsub('.*:|_.*$','',clusterCnts.segs$regionID))
# clusterCnts.segs$chr = gsub(':.*$','',clusterCnts.segs$regionID)
# clusterCnts.segs$altFreq = clusterCnts.segs$altCount / (clusterCnts.segs$altCount + clusterCnts.segs$refCount)
# clusterCnts.segs$totCount = clusterCnts.segs$altCount + clusterCnts.segs$refCount
# clusterCnts.segs$cov = ifelse(clusterCnts.segs$totCount < 5,'<5',
#                               ifelse(clusterCnts.segs$totCount < 10,'<10',
#                                      ifelse(clusterCnts.segs$totCount < 20,'<20','>=20')))
# clusterCnts.segs$cov = factor(clusterCnts.segs$cov,c('<5','<10','<20','>=20'))
# clusterCnts.segs$timePoint = gsub('^.*:','',clusterCnts.segs$cellID)
# clusterCnts.segs$celltype = gsub(':.*$','',clusterCnts.segs$cellID)
# 
# chr_toPlot = paste0('chr',c(5,17,19))
# ct_toPlot = c('Tumour','LE','T_CD4','unsure_Tumour','Tumour_WT','unsureME?',unique(clusterCnts.segs$celltype[grepl('^Tumour',clusterCnts.segs$celltype)]))
# ct_toPlot = c(unique(clusterCnts.segs$celltype[grepl('unsure_LE',clusterCnts.segs$celltype)]))
# ggplot(clusterCnts.segs[clusterCnts.segs$totCount > 5 &
#                           clusterCnts.segs$celltype %in%  ct_toPlot &
#                           clusterCnts.segs$chr %in% chr_toPlot,],aes(pos,altFreq))+
#   geom_point(aes(col=cov),size=0.2,alpha=1)+
#   scale_color_manual(values = c(grey(0.7),col25[c(3,4)]))+
#   facet_grid(celltype + timePoint ~ chr,scales = 'free_x')+
#   theme_classic(base_size = 7)+theme(axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         panel.border = element_rect(fill=F),
#                         axis.line = element_blank())
# 
# 
# ##-------- Subclustering L038 -------------####
# 
# l038 = subset(srat,subset = cellID %in% srat$cellID[srat$donorID == 'L038'])
# l038 = standard_clustering(l038)
# DimPlot(l038,label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# DimPlot(l038,group.by = 'annot_dec23',label = T,repel = T,label.box = T,cols = c(pal34H,col25)) + NoLegend()
# DimPlot(l038,cells.highlight = l038$cellID[grepl('unsure_ME',l038$annot) & l038$timePoint != 'Diagnostic'])
# 
# l038$annot_dec23 = as.character(l038$annot_nov23)
# l038$annot_dec23[l038$annot_nov23 == 'unsure_Tumour'] = 'Tumour'
# l038$annot_dec23[l038$annot_nov23 == 'Tumour_WT'] = 'Tumour'
# l038$annot_dec23[l038$seurat_clusters %in% c(7)] = 'Tumour'
# l038$annot_dec23[l038$seurat_clusters %in% c(10) & l038$annot_nov23 == 'unsure_ME'] = 'Tumour'
# 
# ## Add new annotation to main srat object
# srat$annot_dec23 = as.character(srat$annot_nov23)
# srat$annot_dec23[srat$cellID %in% l038$cellID] = l038$annot_dec23[match(srat$cellID[srat$cellID %in% l038$cellID],l038$cellID)]
# write.csv(srat@meta.data,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_mdat_231221')
# 
# 
# 
# l038$categories = as.character(l038$annot_nov23)
# l038$categories[l038$seurat_clusters %in% c(7) & l038$timePoint != 'Diagnostic'] =paste0(l038$categories[l038$seurat_clusters %in% c(7) & l038$timePoint != 'Diagnostic'],'_7')
# l038$categories[l038$seurat_clusters %in% c(10) & l038$timePoint != 'Diagnostic'] =paste0(l038$categories[l038$seurat_clusters %in% c(10) & l038$timePoint != 'Diagnostic'],'_10')
# 
# l038$categories[l038$seurat_clusters %in% c(20,5) & l038$timePoint == 'Diagnostic' & l038$annot == 'Tumour'] ='TumourD_wTP1'
# l038$categories[l038$seurat_clusters %in% c(4,9,1,29) & l038$timePoint != 'Diagnostic' & l038$annot == 'Tumour'] ='TumourTP1_wD'
# l038$categories[l038$seurat_clusters %in% c(7) & l038$timePoint != 'Diagnostic' & l038$annot == 'Tumour'] ='TumourTP1_wEry'
# l038$categories[l038$seurat_clusters %in% c(7) & l038$timePoint != 'Diagnostic' & l038$annot %in% c('ME','unsure_ME')] ='unsureME?'
# DimPlot(l038,cells.highlight = l038$cellID[l038$seurat_clusters %in% c(10,7) & l038$timePoint != 'Diagnostic' & l038$annot == 'Tumour'])
# 
# 
# 
# 
# ##-------- Subclustering L038 Tumour cells only -------------####
# 
# l038_tumour = subset(srat,subset = cellID %in% srat$cellID[srat$donorID == 'L038' & grepl('Tum|MEP|EE|ME|LE',srat$annot_dec23) & 
#                                                              !grepl('unsure',srat$annot_dec23)])
# l038_tumour = standard_clustering(l038_tumour,clusteringRes = 1.3)
# DimPlot(l038_tumour,label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# DimPlot(l038_tumour,group.by = 'timePoint',label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# 
# # cnClone
# cnClust = c(21,24,3,7,23)
# nonClust = c(8,1,9)
# border = c(19,4,10)
# 
# l038_tumour$group = paste0(as.character(l038_tumour$annot_dec23),':',l038_tumour$seurat_clusters)
# l038_tumour$group[l038_tumour$seurat_clusters %in% cnClust] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% cnClust],'_CN')
# l038_tumour$group[l038_tumour$seurat_clusters %in% nonClust] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% nonClust],'_nonCN')
# l038_tumour$group[l038_tumour$seurat_clusters %in% border] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% border],'_?CN')
# 
# 
# df = l038_tumour@meta.data
# df$umap_1 = l038_tumour@reductions$umap@cell.embeddings[,1]
# df$umap_2 = l038_tumour@reductions$umap@cell.embeddings[,2]
# 
# write.csv(df,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_mdat_231221.csv')
# 
# 
# gCnts.segs = phCnts[phCnts$cellID %in% l038_tumour$cellID,]
# 
# clusterIDs = setNames(l038_tumour@meta.data$cellID,
#                       ifelse(l038_tumour$annot_dec23 == 'Tumour',paste0(l038_tumour@meta.data$annot_dec23,l038_tumour$timePoint),
#                              l038_tumour@meta.data$annot_dec23))
# clusterIDs = setNames(l038_tumour@meta.data$cellID,
#                       ifelse(l038_tumour$annot_dec23 == 'Tumour',paste0(l038_tumour@meta.data$group,':',l038_tumour$timePoint),
#                                                           l038_tumour@meta.data$group))
# 
# clusterIDs = setNames(l038_tumour@meta.data$cellID,paste0(as.character(l038_tumour$seurat_clusters),':',l038_tumour$timePoint))
# 
# gCnts.segs$clusterID = names(clusterIDs)[match(gCnts.segs$cellID,clusterIDs)]
# clusterCnts = aggregateByLists(gCnts.segs, assays = c("altCount", "refCount"), gCnts.segs$clusterID, gCnts.segs$regionID)
# clusterCnts.segs = clusterCnts
# clusterCnts.segs$pos = as.numeric(gsub('.*:|_.*$','',clusterCnts.segs$regionID))
# clusterCnts.segs$chr = gsub(':.*$','',clusterCnts.segs$regionID)
# clusterCnts.segs$altFreq = clusterCnts.segs$altCount / (clusterCnts.segs$altCount + clusterCnts.segs$refCount)
# clusterCnts.segs$totCount = clusterCnts.segs$altCount + clusterCnts.segs$refCount
# clusterCnts.segs$cov = ifelse(clusterCnts.segs$totCount < 5,'<5',
#                               ifelse(clusterCnts.segs$totCount < 10,'<10',
#                                      ifelse(clusterCnts.segs$totCount < 20,'<20','>=20')))
# clusterCnts.segs$cov = factor(clusterCnts.segs$cov,c('<5','<10','<20','>=20'))
# clusterCnts.segs$timePoint = gsub('^.*:','',clusterCnts.segs$cellID)
# clusterCnts.segs$celltype = gsub(':.*$','',clusterCnts.segs$cellID)
# 
# chr_toPlot = paste0('chr',c(5,17,19))
# ct_toPlot = c('Tumour','LE','T_CD4','unsure_Tumour','Tumour_WT','unsureME?',unique(clusterCnts.segs$celltype[grepl('^Tumour',clusterCnts.segs$celltype)]))
# 
# clusterCnts.segs$chr = factor(clusterCnts.segs$chr,paste0('chr',c(1:22,'X')))
# clusterCnts.segs$celltype[clusterCnts.segs$celltype == 'TumourDiagnostic'] = 'Tumour_D'
# clusterCnts.segs$celltype[clusterCnts.segs$celltype == 'TumourTP1'] = 'Tumour_TP1'
# 
# 
# 
# DimPlot(l038_tumour,cells.highlight = l038_tumour$cellID[l038_tumour$annot_dec23 == 'Tumour' & 
#                                                            l038_tumour$seurat_clusters == 5])
# 
# ct_toPlot = unique(clusterCnts.segs$cellID[grepl('Tumour:',clusterCnts.segs$cellID)])
# ct_toPlot = unique(clusterCnts.segs$cellID)
# ct_toPlot = c(18,0,19,
#               9,20,27,8,1,7,5,23,21,2,11,25,24,26)
# ct_toPlot = unique(clusterCnts.segs$cellID)#[! unique(clusterCnts.segs$cellID) %in% c(18,0,19,
#               #9,20,27,8,1,7,5,23,21,2,11,25,24,26)]
# 
# pdf('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/BAF_by_seuratClusters.pdf',width = 10,height = 20)
# p = ggplot(clusterCnts.segs[clusterCnts.segs$totCount > 5 &
#                           clusterCnts.segs$cellID %in%  ct_toPlot &
#                           clusterCnts.segs$chr %in% chr_toPlot &
#                           !grepl('Tum_MK',clusterCnts.segs$celltype),],aes(pos,altFreq))+
#   geom_point(aes(col=cov),size=0.2,alpha=1)+
#   scale_color_manual(values = c(grey(0.7),col25[c(3,4)]))+
#   facet_grid(cellID  ~ chr,scales = 'free_x')+
#   theme_classic(base_size = 10)+theme(axis.text.x = element_blank(),
#                                      axis.ticks.x = element_blank(),
#                                      panel.border = element_rect(fill=F),
#                                      axis.line = element_blank())
# 
# 
# print(p)
# dev.off()
# 
# 
cn_g1 = c(18,19)
cn_g2 = c(2,11)
noCN_g1 = c(1,8,7)
cn_g2.2 = c(13,14,15)
le_1 = c(22,3,4)
le_2 = c(10)

l038_tumour$group = ifelse(l038_tumour$seurat_clusters %in% cn_g1,'cn_g1',
                           ifelse(l038_tumour$seurat_clusters %in% cn_g2,'cn_g2',
                                  ifelse(l038_tumour$seurat_clusters %in% cn_g2.2,'cn_g2.2',
                                         ifelse(l038_tumour$seurat_clusters %in% noCN_g1,'noCN_g1',
                                                ifelse(l038_tumour$seurat_clusters %in% le_1,'le_1',
                                                       ifelse(l038_tumour$seurat_clusters %in% le_2,'le_2','others'))))))
Idents(l038_tumour) = l038_tumour$group

markers = FindMarkers(l038_tumour,ident.1 = 'cn_g1',ident.2 = 'cn_g2')
markers$geneSym = rownames(markers)
markers = markers[markers$p_val_adj < 0.01,]
markers$ensID = geneMap$ensID[match(markers$geneSym,geneMap$geneSym)]
markers$chr = geneMap$chr[match(markers$geneSym,geneMap$geneSym)]
markers$pos = NA
markers$pos[markers$ensID %in% gns$gene_id] = start(gns[match(markers$ensID[markers$ensID %in% gns$gene_id],gns$gene_id)])
View(markers[markers$chr %in% c('5'),])
markers$comparison = 'cnD_vs_cnTP1'
markers$chr = factor(markers$chr,c(1:22,'X'))
ggplot(markers,aes(pos,avg_log2FC))+
  geom_point(size=0.01)+
  geom_hline(yintercept = 0)+
  facet_wrap(vars(chr))+
  theme_classic() + ylim(-1,1)





markers_2 = FindMarkers(l038_tumour,ident.1 = c('cn_g1','cn_g2'),ident.2 = 'noCN_g1')
markers_2$geneSym = rownames(markers_2)
markers_2 = markers_2[markers_2$p_val_adj < 0.01,]
markers_2$ensID = geneMap$ensID[match(markers_2$geneSym,geneMap$geneSym)]
markers_2$chr = geneMap$chr[match(markers_2$geneSym,geneMap$geneSym)]
markers_2$pos = NA
markers_2$pos[markers_2$ensID %in% gns$gene_id] = start(gns[match(markers_2$ensID[markers_2$ensID %in% gns$gene_id],gns$gene_id)])
View(markers_2[markers_2$chr %in% c('5'),])
markers_2$comparison = 'cnD.TP1_vs_noCN.D'

markers_2$chr = factor(markers_2$chr,c(1:22,'X'))
ggplot(markers_2,aes(pos,avg_log2FC))+
  geom_point(size=0.01)+
  geom_hline(yintercept = 0)+
  facet_wrap(vars(chr))+
  theme_classic() + ylim(-1,1)





markers_3 = FindMarkers(l038_tumour,ident.1 = c('cn_g1'),ident.2 = 'noCN_g1')
markers_3$geneSym = rownames(markers_3)
markers_3 = markers_3[markers_3$p_val_adj < 0.01,]
markers_3$ensID = geneMap$ensID[match(markers_3$geneSym,geneMap$geneSym)]
markers_3$chr = geneMap$chr[match(markers_3$geneSym,geneMap$geneSym)]
markers_3$pos = NA
markers_3$pos[markers_3$ensID %in% gns$gene_id] = start(gns[match(markers_3$ensID[markers_3$ensID %in% gns$gene_id],gns$gene_id)])
View(markers_3[markers_3$chr %in% c('5'),])
markers_3$comparison = 'cnD_vs_noCN.D'

markers_3$chr = factor(markers_3$chr,c(1:22,'X'))
ggplot(markers_3,aes(pos,avg_log2FC))+
  geom_point(size=0.01)+
  geom_hline(yintercept = 0)+
  facet_wrap(vars(chr))+
  theme_classic() + ylim(-1,1)


markers_4  = FindMarkers(l038_tumour,ident.1 = c('cn_g2'),ident.2 = 'cn_g2.2')
markers_4$geneSym = rownames(markers_4)
markers_4 = markers_4[markers_4$p_val_adj < 0.01,]
markers_4$ensID = geneMap$ensID[match(markers_4$geneSym,geneMap$geneSym)]
markers_4$chr = geneMap$chr[match(markers_4$geneSym,geneMap$geneSym)]
markers_4$pos = NA
markers_4$pos[markers_4$ensID %in% gns$gene_id] = start(gns[match(markers_4$ensID[markers_4$ensID %in% gns$gene_id],gns$gene_id)])
View(markers_4[markers_4$chr %in% c('5'),])
markers_4$comparison = 'cnTP1_vs_cnTP1ery'

markers_4$chr = factor(markers_4$chr,c(1:22,'X'))
ggplot(markers_4,aes(pos,avg_log2FC))+
  geom_point(size=0.01)+
  geom_hline(yintercept = 0)+
  facet_wrap(vars(chr))+
  theme_classic() + ylim(-1,1)



markers_5  = FindMarkers(l038_tumour,ident.1 = c('le_1'),ident.2 = 'le_2')
markers_5$geneSym = rownames(markers_5)
markers_5 = markers_5[markers_5$p_val_adj < 0.01,]
markers_5$ensID = geneMap$ensID[match(markers_5$geneSym,geneMap$geneSym)]
markers_5$chr = geneMap$chr[match(markers_5$geneSym,geneMap$geneSym)]
markers_5$pos = NA
markers_5$pos[markers_5$ensID %in% gns$gene_id] = start(gns[match(markers_5$ensID[markers_5$ensID %in% gns$gene_id],gns$gene_id)])
View(markers_5[markers_5$chr %in% c('5'),])
markers_5$comparison = 'LE.wTumTP1_vs_LE.withoutTum'

markers_5$chr = factor(markers_5$chr,c(1:22,'X'))
ggplot(markers_5,aes(pos,avg_log2FC))+
  geom_point(size=0.01)+
  geom_hline(yintercept = 0)+
  facet_wrap(vars(chr))+
  theme_classic() + ylim(-1,1)



l038_markers_2 = do.call(rbind,list(markers,markers_2,markers_3,markers_4,markers_5))
write.csv(l038_markers,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_FindMarkers_results.csv')



# Distribution of log2fc
ggplot(l038_markers_2,aes(chr,avg_log2FC))+
  geom_boxplot(outlier.size = 0.01)+
  geom_hline(yintercept = 0)+
  theme_classic() + 
  facet_wrap(vars(comparison),ncol=1)+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank()) 

l038_markers_2$direction = ifelse(l038_markers_2$avg_log2FC > 0,'up','down')
ggplot(l038_markers_2[l038_markers_2$chr %in% c(1:22,'X'),],aes(pos,avg_log2FC,col=direction))+
  geom_point(size=0.01)+
  theme_classic(base_size = 11) + 
  scale_color_manual(values = col25)+
  geom_hline(yintercept = 0,lwd=0.3)+
  facet_grid(comparison ~ chr,scales = 'free_x')+
  theme(panel.border = element_rect(fill=F),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length.x = unit(0,'mm')) + ggtitle('log2FC of significant markers')







l038_markers$pct.diff = l038_markers$pct.1 - l038_markers$pct.2
View(l038_markers[abs(l038_markers$avg_log2FC) > 0.5 & abs(l038_markers$pct.diff) > 0.2 & 
                    l038_markers$comparison == 'cnTP1_vs_cnTP1ery' & 
                    l038_markers$direction == 'up',])
FeaturePlot(l038_tumour,c('DNM3'))












VlnPlot(l038_tumour,group.by = 'group',features = c('nFeature_RNA','nCount_RNA','percent.mt'))




outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/supplementPlots'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
plotDir = outDir
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')

df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_mdat_231221.csv')
#df = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourOnly_subClustering_mdat_240206.csv')
library(RColorBrewer)
df = df[df$annot_jan24 != "Tum_MK?",]
plotFun_celltype = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))

  celltype_cols = c('Tumour'=grey(0.3),
                    'LE' = '#DABE99',
                    'EE' = "#f79083",
                    'ME' = "#EF4E22",
                    'MEP' = '#8870ad')

  plot(df$umap_1,df$umap_2,
       las=1,
       type='n',
       cex.main = 0.85,xaxt='n',yaxt='n',
       xlab='',ylab='',
       main=ifelse(noFrame,'','L038'),
       frame.plot=F)

  if(!noPlot){
    points(df$umap_1,df$umap_2,
           col = celltype_cols[df$annot_jan24],
           pch = 19,
           cex=0.01)
  }
  #legend(x=-8, y=9,legend=unique(df$annot_jan24),fill = celltype_cols[unique(df$annot_jan24)],lwd = 0,cex = 0.7,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
}

saveFig(file.path(plotDir,'SFigxx_L038sub_celltype_UMAP'),plotFun_celltype,rawData=df,width = 3,height = 2.6,res = 500,useDingbats = T)

df$timePoint = factor(df$timePoint,c('TP1','Diagnostic'))
plotFun_timepoint = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))

  tp_cols = c('Diagnostic'='#635547',
                    'TP1' = '#c19f70')
  df = df[order(df$timePoint),]
  plot(df$umap_1,df$umap_2,
       las=1,
       type='n',
       cex.main = 0.85,xaxt='n',yaxt='n',
       xlab='',ylab='',
       main=ifelse(noFrame,'','L038'),
       frame.plot=F)

  if(!noPlot){
    points(df$umap_1,df$umap_2,
           col = tp_cols[as.character(df$timePoint)],
           pch = 19,
           cex=0.01)
  }
  #legend(x=-8, y=9,legend=unique(as.character(df$timePoint)),fill = tp_cols[unique(as.character(df$timePoint))],lwd = 0,cex = 0.4,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
}

saveFig(file.path(plotDir,'SFigxx_L038sub_timepoint_UMAP'),plotFun_timepoint,rawData=df,width = 3,height = 2.6,res = 500,useDingbats = T)







dd=df
dd$GATA1s_status = dd$GATA1s_status2
dd$GATA1s_status[dd$GATA1s_status == 'noGATA1expr'] = 'No GATA1 expression'
dd$GATA1s_status[dd$GATA1s_status == 'WT'] = 'GATA1 wild type'
#dd$GATA1s_status[dd$GATA1s_status == 'noCov'] = 'GATA1 expression, no mutation coverage'
dd$GATA1s_status[dd$GATA1s_status == 'Mut'] = 'GATA1s mutation'
dd$GATA1s_status[dd$GATA1s_status %in% c('unsure','noCov')] = 'Uninformative'
#Mut, noCov, noExpr, unsure, WT
#unsure, WT, Mut, noCov
col_jul23 = c('#c18ed1','#5E90BE',brewer.pal(8,'OrRd')[c(8)],alpha(brewer.pal(8,'OrRd')[c(2)],0.7),grey(0.85),brewer.pal(8,'GnBu')[c(7)])



plotFun_GATA1status = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))



  ccs = c('No GATA1 expression' = grey(0.8),
          'Uninformative' = grey(0.65),
          'GATA1s mutation' = '#A92821',
          'GATA1 wild type' = '#005579')
  plot(dd$umap_1,dd$umap_2,
       las=1,
       type='n',
       #xlim=c(-13,17),
       #ylim=c(-13,17),
       cex.main = 0.85,xaxt='n',yaxt='n',
       xlab='',ylab='',
       main=ifelse(noFrame,'','L038'),
       frame.plot=F)

  if(!noPlot){
    #Add density contours
    #addDensityContours(dd$UMAP_1,dd$UMAP_2,dd$finalAnn,col=colAlpha('black',0.4),nGrid = 2000)
    points(dd$umap_1,dd$umap_2,
           col = ccs[dd$GATA1s_status],
           pch = 19,
           cex=0.01)


  }
  legend(x=-8.5, y=9,legend=unique(dd$GATA1s_status),fill = ccs[unique(dd$GATA1s_status)],lwd = 0,cex = 0.5,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')

}

saveFig(file.path(plotDir,'SFigxx_L038sub_GATA1s_UMAP'),plotFun_GATA1status,rawData=dd,width = 3,height = 2.6,res = 500,useDingbats = F)



# plotFun_seuratClusters = function(noFrame=FALSE,noPlot=FALSE){
#   par(mar=c(0.1,0.1,1,0.1))
#   
#   set.seed(2105)
#   clusters_cols = sample(c(pal34H,col25),n_distinct(df$seurat_clusters))
#   names(clusters_cols) = unique(df$seurat_clusters)
#   
#   
#   plot(df$umap_1,df$umap_2,
#        las=1,
#        type='n',
#        cex.main = 0.85,xaxt='n',yaxt='n',
#        xlab='',ylab='',
#        main=ifelse(noFrame,'','L038'),
#        frame.plot=F)
#   
#   if(!noPlot){
#     points(df$umap_1,df$umap_2,
#            col = clusters_cols[as.character(df$seurat_clusters)],
#            pch = 19,
#            cex=0.01)
#     
#     #Add coloured labels
#     mids = aggregate(cbind(umap_1,umap_2) ~ seurat_clusters,data=df,FUN=mean)
#     
#     mids$label = mids$seurat_clusters
#     library(plotrix)
#     boxed.labels(mids$umap_1,mids$umap_2,
#                  labels=mids$label,cex = 0.3,xpad = 2,ypad = 2,border = T,
#                  bg=clusters_cols[as.character(mids$label)],
#                  col='black')
#   }
#   #legend(x=-8, y=9,legend=unique(as.character(df$timePoint)),fill = tp_cols[unique(as.character(df$timePoint))],lwd = 0,cex = 0.4,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
# }
# 
# saveFig(file.path(plotDir,'SFigxx_L038sub_seuratClusters_UMAP'),plotFun_seuratClusters,rawData=df,width = 3,height = 3,res = 500,useDingbats = F)
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
# df = clusterCnts.segs[!grepl('Tum_MK',df$celltype),]
# chr_toPlot = c('chr5','chr17','chr19')
# ct_toPlot = unique(df$cellID)[grepl('Tumour',df$cellID)]
# df$cov = as.character(df$cov)
# df$cov[df$cov == '<10'] = '5-10'
# df$cov[df$cov == '<20'] = '10-20'
# df$cov = factor(df$cov,c('5-10','10-20','>=20'))
# df$cellID = factor(df$cellID,as.character(c(0,1:27)))
# 
# selected_clusters = c(1,7,8,
#                       0,
#                       18,19,
#                       2,11,21,
#                       13,14,15,
#                       6,10,3)
# ct_toPlot = selected_clusters
# df$cellID = factor(as.character(df$cellID),as.character(selected_clusters))
# 
# 
# 
# 
# plotFun_L038_BAF_byClusters = function(noFrame=FALSE,noPlot=FALSE){
#   par(mar=c(0.1,0.1,1,0.1))
#   
#   p = ggplot(df[df$totCount > 5 &
#                   df$cellID %in%  ct_toPlot &
#                   df$chr %in% chr_toPlot,],aes(pos,altFreq))+
#     geom_hline(yintercept = 0.5,col='black')+
#     geom_point(aes(col=cov),size=0.01,alpha=1)+
#     scale_color_manual(values = c(colAlpha(grey(0.67),0.6),'#65A83E','#532C8A'))+
#     facet_grid(cellID  ~ chr,scales = 'free_x')+
#     theme_classic(base_size = 10)+theme(axis.text.x = element_blank(),
#                                         axis.ticks.x = element_blank(),
#                                         panel.border = element_rect(fill=F),
#                                         axis.line = element_blank())+
#     xlab('Genomic position') + ylab('Alt allele frequency') + 
#     scale_y_continuous(breaks = c(0,0.5,1))
# 
# 
# print(p)
# 
# }
# 
# #saveFig(file.path(plotDir,'SFigxx_L038sub_BAF_bySratClusters_selected'),plotFun_L038_BAF_byClusters,rawData=df,width = 5,height = 10,res = 500,useDingbats = F)
# saveFig(file.path(plotDir,'SFigxx_L038sub_BAF_bySratClusters_selected'),plotFun_L038_BAF_byClusters,rawData=df,width = 5,height = 8,res = 500,useDingbats = F)
# 
# 
# 
# 
# df = clusterCnts.segs
# plotFun_L038_BAF_byClustersTP = function(noFrame=FALSE,noPlot=FALSE){
#   par(mar=c(0.1,0.1,1,0.1))
#   
#   p = ggplot(df[df$totCount > 5 &
#                   #df$cellID %in%  ct_toPlot &
#                   df$celltype %in%  ct_toPlot &
#                   df$chr %in% chr_toPlot,],aes(pos,altFreq))+
#     geom_hline(yintercept = 0.5,col='black')+
#     geom_point(aes(col=cov),size=0.01,alpha=1)+
#     scale_color_manual(values = c(colAlpha(grey(0.67),0.6),'#65A83E','#532C8A'))+
#     facet_grid(celltype+timePoint  ~ chr,scales = 'free_x')+
#     theme_classic(base_size = 10)+theme(axis.text.x = element_blank(),
#                                         axis.ticks.x = element_blank(),
#                                         panel.border = element_rect(fill=F),
#                                         axis.line = element_blank())+
#     xlab('Genomic position') + ylab('Alt allele frequency') + 
#     scale_y_continuous(breaks = c(0,0.5,1))
#   
#   
#   print(p)
#   
#   
# }
# 
# #saveFig(file.path(plotDir,'SFigxx_L038sub_BAF_bySratClusters_selected'),plotFun_L038_BAF_byClusters,rawData=df,width = 5,height = 10,res = 500,useDingbats = F)
# saveFig(file.path(plotDir,'SFigxx_L038sub_BAF_bySratClusters:TP_selected'),plotFun_L038_BAF_byClustersTP,rawData=df,width = 5,height = 20,res = 500,useDingbats = F)
# 
# 
# 
# ggplot(l038_tumour@meta.data,aes(seurat_clusters,fill=timePoint))+
#   geom_bar(position = 'fill')+
#   #facet_wrap(vars(seurat_clusters))+
#   scale_fill_manual(values = col25)+
#   theme_classic(base_size = 13)






























# ## Implement some sort of rolling mean aggregated BAF across all SNPs along chromosomes
# binSize = 1e7
# tgtChrs = paste0('chr',c(1:22,'X'))
# genomicBins = GRanges(rep(tgtChrs,1e9/binSize),
#                       IRanges(sapply(seq(1,1e9,binSize),function(i){rep(i,length(tgtChrs))}),
#                               sapply(seq(binSize,1e9,binSize),function(i){rep(i,length(tgtChrs))})))
# genomicBins$binID = paste0(seqnames(genomicBins),':',seq(1:length(genomicBins)))
# # genomicBins$chr = as.character(seqnames(genomicBins))
# # genomicBins$pos = start(genomicBins)
# # View(as.data.frame(mcols(genomicBins[seqnames(genomicBins) == 'chr21',])))
# clusterCnts = clusterCnts.segs
# clusterCnts_gr = GRanges(clusterCnts$chr,IRanges(clusterCnts$pos,clusterCnts$pos),
#                          clusterID = clusterCnts$cellID,altCount = clusterCnts$altCount,refCount = clusterCnts$refCount,
#                          regionID = clusterCnts$regionID,cov=clusterCnts$cov)
# 
# clusterCnts_gr$varID = seq(1:length(clusterCnts_gr))
# tmp = mergeByOverlaps(clusterCnts_gr,genomicBins)
# clusterCnts_gr$binID = tmp$binID[match(clusterCnts_gr$varID,tmp$varID)]
# clusterCnts_gr = clusterCnts_gr[clusterCnts_gr$cov != '<5']
# 
# gCnts.segs.byChrBins = aggregateByLists(gCnts = clusterCnts_gr, assays = c("altCount", "refCount"), cellList = clusterCnts_gr$clusterID, regionList = clusterCnts_gr$binID)
# gCnts.segs.byChrBins$chr = as.character(seqnames(genomicBins))[match(gCnts.segs.byChrBins$regionID,genomicBins$binID)]
# gCnts.segs.byChrBins$pos = as.numeric(start(genomicBins))[match(gCnts.segs.byChrBins$regionID,genomicBins$binID)]
# 
# gCnts.segs.byChrBins$totCount = gCnts.segs.byChrBins$altCount + gCnts.segs.byChrBins$refCount
# gCnts.segs.byChrBins$majorFreq = ifelse(gCnts.segs.byChrBins$altCount > gCnts.segs.byChrBins$refCount,
#                                         gCnts.segs.byChrBins$altCount/gCnts.segs.byChrBins$totCount,
#                                         gCnts.segs.byChrBins$refCount/gCnts.segs.byChrBins$totCount)
# gCnts.segs.byChrBins$altFreq = gCnts.segs.byChrBins$altCount / (gCnts.segs.byChrBins$altCount + gCnts.segs.byChrBins$refCount)
# 
# gCnts.segs.byChrBins$cov = ifelse(gCnts.segs.byChrBins$totCount < 5,'<5',
#                                   ifelse(gCnts.segs.byChrBins$totCount < 10,'<10',
#                                          ifelse(gCnts.segs.byChrBins$totCount < 20,'<20','>=20')))
# gCnts.segs.byChrBins$cov = factor(gCnts.segs.byChrBins$cov,c('<5','<10','<20','>=20'))
# gCnts.segs.byChrBins$tissue = gsub('^.*:','',gCnts.segs.byChrBins$cellID)
# gCnts.segs.byChrBins$celltype = gsub(':.*$','',gCnts.segs.byChrBins$cellID)
# 
# 
# ggplot(gCnts.segs.byChrBins[gCnts.segs.byChrBins$totCount > 5 &
#                               gCnts.segs.byChrBins$celltype %in%  ct_toPlot &
#                               gCnts.segs.byChrBins$chr %in% chr_toPlot,],aes(pos,majorFreq))+
#   geom_point(aes(col=cov),size=0.8,alpha=1)+
#   scale_color_manual(values = c(grey(0.7),col25[c(3,4)]))+
#   geom_hline(yintercept = 0.5,lwd=0.4)+
#   facet_grid(celltype + tissue ~ chr,scales = 'free_x')+
#   theme_classic()+theme(axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         panel.border = element_rect(fill=F),
#                         axis.line = element_blank())+ylim(0,1)
# 
# 
# ggplot(clusterCnts.segs[clusterCnts.segs$totCount > 5 &
#                           clusterCnts.segs$celltype %in%  ct_toPlot &
#                           clusterCnts.segs$chr %in% chr_toPlot,],aes(pos,altFreq))+
#   geom_point(aes(col=cov),size=0.2,alpha=1)+
#   geom_point(data = gCnts.segs.byChrBins[gCnts.segs.byChrBins$totCount > 5 &
#                                            gCnts.segs.byChrBins$celltype %in%  ct_toPlot &
#                                            gCnts.segs.byChrBins$chr %in% chr_toPlot,],
#              aes(pos,majorFreq,col=cov),size=0.8,alpha=1,col='red')+
#   scale_color_manual(values = c(grey(0.7),col25[c(3,4)]))+
#   geom_hline(yintercept = 0.5,lwd=0.4)+
#   geom_hline(yintercept = c(1/3,2/3),lwd=0.2,col='red')+
#   facet_grid(celltype + tissue ~ chr,scales = 'free_x')+
#   theme_classic()+theme(axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         panel.border = element_rect(fill=F),
#                         axis.line = element_blank())
# 
# 
# 
# 
# 
# # 
# # df = clusterCnts.segs[clusterCnts.segs$chr == 'chr5' & clusterCnts.segs$cellID %in% c('Mono_CD14','Tumour'),]
# 
# 
# 
# # 
# # 
# # ##################
# # #Calibrate model
# # #Specify the error rate
# # gCnts$errRate = c('Exonic'=0.01,'Intronic'=0.05,'Intergenic'=0.15)[gCnts$regionType]
# # 
# # #Detect allele specific expression
# # gCnts = calcASE(gCnts)
# # 
# # 
# # #Get over-dispersion
# # od = calcOverDispersion(gCnts)
# # 
# # 
# # ############
# # # Inference
# # #gCnts@metadata$segs$tumFrac = gCnts@metadata$segs$matNum/(gCnts@metadata$segs$matNum + gCnts@metadata$segs$patNum)
# # pp = abbSegProb(gCnts,od)
# # 
# # 
# # #############
# # # Validation
# # if(normREF){
# #   pdf(file.path(outDir,paste0(PDID,'_rawAI_output.pdf')))
# # }else{
# #   pdf(file.path(outDir,paste0(PDID,'_rawAI_output_noNormREF.pdf')))
# # }
# # names(gCnts)
# # gCnts$regionID = gsub(':.*$','',gCnts$regionID)
# # dat = plotRawData(gCnts,segs = segs,returnData=TRUE)
# # p = plotPosteriorHeatmap(pp,'nLL')
# # p = plotPosteriorHeatmap(pp,'nLL',split='',km = 4)
# # print(p)
# # dev.off()
# 
