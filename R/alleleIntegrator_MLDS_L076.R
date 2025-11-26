# Run alleleIntegrator on tumour samples

# sudo apt-get install bcftools
# sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
# sudo chmod +x /usr/local/bin/alleleCounter


outDir <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L076/L076_alleleIntegrator"
if(!dir.exists(outDir)){
  dir.create(outDir)
}
setwd(outDir)

#############
# Libraries #
#############
# Load libraries
library(alleleIntegrator)
library(ggplot2)
source("~/ML-DS/utils/misc.R")
source("~/ML-DS/utils/sc_utils.R")
source("~/lustre_mt22/alleleIntegrator_mt22/R/alleleIntegrator_helperFunctions.R")

#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
skipIfExists = T
normREF = T
mainDir = outDir

refGenome = '/lustre/scratch125/cellgen/behjati/reference_files/hg38/GRCh38_full_analysis_set_plus_decoy_hla.fa'
refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa' # hg38  
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf' # hg38
liftChain = '/lustre/scratch125/casm/team274sb/mt22/hg19ToHg38_noChr.over.chain'
nParallel=48

############
# PD62331a #
############
#########################
# Sample specific params
tumourDNA = '/nfs/cancer_ref01/nst_links/live/3484/PD62331a/PD62331a.sample.dupmarked.bam'
patientDNA = tumourDNA
#patientDNA = '/nfs/cancer_ref01/nst_links/live/3484/PD64665a/PD64665a.sample.dupmarked.bam'
bams10X = c('~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_47580_SB_Leuk13760338_GRCh38-2020-A/possorted_genome_bam.bam',
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_47260_SB_Leuk13697519_GRCh38-2020-A/possorted_genome_bam.bam',
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_48665_MY_200531_14635833_GRCh38-2020-A/possorted_genome_bam.bam',
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_49031_ALeuk_RNA14832000_GRCh38-2020-A/possorted_genome_bam.bam',
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_49031_ALeuk_RNA14832001_GRCh38-2020-A/possorted_genome_bam.bam')
bams10X = setNames(bams10X,gsub('_','.',gsub('.*MY_','MY_',gsub('.*SB_|.*ALeuk_|_GRCh38-2020-A','',basename(dirname(bams10X))))))


#Define CN segments roughly
altChrs = c('chr5','chr8','chr13',"chr21")
segs = GRanges(altChrs,IRanges(c(3.067e6, 1,   3.321e7,1),
                               c(4.518e6, 1e9, 5.28e7,1e9)))
segs$matNum = c(1,2,1,2)
segs$patNum = c(0,1,0,1)
segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
segs$idx = seq(1:length(segs))
names(segs) = altChrs

PDID = 'PD62331a'


if(skipIfExists & file.exists(file.path(outDir,paste0(PDID,'_phCnts.RDS')))){
  phCnts = readRDS(file.path(outDir,paste0(PDID,'_phCnts.RDS')))
  
}else{
  ######################
  # Call and phase SNPs
  
  hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_tumour_hetSNPs.vcf')),
                      nParallel=nParallel, minDeviation = 0.2)
  
  #Expectation is that we'll find ~ 3 million of them
  message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
  #Found 1,490,566 heterozygous SNPs
  #Found 2,250,268 heterozygous SNPs in PD64665a
  
  ## Plot B-allele frequency of hSNPs on CN segments
  df = hSNPs[seqnames(hSNPs) %in% c('chr8','chr13','chr21','chr22')]
  df = as.data.frame(mcols(df))
  df$pos = as.numeric(gsub('.*:|_.*$','',rownames(df)))
  df$chr = gsub(':.*$','',rownames(df))

  
  p = ggplot(df,aes(pos,BAF))+
    geom_point(size=0.01,alpha=0.2)+
    geom_hline(yintercept = 0.5)+
    facet_wrap(vars(chr))+
    geom_vline(xintercept = c(3.4e7,5.3e7))+
    theme_classic() #+ theme(axis.text.x = element_blank())

  print(p)
  
  
  
  #Use tumour DNA to phase them.
  # As we are using tumour as patientDNA (i.e. no matched normal), we will switch off EM, and set minPhasable to be very low.
  # This is because if nPhased/nSNPs < minPhasable, it will set phSNPs$passSanity to be False, as most likely there's no real CN changes
  phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tBAM = tumourDNA,
                           useEM = T,#FDR=0.1,minPhasable = 1e-20,
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
  
  # Save this object so that we can process it faster next time!
  saveRDS(phCnts,file.path(outDir,paste0(PDID,'_phCnts.RDS')))
}


##############
# Add cell type annotation
l076_srat = readRDS(file.path("~/ML-DS/Results/06_MLDS_refractory_relapse/L076", "L076_sratObj.RDS"))
l076_srat@meta.data$cellID = rownames(l076_srat@meta.data)
l076_srat$annot = l076_srat$finalAnn_broad

checkmate::assert_true(all(l076_srat$cellID %in% phCnts$cellID))
table(l076_srat$tissue[l076_srat$cellID %in% phCnts$cellID])

## plot BAF in each celltype
l076_srat$group = paste0(l076_srat$broadLineage,':',l076_srat$tissue)
df = plot_BAF_byCellClusters(mDat=l076_srat@meta.data, cellID_column='cellID',
                             group = 'group',
                             normalGroups = c('Monocyte/Macrophage:Blood','Monocyte/Macrophage:BM'),
                             outDir=outDir,patientID=PDID,PDID=PDID,
                             phCnts_fp=file.path(outDir,paste0(PDID,'_phCnts.RDS')),
                             tgtChrs=tgtChrs)


## Calculate single-cell probability of having each CN segments ####
segs = segs[seqnames(segs) != 'chr21']
pp = calculate_scCNA_prob(patientID=PDID,outDir=outDir,segs = segs,
                          phCnts_fp=file.path(outDir,paste0(PDID,'_phCnts.RDS')),
                          mDat=l076_srat@meta.data,cellID_column='cellID',mode='tumourDNA_only',
                          group = 'group', normalGroups =c('Monocyte/Macrophage:Blood','Monocyte/Macrophage:BM'),
                          normREF = normREF,aggToClust=FALSE,skipIfExists=skipIfExists)


## Add single-cell call to srat object ##
m = match(pp[seqnames(pp) == 'genomeWide',]$cellID,l076_srat@meta.data$cellID)
l076_srat$AI_output = '?'
l076_srat$AI_output[m] = ifelse(pp[seqnames(pp) == 'genomeWide',]$maxPostProb>0.95,pp[seqnames(pp) == 'genomeWide',]$mostLikelyState,'Uncalled')
m = match(l076_srat@meta.data$cellID,pp[seqnames(pp) == 'genomeWide',]$cellID)
l076_srat$AI_output_pp = pp[seqnames(pp) == 'genomeWide',]$maxPostProb[m]
FeaturePlot(l076_srat,'AI_output_pp')

m = match(pp[seqnames(pp) == 'chr5',]$cellID,l076_srat$cellID)
l076_srat$AI_output_chr5 = '?'
l076_srat$AI_output_chr5[m] = ifelse(pp[seqnames(pp) == 'chr5' ,]$maxPostProb>0.9,pp[seqnames(pp) == 'chr5',]$mostLikelyState,'Uncalled')
m = match(l076_srat$cellID,pp[seqnames(pp) == 'chr5',]$cellID)
l076_srat$AI_output_chr5_pp = pp[seqnames(pp) == 'chr5',]$maxPostProb[m]
FeaturePlot(l076_srat,'AI_output_chr21_pp')

m = match(pp[seqnames(pp) == 'chr8' ,]$cellID,l076_srat$cellID)
l076_srat$AI_output_chr8 = '?'
l076_srat$AI_output_chr8[m] = ifelse(pp[seqnames(pp) == 'chr8',]$maxPostProb>0.9,pp[seqnames(pp) == 'chr8',]$mostLikelyState,'Uncalled')
m = match(l076_srat$cellID,pp[seqnames(pp) == 'chr8',]$cellID)
l076_srat$AI_output_chr8_pp = pp[seqnames(pp) == 'chr8',]$maxPostProb[m]

m = match(pp[seqnames(pp) == 'chr13' ,]$cellID,l076_srat$cellID)
l076_srat$AI_output_chr13 = '?'
l076_srat$AI_output_chr13[m] = ifelse(pp[seqnames(pp) == 'chr13',]$maxPostProb>0.9,pp[seqnames(pp) == 'chr13',]$mostLikelyState,'Uncalled')
m = match(l076_srat$cellID,pp[seqnames(pp) == 'chr13',]$cellID)
l076_srat$AI_output_chr13_pp = pp[seqnames(pp) == 'chr13',]$maxPostProb[m]

m = match(pp[seqnames(pp) == 'chr21' ,]$cellID,l076_srat$cellID)
l076_srat$AI_output_chr21 = '?'
l076_srat$AI_output_chr21[m] = ifelse(pp[seqnames(pp) == 'chr21',]$maxPostProb>0.9,pp[seqnames(pp) == 'chr21',]$mostLikelyState,'Uncalled')
m = match(l076_srat$cellID,pp[seqnames(pp) == 'chr21',]$cellID)
l076_srat$AI_output_chr21_pp = pp[seqnames(pp) == 'chr21',]$maxPostProb[m]

write.csv(l076_srat@meta.data,file.path(outDir,'L076_AI_res.csv'))

# ---- OLD - WITHOUT WRAPPER -------
#If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
passCellIDs = rownames(l076_srat@meta.data)
clusterIDs = setNames(l076_srat@meta.data$cellID,l076_srat$group)
normIDs = setNames(l076_srat@meta.data$cellID[l076_srat$annot %in% c('Mono_CD14')],l076_srat@meta.data$group[l076_srat$annot %in% c('Mono_CD14')])

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
  #gCnts = filterCells2(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs,verbose = 2,segs=segs)
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
gCnts@metadata$segs = gCnts@metadata$segs[seqnames(gCnts@metadata$segs) != 'chr21']
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
