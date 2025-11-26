##---   MDS - L067:  single-cell genotyping RNAseq data   -----##
##---   1. run alleleIntegrator to detect CN changes (3p CN-neutral LOH), and determine status of chr21 (if there's enough coverage) ?  ##
##---   2. detect GATA1s mutation (11bp insertion)


outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/aug24/L067/'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)

#############
# Libraries #
#############
# Load libraries
library(alleleIntegrator)
library(GenomicFeatures)
library(tidyverse)
#source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
#source("~/lustre_mt22/generalScripts/utils/misc.R")

#########################
# Set Global parameters #
#########################
tgtChrs=c(1:22) 
skipIfExists = T
normREF = T
mainDir = outDir

refGenome = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa'
refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
liftChain = '~/lustre_mt22/hg19ToHg38_noChr.over.chain'
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
nParallel=24





############
# PD61857a #
############
#########################
# Sample specific params
tumourDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD61857a/PD61857a.sample.dupmarked.bam' # D0 sample
patientDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD61857a/PD61857a.sample.dupmarked.bam' # Diagnostic sample, but only 1.8% blasts --> pretty much just normal
bams10X = c('/lustre/scratch126/casm/team274sb/project_folders/GOSH_Leuk/sc_raw_data/cellranger700_count_47089_SB_Leuk13645530_GRCh38-2020-A/possorted_genome_bam.bam') # D0 sample
bams10X = setNames(bams10X,gsub('.*SB_|_GRCh38-2020-A','',basename(dirname(bams10X))))


#Define CN segments roughly
altChrs = c('chr3','chr21')
segs = GRanges(altChrs,IRanges(c(115e6,1),c(1e9,1e9)))
segs$matNum = c(2,2)
segs$patNum = c(0,1)
segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
names(segs) = c('chr3','chr21')

PDID = 'PD61857a'





if(skipIfExists & file.exists(file.path(outDir,paste0(PDID,'_phCnts.RDS')))){
  phCnts = readRDS(file.path(outDir,paste0(PDID,'_phCnts.RDS')))
  phCnts$passSanity = T
  
}else{
  ######################
  # Call and phase SNPs
  
  hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
  
  #Expectation is that we'll find ~ 3 million of them
  message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
  #Found 1,821,916 heterozygous SNPs
  
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
  # Use tumour DNA to phase them.
  # As we are using tumour as patientDNA (i.e. no matched normal), we will switch off EM, and set minPhasable to be very low.
  # This is because if nPhased/nSNPs < minPhasable, it will set phSNPs$passSanity to be False, as most likely there's no real CN changes
  phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tBAM = tumourDNA,
                           useEM = F,FDR = 0.2,
                           outPath=file.path(outDir,paste0(PDID,'_tumour_countAtHetSNPs.tsv')),nParallel=nParallel)
  
  table(phSNPs$informative,seqnames(phSNPs))
  table(phSNPs$altIsMum,seqnames(phSNPs))
  #Liftover to GRCh38
  #phSNPs38 = changeGenomeVersion(phSNPs,liftChain)
  
  phSNPs38 = phSNPs
  # ## Plot B-allele frequency of hSNPs on chr5 and chr7
  # df = phSNPs38[seqnames(phSNPs38) %in% c('chr5','chr17','chr21','chr22')]
  # df = as.data.frame(mcols(df))
  # df$pos = as.numeric(gsub('.*:|_.*$','',rownames(df)))
  # df$chr = gsub(':.*$','',rownames(df))
  # df$tumBAF = df$altCountTum / df$totCountTum
  # 
  # p = ggplot(df[df$totCountTum > 50,],aes(pos,tumBAF))+
  #   geom_point(size=0.01,alpha=0.2)+
  #   geom_hline(yintercept = 0.5)+
  #   facet_wrap(vars(chr))+
  #   geom_vline(xintercept = 1000000000,col='red')+
  #   theme_classic() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  # 
  # print(p)
  
  
  #Annotate SNPs using GTF
  phSNPs38 = annotateSNPs(phSNPs38,gtf)
  
  ########################
  # Integrate with 10X.
  # If the majority of the high coverage SNPs don't look heterozygous, something has gone wrong...
  phCnts = getAllelicExpression(loci=phSNPs38,refGenome = refGenome10X,bams = bams10X,
                                outputs=file.path(outDir,paste0(PDID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
                                nParallel=nParallel)
  
  #dev.off()
  
  # phCnts$pos = as.numeric(gsub('.*:|_.*$','',phCnts$regionID))
  # View(as.data.frame(mcols(phCnts[seqnames(phCnts) == 'chr5' & phCnts$pos > 35000000 & phCnts$pos < 48800000])))
  # table(phCnts[seqnames(phCnts) == 'chr5' & phCnts$pos > 35000000 & phCnts$pos < 48800000]$informative)
  
  # Save this object so that we can process it faster next time!
  saveRDS(phCnts,file.path(outDir,paste0(PDID,'_phCnts.RDS')))
}


##------------------------------##
## Add single-cell annotation ####
##------------------------------##
otherLeuk = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2503.RDS')
mds = subset(otherLeuk,subset = cellID %in% otherLeuk$cellID[otherLeuk$disease == 'MDS' & otherLeuk$donorID == 'L067'])
mds$finalAnn = mds$annot_aug24
# mds_old = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MDS/MDS_clean_annotated_tmp.RDS')
mdat = mds@meta.data
table(mdat$finalAnn)

normREF = T
passCellIDs = mdat$cellID
clusterIDs = setNames(mdat$cellID,mdat$finalAnn)
# Group all normal cell clusters into "Normal"
clusterIDs = setNames(mdat$cellID,ifelse(grepl('Tumour',mdat$finalAnn),mdat$finalAnn,'Normal'))
normIDs = setNames(mdat$cellID[grepl('NK',mdat$finalAnn)],mdat$finalAnn[grepl('NK',mdat$finalAnn)])

#If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
aggToClust=FALSE
if(aggToClust & !is.null(clusterIDs)){
  gCnts = aggregateByClusters(phCnts,clusterIDs)
  gCnts = filterCells(gCnts,passCellIDs=levels(clusterIDs),normIDs=normIDs)
}else if (!aggToClust & !is.null(clusterIDs)){
  # Not aggToClust but using clusterInfo, including normCells being Leukocytes
  gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=normIDs)
}else if (!aggToClust & is.null(clusterIDs)){
  # No annotation info is passed
  gCnts = filterCells(phCnts,clusterIDs=NULL,passCellIDs = passCellIDs)
}



##################
#Calibrate model
#Specify the error rate
gCnts$errRate = c('Exonic'=0.01,'Intronic'=0.05,'Intergenic'=0.15)[gCnts$regionType]

#Detect allele specific expression
gCnts = calcASE(gCnts,priorKappa=40)

# #Get over-dispersion
# od = calcOverDispersion(gCnts)
# 
# 
# ##----------------##
# ##  Inference   ####
# ##----------------##
# pp = abbSegProb(gCnts,od)
# 
# ##-----------------##
# ##  Validation   ####
# ##-----------------##
# dat = plotRawData(gCnts,segs = gCnts@metadata$segs,returnData=TRUE)
# p = plotPosteriorHeatmap(pp,'nLL')
# print(p)
# dev.off()




##------ Generate BAF Copy Number plots  ----####
#mds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MDS/MDS_clean_annotated_tmp.RDS')
otherLeuk = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2503.RDS')
mds = subset(otherLeuk,subset = cellID %in% otherLeuk$cellID[otherLeuk$disease == 'MDS' & otherLeuk$donorID == 'L067'])
mds$finalAnn = mds$annot_aug24
mdat = mds@meta.data
table(mdat$finalAnn)

normREF = T
passCellIDs = mdat$cellID
#clusterIDs = setNames(mdat$cellID,mdat$finalAnn)
clusterIDs = setNames(mdat$cellID,ifelse(grepl('Tumour',mdat$finalAnn),mdat$finalAnn,'Normal'))

## Keep uninformative SNPs too
table(phCnts$informative)
phCnts$matCount[is.na(phCnts$matCount)] = ifelse(phCnts$altCountTum[is.na(phCnts$matCount)]>phCnts$refCountTum[is.na(phCnts$matCount)],phCnts$altCount[is.na(phCnts$matCount)],phCnts$refCount[is.na(phCnts$matCount)])
phCnts$patCount[is.na(phCnts$patCount)] = ifelse(phCnts$altCountTum[is.na(phCnts$patCount)]<phCnts$refCountTum[is.na(phCnts$patCount)],phCnts$refCount[is.na(phCnts$patCount)],phCnts$refCount[is.na(phCnts$patCount)])

table(phCnts$altIsMum)
phCnts$altIsMum[is.na(phCnts$altIsMum)] = ifelse(phCnts$altCountTum[is.na(phCnts$altIsMum)]>phCnts$refCountTum[is.na(phCnts$altIsMum)],T,F)

View(as.data.frame(mcols(phCnts[seqnames(phCnts) == 'chr21'])))

## Filter SNPs from imprinting regions, but keep uninformative SNPs
normIDs = setNames(mdat$cellID[grepl('T_CD4|T_CD8',mdat$finalAnn)],mdat$finalAnn[grepl('T_CD4|T_CD8',mdat$finalAnn)])
gCnts = filterCells(phCnts,dropUninformative = F,dropInsane = F,clusterIDs=clusterIDs,normIDs=normIDs,regionsToKeep = c('Exonic','Genic','Intronic'))
View(as.data.frame(mcols(gCnts[seqnames(gCnts) == 'chr21'])))

gCnts.sub = gCnts[gCnts$cellID %in% mdat$cellID]

gCnts.segs = gCnts.sub
gCnts.segs$regionID = names(gCnts.segs)
gCnts.segs$clusterID = names(clusterIDs)[match(gCnts.segs$cellID,clusterIDs)]

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


# ## Investigate variants with weirdly low VAF in scRNAseq but good BAF in WGS
# var = unique(clusterCnts.segs$regionID[clusterCnts.segs$altFreq < 0.2 & clusterCnts.segs$cellID == 'Normal' & clusterCnts.segs$cov == '>=20'])
# var = var[var %in% names(phCnts[phCnts$BAF > 0.4 & phCnts$BAF < 0.6])]
# table(gsub(':.*$','',var))
# table(phCnts[names(phCnts) %in% var]$geneName)
# View(clusterCnts.segs[clusterCnts.segs$regionID %in% var & clusterCnts.segs$cellID == 'Normal',])
# View(as.data.frame(mcols(phCnts[names(phCnts) %in% var])))



chr_toPlot = paste0('chr',c(1:22,'X'))
chr_toPlot = paste0('chr',c(3,21,2))
chrom = chr_toPlot



plotDir='~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots/'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')

allSNP = T
figSupXX_alleleIntegrator_L067 = function(){
  if(file.exists(file.path(plotDir,paste0('SFigxx_L038Tum_rawBAF_',chrom,'_rawData.tsv')))){
    #df = read.delim(file.path(plotDir,paste0('SFigxx_L038Tum_rawBAF_',chrom,'_rawData.tsv')),sep = '\t')
    df = read.delim(file.path(plotDir,paste0('SFigxx_L067.MDS_Tum_rawBAF_allSNPs_red_rawData.tsv')),sep = '\t')
  }

  if(!allSNP){
    df_tum = clusterCnts.segs[clusterCnts.segs$totCount >= 5 & clusterCnts.segs$chr %in% chrom & clusterCnts.segs$cellID == 'Tumour',]
    df_norm = clusterCnts.segs[clusterCnts.segs$totCount >= 5 & clusterCnts.segs$chr %in% chrom & clusterCnts.segs$cellID != 'Tumour' & clusterCnts.segs$regionID %in% df_tum$regionID,]
    df = rbind(df_norm,df_tum)
  }else{
    df = clusterCnts.segs[clusterCnts.segs$totCount >= 5 & clusterCnts.segs$chr %in% chrom,]
  }

  df$chr = factor(df$chr,chr_toPlot)

  plotFun_rawBAF_L067Tum = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    if(noPlot & !noFrame){
      tmp = rbind(df[df$altFreq %in% min(df$altFreq),],
                  df[df$altFreq == max(df$altFreq),],
                  df[df$pos == min(df$pos),],
                  df[df$pos == max(df$pos),])
      p = ggplot(tmp,aes(pos/1e6,altFreq))+
        #geom_point(aes(col=cov),size=0.3,alpha=1)+
        #geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+
        geom_hline(yintercept = 0.5,lty=1,lwd=0.3,col='black')+
        geom_point(aes(col=cov),size=0.25,alpha=1)+
        #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
        facet_grid(cellID  ~ chr,scales = 'free_x')+

        #scale_color_manual(values = c('black','red',grey(0.5)))+
        #scale_color_manual(values = c(grey(0.7),col25[1:3]))+
        #scale_color_manual(values = c(grey(0.7),col25[1]))+
        #scale_color_manual(values = c(grey(0.6),grey(0.3),grey(0.05)))+
        scale_color_manual(values = rev(c('#DBB2B2','#C68484','#890000')))+
        #geom_vline(xintercept = 115,col='red')+
        #geom_vline(xintercept = 35e6/1e6,col='red')+
        scale_y_continuous(breaks = c(0.0,0.5,1.0))+
        theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(),
          panel.border = element_rect(fill=F,colour = 'black'),
          strip.background = element_blank(),
          axis.line = element_blank(),axis.text = element_text(size=8,color = 'black'))+
        xlab('Genomic position') + ylab('Aggregated Alt allele frequency')

    }

    if(noFrame & !noPlot){
      p = ggplot(df,aes(pos/1e6,altFreq))+
        #geom_point(aes(col=cov),size=0.3,alpha=1)+
        #geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+

        geom_point(aes(col=cov),size=0.25,alpha=1)+
        #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
        facet_grid(cellID  ~ chr,scales = 'free_x')+

        #scale_color_manual(values = c('black','red',grey(0.5)))+
        #scale_color_manual(values = c(grey(0.7),col25[1:3]))+
        #scale_color_manual(values = c(grey(0.7),col25[1]))+
        #scale_color_manual(values = c(grey(0.6),grey(0.3),grey(0.05)))+
        scale_color_manual(values = rev(c('#DBB2B2','#C68484','#890000')))+
        #geom_vline(xintercept = 115,col='red')+
        #geom_vline(xintercept = 35e6/1e6,col='red')+
        scale_y_continuous(breaks = c(0.0,0.5,1.0))+
        theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(),
          panel.border = element_rect(fill=F,colour = 'white'),
          strip.background = element_blank(),
          axis.line = element_blank(),axis.text = element_text(size=8,color = 'black'))+
        xlab('Genomic position') + ylab('Aggregated Alt allele frequency')

    }

    if(!noFrame & !noPlot){
      p = ggplot(df,aes(pos/1e6,altFreq))+
        #geom_point(aes(col=cov),size=0.3,alpha=1)+
        #geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+
        geom_hline(yintercept = 0.5,lty=1,lwd=0.3,col='black')+
        geom_point(aes(col=cov),size=0.25,alpha=1)+
        #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
        facet_grid(cellID  ~ chr,scales = 'free_x')+

        #scale_color_manual(values = c('black','red',grey(0.5)))+
        #scale_color_manual(values = c(grey(0.7),col25[1:3]))+
        #scale_color_manual(values = c(grey(0.7),col25[1]))+
        #scale_color_manual(values = c(grey(0.6),grey(0.3),grey(0.05)))+
        scale_color_manual(values = rev(c('#DBB2B2','#C68484','#890000')))+
        #geom_vline(xintercept = 115,col='red')+
        #geom_vline(xintercept = 35e6/1e6,col='red')+
        scale_y_continuous(breaks = c(0.0,0.5,1.0))+
        theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(),
          panel.border = element_rect(fill=F,colour = 'black'),
          strip.background = element_blank(),
          axis.line = element_blank(),axis.text = element_text(size=8,color = 'black'))+
        xlab('Genomic position') + ylab('Aggregated Alt allele frequency')

    }
    print(p)


  }

  if(!allSNP){
    saveFig(file.path(plotDir,'SFigxx_L067.MDS_Tum_rawBAF_sameSNPs_red'),plotFun_rawBAF_L067Tum,rawData=df,width = 8,height = 3,res = 500)
  }else{
    saveFig(file.path(plotDir,'SFigxx_L067.MDS_Tum_rawBAF_allSNPs_allChr'),plotFun_rawBAF_L067Tum,rawData=df,width = 50,height = 3,res = 500)
    saveFig(file.path(plotDir,'SFigxx_L067.MDS_Tum_rawBAF_allSNPs'),plotFun_rawBAF_L067Tum,rawData=df,width = 7,height = 3,res = 500)
    saveFig(file.path(plotDir,'SFigxx_L067.MDS_Tum_rawBAF_allSNPs_red'),plotFun_rawBAF_L067Tum,rawData=df,width = 9,height = 3,res = 500)
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

## Plot density plot of BAF. There might be some hint of T21 in tumour population, but there's too few SNP coverage to say anything....
mdat$group = ifelse(mdat$finalAnn == 'Tumour','Tumour','Normal')
patientID = 'L067'
source('~/lustre_mt22/alleleIntegrator_BehjatiLab/R/alleleIntegrator_helperFunctions.R')
df = plot_BAF_byCellClusters(mDat=mdat, cellID_column='cellID',group = 'group', normalGroups = 'Normal',
                             outDir,patientID,PDID,
                             phCnts_fp=file.path(outDir,paste0(PDID,'_phCnts.RDS')),
                             tgtChrs=tgtChrs)

## Plot average gene expression between Tumour and other clusters
mds = standard_clustering(mds)

txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

geneMap = read.delim('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data/Donor15680/liver/fLiver_MY_200531_10043298/filtered_feature_bc_matrix/features.tsv.gz',header = F)
colnames(geneMap) = c('ensID','geneSym','GEX')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = NA
geneMap$chr[geneMap$ensID %in% gns$gene_id] = as.character(seqnames(gns[match(geneMap$ensID[geneMap$ensID %in% gns$gene_id],gns$gene_id)]))

avgExpr = AverageExpression(mds,features = geneMap$geneSym[geneMap$chr == 'chr21'],
                            group.by = 'finalAnn')
avgExpr_mtx = as.data.frame(avgExpr$RNA)
library(ComplexHeatmap)
Heatmap(avgExpr_mtx[rowSums(avgExpr_mtx) > 0,],cluster_columns = F,cluster_rows = T,show_row_dend = F,show_column_dend = F,show_row_names = F)

avgExpr_df = avgExpr_mtx
avgExpr_df$geneSym = rownames(avgExpr_df)
avgExpr_df=pivot_longer(avgExpr_df,cols = -geneSym,names_to = 'finalAnn',values_to = 'exp')

ggplot(avgExpr_df,aes(finalAnn,exp))+
  geom_boxplot()+
  scale_y_log10()
  


##----  Plot aggregated major allele frequency across all SNP on copy number region, per cell, and per cluster  ----##
##  Here, since we couldn't confidently phase SNPs, the result is very noisy as we might be including germline SNP at single cell level, where coverage is already quite sparse
##  Filter bad SNPs per cell, only keep SNPs with >5 reads within that cell
##  Overall, I think it's a bad idea...
gCnts$totCount = gCnts$matCount + gCnts$patCount
dat = plotRawData(gCnts[gCnts$totCount > 5,],segs = gCnts@metadata$segs,returnData=TRUE)
dd = as.data.frame(mcols(dat$cellCnts))

ggplot(dd,aes(MAF))+
  geom_histogram(binwidth = 0.01,fill='white',col='black')+
  facet_grid(clusterID~regionID,scales = 'free_y')+
  #geom_vline(xintercept = c(0.5,0.66),col='red')+
  theme_classic()

dd$covBins = cut(dd$totCount, breaks = c(0, 10, 20, Inf))
ggplot(dd, aes(clusterID, MAF)) +
  #geom_hline(data = tmp, aes(yintercept = tumFrac), colour = "red", linetype = "dashed") +
  geom_jitter(aes(colour = covBins), size = 0.2, alpha = 1/1, height = 0) +
  #geom_boxplot(outlier.shape = NA, alpha = 1/100, lwd = 1/4) +
  facet_wrap(~regionID) +
  ylim(0, 1) + ylab("Maternal allele frequency") + xlab("Cluster") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(colour = guide_legend(override.aes = list(size = 2), title = "Coverage")) +
  #geom_point(data = clCnts, colour = "red") +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))


















##--------------------------##
##    GATA1s genotyping   ####
##--------------------------##
mds$donorID = 'L067'
out = mds@meta.data[,c('cellID','donorID','finalAnn')]
out$GATA1s_cov = '?'   # Number of reads mapped to GATA1 gene - which overlaps with the mutation site? (somewhere on exon 2)
out$n_GATA1s_mutant = '?'
out$n_GATA1s_wt = '?'
for(donorID in unique(mds$donorID)){
  print(donorID)
  if(donorID %in% c('L041','CC1')){next()}
  channels = unique(as.character(mds$orig.ident[mds$donorID == donorID]))

  for(channel in channels){
    print(channel)

    ##-- Which cells have GATA1 reads covering that region?

    # All reads covering the region of interest (mutation site) in GATA1 gene
    reads = read.delim(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/',donorID,
                                 paste0(donorID,'_',channel,'_GATA1_reads.txt')),header = F)
    # Remove rows which are not real reads (just another column in the reads)
    reads = reads[!grepl('^UB:Z:',reads$V1),]
    # Extract cellID from read
    reads$cellID = sapply(1:nrow(reads),function(i){reads[i,which(grepl('^CB\\:Z',reads[i,]))]})
    # if(grepl('cellranger',channel)){
    #   reads$cellID = paste0(gsub('_','.',gsub('.*_MY_','MY_',channel)),'_',gsub('^CB:Z:','',reads$cellID))
    # }else{
    #   reads$cellID = paste0(channel,'_',gsub('^CB:Z:','',reads$cellID))
    # }
    reads$cellID = paste0(channel,'_',gsub('^CB:Z:','',reads$cellID))

    message(sprintf('[%s - %s] %d / %d cells found with mutation-site GATA1 reads',donorID,channel,sum(unique(reads$cellID) %in% mds$cellID),length(mds$cellID[mds$donorID == donorID & mds$orig.ident == channel])))
    #mds$GATA1s_cov[mds$cellID %in% reads$cellID] = TRUE
    reads_summary = reads %>% group_by(cellID) %>% summarise(nRead= n())

    out$GATA1s_cov[out$cellID %in% reads_summary$cellID] = reads_summary$nRead[match(out$cellID[out$cellID %in% reads_summary$cellID],reads_summary$cellID)]



    ##-- GATA1 reads with actual mutation
    fp = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads',donorID,
                   paste0(donorID,'_',channel,'_GATA1mut_reads.txt'))

    if(file.size(fp) == 0){
      mut = data.frame(readID = c(),
                       V25 = c(),
                       cellID = c(),
                       containsMUT=c())
    }else{
      mut = read.delim(fp,sep = '\t',header = F,strip.white = T,fill = T,col.names = paste0('V',seq(1:29)))
      mut = mut[!grepl('^UB:Z',mut$V1),]

      mut$cellID = sapply(1:nrow(mut),function(i){mut[i,which(grepl('^CB\\:Z',mut[i,]))]})
      mut$cellID = paste0(channel,'_',gsub('^CB:Z:','',mut$cellID))
      colnames(mut)[1] = c('readID')
      mut = mut[,c('readID','cellID')]
      mut$containsMUT = T
      # check how many cellIDs are included in the seurat object
      message(sprintf('[%s - %s] %d / %d cells found with mutant GATA1 reads',donorID,channel,sum(unique(mut$cellID) %in% mds$cellID),length(mds$cellID[mds$donorID == donorID & mds$orig.ident == channel])))

      # add to output
      mut_summary = mut %>% group_by(cellID) %>% summarise(nRead = n())
      out$n_GATA1s_mutant[out$cellID %in% mut_summary$cellID] = mut_summary$nRead[match(out$cellID[out$cellID %in% mut_summary$cellID],mut_summary$cellID)]
    }


    ##-- GATA1 reads withOUT mutation
    if(file.size(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads',donorID,
                           paste0(donorID,'_',channel,'_GATA1_WT_reads.txt'))) == 0){
      out$n_GATA1s_wt[out$cellID %in% reads_summary$cellID] = NA
    }else{
      wt = read.delim(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads',donorID,
                                paste0(donorID,'_',channel,'_GATA1_WT_reads.txt')),sep = '\t',header = F)
      wt = wt[!grepl('^UB:Z',wt$V1),]
      wt$cellID = sapply(1:nrow(wt),function(i){wt[i,which(grepl('^CB\\:Z',wt[i,]))]})
      wt$cellID = paste0(channel,'_',gsub('^CB:Z:','',wt$cellID))

      colnames(wt)[1] = c('readID')
      wt = wt[,c('readID','cellID')]
      wt$containsMUT = F

      # check how many cellIDs are included in the seurat object
      message(sprintf('[%s - %s] %d / %d cells found with wild-type GATA1 reads',donorID,channel,sum(unique(wt$cellID) %in% mds$cellID),length(mds$cellID[mds$donorID == donorID & mds$orig.ident == channel])))

      ## check if there are duplicated readID being called between MUT vs WT
      if(sum(mut$readID %in% wt$readID) > 0){
        stop(sprintf('Duplicated readID for channel %s from patient %s. Please check!',channel,donorID))
        View(mut[mut$readID %in% wt$readID,])
        View(wt[wt$readID %in% mut$readID,])
      }

      # add to output
      wt_summary = wt %>% group_by(cellID) %>% summarise(nRead = n())
      out$n_GATA1s_wt[out$cellID %in% wt_summary$cellID] = wt_summary$nRead[match(out$cellID[out$cellID %in% wt_summary$cellID],wt_summary$cellID)]

      unsureCells = sum(as.numeric(out$n_GATA1s_mutant[out$n_GATA1s_mutant !='?' & out$n_GATA1s_wt !='?']) > 0 & as.numeric(out$n_GATA1s_wt[out$n_GATA1s_mutant !='?' & out$n_GATA1s_wt !='?']) > 0)
      message(sprintf('[%s - %s] %d cells found with both mutation and WT GATA1 reads',donorID,channel,unsureCells))

    }

    # ## summarise WT / MUT by cellID
    # df = rbind(wt,mut)
    # df = df %>% group_by(cellID,containsMUT) %>% summarise(nReads = n_distinct(readID)) %>%
    #   group_by(cellID) %>% mutate(cat = n_distinct(containsMUT),
    #                               annot = ifelse(cat == 1 & containsMUT,'Tum',
    #                                              ifelse((cat == 1 & !containsMUT),'Norm','?')))
    #
    # # Cells with '?' are most likely those with both normal and mutated reads
    # df.unsure = df[df$cat == 2,]
    #
    #
    # if(nrow(df.unsure) >= 1){
    #   df.unsure$containsMUT = ifelse(df.unsure$containsMUT,'Mut','WT')
    #   df.unsure = pivot_wider(df.unsure,id_cols = c('cellID','cat','annot'),names_from = 'containsMUT',values_from = 'nReads')
    #   df.unsure$annot = ifelse(df.unsure$WT <= 5 & df.unsure$Mut > 2*df.unsure$WT,'Tum',
    #                            ifelse(df.unsure$Mut <= 5 & df.unsure$WT > 2*df.unsure$Mut,'WT','unsure'))
    #
    #   df.unsure$channelID = channel
    #   df.unsure$donorID = donorID
    #   unsure.GATAstatus.all = rbind(unsure.GATAstatus.all,df.unsure)
    #   df$annot[df$cellID %in% df.unsure$cellID] = df.unsure$annot[match(df$cellID[df$cellID %in% df.unsure$cellID],df.unsure$cellID)]
    # }
    #
    #
    # df$channelID = channel
    # df$donorID = donorID
    #
    # tumCells_annot = rbind(tumCells_annot,df)
  }
}

# Cells with GATA1s-site coverage but no WT reads --> WT = 0
out$n_GATA1s_wt[out$GATA1s_cov != '?' & out$n_GATA1s_wt == '?'] = 0
out$n_GATA1s_mutant[out$GATA1s_cov != '?' & out$n_GATA1s_mutant == '?'] = 0

##--------------  GATA expression in each cell ---------------####
# Calculate number of reads mapped to GATA1 gene - does the cell express GATA1 or not?
gata1_expr = mds@assays$RNA@counts['GATA1',]
out$GATA1_UMI_soupXed_count = gata1_expr[match(out$cellID,names(gata1_expr))]

View(out[is.na(out$n_GATA1s_wt) & out$GATA1s_cov > 0,])
out$n_GATA1s_wt[is.na(out$n_GATA1s_wt)] = 0

out$GATA1_status = '?'
out$GATA1_status[out$GATA1s_cov != '?' & out$n_GATA1s_mutant > 0 & out$n_GATA1s_wt == 0] = 'GATA1s_mutant'
out$GATA1_status[out$GATA1s_cov != '?' & out$n_GATA1s_mutant == 0 & out$n_GATA1s_wt > 0] = 'GATA1s_WT'
out$GATA1_status[out$GATA1s_cov != '?' & out$n_GATA1s_mutant > 0 & out$n_GATA1s_wt > 0] = 'GATA1s_unsure'
out$GATA1_status[out$GATA1_status == '?' & out$GATA1_UMI_soupXed_count > 0 ] = 'uninformative'
out$GATA1_status[out$GATA1_status == '?' & out$GATA1_UMI_soupXed_count == 0 ] = 'no_GATA1_expr'

table(out$GATA1_status,out$donorID)

# Save the output
write.csv(out, '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/GATA1mut_reads_summary_perCell_L067_mar24.csv')




##-------    Add results back to MDS   ------####
match(mds$cellID,out$cellID)
colnames(out)[colnames(out) == 'mut_vs_WT_diff'] = 'GATA1_mut.vs.WT_diff'
mds@meta.data = cbind(mds@meta.data,out[match(mds$cellID,out$cellID),!colnames(out) %in% c('cellID','donorID','finalAnn','current_GATA1s_status')])

## Import big.srat (combined leukaemia) annotation
big.srat.mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/combinedLeuk_mar24_mdat.csv')
table(big.srat.mdat$cellID %in% mdat$cellID)
table(big.srat.mdat$finalAnn_broad[big.srat.mdat$cellID %in% mdat$cellID])

mds$finalAnn_v2 = big.srat.mdat$finalAnn_broad[match(mds$cellID,big.srat.mdat$cellID)]
mds$finalAnn_broad = mds$finalAnn_v2
mds$broadLineage = big.srat.mdat$broadLineage[match(mds$cellID,big.srat.mdat$cellID)]

DimPlot(mds,cells.highlight = mds$cellID[mds$GATA1_status == 'GATA1s_mutant'])
DimPlot(mds,group.by = 'broadLineage',cols = col25,label = T,repel = T)

## Create final metadata object
mdat = cbind(mds@meta.data,mds@reductions$umap@cell.embeddings)

# Save the output
write.csv(mdat, '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MDS/MDS_clean_annotated_2404_mdat.csv')






