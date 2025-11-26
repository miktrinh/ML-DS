##---   alleleIntegrator - single-cell genotyping RNAseq data   -----##
##      MLDS - L038, using TP1 sample against diagnostic sample, 
##      as we are interested in detecting the relapse clone ##


# Run alleleIntegrator on L038

# As no matched normal is available for L038, we will try the following:
# 1. Call hSNPs using Diagnostic Tumour DNA (the predominant clone at Diagnosis is one without CN changes)
# 2. Phase hSNPs using TP1 Tumour DNA (the predominant clone at TP1 is one WITH CN changes)

# sudo apt-get install bcftools
# sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
# sudo chmod +x /usr/local/bin/alleleCounter



outDir <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L038/L038_alleleIntegrator"
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
# source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')


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
# PD60301a #
############
#########################
# Sample specific params
tumourDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD61846a/PD61846a.sample.dupmarked.bam' # TP1 sample
patientDNA = '/nfs/cancer_ref01/nst_links/live/3030/PD60301a/PD60301a.sample.dupmarked.bam' # Diagnostic sample
bams10X = c('~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_46351_SB_Leuk13234191_GRCh38-2020-A/possorted_genome_bam.bam', # Diagnostic
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_46351_SB_Leuk13234192_GRCh38-2020-A/possorted_genome_bam.bam', # Diagnostic
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_46351_SB_Leuk13234193_GRCh38-2020-A/possorted_genome_bam.bam', # Diagnostic
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_46682_SB_Leuk13415922_GRCh38-2020-A/possorted_genome_bam.bam', # TP1
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_46682_SB_Leuk13415923_GRCh38-2020-A/possorted_genome_bam.bam', # TP1
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_46682_SB_Leuk13415924_GRCh38-2020-A/possorted_genome_bam.bam') # TP1
bams10X = setNames(bams10X,gsub('.*SB_|_GRCh38-2020-A','',basename(dirname(bams10X))))


#Define CN segments roughly
altChrs = c('chr5','chr17','chr21','chr5')
segs = GRanges(altChrs,IRanges(c(rep(1,3),48.8e6),c(35e6,24e6,1e9,1e9)))

segs$matNum = c(1,1,2,1)
segs$patNum = c(0,0,1,0)
segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
names(segs) = c('chr5','chr17','chr21','chr5.1')

PDID = 'PD60301a'

if(skipIfExists & file.exists(file.path(outDir,paste0(PDID,'.PD61846a_phCnts.RDS')))){
  phCnts = readRDS(file.path(outDir,paste0(PDID,'.PD61846a_phCnts.RDS')))
  
  
}else{
  ######################
  # Call and phase SNPs
  
  hSNPs = findHetSNPs(patientDNA,refGenome,file.path(outDir,paste0(PDID,'_patient_hetSNPs.vcf')),nParallel=nParallel)
  
  #Expectation is that we'll find ~ 3 million of them
  message(sprintf("Found %s heterozygous SNPs",prettyNum(length(hSNPs),big.mark=',')))
  #Found 1,918,678 heterozygous SNPs
  
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
  
  
  pdf(file.path(outDir,paste0(PDID,'_AIplots.pdf')))
  # Use tumour DNA to phase them.
  # As we are using tumour as patientDNA (i.e. no matched normal), we will switch off EM, and set minPhasable to be very low.
  # This is because if nPhased/nSNPs < minPhasable, it will set phSNPs$passSanity to be False, as most likely there's no real CN changes
  phSNPs = phaseSNPsFromCN(hSNPs,segs,refGenome,tBAM = tumourDNA, FDR = 0.2,
                           #useEM = F,
                           outPath=file.path(outDir,paste0(PDID,'PD61846a_tumour_countAtHetSNPs.tsv')),nParallel=nParallel)
  dev.off()
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
  
  dev.off()

  # phCnts$pos = as.numeric(gsub('.*:|_.*$','',phCnts$regionID))
  # View(as.data.frame(mcols(phCnts[seqnames(phCnts) == 'chr5' & phCnts$pos > 35000000 & phCnts$pos < 48800000])))
  # table(phCnts[seqnames(phCnts) == 'chr5' & phCnts$pos > 35000000 & phCnts$pos < 48800000]$informative)
  
  # Save this object so that we can process it faster next time!
  saveRDS(phCnts,file.path(outDir,paste0(PDID,'.PD61846a_phCnts.RDS')))
}


table(phCnts$informative,seqnames(phCnts))

##############
# Add cell type annotation
l038_srat = readRDS(file.path("~/ML-DS/Results/06_MLDS_refractory_relapse/L038", "L038_sratObj.RDS"))
l038_srat@meta.data$cellID = rownames(l038_srat@meta.data)

checkmate::assert_true(all(l038_srat$cellID %in% phCnts$cellID))

## plot BAF in each celltype
l038_srat$group = paste0(l038_srat$broadLineage,':',l038_srat$timePoint)
df = plot_BAF_byCellClusters(mDat=l038_srat@meta.data, cellID_column='cellID',
                             group = 'group', 
                             normalGroups = c('Monocyte/Macrophage:Diagnostic','Monocyte/Macrophage:TP1'),
                             outDir=outDir,patientID=PDID,PDID=PDID,
                             phCnts_fp=file.path(outDir,paste0(PDID,'.PD61846a_phCnts.RDS')),
                             tgtChrs=tgtChrs)

## Calculate single-cell probability of having each CN segments ####
normREF = T
passCellIDs = rownames(l038_srat@meta.data)
clusterIDs = setNames(l038_srat@meta.data$cellID,ifelse(l038_srat@meta.data$annot == 'Tumour',paste0(l038_srat$annot,':',l038_srat$timePoint),l038_srat$annot))
normIDs = setNames(l038_srat@meta.data$cellID[l038_srat$annot %in% c('Mono_CD14')],l038_srat@meta.data$annot[l038_srat$annot %in% c('Mono_CD14')])


#If you don't have much power at the individual cell level, you could consider collapsing all counts into small clusters and then treating each cluster as a cell.  The code below demonstrates how to do this.  After running aggregateByClusters you can proceed in the same way as if there had been no aggregation.
aggToClust=FALSE
if(aggToClust & !is.null(clusterIDs)){
  gCnts = aggregateByClusters(phCnts,clusterIDs)
  gCnts = filterCells(gCnts,passCellIDs=levels(clusterIDs),normIDs=normIDs)
}else if (!aggToClust & !is.null(clusterIDs)){
  # Not aggToClust but using clusterInfo, including normCells being Leukocytes
  gCnts = filterCells(phCnts,clusterIDs=clusterIDs,normIDs=normIDs,regionsToKeep = c('Exonic','Genic','Intronic'))
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

# Re define CN segments
altChrs = c('chr5','chr17')
segs = GRanges(altChrs,IRanges(c(rep(1,2)),c(1e9,24e6)))

segs$matNum = c(1,1)
segs$patNum = c(0,0)
segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
names(segs) = altChrs

gCnts@metadata$segs = segs


pp = abbSegProb(gCnts,od)

saveRDS(pp,file.path(outDir,'L038_PD60301a.PD61846a_pp.RDS'))
pp = readRDS(file.path(outDir,'L038_PD60301a.PD61846a_pp.RDS'))


#############
# Validation
if(normREF){
  pdf(file.path(outDir,paste0(PDID,'_rawAI_output.pdf')))
}else{
  pdf(file.path(outDir,paste0(PDID,'_rawAI_output_noNormREF.pdf')))
}

gCnts$regionID = gsub(':.*$','',gCnts$regionID)

dat = plotRawData(gCnts,segs = segs,returnData=TRUE)
p = plotPosteriorHeatmap(pp,'nLL')
print(p)

pTum = plotPosteriorHeatmap(pp[grepl('^Tumour',pp$clusterID)],'nLL',row_gap = unit(0.5,'cm'))
print(pTum)
dev.off()

## Add single-cell call to srat object ##
m = match(pp[seqnames(pp) == 'genomeWide',]$cellID,l038_srat@meta.data$cellID)
l038_srat$AI_output = '?'
l038_srat$AI_output[m] = ifelse(pp[seqnames(pp) == 'genomeWide',]$maxPostProb>0.95,pp[seqnames(pp) == 'genomeWide',]$mostLikelyState,'Uncalled')
m = match(l038_srat@meta.data$cellID,pp[seqnames(pp) == 'genomeWide',]$cellID)
l038_srat$AI_output_pp = pp[seqnames(pp) == 'genomeWide',]$maxPostProb[m]
FeaturePlot(l038_srat,'AI_output_pp')

m = match(pp[seqnames(pp) == 'chr5',]$cellID,l038_srat$cellID)
l038_srat$AI_output_chr5 = '?'
l038_srat$AI_output_chr5[m] = ifelse(pp[seqnames(pp) == 'chr5' ,]$maxPostProb>0.9,pp[seqnames(pp) == 'chr5',]$mostLikelyState,'Uncalled')
m = match(l038_srat$cellID,pp[seqnames(pp) == 'chr5',]$cellID)
l038_srat$AI_output_chr5_pp = pp[seqnames(pp) == 'chr5',]$maxPostProb[m]

m = match(pp[seqnames(pp) == 'chr17' ,]$cellID,l038_srat$cellID)
l038_srat$AI_output_chr17 = '?'
l038_srat$AI_output_chr17[m] = ifelse(pp[seqnames(pp) == 'chr17',]$maxPostProb>0.9,pp[seqnames(pp) == 'chr17',]$mostLikelyState,'Uncalled')
m = match(l038_srat$cellID,pp[seqnames(pp) == 'chr17',]$cellID)
l038_srat$AI_output_chr17_pp = pp[seqnames(pp) == 'chr17',]$maxPostProb[m]

write.csv(l038_srat@meta.data,file.path(outDir,'L038_AI_res.csv'))








# # 
# # 
# # 
# # ## Add to L038 sub seurat object ## 
# # l038_Tum_mdat = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_mdat_231221.csv'
# # if(file.exists(l038_Tum_mdat)){
# #   df = read.csv(l038_Tum_mdat)
# #   ggplot(df,aes(umap_1,umap_2))+
# #     geom_point(size=0.001,aes(col=as.factor(seurat_clusters)))+
# #     scale_color_manual(values = c(col25,pal34H)) +
# #     theme_classic()
# # }
# # 
# # pp2 = pp[pp$cellID %in% df$cellID]
# # m = match(pp2[seqnames(pp2) == 'genomeWide',]$cellID,df$cellID)
# # df$AI_output = '?'
# # df$AI_output[m] = ifelse(pp2[seqnames(pp2) == 'genomeWide',]$maxPostProb>0.95,pp2[seqnames(pp2) == 'genomeWide',]$mostLikelyState,'Uncalled')
# # ggplot(df,aes(umap_1,umap_2))+
# #   geom_point(size=0.001,aes(col=AI_output))+
# #   scale_color_manual(values = c(col25,pal34H)) +
# #   theme_classic()
# # ggplot(df,aes(umap_1,umap_2))+
# #   geom_point(size=0.001,aes(col=GATA1s_status2))+
# #   scale_color_manual(values = c(col25,pal34H)) +
# #   theme_classic()
# # table(df$AI_output,df$timePoint)
# # 
# # 
# # m = match(pp2[seqnames(pp2) == 'chr5',]$cellID,df$cellID)
# # df$AI_output_chr5 = '?'
# # df$AI_output_chr5[m] = ifelse(pp2[seqnames(pp2) == 'chr5' ,]$maxPostProb>0.9,pp2[seqnames(pp2) == 'chr5',]$mostLikelyState,'Uncalled')
# # 
# # m = match(pp2[seqnames(pp2) == 'chr17' ,]$cellID,df$cellID)
# # df$AI_output_chr17 = '?'
# # df$AI_output_chr17[m] = ifelse(pp2[seqnames(pp2) == 'chr17',]$maxPostProb>0.9,pp2[seqnames(pp2) == 'chr17',]$mostLikelyState,'Uncalled')
# # ggplot(df,aes(umap_1,umap_2))+
# #   geom_point(size=0.001,aes(col=AI_output_chr17))+
# #   scale_color_manual(values = c(col25,pal34H)) +
# #   theme_classic()
# # 
# # 
# # write.csv(df,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_AIRes_240131.csv')
# # 
# # ggplot(df,aes(umap_1,umap_2))+
# #   geom_point(size=0.001,aes(col=AI_output))+
# #   scale_color_manual(values = c(grey(0.9),col25[c(2,1)],grey(0.3),pal34H),name='alleleIntegrator result') +
# #   xlab('UMAP 1') + ylab('UMAP 2') + 
# #   #scale_color_manual(values = c('grey','grey','grey','grey','grey','grey','red','grey','grey')) +
# #   theme_classic(base_size = 13)+
# #   theme(panel.border = element_rect(fill=F),axis.line = element_blank())
# # 
# # ## bar Plot quantify AI call per cell type
# # d = as.data.frame(table(df$annot_dec23,df$timePoint,df$AI_output))
# # d = as.data.frame(table(srat$annot_dec23[srat$donorID == 'L038'],srat$timePoint[srat$donorID == 'L038'],srat$AI_output[srat$donorID == 'L038']))
# # colnames(d) = c('celltype','timePoint','AI_call','nCell')
# # d = d %>% group_by(celltype,timePoint) %>% mutate(totalCell = sum(nCell),
# #                                                   frac = nCell/totalCell)
# # d$AI_call = as.character(d$AI_call)
# # d$AI_call[d$AI_call == '?'] = 'unInformative'
# # d$AI_call = factor(d$AI_call,rev(c('abbFrac','normFrac','Uncalled','unInformative')))
# # d$celltype = paste0(d$celltype,' (n=',d$totalCell,')')
# # ggplot(d,aes(celltype,frac,fill = AI_call))+
# #   geom_col(position = position_fill())+
# #   scale_fill_manual(values = c(grey(0.9),grey(0.6),'#3477B6','#BE3B4A',pal34H),name='alleleIntegrator result')+
# #   facet_wrap(vars(timePoint),scales = 'free_x')+
# #   theme_classic(base_size = 13) + xlab('') + ylab('Fraction of cells')+
# #   theme(panel.border = element_rect(fill=F,linewidth = 1),axis.line = element_blank(),
# #                                      axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) 
# # 
# # 
# # ## bar Plot quantify AI call per cell type per seurat clusters ##
# # d = as.data.frame(table(df$annot_dec23,df$seurat_clusters, df$timePoint,df$AI_output))
# # colnames(d) = c('celltype','seurat_clusters','timePoint','AI_call','nCell')
# # d = d %>% group_by(celltype,seurat_clusters,timePoint) %>% mutate(totalCell = sum(nCell),
# #                                                   frac = nCell/totalCell)
# # d = d[d$totalCell > 0,]
# # 
# # 
# # d$AI_call = as.character(d$AI_call)
# # d$AI_call[d$AI_call == '?'] = 'unInformative'
# # d$AI_call = factor(d$AI_call,rev(c('abbFrac','normFrac','Uncalled','unInformative')))
# # d$celltype = paste0(d$celltype,':',d$seurat_clusters,' (n=',d$totalCell,')')
# # ggplot(d[grepl('Tumour:2 ',d$celltype),],aes(celltype,frac,fill = AI_call))+
# #   geom_col(position = position_fill())+
# #   scale_fill_manual(values = c(grey(0.9),grey(0.6),'#3477B6','#BE3B4A',pal34H),name='alleleIntegrator result')+
# #   facet_wrap(vars(timePoint),scales = 'free_x')+
# #   theme_classic(base_size = 13) + xlab('') + ylab('Fraction of cells')+
# #   theme(panel.border = element_rect(fill=F,linewidth = 1),axis.line = element_blank(),
# #         axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) 
# # 
# # 
# # 
# # only5 = df[df$AI_output_chr5 == 'abbFrac' & df$AI_output_chr17 == 'normFrac',]
# # only17 = df[df$AI_output_chr17 == 'abbFrac' & df$AI_output_chr5 == 'normFrac',]
# # 
# # 
# # ##-------- Subclustering L038 -------------####
# # l038_mdat = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_mdat_231221'
# # if(file.exists(l038_mdat)){
# #   df = read.csv(l038_mdat)
# # }else{
# #   l038 = subset(srat,subset = cellID %in% srat$cellID[srat$donorID == 'L038'])
# #   l038 = standard_clustering(l038)
# #   DimPlot(l038,label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# #   DimPlot(l038,group.by = 'annot_dec23',label = T,repel = T,label.box = T,cols = c(pal34H,col25)) + NoLegend()
# #   DimPlot(l038,cells.highlight = l038$cellID[grepl('unsure_ME',l038$annot) & l038$timePoint != 'Diagnostic'])
# # 
# #   l038$annot_dec23 = as.character(l038$annot_nov23)
# #   l038$annot_dec23[l038$annot_nov23 == 'unsure_Tumour'] = 'Tumour'
# #   l038$annot_dec23[l038$annot_nov23 == 'Tumour_WT'] = 'Tumour'
# #   l038$annot_dec23[l038$seurat_clusters %in% c(7)] = 'Tumour'
# #   l038$annot_dec23[l038$seurat_clusters %in% c(10) & l038$annot_nov23 == 'unsure_ME'] = 'Tumour'
# # 
# #   ## Add new annotation to main srat object
# #   srat$annot_dec23 = as.character(srat$annot_nov23)
# #   srat$annot_dec23[srat$cellID %in% l038$cellID] = l038$annot_dec23[match(srat$cellID[srat$cellID %in% l038$cellID],l038$cellID)]
# #   write.csv(srat@meta.data,l038_mdat)
# # 
# # 
# # 
# #   l038$categories = as.character(l038$annot_nov23)
# #   l038$categories[l038$seurat_clusters %in% c(7) & l038$timePoint != 'Diagnostic'] =paste0(l038$categories[l038$seurat_clusters %in% c(7) & l038$timePoint != 'Diagnostic'],'_7')
# #   l038$categories[l038$seurat_clusters %in% c(10) & l038$timePoint != 'Diagnostic'] =paste0(l038$categories[l038$seurat_clusters %in% c(10) & l038$timePoint != 'Diagnostic'],'_10')
# # 
# #   l038$categories[l038$seurat_clusters %in% c(20,5) & l038$timePoint == 'Diagnostic' & l038$annot == 'Tumour'] ='TumourD_wTP1'
# #   l038$categories[l038$seurat_clusters %in% c(4,9,1,29) & l038$timePoint != 'Diagnostic' & l038$annot == 'Tumour'] ='TumourTP1_wD'
# #   l038$categories[l038$seurat_clusters %in% c(7) & l038$timePoint != 'Diagnostic' & l038$annot == 'Tumour'] ='TumourTP1_wEry'
# #   l038$categories[l038$seurat_clusters %in% c(7) & l038$timePoint != 'Diagnostic' & l038$annot %in% c('ME','unsure_ME')] ='unsureME?'
# #   DimPlot(l038,cells.highlight = l038$cellID[l038$seurat_clusters %in% c(10,7) & l038$timePoint != 'Diagnostic' & l038$annot == 'Tumour'])
# # 
# # }
# # 
# # 
# # 
# # ##-------- Subclustering L038 Tumour cells + EE/ME/LE only -------------####
# # l038_Tum_mdat = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_TumourEry_subClustering_mdat_231221.csv'
# # if(file.exists(l038_Tum_mdat)){
# #   df = read.csv(l038_Tum_mdat)
# #   ggplot(df,aes(umap_1,umap_2))+
# #     geom_point(size=0.001,aes(col=as.character(seurat_clusters)))+
# #     scale_color_manual(values = c(col25,pal34H)) +
# #     theme_classic()
# # 
# # }else{
# #   l038_tumour = subset(srat,subset = cellID %in% srat$cellID[srat$donorID == 'L038' & grepl('Tum|MEP|EE|ME|LE',srat$annot_dec23) &
# #                                                                !grepl('unsure',srat$annot_dec23)])
# #   l038_tumour = standard_clustering(l038_tumour,clusteringRes = 1.3)
# #   DimPlot(l038_tumour,label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# #   DimPlot(l038_tumour,group.by = 'timePoint',label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# # 
# #   # cnClone
# #   cnClust = c(21,24,3,7,23)
# #   nonClust = c(8,1,9)
# #   border = c(19,4,10)
# # 
# #   l038_tumour$group = paste0(as.character(l038_tumour$annot_dec23),':',l038_tumour$seurat_clusters)
# #   l038_tumour$group[l038_tumour$seurat_clusters %in% cnClust] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% cnClust],'_CN')
# #   l038_tumour$group[l038_tumour$seurat_clusters %in% nonClust] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% nonClust],'_nonCN')
# #   l038_tumour$group[l038_tumour$seurat_clusters %in% border] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% border],'_?CN')
# # 
# # 
# #   df = l038_tumour@meta.data
# #   df$umap_1 = l038_tumour@reductions$umap@cell.embeddings[,1]
# #   df$umap_2 = l038_tumour@reductions$umap@cell.embeddings[,2]
# # 
# #   write.csv(df,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_mdat_231221.csv')
# # 
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ##-------- Subclustering L038 Tumour cells + EE/ME/LE only -------------####
# # l038_Tum_mdat = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_mdat_231221.csv'
# # if(file.exists(l038_Tum_mdat)){
# #   df = read.csv(l038_Tum_mdat)
# #   ggplot(df,aes(umap_1,umap_2))+
# #     geom_point(size=0.001,aes(col=as.character(seurat_clusters)))+
# #     scale_color_manual(values = c(col25,pal34H)) +
# #     theme_classic()
# #   
# # }else{
# #   l038_tumour = subset(srat,subset = cellID %in% srat$cellID[srat$donorID == 'L038' & grepl('Tum|MEP|EE|ME|LE',srat$annot_dec23) &
# #                                                                !grepl('unsure',srat$annot_dec23)])
# #   l038_tumour = standard_clustering(l038_tumour,clusteringRes = 1.3)
# #   DimPlot(l038_tumour,label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# #   DimPlot(l038_tumour,group.by = 'timePoint',label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# #   
# #   # cnClone
# #   cnClust = c(21,24,3,7,23)
# #   nonClust = c(8,1,9)
# #   border = c(19,4,10)
# #   
# #   l038_tumour$group = paste0(as.character(l038_tumour$annot_dec23),':',l038_tumour$seurat_clusters)
# #   l038_tumour$group[l038_tumour$seurat_clusters %in% cnClust] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% cnClust],'_CN')
# #   l038_tumour$group[l038_tumour$seurat_clusters %in% nonClust] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% nonClust],'_nonCN')
# #   l038_tumour$group[l038_tumour$seurat_clusters %in% border] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% border],'_?CN')
# #   
# #   
# #   df = l038_tumour@meta.data
# #   df$umap_1 = l038_tumour@reductions$umap@cell.embeddings[,1]
# #   df$umap_2 = l038_tumour@reductions$umap@cell.embeddings[,2]
# #   
# #   write.csv(df,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_mdat_231221.csv')
# #   
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # ##-------- Subclustering L038 Tumour cells only -------------####
# # l038_Tum_mdat = file.path(outDir,'L038_TumourOnly_subClustering_mdat_240206.csv')
# # if(file.exists(l038_Tum_mdat)){
# #   df = read.csv(l038_Tum_mdat)
# #   ggplot(df,aes(umap_1,umap_2))+
# #     geom_point(size=0.001,aes(col=as.character(seurat_clusters)))+
# #     scale_color_manual(values = c(col25,pal34H)) +
# #     theme_classic()
# #   l038_tumour = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'L038' & grepl('Tumour',mlds$annot_aug24) &
# #                                                                !grepl('unsure',mlds$annot_aug24)])
# # }else{
# #   l038_tumour = subset(srat,subset = cellID %in% srat$cellID[srat$donorID == 'L038' & grepl('Tumour',srat$annot_jan24) &
# #                                                                !grepl('unsure',srat$annot_jan24)])
# #   l038_tumour = standard_clustering(l038_tumour)
# #   
# #   aiRes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_AIRes_240131.csv')
# #   l038_tumour$AIres = '?'
# #   l038_tumour$AIres[l038_tumour$cellID %in% aiRes$cellID] = aiRes$AI_output[match(l038_tumour$cellID[l038_tumour$cellID %in% aiRes$cellID],aiRes$cellID)]
# #   
# #   DimPlot(l038_tumour,label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# #   DimPlot(l038_tumour,group.by = 'timePoint',label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
# #   FeaturePlot(l038_tumour,'GATA1')
# #   DimPlot(l038_tumour,cells.highlight = l038_tumour$cellID[l038_tumour$AIres == 'abbFrac']) + NoLegend()
# # 
# #   # # cnClone
# #   # cnClust = c(21,24,3,7,23)
# #   # nonClust = c(8,1,9)
# #   # border = c(19,4,10)
# #   # 
# #   # l038_tumour$group = paste0(as.character(l038_tumour$annot_dec23),':',l038_tumour$seurat_clusters)
# #   # l038_tumour$group[l038_tumour$seurat_clusters %in% cnClust] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% cnClust],'_CN')
# #   # l038_tumour$group[l038_tumour$seurat_clusters %in% nonClust] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% nonClust],'_nonCN')
# #   # l038_tumour$group[l038_tumour$seurat_clusters %in% border] = paste0(l038_tumour$group[l038_tumour$seurat_clusters %in% border],'_?CN')
# # 
# # 
# #   df = l038_tumour@meta.data
# #   df$umap_1 = l038_tumour@reductions$umap@cell.embeddings[,1]
# #   df$umap_2 = l038_tumour@reductions$umap@cell.embeddings[,2]
# # 
# #   write.csv(df,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourOnly_subClustering_mdat_240206.csv')
# # 
# # }
# # 
# # 
# ## Plot L038 Tumour only UMAP
# plotDir = '~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/'
# 
# source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
# library(RColorBrewer)
# 
# ## Import L038 Cancer cells only
# df = read.csv(file.path(outDir,'../L038_TumourOnly_subClustering_mdat_240206.csv'))
# table(df$GATA1s_status2,df$AIres,df$timePoint)
# # subset to remove cells with WT GATA1s
# df = df[df$GATA1s_status2 != 'WT',]
# 
# 
# # plotFun_celltype = function(noFrame=FALSE,noPlot=FALSE){
# #   par(mar=c(0.1,0.1,1,0.1))
# #
# #   celltype_cols = c('Tumour'=grey(0.3),
# #                     'LE' = '#DABE99',
# #                     'EE' = "#f79083",
# #                     'ME' = "#EF4E22",
# #                     'MEP' = '#8870ad')
# #
# #   plot(df$umap_1,df$umap_2,
# #        las=1,
# #        type='n',
# #        cex.main = 0.85,xaxt='n',yaxt='n',
# #        xlab='',ylab='',
# #        main=ifelse(noFrame,'','L038'),
# #        frame.plot=F)
# #
# #   if(!noPlot){
# #     points(df$umap_1,df$umap_2,
# #            col = celltype_cols[df$annot_jan24],
# #            pch = 19,
# #            cex=0.01)
# #   }
# #   #legend(x=-8, y=9,legend=unique(df$annot_jan24),fill = celltype_cols[unique(df$annot_jan24)],lwd = 0,cex = 0.7,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
# # }
# #
# # saveFig(file.path(plotDir,'SFigxx_L038sub_celltype_UMAP'),plotFun_celltype,rawData=df,width = 3,height = 2.6,res = 500,useDingbats = T)
# 
# df$timePoint = factor(df$timePoint,c('TP1','Diagnostic'))
# plotFun_timepoint = function(noFrame=FALSE,noPlot=FALSE){
#   par(mar=c(0.1,0.1,1,0.1))
# 
#   # Pink / purple colors
#   # tp_cols = c('Diagnostic'='#f03aea',#'#84a8d1',
#   #             'TP1' = '#5d07a8')#'#3C5179'
#   # blue colors
#   tp_cols = c('Diagnostic'='#84a8d1',
#               'TP1' = '#3C5179')
#   df = df[order(df$timePoint),]
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
#            col = tp_cols[as.character(df$timePoint)],
#            pch = 19,
#            cex=0.01)
#   }
#   #legend(x=-8, y=9,legend=unique(as.character(df$timePoint)),fill = tp_cols[unique(as.character(df$timePoint))],lwd = 0,cex = 0.4,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
# }
# 
# saveFig(file.path(plotDir,'SFigxx_L038Tum_timepoint_UMAP'),plotFun_timepoint,rawData=df,width = 3,height = 2.6,res = 500,useDingbats = T)
# 
# 
# 
# 
# df$AIres[df$AIres %in% c('?','Uncalled')] = 'Uninformative'
# df$AIres[df$AIres %in% c('normFrac')] = 'without CNA'
# df$AIres[df$AIres %in% c('abbFrac')] = 'with CNA'
# df$AIres = ifelse(df$group %in% c('D_noCNA','TP1_noCNA'),'without CNA',
#                   ifelse(df$group %in% c('D_wCNA','TP1_wCNA'),'with CNA','Uninformative'))
# df$AIres = factor(df$AIres,c('Uninformative','with CNA','without CNA'))
# 
# plotFun_AIresult = function(noFrame=FALSE,noPlot=FALSE){
#   par(mar=c(0.1,0.1,1,0.1))
# 
#   group_cols = c('Uninformative'=grey(0.8),
#                  'without CNA' = grey(0.3),
#                  'with CNA' = '#A92821')
#   df = df[order(df$AIres),]
#   plot(df$umap_1,df$umap_2,
#        las=1,
#        type='n',
#        cex.main = 0.85,xaxt='n',yaxt='n',
#        xlab='',ylab='',
#        main=ifelse(noFrame,'','L038'),
#        frame.plot=F)
# 
#   if(!noPlot){
#     for(g in c('Uninformative','with CNA','without CNA')){
#       points(df$umap_1[df$AIres == g],df$umap_2[df$AIres == g],
#              col = group_cols[as.character(df$AIres[df$AIres == g])],
#              pch = 19,
#              cex=0.1)
#     }
#   }
#   #legend(x=-8, y=9,legend=unique(as.character(df$timePoint)),fill = tp_cols[unique(as.character(df$timePoint))],lwd = 0,cex = 0.4,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
# }
# 
# saveFig(file.path(plotDir,'SFigxx_L038Tum_AIresult_UMAP'),plotFun_AIresult,rawData=df,width = 3,height = 2.6,res = 500,useDingbats = T)
# 
# 
# 
# 
# 
# dd=df
# dd$GATA1s_status = dd$GATA1s_status2
# dd$GATA1s_status[dd$GATA1s_status == 'noGATA1expr'] = 'No GATA1 expression'
# dd$GATA1s_status[dd$GATA1s_status == 'WT'] = 'GATA1 wild type'
# dd$GATA1s_status[dd$GATA1s_status == 'Mut'] = 'GATA1s mutation'
# dd$GATA1s_status[dd$GATA1s_status %in% c('unsure','noCov')] = 'Uninformative'
# #Mut, noCov, noExpr, unsure, WT
# #unsure, WT, Mut, noCov
# col_jul23 = c('#c18ed1','#5E90BE',brewer.pal(8,'OrRd')[c(8)],alpha(brewer.pal(8,'OrRd')[c(2)],0.7),grey(0.85),brewer.pal(8,'GnBu')[c(7)])
# 
# 
# plotFun_GATA1status = function(noFrame=FALSE,noPlot=FALSE){
#   par(mar=c(0.1,0.1,1,0.1))
# 
# 
# 
#   ccs = c('No GATA1 expression' = grey(0.9),
#           'Uninformative' = grey(0.55),
#           'GATA1s mutation' = '#A92821',
#           'GATA1 wild type' = '#005579')
#   plot(dd$umap_1,dd$umap_2,
#        las=1,
#        type='n',
#        #xlim=c(-13,17),
#        #ylim=c(-13,17),
#        cex.main = 0.85,xaxt='n',yaxt='n',
#        xlab='',ylab='',
#        main=ifelse(noFrame,'','L038'),
#        frame.plot=F)
# 
#   if(!noPlot){
#     #Add density contours
#     #addDensityContours(dd$UMAP_1,dd$UMAP_2,dd$finalAnn,col=colAlpha('black',0.4),nGrid = 2000)
#     points(dd$umap_1,dd$umap_2,
#            col = ccs[dd$GATA1s_status],
#            pch = 19,
#            cex=0.01)
# 
# 
#   }
#   #legend(x=-8.5, y=9,legend=unique(dd$GATA1s_status),fill = ccs[unique(dd$GATA1s_status)],lwd = 0,cex = 0.5,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
# 
# }
# 
# saveFig(file.path(plotDir,'SFigxx_L038sub_GATA1s_UMAP'),plotFun_GATA1status,rawData=dd,width = 3,height = 2.6,res = 500,useDingbats = F)
# 
# 
# 
# d = df[df$AIres == 'without CNA' & df$timePoint == 'Diagnostic',]
# df$group = ifelse(df$AIres == 'without CNA' & df$timePoint == 'Diagnostic','D_noCNA',
#                   ifelse(df$AIres == 'with CNA' & df$timePoint == 'Diagnostic','D_wCNA',
#                          ifelse(df$AIres == 'with CNA' & df$timePoint == 'TP1','TP1_wCNA',
#                                 ifelse(df$AIres == 'without CNA' & df$timePoint == 'TP1','TP1_noCNA','others'))))
# df$group = l038_tumour$group_2[match(df$cellID,l038_tumour$cellID)]
# 
# for(g in c('D_noCNA','D_wCNA','TP1_wCNA','TP1_noCNA')){
#   #group_cols=c(grey(0.2),rep(grey(0.8),4))
#   if(grepl('noCNA$',g)){
#     group_cols=c('black',rep(grey(0.8),6))
#   }else if(grepl('wCNA',g)){
#     group_cols=c('#9a6699',rep(grey(0.8),6))
#   }
# 
#   names(group_cols) = c(g,unique(df$group[df$group != g]))
# 
#   plotFun_diagnostic_noCNA = function(noFrame=FALSE,noPlot=FALSE){
#     par(mar=c(0.1,0.1,1,0.1))
# 
#     plot(df$umap_1,df$umap_2,
#          las=1,
#          type='n',
#          cex.main = 0.85,xaxt='n',yaxt='n',
#          xlab='',ylab='',
#          main=ifelse(noFrame,'','L038'),
#          frame.plot=F)
# 
#     if(!noPlot){
#       points(df$umap_1[df$group != g],df$umap_2[df$group != g],
#              col = grey(0.8),
#              pch = 19,
#              cex=0.1)
#       points(df$umap_1[df$group == g],df$umap_2[df$group == g],
#              col = ifelse(grepl('noCNA$',g),'black',
#                           ifelse(grepl('wCNA',g),'#9a6699','red')),
#              pch = 19,
#              cex=0.1)
#     }
#   }
# 
#   saveFig(file.path(plotDir,paste0('SFigxx_L038Tum_highlight_',g,'_UMAP')),plotFun_diagnostic_noCNA,rawData=df,width = 3,height = 2.6,res = 500,useDingbats = T)
# }
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
