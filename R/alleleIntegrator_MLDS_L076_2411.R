## Unmatched alleleIntegrator - MLDS sample L076

## Set (and create if not exist) output directory
outDir <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L076/L076_alleleIntegrator"
if(!dir.exists(outDir)){
  dir.create(outDir)
}
setwd(outDir)


##-------------------##
##   Libraries     ####
##-------------------##
# Load libraries
library(alleleIntegrator)
library(ggplot2)
source("~/ML-DS/utils/misc.R")
source("~/ML-DS/utils/sc_utils.R")
source("~/lustre_mt22/alleleIntegrator_mt22/R/alleleIntegrator_helperFunctions.R")


##----------------------------##
##   Set Global parameters  ####
##----------------------------##
tgtChrs=c(1:22) 
wgs_ref_version = 'hg38' # or 'hg19'
sc_ref_version = 'hg38' # or 'hg19'
nParallel=48
skipIfExists = T
importOutput_ifExist=T

##--------------------------------------##
##   Set patient specific parameters  ####
##--------------------------------------##
patientID = 'L076' 

# Normal DNA (full path)
normDNA = NA

# Tumour DNA (full path)
# Uppon investigation, it looks like PD64665a is the most promising for calling hetSNPs and phasing hetSNPs
tumourDNA = c('/nfs/cancer_ref01/nst_links/live/3484/PD62331c/PD62331c.sample.dupmarked.bam',
              '/nfs/cancer_ref01/nst_links/live/3484/PD62331a/PD62331a.sample.dupmarked.bam',
              '/nfs/cancer_ref01/nst_links/live/3484/PD64665a/PD64665a.sample.dupmarked.bam',
              '/nfs/cancer_ref01/nst_links/live/3484/PD66167a/PD66167a.sample.dupmarked.bam')
names(tumourDNA) = basename(dirname(tumourDNA))
tumourDNA = tumourDNA[names(tumourDNA) == 'PD62331a']

# vector of scRNA-seq bam files
bams10X = c('~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_47580_SB_Leuk13760338_GRCh38-2020-A/possorted_genome_bam.bam',
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_47260_SB_Leuk13697519_GRCh38-2020-A/possorted_genome_bam.bam',
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_48665_MY_200531_14635833_GRCh38-2020-A/possorted_genome_bam.bam',
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_49031_ALeuk_RNA14832000_GRCh38-2020-A/possorted_genome_bam.bam',
            '~/ML-DS/Data/MLDS_scRNAseq/cellranger700_count_49031_ALeuk_RNA14832001_GRCh38-2020-A/possorted_genome_bam.bam')
bams10X = setNames(bams10X,gsub('_','.',gsub('.*MY_','MY_',gsub('.*SB_|.*ALeuk_|_GRCh38-2020-A','',basename(dirname(bams10X))))))

## Specify file path to either seurat object or metadata table containing the following columns (for step 2):
# - cellID: individual cell barcodes in the format sampleID_barcode. sampleID part must match the name of the corresponding bam file as specified in bams10X, eg. sample1_ACTG...-1 
# - group: categories to group cells by. Eg. normal, tum_1, tum_2 (different tumour clones with different copy number profile as called with inferCNV)
srat_fp = "~/ML-DS/Results/06_MLDS_refractory_relapse/L076/L076_sratObj.RDS" 
mdat_fp = NULL


mode = 'tumourDNA_only'
WGS_CNV_method = NULL
purple_fp = NULL


if(mode %in% c('unmatched','tumourDNA_only')){
  normDNA = '/nfs/cancer_ref01/nst_links/live/3484/PD64665a/PD64665a.sample.dupmarked.bam'
}


##--- Analysed PURPLE outputs
# library(GenomicRanges)
# minSegLen = 1e6
# tgtChrs = paste0(c(1:22))
# segs_PD62331a = processPURPLE(purple_fp='~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_PURPLE/L076/PD62331a/purple_output/without_GRIDDS/PD62331a.purple.cnv.somatic.tsv',
#                               minSegLen = minSegLen,PDID='PD62331a',tgtChrs=tgtChrs,removeBalancedSegs=T,method='allelicRatio')
# segs_PD62331c = processPURPLE(purple_fp='~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_PURPLE/L076/PD62331c/purple_output/without_GRIDDS/PD62331c.purple.cnv.somatic.tsv',
#                               minSegLen = minSegLen,PDID='PD62331c',tgtChrs=tgtChrs,removeBalancedSegs=T,method='allelicRatio')
# 
# segs_PD64665a = processPURPLE(purple_fp='~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_PURPLE/L076/PD64665a/purple_output/without_GRIDDS/PD64665a.purple.cnv.somatic.tsv',
#                               minSegLen = minSegLen,PDID='PD64665a',tgtChrs=tgtChrs,removeBalancedSegs=T,method='allelicRatio')
# segs_PD66167a = processPURPLE(purple_fp='~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_PURPLE/L076/PD66167a/purple_output/without_GRIDDS/PD66167a.purple.cnv.somatic.tsv',
#                               minSegLen = minSegLen,PDID='PD66167a',tgtChrs=tgtChrs,removeBalancedSegs=T,method='allelicRatio')
# 
# segs_combined = do.call(c,list(segs_PD62331a,segs_PD62331c,segs_PD64665a,segs_PD66167a))
# segs_combined$segLen = end(segs_combined) - start(segs_combined) + 1
# segs_combined$segLen_mb = segs_combined$segLen/1e6
# View(as.data.frame(mcols(segs_combined)))
# 
# ## chr22 does not really starts until ~ 17.1mb 
# w = (segs_combined$chr == '22' & end(segs_combined) < 17.2*1e6)
# segs_combined = segs_combined[!w]
# ## I manually inspected chr3, chr4, chr5(46mb-48mb), chr7, chr10, chr12, chr16, chr19 segments on jbrowse. These are low quality regions --> remove them from the list
# ## Purple also called the wrong chr9 loss... manually defined these from jbrowse
# w = (segs_combined$chr %in% c(3,4,7,10,12,16,19))
# segs_combined = segs_combined[!w]




# #Define CN segments roughly
# altChrs = c('chr1','chr2','chr5','chr6','chr8','chr9','chr13','chr14','chr15','chr21')
# segs = GRanges(altChrs,IRanges(c(119.733e6, 62.621e6, 3.067e6, 121.621e6, 1,   21.238e6, 3.321e7, 106.027e6, 83.674e6, 1),
#                                c(123.605e6, 64.278e6, 4.518e6, 122.704e6, 1e9, 22.344e6, 5.28e7,  107.043e6, 84.794e6, 1e9)))
# 
# 
# segs$matNum = c(1,1,1,1,2,1,1,1,1,2)
# segs$patNum = c(0,0,0,0,1,0,0,0,0,1)
# segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
# segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
# segs$idx = seq(1:length(segs))

altChrs = c('chr5','chr8','chr13')
segs = GRanges(altChrs,IRanges(c(3.067e6, 1,   3.321e7),
                               c(4.518e6, 1e9, 5.28e7)))


segs$matNum = c(1,2,1)
segs$patNum = c(0,1,0)
segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
segs$idx = seq(1:length(segs))
names(segs) = altChrs





## Other options to change if needed
# returnData
# PDID: if NULL, PDID = names(tumourDNA)
# refGenome
# refGenome10X
# liftChain
# 
# # Minimum CN segment length to consider. Any segments with shorter length will be ignored
# minSegLen=1e6
# # Minimum sub-clonal CN segment length to consider. Any segments with shorter length will be ignored
# subCl.minSegLen=5e6
# 
# autoSegsFilter = T # Should automated filter to remove segments with low EM confidence be on?
# minVarQual_lowLimit
# min_hSNPs
# 
# autoSegsFilter
# useEM
if(!is.null(WGS_CNV_method)){
  outDir=file.path(outDir,patientID,paste0(mode,'Analysis_w',WGS_CNV_method))  
}else{
  outDir=file.path(outDir,patientID,paste0(mode,'Analysis'))
}

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)

alleleIntegrator_wrapper(WGS_CNV_method = WGS_CNV_method,mode = mode,purple_fp = purple_fp,
                         outDir=outDir,
                         tgtChrs=tgtChrs,nParallel=nParallel, # General parameters
                         patientID=patientID,patientDNA=normDNA,tumourDNA=tumourDNA,#PDID=NULL,
                         wgs_ref_version=wgs_ref_version,sc_ref_version=sc_ref_version,
                         #refGenome=NULL,refGenome10X=NULL,
                         # lifting coordinates from hg19 to hg38 - Only needed if WGS is in hg19. scRNAseq is always mapped to hg38
                         #liftChain = '/lustre/scratch125/casm/team274sb/mt22/hg19ToHg38_noChr.over.chain',
                         # processBTB parameters
                         #minSegLen = 1e6,subCl.minSegLen = 5e6,
                         # findHetSNPs parameters
                         #minVarQual_lowLimit=100,min_hSNPs=900000,skipIfExists=T,
                         # phaseSNPsFromCN_v2 parameters
                         #autoSegsFilter=T,useEM=T,
                         segs=segs,
                         # getAllelicExpression parameters
                         bams10X=bams10X,
                         importOutput_ifExist=importOutput_ifExist
)

##############
# Add cell type annotation
l076_srat = readRDS(file.path("~/ML-DS/Results/06_MLDS_refractory_relapse/L076", "L076_sratObj.RDS"))
l076_srat@meta.data$cellID = rownames(l076_srat@meta.data)

checkmate::assert_true(all(l076_srat$cellID %in% phCnts$cellID))

## plot BAF in each celltype
l076_srat$group = ifelse(l076_srat$broadLineage == 'Tumour',
                         paste0(l076_srat$broadLineage,':',l076_srat$timePoint,':',l076_srat$tissue),
                         l076_srat$broadLineage)
df = plot_BAF_byCellClusters(mDat=l076_srat@meta.data, cellID_column='cellID',
                             group = 'group', 
                             normalGroups = c('Monocyte/Macrophage',"B lineage"),
                             outDir=outDir,patientID=PDID,PDID=PDID,
                             phCnts_fp=file.path(outDir,paste0(PDID,'_phCnts.RDS')),
                             tgtChrs=tgtChrs)


pp = calculate_scCNA_prob(patientID=patientID,outDir,
                          phCnts_fp=file.path(outDir,paste0(PDID,'_phCnts.RDS')),segs = segs['chr8'],
                          mDat=l076_srat@meta.data,cellID_column='cellID',group = 'group', normalGroups = c('Monocyte/Macrophage'),
                          normREF = T,aggToClust=FALSE,skipIfExists=skipIfExists)


# refGenome = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa' # hg38  
# refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa' # hg38  
# gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf' # hg38
# BAF_lim=c(0,1)
# minVarQual=200
# # hSNPs_consensus = c()
# # for(i in 1:length(tumourDNA)){
# #   hSNPs = findHetSNPs_fromTumWGS(bam = tumourDNA[i],refGenome,BAF_lim=BAF_lim,minVarQual=minVarQual,
# #                                  #outVCF = file.path(outDir,paste0(PDID,sprintf('_tumour_BAFlim%.1fto%.1f_hetSNPs.vcf',min(BAF_lim),max(BAF_lim)))),
# #                                  outVCF = file.path(outDir,paste0(names(tumourDNA)[i],'_tumour_hetSNPs.vcf')),
# #                                  nParallel=nParallel)  
# #   if(length(hSNPs_consensus) == 0){
# #     hSNPs_consensus = names(hSNPs)
# #   }else{
# #     hSNPs_consensus = intersect(hSNPs_consensus,names(hSNPs))
# #   }
# # }
# # 
# # 
# # 
# # hSNPs = hSNPs[names(hSNPs) %in% hSNPs_consensus]
# # autoSegsFilter=T
# # phSNPs = phaseSNPsFromCN_v2(hSNPs,segs,refGenome,tBAM = '/nfs/cancer_ref01/nst_links/live/3484/PD66167a/PD66167a.sample.dupmarked.bam',useEM = T,mode=mode,
# #                             outPath=file.path(outDir,paste0('PD66167a_tumour_countAtHetSNPs.tsv')),
# #                             max_taus = 0.8,minPhasable = 0,
# #                             nParallel=nParallel,autoSegsFilter=autoSegsFilter)
# # 
# # table(phSNPs$altIsMum,seqnames(phSNPs))
# # 
# # 
# # # If hSNPs from WGS is in hg19, need to liftover to hg38
# # if(wgs_ref_version == 'hg19' & sc_ref_version == 'hg38'){
# #   phSNPs = changeGenomeVersion(phSNPs_parent,liftChain)
# # }
# # 
# # # Annotate SNPs using GTF
# # phSNPs = annotateSNPs(phSNPs,gtf)
# # 
# # 
# # # Make diagnostic hSNPs per chromosome plot
# # chromInfo = read.delim('/lustre/scratch125/casm/team274sb/mt22/generalResources/chrom_abspos_kb.txt',sep = '\t')
# # chromInfo$chr = paste0('chr',chromInfo$chrom)
# # chromInfo$chr[chromInfo$chr == 'chr23'] = 'chrX'
# # chromInfo$chr = factor(chromInfo$chr,paste0('chr',c(1:22,'X')))
# # chromInfo = chromInfo[chromInfo$arm == 'q',]
# # 
# # 
# # 
# # df = phSNPs
# # df$pos = start(phSNPs)
# # df$chr = seqnames(df)
# # df = as.data.frame(mcols(df))
# # df$chrMax = chromInfo$end[match(as.character(df$chr),chromInfo$chr)]
# # df$tumBAF = df$altCountTum/df$totCountTum
# # df$altIsMum[!df$informative] = 'uninformative'
# # df$chr = factor(df$chr,paste0('chr',c(1:22,'X')))
# # #df$inHSNPs = (paste0(df$chr,':',df$pos) %in% paste0(commonSNPs$Chromosome,':',commonSNPs$Position))
# # PDID = 'PD66167a'
# # png(file.path(outDir,paste0(PDID,'_tumour_countAtHetSNPs.png')),width=1200,height=900)
# # p = ggplot(df,aes(pos/1e6,tumBAF,col=altIsMum))+
# #   #geom_point(size=0.0001,alpha=0.1,aes(col=inHSNPs))+
# #   geom_point(size=0.0001,alpha=0.1)+
# #   geom_hline(yintercept = 0.5)+
# #   scale_color_manual(values = c('black','orange',grey(0.8)))+
# #   scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.0),labels = c(0,0.25,0.5,0.75,1),limits = c(0,1))+
# #   geom_vline(data = chromInfo[chromInfo$chr %in% unique(df$chr),],aes(xintercept = end/1e3),lty=2,col='red')+
# #   theme_classic()+facet_wrap(vars(chr),scales = 'free_x')+
# #   theme(axis.text.x = element_blank())
# # print(p)
# # dev.off()    
# 
# 
# 
# 
# ##----------------------------##
# ##    Integrate with 10X.   ####
# ##----------------------------##
# phCnts_fp = file.path(outDir,paste0(patientID,'_phCnts.RDS'))
# #If the majority of the high coverage SNPs don't look heterozygous, something has gone wrong...
# # phCnts = getAllelicExpression(loci=phSNPs,refGenome = refGenome10X,bams = bams10X,segs=NULL,
# #                               outputs=file.path(outDir,paste0(patientID,'_',names(bams10X),'_scRNA_alleleCounts.tsv')),
# #                               nParallel=nParallel)
# # if(mode != 'normalDNA_only'){
# #   phCnts@metadata$segs = segs  
# #   names(phCnts@metadata$segs) = make.unique(names(phCnts@metadata$segs))
# # }
# # 
# # # set segs as null above to retain non-CNA regions in phCnts object, but then should add segs back in the metadata
# # # Save this object so that we can process it faster next time!
# # saveRDS(phCnts,phCnts_fp)  
# 
# phCnts = readRDS(phCnts_fp)
# 
# 
# 
# 
# 
# #---------------------------------##
# # Import single-cell annotation ####
# #---------------------------------##
# ## Generate L076-only seurat object
# # mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS')
# # l076 = subset(mlds,donorID == 'L076')
# # l076 = standard_clustering(l076)
# # saveRDS(l076,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/L076_sratObj.RDS')
# #
# # DimPlot(l076,group.by = 'timePoint',repel = T,cols = col25)
# # FeaturePlot(l076,c('GATA1','TP53'))
# 
# 
# library(Seurat)
# cellID_column = 'cellID'
# group = 'annot'
# #normalGroups = c(13,7,24,3,11,19)
# normalGroups = c('naive.B','T_cells','T_CD4','T_CD8','Neutrophil','Mono_CD14')
# #normalGroups = NULL
# if(is.null(normalGroups)){
#   normREF=F
# }else{
#   normREF=T
# }
# 
# 
# ## Import scRNAseq metadata
# print('Reading in scRNAseq metadata')
# if(!is.null(srat_fp)){
#   srat = readRDS(srat_fp)
#   srat$annot = ifelse(srat@meta.data[['annot_aug24']] %in% normalGroups,'Normal',
#                       ifelse(srat@meta.data[['annot_aug24']] == 'Tumour',paste0(srat$timePoint,'_',srat$tissue),'others'))
#   srat$cellID = srat[[cellID_column]]
#   # Subset to just the relevant cells if needed
#   mDat = srat@meta.data
# }
# 
# normalGroups = c('Normal')
# group = 'annot'
# 
# ## Import scRNA-seq metadata
# if(!is.null(mdat_fp)){
#   if(grepl('.csv',mdat_fp)){
#     mDat = read.csv(mdat_fp)
#   }else if(grepl('.tsv',mdat_fp)){
#     mDat = read.delim(mdat_fp,sep = '\t')
#   }
# }
# 
# 
# 
# 
# 
# ##-----------------------------------------------------##
# ## STEP2: Filter phCnts using single-cell annotation ####
# ##-----------------------------------------------------##
# # if(mode == 'normalDNA_only'){
# #   PDID = patientID
# #   dropImprinted=F
# #   dropUninformative=F
# #   # User provided CN-segments
# #   altChrs = c('chr8','chr13','chr21')
# #   segs = GRanges(altChrs,IRanges(c(1,3.4e7,1),c(1e9,5.3e7,1e9)))
# #   
# #   segs$matNum = c(2,1,2)
# #   segs$patNum = c(1,0,1)
# #   segs$tumFrac = segs$matNum / (segs$matNum + segs$patNum)
# #   segs$normFrac = segs$patNum / (segs$matNum + segs$patNum)
# #   names(segs) = altChrs
# #   
# # 
# # 
# #   #
# #   # chromInfo = read.delim('/lustre/scratch125/casm/team274sb/mt22/generalResources/chrom_abspos_kb.txt',sep = '\t')
# #   # chromInfo$chr = paste0('chr',chromInfo$chrom)
# #   # chromInfo$chr[chromInfo$chr == 'chr23'] = 'chrX'
# #   # chromInfo$chr = factor(chromInfo$chr,paste0('chr',c(1:22,'X')))
# #   # 
# #   # # User provided CN-segments
# #   # altChrs = chromInfo$chr
# #   # segs = GRanges(altChrs,IRanges(chromInfo$start*1e3,chromInfo$end*1e3),
# #   #                segID = paste0(chromInfo$chr,':',chromInfo$arm))
# #   # names(segs) = segs$segID
# #   # 
# #   # altChrs = chromInfo$chr
# #   # segs = GRanges(altChrs,IRanges(chromInfo$start*1e3,chromInfo$end*1e3),
# #   #                segID = paste0(chromInfo$chr,':',chromInfo$arm))
# #   # names(segs) = segs$segID
# #   # 
# #   # segs = GRanges(do.call(c,lapply(paste0('chr',c(1:22,'X')),function(i){rep(i,10)})),IRanges(rep(seq(1,9e8+1,1e8),23),rep(seq(1e8,1e9,1e8),23)))
# #   # segs$segID = paste0(as.character(seqnames(segs)),'.',rep(1:10,length(unique(seqnames(segs)))))
# # }else{
# #   PDID = names(tumourDNA)
# #   dropImprinted=T
# #   dropUninformative=T
# #   minSegLen=1e6
# #   if(is.null(segs)){
# #     if(WGS_CNV_method == 'PURPLE'){
# #       segs = processPURPLE(purple_fp=purple_fp,minSegLen = minSegLen,PDID=PDID,tgtChrs=tgtChrs,removeBalancedSegs=T,method='allelicRatio')
# #     }else if(WGS_CNV_method == 'BTB'){
# #       btb.fp = gsub('sample.dupmarked.bam$','battenberg.subclones.txt.gz',tumourDNA)
# #       if(wgs_ref_version=='hg38' & !all(grepl('^chr',tgtChrs))){
# #         tgtChrs = paste0('chr',tgtChrs)
# #       }
# #       segs = processBTB(btb.fp = btb.fp,minSegLen = minSegLen,subCl.minSegLen = subCl.minSegLen,PDID = PDID,tgtChrs=tgtChrs,
# #                         removeBalancedSegs=T,longFormat = F,keepSubClonalSegs = T,method = 'allelicRatio')
# #       
# #     }
# #   }
# #   
# # }
# 
# dropImprinted=T
# dropUninformative=T
# print(segs)
# min_normCells=100
# if(sum(mDat[[group]] %in% normalGroups) < min_normCells){
#   normalGroups = NULL
#   dropImprinted=F
# }
# 
# PDID = 'L076'
# # df = plot_BAF_byCellClusters(mDat=mDat, cellID_column=cellID_column,group = group, normalGroups = normalGroups,
# #                              outDir=outDir,minRead = 5,
# #                              patientID = patientID,PDID = PDID,
# #                              phCnts_fp=NULL,dropImprinted=dropImprinted,segs=segs,
# #                              tgtChrs=tgtChrs)
# 
# 
# 
# 
# 
# #------------------------------------------------------##
# ## Step3: Filter phCnts using single-cell annotation ####
# ##-----------------------------------------------------##
# 
# pp = calculate_scCNA_prob(patientID=patientID,outDir,PDID = PDID,segs=segs,mode=mode,
#                           mDat=mDat,cellID_column=cellID_column,group = group, normalGroups = normalGroups,
#                           dropImprinted=dropImprinted,dropUninformative=dropUninformative,
#                           normREF = normREF,aggToClust=FALSE,skipIfExists=skipIfExists)
# 
# 
# 
# 
# 
# 
# # ## Plot aggregate counts by cluster
# # gCnts_fp = file.path(outDir,paste0(PDID,'_gCnts.RDS'))
# # if(skipIfExists & file.exists(gCnts_fp)){
# #   gCnts = readRDS(gCnts_fp)
# # }
# # dat = plotRawData(gCnts,returnData=TRUE)
# #
# #
# # ## Plot heatmap of postProb
# # fp = file.path(outDir,paste0(PDID,'_pp.RDS'))
# # pp = readRDS(fp)
# # p = plotPosteriorHeatmap(pp,'postProb',plotStates = c('normFrac'))
# # print(p)
# # dev.off()
# #
# #
# # ##----  Plott heatmap of probability of each segment -------####
# # # CN states: c(pLoss = 1, mGain = 2/3, diploid = 1/2, pGain = 1/3, mLoss = 0)
# # fp_multipleStates = file.path(outDir,paste0(PDID,'_pp_multipleStates.RDS'))
# # pp_mutiple_CNstates = readRDS(fp_multipleStates)
# # pp_tmp = pp_mutiple_CNstates_manual[names(pp_mutiple_CNstates_manual) == 'chr1.2']
# # pp_tmp = pp_mutiple_CNstates_manual[names(pp_mutiple_CNstates_manual) == 'chr1.2' & pp_mutiple_CNstates_manual$clusterID == 'Tumour']
# # mtx = as.matrix(mcols(pp_tmp)[grepl('postProb',colnames(mcols(pp_tmp)))])
# # rownames(mtx) = pp_tmp$cellID
# #
# # Heatmap(mtx,
# #         show_column_dend = F,show_row_dend = F,cluster_rows = T,cluster_columns = F,
# #         show_row_names = F,show_column_names = T,km=8,
# #         #column_split = gsub('\\..*$','',colnames(mtx)),
# #         column_names_gp = gpar(fontsize=7),row_title_rot = 0,row_title_gp = gpar(fontsize=7),
# #         #row_split = pp_tmp$clusterID[match(rownames(mtx),pp_tmp$cellID)]
# # )
# #
# 
# 
# 
# 
# 
# #
# # pp_og = pp
# # pp = pp_og
# # p = plotPosteriorHeatmap(pp = pp[grepl('chr7$',names(pp))],'nLL',na_col = 'grey',cluster_rows=T)
# # print(p)
# # dev.off()
# #
# # tmp = as.data.frame(mcols(pp[grepl('chr1$|chr1\\.',names(pp))]))
# # ## For each chromosome, determine the main CN segments, then extract cells by postProb_abbFrac
# # chr = 'chr1'
# 
# ##---------------------------------------------##
# ##    Compare hSNPs matched and unmatched   ####
# ##---------------------------------------------##
# 
# # matched_hSNPs = findHetSNPs(bam = '/nfs/cancer_ref01/nst_links/live/3030/PD59701d/PD59701d.sample.dupmarked.bam',
# #                             refGenome='/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa',
# #                             outVCF = '~/lustre_mt22/alleleIntegrator_test/PD59701a/matchedAnalysis_wBTB/PD59701a_patient_hetSNPs.vcf',nParallel=nParallel)
# # BAF_lim = c(0,1)
# # minVarQual=200
# # unmatched_hSNPs = findHetSNPs_fromTumWGS(bam = tumourDNA,refGenome='/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa',
# #                                          BAF_lim=BAF_lim,minVarQual=minVarQual,
# #                                          outVCF = '~/lustre_mt22/alleleIntegrator_test/PD59701a/unmatchedAnalysis_wBTB/PD59701a_tumour_hetSNPs.vcf',
# #                                          nParallel=nParallel)
# #
# # library(UpSetR)
# # upset(fromList(list('matched_hSNPs' = names(matched_hSNPs),
# #                     'unmatched_hSNPs'= names(unmatched_hSNPs))),text.scale = 2)
# 
# 
# ## Compare phSNPs matched and unmatched
# 
# # matched_hSNPs = findHetSNPs(bam = '/nfs/cancer_ref01/nst_links/live/3030/PD59701d/PD59701d.sample.dupmarked.bam',
# #                             refGenome='/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa',
# #                             outVCF = '~/lustre_mt22/alleleIntegrator_test/PD59701a/matchedAnalysis_wBTB/PD59701a_patient_hetSNPs.vcf',nParallel=nParallel)
# # BAF_lim = c(0,1)
# # minVarQual=200
# # unmatched_hSNPs = findHetSNPs_fromTumWGS(bam = tumourDNA,refGenome='/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa',
# #                                          BAF_lim=BAF_lim,minVarQual=minVarQual,
# #                                          outVCF = '~/lustre_mt22/alleleIntegrator_test/PD59701a/unmatchedAnalysis_wBTB/PD59701a_tumour_hetSNPs.vcf',
# #                                          nParallel=nParallel)
# #
# # library(UpSetR)
# # upset(fromList(list('matched_hSNPs' = names(matched_hSNPs),
# #                     'unmatched_hSNPs'= names(unmatched_hSNPs))),text.scale = 2)
# #
# #
# #
# # ## Comparing phSNPs
# # phSNPs_matched_EM = phSNPs
# #
# # phSNPs = phaseSNPsFromCN_v2(hSNPs,segs,refGenome,tBAM = tumourDNA,useEM = F,
# #                             outPath=file.path(outDir,paste0(patientID,'_tumour_countAtHetSNPs.tsv')),
# #                             nParallel=nParallel,autoSegsFilter=autoSegsFilter)
# # phSNPs_matched_binom = phSNPs
# # phSNPs_matched_binom = phSNPs_matched_binom[match(as.character(phSNPs_matched_EM),as.character(phSNPs_matched_EM))]
# # table(phSNPs_matched_EM$informative,phSNPs_matched_binom$informative)
# #
# # upset(fromList(list('phSNPs_EM_informative' = as.character(phSNPs_matched_EM[phSNPs_matched_EM$informative]),
# #                     'phSNPs_EM_UNinformative' = as.character(phSNPs_matched_EM[!phSNPs_matched_EM$informative]),
# #                     'phSNPs_binom_informative' = as.character(phSNPs_matched_binom[phSNPs_matched_binom$informative]),
# #                     'phSNPs_binom_UNinformative' = as.character(phSNPs_matched_binom[!phSNPs_matched_binom$informative]))),text.scale = 2,nsets = 20)
# #
# # phSNPs_matched_EM$altIsMum[!phSNPs_matched_EM$informative] = 'uninformative'
# # phSNPs_matched_binom$altIsMum[!phSNPs_matched_binom$informative] = 'uninformative'
# # df = as.data.frame(table(seqnames(phSNPs_matched_EM),phSNPs_matched_EM$informative,phSNPs_matched_EM$altIsMum))
# # colnames(df) = c('chr','informative','altIsMum','nSNP')
# # df$mode = 'EM'
# #
# # df2 = as.data.frame(table(seqnames(phSNPs_matched_binom),phSNPs_matched_binom$informative,phSNPs_matched_binom$altIsMum))
# # colnames(df2) = c('chr','informative','altIsMum','nSNP')
# # df2$mode = 'binom'
# #
# #
# # df = rbind(df,df2)
# # ggplot(df,aes(mode,nSNP,fill=altIsMum))+
# #   geom_col(position = 'stack')+
# #   facet_wrap(vars(chr))+
# #   #facet_grid(mode~chr)+
# #   theme_bw()+
# #   theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
# #
# 
