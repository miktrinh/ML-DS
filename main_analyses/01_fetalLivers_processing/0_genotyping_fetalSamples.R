# Genotyping WGS + scRNAseq of:
# 1. foetal samples (2n + AK)
# 2. TAM / ML-DS samples
# 3. other leukaemias

# Install relevant packages
#install.packages('/nfs/users/nfs_m/my4/alleleIntegrator_0.7.3.tar.gz',repos = NULL,type='source')
#BiocManager::install("VariantAnnotation")
#BiocManager::install('SNPRelate')

#sudo apt-get install bcftools
#sudo cp /nfs/users/nfs_m/my4/bin/alleleCounter /usr/local/bin/alleleCounter
#sudo chmod +x /usr/local/bin/alleleCounter

outDir = "~/lustre_mt22/Aneuploidy/genotypeCheck_sep24/GenotypingResults"
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)


##-------------------##
##   Libraries     ####
##-------------------##
# Load libraries
library(alleleIntegrator)
library(GenomicFeatures)
library(tidyverse)
library(readxl)
source('/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/misc.R')


matchBAMs_notsureWhy = function (BAMs, refGenomes, snps = ExAcSNPs, outputs = NULL, 
          liftOvers = NULL, is10X = FALSE, usesChr = FALSE, nParallel = 1, 
          doPlot = TRUE, colPal = "Greens", bamGrouping = NULL, ...) 
{
  if (!is.null(bamGrouping) && length(bamGrouping) != length(BAMs)) 
    stop("Length of bamGrouping must match BAMs")
  if (is.null(names(BAMs))) 
    names(BAMs) = paste0("BAM_number", seq_along(BAMs))
  w = which(names(BAMs) == "")
  if (length(w) > 0) 
    names(BAMs)[w] = paste0("BAM_number", seq(length(w)))
  if (any(duplicated(names(BAMs)))) 
    names(BAMs) = make.unique(names(BAMs))
  if (is.null(liftOvers)) 
    liftOvers = NA
  if (length(liftOvers) == 1) 
    liftOvers = rep(liftOvers, length(BAMs))
  if (length(is10X) == 1) 
    is10X = rep(is10X, length(BAMs))
  if (length(refGenomes) == 1) 
    refGenomes = rep(refGenomes, length(BAMs))
  if (length(usesChr) == 1) 
    usesChr = rep(usesChr, length(BAMs))
  if (!is.null(outputs) && length(outputs) != length(BAMs)) 
    stop("outputs must have same length as BAMs")
  if (all(c("REF", "ALT") %in% colnames(snps))) 
    stop("snps must have REF and ALT")
  if (!is.character(snps$REF)) 
    stop("REF and ALT must be a character vector")
  if (!is.character(snps$ALT)) 
    stop("REF and ALT must be a character vector")
  refSNPs = snps[nchar(snps$REF) == 1 & nchar(snps$ALT) == 
                   1]
  refSNPs = renameSeqlevels(refSNPs, setNames(gsub("^chr", 
                                                   "", seqlevels(refSNPs)), seqlevels(refSNPs)))
  refSNPs$snpID = paste0("SNP", seq_along(refSNPs))
  if (!"AF" %in% colnames(refSNPs)) 
    refSNPs$AF = 0.5
  liftOvers[is.na(liftOvers)] = "NA"
  snpCnts = list()
  for (tgtLiftOver in unique(liftOvers)) {
    for (tgt10X in c(TRUE, FALSE)) {
      for (tgtUsesChr in c(TRUE, FALSE)) {
        for (refGenome in unique(refGenomes)) {
          w = liftOvers == tgtLiftOver & is10X == tgt10X & 
            usesChr == tgtUsesChr & refGenomes == refGenome
          if (any(w)) {
            message(sprintf("Processing %d samples", 
                            sum(w)))
            if (is.null(outputs)) {
              tgtOuts = NULL
            }
            else {
              tgtOuts = outputs[w]
            }
            if (tgtLiftOver != "NA") {
              ch = rtracklayer::import.chain(tgtLiftOver)
              tgtSNPs = unlist(rtracklayer::liftOver(refSNPs, 
                                                     ch))
            }
            else {
              tgtSNPs = refSNPs
            }
            if (tgtUsesChr) {
              tgtSNPs = renameSeqlevels(tgtSNPs, setNames(paste0("chr", 
                                                                 seqlevels(tgtSNPs)), gsub("chr", "", 
                                                                                           seqlevels(tgtSNPs))))
            }
            if (tgt10X) {
              out = alleleCounter(BAMs[w], refGenome, autoChr = F,
                                  tgtSNPs, outputs = tgtOuts, nParallel = nParallel, 
                                  ...)
            }
            else {
              out = alleleCounter(BAMs[w], refGenome, 
                                  tgtSNPs, outputs = tgtOuts, x = FALSE, 
                                  f = 3, F = 3852, m = 20, q = 35, nParallel = nParallel, 
                                  ...)
            }
            if (tgt10X) {
              out = lapply(out, function(e) {
                tmp = mcols(e)
                tmp$mark = as.character(e)
                tmp = aggregate(cbind(A, C, G, T, Tot) ~ 
                                  mark, data = tmp, FUN = sum)
                m = match(tmp$mark, as.character(e))
                e = e[m]
                e$barcode = NULL
                e$A = tmp$A
                e$C = tmp$C
                e$G = tmp$G
                e$T = tmp$T
                e$Tot = tmp$Tot
                e
              })
            }
            for (i in seq_along(out)) {
              snpCnts[[(names(BAMs)[w])[i]]] = out[[i]]
            }
          }
        }
      }
    }
  }
  cnts = matrix(NA, nrow = length(refSNPs), ncol = length(BAMs), 
                dimnames = list(refSNPs$snpID, names(BAMs)))
  altCnts = cnts
  for (i in seq_along(snpCnts)) {
    jj = match(names(snpCnts)[i], colnames(cnts))
    ii = match(snpCnts[[i]]$snpID, rownames(cnts))
    bases = c("A", "C", "G", "T")
    tmp = mcols(snpCnts[[i]])
    tmp = as.matrix(tmp[, bases])
    cnts[ii, jj] = rowSums(tmp)
    altCnts[ii, jj] = tmp[cbind(seq_along(ii), match(snpCnts[[i]]$ALT, 
                                                     bases))]
  }
  pAA = dbinom(altCnts, cnts, 0.1, log = TRUE)
  pAa = dbinom(altCnts, cnts, 0.5, log = TRUE)
  paa = dbinom(altCnts, cnts, 0.9, log = TRUE)
  pp = data.frame(pAA = as.vector(pAA), pAa = as.vector(pAa), 
                  paa = as.vector(paa))
  pp = exp(pp)/rowSums(exp(pp), na.rm = TRUE)
  pp$sampleID = rep(colnames(cnts), each = nrow(cnts))
  pp$snpID = rep(rownames(cnts), ncol(cnts))
  pp$snpAF = rep(refSNPs$AF, ncol(cnts))
  pp$numAlt = rep(3, nrow(pp))
  pp$numAlt[pp$paa > 0.9] = 2
  pp$numAlt[pp$pAa > 0.9] = 1
  pp$numAlt[pp$pAA > 0.9] = 0
  gtMat = matrix(3, nrow = length(refSNPs), ncol = length(BAMs), 
                 dimnames = list(refSNPs$snpID, names(BAMs)))
  gtMat[cbind(match(pp$snpID, rownames(gtMat)), match(pp$sampleID, 
                                                      colnames(gtMat)))] = pp$numAlt
  gdsFile = tempfile()
  snpgdsCreateGeno(gdsFile, genmat = gtMat, sample.id = colnames(gtMat), 
                   snp.id = rownames(gtMat), snp.chromosome = as.character(seqnames(refSNPs)), 
                   snp.position = start(refSNPs), snp.allele = paste0(refSNPs$ALT, 
                                                                      "/", refSNPs$REF), snpfirstdim = TRUE, )
  gds = snpgdsOpen(gdsFile)
  gdsSub = snpgdsLDpruning(gds, ld.threshold = 0.2)
  ibs = snpgdsIBS(gds, num.thread = nParallel)
  rownames(ibs$ibs) = colnames(ibs$ibs) = names(BAMs)
  if (doPlot) {
    colFun = suppressWarnings(brewer.pal(100, colPal))
    colFun = colorRamp2(seq(0.5, 1, length.out = length(colFun)), 
                        colFun)
    hm = Heatmap(ibs$ibs, col = colFun, name = "IBS", show_row_names = TRUE, 
                 show_column_names = TRUE, show_row_dend = FALSE, 
                 show_column_dend = FALSE, row_title_rot = 0, column_split = bamGrouping, 
                 row_split = bamGrouping)
    draw(hm)
  }
  return(list(ibs = ibs, pp = pp, snpCnts = snpCnts, refSNPs))
}



##----------------------------##
##   Set Global parameters  ####
##----------------------------##

refGenome_hg38 = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa'
refGenome_hg19 = '/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh37d5/genome.fa'
refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
liftChain = '/lustre/scratch125/casm/team274sb/mt22/hg19ToHg38_noChr.over.chain'

nParallel=38
skipIfExists = TRUE
#nMaxSim = 60 # Max number of samples being processed at once
plotDir = outDir



# ##-------------------------------##
# ##   Import list of bam files  ####
# ##-------------------------------##
# 
# ##---   Foetal samples
# projMani = read_excel("/lustre/scratch125/casm/team274sb/mt22/projectManifest.xlsx",sheet = "Aneuploidy_mani")
# projMani = projMani[!is.na(projMani$donorID),] %>% filter(assay != 'Hi-C')
# 
# #----- RNA BAMs
# projMani.rna = projMani[!projMani$assay %in% c('WGS','t-NanoSeq'),]
# cellrangerID = basename(projMani.rna$`irods location`)
# bams10X = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data',cellrangerID,'possorted_genome_bam.bam')
# names(bams10X) = paste0('f',projMani.rna$Tissue,'_',projMani.rna$sangerSampleID)
# length(bams10X) == n_distinct(bams10X)
# # Check that each bams10X file has a unique name
# if(length(unique(names(bams10X))) != length(bams10X)){
#   stop(sprintf('Duplicated bams10X names detected: %s',names(bams10X)[duplicated(names(bams10X))]))
# }
# 
# bams10X_out = file.path(outDir,paste0(names(bams10X),'_genotypeCheck.tsv'))
# w = file.exists(bams10X_out)
# bams10X = bams10X[!w]
# table(file.exists(bams10X))
# 
# 
# # Download bam files from irods
# existing_bams = list.files('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data',pattern = 'bam$',recursive = T,full.names = T)
# names(existing_bams) = gsub('_GRCh38-2020-A|_GRCh38-1_2_0','',gsub('.*_MY_','MY_',basename(dirname(existing_bams))))
# names(existing_bams) = paste0('f',projMani.rna$Tissue[match(names(existing_bams),projMani.rna$sangerSampleID)],'_',names(existing_bams))
# table(names(bams10X) %in% names(existing_bams))
# 
# # for(b in bams10X){
# #   existingResult_fp = file.path(outDir,paste0(c(names(bams10X)[bams10X == b]),'_genotypeCheck.tsv'))
# #   bam_fp = b
# #   bai_fp = paste0(b,'.bai')
# #   system(sprintf('rm %s',bai_fp))
# #   if((!file.exists(bam_fp)) | (!file.exists(bai_fp))){
# #     print(b)  
# #     
# #     sangerID = gsub('f.*_MY','MY',names(bams10X)[bams10X == b])
# #     irod_path = projMani.rna$`irods location`[projMani.rna$sangerSampleID == sangerID]
# #     cellrangerID = basename(irod_path)
# #     d = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data',cellrangerID)
# #     if(!dir.exists(d)){
# #       dir.create(d,recursive = T)
# #     }
# #     
# #     o1 = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data',cellrangerID)
# #     o2 = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data',cellrangerID)
# #     if(!file.exists(bam_fp)){
# #       print(bam_fp)
# #       system(sprintf('irods iget %s/possorted_genome_bam.bam %s/possorted_genome_bam.bam',irod_path,o1))  
# #     }
# #     if(!file.exists(bai_fp)){
# #       print(bai_fp)
# #       #print(sprintf('irods iget %s/possorted_genome_bam.bam.bai %s/possorted_genome_bam.bam.bai',irod_path,o2))
# #       system(sprintf('irods iget %s/possorted_genome_bam.bam.bai %s/possorted_genome_bam.bam.bai',irod_path,o2))
# #     }
# #   #}
# #   # if((!file.exists(existingResult_fp)) & (!file.exists(b))){
# #   #   print(b)  
# #   #   
# #   #   sangerID = gsub('f.*_MY','MY',names(bams10X)[bams10X == b])
# #   #   irod_path = projMani.rna$`irods location`[projMani.rna$sangerSampleID == sangerID]
# #   #   cellrangerID = basename(irod_path)
# #   #   d = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data',cellrangerID)
# #   #   if(!dir.exists(d)){
# #   #     dir.create(d,recursive = T)
# #   #   }
# #   # 
# #   #   o1 = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data',cellrangerID)
# #   #   o2 = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data',cellrangerID)
# #   #   system(sprintf('irods iget %s/possorted_genome_bam.bam %s/possorted_genome_bam.bam',irod_path,o1))
# #   #   system(sprintf('irods iget %s/possorted_genome_bam.bam %s/possorted_genome_bam.bam.bai',irod_path,o2))
# #   #   print('done downloading')
# #   }else{
# #     print('Next')
# #   }
# #   
# # }
# 
# #----- DNA BAMs 
# dnaBAMs = c(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2623',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T),
#             list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2810',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T),
#             list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3174',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T),
#             list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3374',recursive = T,pattern = 'sample.merged.bam$',full.names = T)) # targeted nanoseq
# names(dnaBAMs) = basename(dirname(dnaBAMs))
# 
# ## Cross check with sample manifest
# projMani$PDID[projMani$assay == 't-NanoSeq'] = paste0(projMani$PDID[projMani$assay == 't-NanoSeq'],'_tds0001')
# projMani.dna = projMani[projMani$assay %in% c("WGS",'t-NanoSeq') & !projMani$PDID %in% c("PD51325d", "PD51609b",'PD51610d','PD57607c','PD57607b','PD52544b'),]
# #projMani.dna$PDID[projMani.dna$assay == 't-NanoSeq'] = paste0(projMani.dna$PDID[projMani.dna$assay == 't-NanoSeq'],'_tds0001')
# table(projMani.dna$PDID %in% names(dnaBAMs))
# table(projMani.dna$PDID[!projMani.dna$PDID %in% names(dnaBAMs)])
# table(names(dnaBAMs) %in% projMani.dna$PDID)
# names(dnaBAMs)[!names(dnaBAMs) %in% projMani.dna$PDID]
# 
# dna_refGenomes =  rep(c(refGenome_hg19,refGenome_hg38,refGenome_hg38,refGenome_hg19),
#                       c(length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2623',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                         length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2810',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                         length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3174',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                         length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3374',recursive = T,pattern = 'sample.merged.bam$',full.names = T))))
# 
# dna_liftChain = rep(c(NA,liftChain,liftChain,NA),
#                     c(length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2623',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                       length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2810',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                       length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3174',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                       length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3374',recursive = T,pattern = 'sample.merged.bam$',full.names = T))))
# 
# 
# # dnaBAMs_out = file.path(outDir,paste0(names(dnaBAMs),'_genotypeCheck.tsv'))
# # table(file.exists(dnaBAMs_out))
# # dnaBAMs_out[!file.exists(dnaBAMs_out)]
# # w = file.exists(dnaBAMs_out) & file.exists(dnaBAMs)
# # dnaBAMs = dnaBAMs[w]
# # dna_refGenomes = dna_refGenomes[w]
# 
# # Check that each dnaBAM file has a unique name
# if(length(unique(names(dnaBAMs))) != length(dnaBAMs)){
#   stop(sprintf('Duplicated bams10X names detected: %s',names(dnaBAMs)[duplicated(names(dnaBAMs))]))
# }
# 
# 
# # Specify useChr based on the following rules:
# # hg19: FALSE (no 'chr' prefix)
# # hg38: TRUE (yes 'chr' prefix)
# # EXCEPTION: samples with existing alleleCounter results - FALSE. This is because alleleCounter does not print 'chr' prefix regardless of the reference genome
# #            hg38 without autoChr check c
# # usesChr = c(rep(FALSE,length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2623',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T))), #hg19 WGS
# #             rep(TRUE,length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2810',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T))), #hg38 WGS
# #             rep(TRUE,length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3174',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T))), #hg38 WGS
# #             rep(FALSE,length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3374',recursive = T,pattern = 'sample.merged.bam$',full.names = T))), #hg19 t-NS
# #             ifelse(file.exists(file.path(outDir,paste0(names(bams10X),'_genotypeCheck.tsv'))),FALSE,TRUE))
# 
# #############################
# # Check genotype consistency
# #Are all the BAMs you're going to use from the same individual?  Check before you start
# genoCheck = matchBAMs(BAMs = c(dnaBAMs,bams10X),doPlot = F,
#                       refGenomes = c(dna_refGenomes,rep(refGenome10X,length(bams10X))),
#                       outputs = file.path(outDir,paste0(c(names(dnaBAMs),names(bams10X)),'_genotypeCheck.tsv')),
#                       liftOvers=c(dna_liftChain,rep(liftChain,length(bams10X))),
#                       is10X=rep(c(FALSE,TRUE),c(length(dnaBAMs),length(bams10X))),
#                       usesChr = rep(c(FALSE,TRUE,TRUE,FALSE,TRUE),
#                                     c(length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2623',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                                       length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/2810',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                                       length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3174',recursive = T,pattern = 'sample.dupmarked.bam$',full.names = T)),
#                                       length(list.files('/lustre/scratch124/casm/team78pipelines/nst_links/live/3374',recursive = T,pattern = 'sample.merged.bam$',full.names = T)),
#                                       length(bams10X))),
#                       nParallel=nParallel,skipIfExists=skipIfExists,nMaxSim=12,nChunks=6)
# 
# 
# 
# #If anything is less than 0.8 and you should be concerned...
# message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))
# if(!is.null(plotDir)){
#   library(RColorBrewer)
#   library(circlize)
#   library(ComplexHeatmap)
#   bamGrouping = NULL
#   colPal = 'Greens'
#   colFun = suppressWarnings(brewer.pal(100, colPal))
#   colFun = colorRamp2(seq(0.5, 1, length.out = length(colFun)),
#                       colFun)
#   hm = Heatmap(genoCheck$ibs$ibs, col = colFun, name = "IBS", show_row_names = TRUE,
#                show_column_names = TRUE, show_row_dend = FALSE,
#                show_column_dend = FALSE, row_title_rot = 0, column_split = bamGrouping,
#                row_split = bamGrouping)
#   pdf(paste0(plotDir,'/genotypeCheck_sep24.pdf'),width = 25,height = 20)
#   draw(hm)
#   dev.off()
# 
# }


##----------------------------------##
##          ML-DS dataset         ####
##----------------------------------##
###### GENOTYPING ######

#############
# Libraries #
#############
# Load libraries
library(tidyverse)
library(alleleIntegrator)
library(readxl)

refGenome10X = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa'

#outDir = file.path('~/lustre_mt22/Aneuploidy/genotypeCheck_sep24/GenotypingResults')
outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/genotypeCheck/mar23_v2'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

plotDir = outDir
nParallel=1
skipIfExists = TRUE

#### Get list of RNA and DNA BAM files ####
library(readxl)
projMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'MLDS_GOSH')
projMani = projMani[!is.na(projMani$donorID),]
projMani = projMani[,c('donorID','physical_sampleID','PDID','sangerSampleID','assay','Disease','Tissue','Point in treatment')]

projMani_otherLeuk = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'GOSH_others_included')
projMani_otherLeuk = projMani_otherLeuk[!is.na(projMani_otherLeuk$donorID),]
projMani_otherLeuk = projMani_otherLeuk[,c('donorID','physical_sampleID','PDID','sangerSampleID','assay','Disease','Tissue','Point in treatment')]

projMani = rbind(projMani,projMani_otherLeuk)
projMani = projMani[!(is.na(projMani$PDID) & is.na(projMani$sangerSampleID)) & 
                      projMani$assay %in% c("5' V2 Dual Index",'GEX','t-NanoSeq','WGS'),]
projMani$sangerSampleID[(projMani$sangerSampleID == 'NA' | is.na(projMani$sangerSampleID)) & projMani$assay %in% c('t-NanoSeq','WGS')] = projMani$PDID[(projMani$sangerSampleID == 'NA' | is.na(projMani$sangerSampleID)) & projMani$assay %in% c('t-NanoSeq','WGS')]
projMani = projMani[projMani$sangerSampleID != '???',]

#----- RNA BAMs
existing_results = list.files('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/genotypeCheck/mar23_v2',pattern = '_genotypeCheck.tsv')
names(existing_results) = gsub('_genotypeCheck.tsv','',existing_results)
names(existing_results)[grepl('^MY\\d+',names(existing_results))] = gsub('^MY','MY.',names(existing_results)[grepl('^MY\\d+',names(existing_results))])
projMani$sangerSampleID_v2 = gsub('_','.',gsub('^.*SB_|^.*ALeuk_','',projMani$sangerSampleID))
projMani$sangerSampleID_v2[grepl('tds',projMani$sangerSampleID_v2)] = gsub('\\.','_',projMani$sangerSampleID_v2[grepl('tds',projMani$sangerSampleID_v2)])
table(projMani$sangerSampleID_v2[!is.na(projMani$sangerSampleID_v2)] %in% names(existing_results))
#table(projMani$sangerSampleID_v2[!is.na(projMani$sangerSampleID_v2) & !projMani$sangerSampleID_v2 %in% names(existing_results)])

# RNA BAMs - Look for these bam files on irods and download them...

rnaBams_toDownload = projMani[!is.na(projMani$sangerSampleID_v2) & !projMani$sangerSampleID_v2 %in% names(existing_results) & 
                                projMani$assay %in% c("5' V2 Dual Index",'GEX'),]
bams10X = list.files('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data',recursive = T,pattern = 'bam$',full.names = T)
names(bams10X) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_|^.*_ALeuk_','',bams10X))
names(bams10X) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_MY_','MY_',names(bams10X)))
names(bams10X) = gsub('_','.',names(bams10X))
# bams10X = bams10X[names(bams10X) %in% rnaBams_toDownload$sangerSampleID_v2]
# table(names(bams10X) %in% rnaBams_toDownload$sangerSampleID_v2)
# table(rnaBams_toDownload$sangerSampleID_v2 %in% names(bams10X))
# rnaBams_toDownload$sangerSampleID_v2[!rnaBams_toDownload$sangerSampleID_v2 %in% names(bams10X)]
# names(bams10X)[!names(bams10X) %in% rnaBams_toDownload$sangerSampleID_v2]
# bams10X=c()
# for(s in unique(rnaBams_toDownload$sangerSampleID)){
#   irod_path = system(sprintf('irods imeta qu -z seq -C sample = %s',s),intern = T)
#   irod_path = gsub('collection: ','',irod_path[grepl('cellranger700_count',irod_path)])
#   if(length(irod_path) == 0){
#     stop(sprintf('Cannot find irods path for sample %s',s))
#   }
#   
#   cellrangerID = basename(irod_path)
#   d = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data',cellrangerID)
#   #check if bam file exist
#   bam_fp = file.path(d,'possorted_genome_bam.bam')
#   
#   if(!file.exists(bam_fp)){
#     bams10X = c(bams10X,bam_fp)
#     if(!dir.exists(d)){
#       dir.create(d,recursive = T)
#     }
#     
#     o1 = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data',cellrangerID)
#     #o2 = file.path('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data',cellrangerID)
#     system(sprintf('irods iget %s/possorted_genome_bam.bam %s/possorted_genome_bam.bam',irod_path,o1))
#     system(sprintf('irods iget %s/possorted_genome_bam.bam.bai %s/possorted_genome_bam.bam.bai',irod_path,o1))
#     print('done downloading')  
#   }
#   
# }


# names(bams10X) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_SB_|^.*_ALeuk_','',bams10X))
# names(bams10X) = gsub('_GRCh38-2020-A.*$','',gsub('^.*_MY_','MY_',names(bams10X)))
bams10X = bams10X[file.exists(bams10X)]
print(table(file.exists(bams10X)))



dnaBams = projMani[!is.na(projMani$sangerSampleID_v2) & #!projMani$sangerSampleID_v2 %in% names(existing_results) & 
                     projMani$assay %in% c('WGS','t-NanoSeq'),]
dnaBams$sangerSampleID[dnaBams$sangerSampleID == 'NA'] = dnaBams$PDID[dnaBams$sangerSampleID == 'NA']

dnaBAMs = c(file.path('/nfs/cancer_ref01/nst_links/live/3030',dnaBams$sangerSampleID[dnaBams$assay == 'WGS' & !dnaBams$sangerSampleID %in% c('PD64665a','PD64666a')],paste0(dnaBams$sangerSampleID[dnaBams$assay == 'WGS' & !dnaBams$sangerSampleID %in% c('PD64665a','PD64666a')],'.sample.dupmarked.bam')),
            '/nfs/cancer_ref01/nst_links/live/3484/PD64665a/PD64665a.sample.dupmarked.bam',
            '/nfs/cancer_ref01/nst_links/live/3484/PD64666a/PD64666a.sample.dupmarked.bam',
            file.path('/nfs/cancer_ref01/nst_links/live/3327',dnaBams$sangerSampleID[dnaBams$assay == 't-NanoSeq'],paste0(dnaBams$sangerSampleID[dnaBams$assay == 't-NanoSeq'],'.sample.merged.bam')),
            '/nfs/cancer_ref01/nst_links/live/3309/PD63417a/PD63417a.sample.dupmarked.bam')
table(file.exists(dnaBAMs))
dnaBAMs = dnaBAMs[file.exists(dnaBAMs)]
names(dnaBAMs) = basename(dirname(dnaBAMs))

genoCheck = matchBAMs(BAMs = c(dnaBAMs,bams10X),doPlot = F,
                      refGenomes = rep(c(refGenome_hg38,refGenome10X),c(length(dnaBAMs),length(bams10X))),
                      outputs = file.path(outDir,paste0(c(names(dnaBAMs),names(bams10X)),'_genotypeCheck.tsv')),
                      liftOvers=rep(c(liftChain,liftChain),c(length(dnaBAMs),length(bams10X))),
                      is10X=rep(c(FALSE,TRUE),c(length(dnaBAMs),length(bams10X))),
                      nParallel=nParallel,skipIfExists=skipIfExists,nChunks = 24, nMaxSim = 3)


#If anything is less than 0.8 and you should be concerned...
message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))
if(!is.null(plotDir)){
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  bamGrouping = NULL
  colPal = 'Greens'
  colFun = suppressWarnings(brewer.pal(100, colPal))
  colFun = colorRamp2(seq(0.5, 1, length.out = length(colFun)),
                      colFun)
  hm = Heatmap(genoCheck$ibs$ibs, col = colFun, name = "IBS", show_row_names = TRUE,
               show_column_names = TRUE, show_row_dend = FALSE,
               show_column_dend = FALSE, row_title_rot = 0, column_split = bamGrouping,
               row_split = bamGrouping)
  pdf(paste0(plotDir,'/genotypeCheck_Leuk_oct24.pdf'),width = 45,height = 45)
  draw(hm)
  dev.off()


}

# #############################
# # Check genotype consistency
# #Are all the BAMs you're going to use from the same individual?  Check before you start
# snps = rowRanges(vcf.pass)
# # convert ALT from DNAStringSetList to Character
# snps$ALT = unlist(unstrsplit(CharacterList(snps$ALT), sep = ','))
# snps$REF = as.character(snps$REF)
# genoCheck = matchBAMs(BAMs = c(bams10X),snps = snps,
#                       refGenomes = rep(c(refGenome10X),c(length(bams10X))),
#                       outputs = file.path(outDir,paste0(c(names(bams10X)),'_cavemanSNV.tsv')),
#                       liftOvers=rep(c(NA),c(length(bams10X))),
#                       is10X=rep(c(TRUE),c(length(bams10X))),
#                       nParallel=nParallel,nMaxSim=20,nChunks=4,skipIfExists=skipIfExists)
#
#
#
# #If anything is less than 0.8 and you should be concerned...
# message(sprintf("The minimum similarity found was %g",min(genoCheck$ibs$ibs)))
# if(!is.null(plotDir)){
#   library(RColorBrewer)
#   library(circlize)
#   library(ComplexHeatmap)
#   bamGrouping = NULL
#   colPal = 'Greens'
#   colFun = suppressWarnings(brewer.pal(100, colPal))
#   colFun = colorRamp2(seq(0.5, 1, length.out = length(colFun)),
#                       colFun)
#   hm = Heatmap(genoCheck$ibs$ibs, col = colFun, name = "IBS", show_row_names = TRUE,
#                show_column_names = TRUE, show_row_dend = FALSE,
#                show_column_dend = FALSE, row_title_rot = 0, column_split = bamGrouping,
#                row_split = bamGrouping)
#   pdf(paste0(plotDir,'/cavemanSNV.pdf'),width = 15,height = 12)
#   draw(hm)
#   dev.off()
#
# }
