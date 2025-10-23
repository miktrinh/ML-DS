## Helper functions for processing and gene_module scoring in bulk RNA-seq dataset

#' Import and preprocess bulk RNA-seq data
#'
#' Loads raw count and TPM data for bulk RNA-seq samples, optionally adds Oxford data, 
#' removes specified samples, applies expression filtering, and returns processed matrices.
#'
#' @param gns Optional gene annotation (GRanges). If NULL, uses the provided GTF path.
#' @param rm_henning_healthy Logical, whether to exclude healthy control samples.
#' @param oxfordBulk Logical, whether to include Oxford cohort data.
#' @param filter_lowExpr_gene Logical, whether to filter out lowly expressed genes.
#' @param plot_filter Logical, whether to show histogram of TPM values.
#' @param tpm_path Character vector of file paths to pre-computed TPM CSVs (if available).
#' @param gtf_path Path to GTF file used to generate gene annotation.
#' @param raw_counts_xlsx Path to Excel file containing raw count data.
#' @param annotation_xlsx Path to Excel file with sample metadata.
#' @param survival_xlsx Path to Excel file with survival annotations.
#' @param oxford_meta_fp Path to Oxford sample metadata file.
#' @param oxford_tpm_path Path to Oxford TPM CSV.
#'
#' @return A list containing:
#' \describe{
#'   \item{tpm}{Data frame of TPM values (merged if Oxford data included)}
#'   \item{raw}{Matrix of raw counts}
#'   \item{cpm}{Matrix of CPM values}
#'   \item{metadata}{Data frame of sample metadata}
#'   \item{gene_map}{Data frame with gene metadata}
#' }
#' @export

import_bulk_data <- function(
  gns = NULL,
  rm_henning_healthy = FALSE,
  oxfordBulk = TRUE,
  filter_lowExpr_gene = TRUE,
  plot_filter = FALSE,
  tpm_fp = c("~/Dropbox (Partners HealthCare)/project__T21TAMMLDS/rnaseq/bulk_rnaseq/rawdata/StarSalmon_tpm.csv"), 
  gtf_path = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf',
  
  raw_counts_xlsx = '~/lustre_mt22/MLDS_scRNAseq/bulkRNA_Data/Down_Leukemia_Klusmann_Apr23.xlsx',
  annotation_xlsx = '~/lustre_mt22/MLDS_scRNAseq/bulkRNA_Data/Down_Leukemia_Klusmann_Apr23.xlsx',
  survival_xlsx = '~/lustre_mt22/MLDS_scRNAseq/bulkRNA_Data/Hennings_bulkRNAseq_metadata.xlsx',
  oxford_meta_fp = '~/lustre_mt22/Down_Leukemia/Oxford_flowSorted_Blineage_bulkRNA/GSE122982_series_matrix.txt',
  oxford_tpm_fp = '~/lustre_mt22/Down_Leukemia/Oxford_flowSorted_Blineage_bulkRNA/GSE122982_FBM_TPM.csv'
) {
  
  ## ---------- 1. Load and check inputs ----------
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) stop("Please install the GenomicFeatures package.")
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Please install the readxl package.")
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Please install the edgeR package.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Please install the tibble package.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install the dplyr package.")
  
  # GTF fallback
  if (is.null(gns)) {
    if (!file.exists(gtf_path)) stop("GTF file does not exist at given path.")
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path)
    gns <- GenomicFeatures::genes(txdb)
  }
  
  ## ---------- 2. Read gene counts + geneMap ----------
  rawCnt <- readxl::read_excel(raw_counts_xlsx, sheet = 'StarSalmon_Counts')
  bulk_geneMap <- rawCnt[, c('gene_id', 'gene_name', 'EffectiveLength')]
  m <- match(bulk_geneMap$gene_id, gns$gene_id)
  bulk_geneMap$chr <- 'unknown'
  bulk_geneMap$chr[!is.na(m)] <- as.character(GenomicRanges::seqnames(gns[m[!is.na(m)]]))
  
  ## ---------- 3. Sample metadata ----------
  message('Importing bulk samples metadata....')
  mdat <- readxl::read_excel(annotation_xlsx, sheet = 'Annotation')
  mdat <- mdat[, c('Sample', 'Source', 'Subgroup')]
  mdat$Sample <- gsub('-', '.', mdat$Sample)
  
  # Remove healthy control from Hennings data
  if(rm_henning_healthy){
    mdat <- mdat[!grepl('healthy', mdat$Sample), ]
  }
  
  # Exclude specific problematic samples
  bad_samples <- c('Patient_TMD_1020','Patient_TMD_1037','Patient_TMD_1303','Patient_TMD_BCKR',
                   'Patient_TMD_KHLR','Patient_MLL_1569','Patient_MLL_1569_2','Patient_MLL_615',
                   'Patient_MLL_DEE7','Patient_MLL_DEEF3','Patient_MLL_DEHS1','Patient_MLL_DEMK11','Patient_MLL_DES10',
                   'Patient_MLL_642')
  mdat <- mdat[!mdat$Sample %in% bad_samples, ]
  
  # Read survival data
  survival_mdat <- readxl::read_excel(survival_xlsx)
  survival_mdat$Sample <- gsub('-', '.', survival_mdat$Sample)
  table(survival_mdat$Sample %in% mdat$Sample)
  survival_mdat <- survival_mdat[survival_mdat$Sample %in% mdat$Sample, ]
  # # Remove samples with high duplication rate
  # mdat = mdat[mdat$Sample %in% survival_mdat$Sample[as.numeric(survival_mdat$Duplication_rate) < 0.65 | survival_mdat$Duplication_rate == '-'],]
  
  # Merge metadata
  mdat <- merge(mdat, survival_mdat[, !grepl('\\.\\.', colnames(survival_mdat))], by = c('Sample', 'Source', 'Subgroup'), all = TRUE)
  
  ## ---------- 4. Process raw counts + Compute bulkRNAdata CPM counts ----------
  rawCnt <- tibble::column_to_rownames(rawCnt, var = 'gene_id')
  rawCnt <- rawCnt[, colnames(rawCnt) %in% mdat$Sample]
  ## Filter out lowly expressed genes
  bulk_dge <- edgeR::DGEList(counts = rawCnt, genes = rownames(rawCnt))
  cpmCnt <- edgeR::cpm(bulk_dge)
  
  
  ## ---------- 5. Compute or import TPM ----------
  message('Importing TPM count....')
  
  if (filter_lowExpr_gene && !oxfordBulk) {
    prop_expressed <- rowMeans(edgeR::cpm(bulk_dge) > 1)
    keep <- prop_expressed > 0.5
    
    if (plot_filter) {
      op <- par(no.readonly = TRUE)
      par(mfrow = c(1, 2))
      hist(edgeR::cpm(bulk_dge, log = TRUE), main = 'Unfiltered', xlab = 'logCPM')
      abline(v = log(1), lty = 2, col = 2)
      hist(edgeR::cpm(bulk_dge[keep, ], log = TRUE), main = 'Filtered', xlab = 'logCPM')
      abline(v = log(1), lty = 2, col = 2)
      par(op)
    }
    
    ## Subset the count matrix
    bulk_dge <- bulk_dge[keep, , keep.lib.sizes = FALSE]
    rawCnt <- rawCnt[keep, ]
    
    ##---- Calculate TPM ----##
    
    geneLen <- as.vector(bulk_geneMap$EffectiveLength[match(rownames(rawCnt), bulk_geneMap$gene_id)])
    rpk <- apply(rawCnt[,!colnames(rawCnt) %in% c('Length','gene_name','Geneid')], 2, function(x) x / (as.vector(geneLen) / 1000))
    #normalize by the sample size using rpk values
    tpm <- apply(rpk, 2, function(x) x / (sum(as.numeric(x)) / 1e6)) %>% as.data.frame()
    tpm$ensID <- rownames(rawCnt)
    tpm$geneLength <- bulk_geneMap$EffectiveLength[match(tpm$ensID,bulk_geneMap$gene_id)]
    tpm_count <- tpm[, c('ensID', 'geneLength', setdiff(colnames(tpm), c('ensID', 'geneLength')))]
  } else {
    ## Import Hennings TPM count
    tpm_count <- data.frame()
    for (file in tpm_fp) {
      if (!file.exists(file)) stop(paste("TPM file not found:", file))
      tpm <- read.csv(file, row.names = 1)
      colnames(tpm) <- gsub('-', '.', colnames(tpm))
      colnames(tpm) <- gsub('\\.markdup\\.sorted\\.bam', '', colnames(tpm))
      rownames(tpm) <- tpm$ensID
      if (ncol(tpm_count) == 0) {
        tpm_count <- tpm
      } else {
        m <- match(colnames(tpm_count), colnames(tpm))
        print(sum(is.na(m)))
        tpm_count <- rbind(tpm_count, tpm[, m])
      }
    }
    tpm_count <- tpm_count[, !grepl('^X\\.', colnames(tpm_count))]
  }
  
  # Validate sample match
  if (!all(mdat$Sample %in% colnames(tpm_count))) {
    stop("Some samples in metadata are missing from TPM count matrix.")
  }
  
  # subset TMP count to include only samples in metadata
  tpm_count <- tpm_count[, c('ensID', mdat$Sample)]
  
  
  
  ## ---------- 6. Merge Oxford Bulk ----------
  if (oxfordBulk) {
    if (!file.exists(oxford_meta_fp) || !file.exists(oxford_tpm_fp)) {
      stop("Oxford metadata or TPM file not found.")
    }
    ## 1. Meta data
    mdat_ox <- read.delim(oxford_meta_fp, skip = 28, sep = '\t')
    mdat_ox <- as.data.frame(t(mdat_ox[c(10,11), -1]))
    colnames(mdat_ox) <- c('celltype', 'uniqueID')
    mdat_ox$uniqueID = gsub('unique identifier: ','',mdat_ox$uniqueID)
    mdat_ox$Source = gsub('population: ','',mdat_ox$celltype)
    mdat_ox$Subgroup = mdat_ox$Source
    mdat_ox$Sample = paste0(mdat_ox$Source,'_',mdat_ox$uniqueID)
    ## Merge metadata                                    
    mdat_ox[, setdiff(names(mdat), names(mdat_ox))] <- NA
    mdat <- rbind(mdat, mdat_ox[, colnames(mdat)])
    
    ## 2. TPM count
    
    ## Import Oxford bulk
    oxford_bulk <- read.csv(oxford_tpm_fp)
    colnames(oxford_bulk)[1] <- 'ensID'
    colnames(oxford_bulk)[-1] <- mdat_ox$Sample[match(colnames(oxford_bulk)[-1], mdat_ox$uniqueID)]
    tpm_count <- merge(tpm_count, oxford_bulk, by = 'ensID')
    tpm_count <- tibble::column_to_rownames(tpm_count, 'ensID')
  }
  
  
  
  return(list(
    tpm_count = tpm_count,
    raw_count = rawCnt,
    cpm_count = cpmCnt,
    mdat = mdat,
    bulk_geneMap = bulk_geneMap
  ))
}











#' Compute module scores for bulk RNA-seq using singscore
#'
#' Takes gene modules and scores TPM data (imported via import_bulk_data),
#' adding metadata, OS grouping, and moduleType annotation.
#'
#' @param module_list Named list of gene sets (vectors). Should include "all" for combined sets.
#' @param module_type Character string: one of "GATA1s_topGenes", "MLDS_topGenes", or "all".
#' @param gns GRanges gene annotation; passed to import_bulk_data.
#' @param rm_henning_healthy Logical; remove healthy control samples.
#' @param oxfordBulk Logical; include Oxford bulk data.
#' @param filter_lowExpr_gene Logical; whether to filter low-expressed genes.
#' @param tpm_fp Character vector of TPM file paths.
#' @param oxford_meta_fp Oxford metadata file path.
#' @param oxford_tpm_fp Oxford TPM counts file path.
#' @param plot_filter Logical; enable filtering histogram plotting.
#' @param min_genes Minimum number of genes required in each gene module (default = 5)
#' @return Data frame with module scores and sample metadata.
#' @importFrom singscore rankGenes simpleScore
#' @export
compute_bulk_module_scores <- function(
  moduleList,
  module_type = c("GATA1s_topGenes", "MLDS_topGenes","MLDS_imprints_T21FL", "all"),
  gns = NULL,
  rm_henning_healthy = FALSE,
  oxfordBulk = FALSE,
  filter_lowExpr_gene = TRUE,
  tpm_fp = NULL,
  oxford_meta_fp = NULL,
  oxford_tpm_fp = NULL,
  plot_filter = TRUE,
  min_genes=5
) {
  
  # Load required packages
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' required but not installed.")
  if (!requireNamespace("singscore", quietly = TRUE)) stop("Package 'singscore' required but not installed.")
  
  ## ----- Helper function: get gene sets based on module_type -----
  get_module_gene_sets <- function(module_list, module_type) {
    if (module_type == "GATA1s_topGenes") {
      if (!all(c("TAM.MLDS.up", "TAM.MLDS.down") %in% names(module_list))) {
        stop("Module list must contain 'TAM.MLDS.up' and 'TAM.MLDS.down' for GATA1s_topGenes.")
      }
      return(list(
        up = module_list[["TAM.MLDS.up"]],
        down = module_list[["TAM.MLDS.down"]]
      ))
    } else if (module_type == "MLDS_topGenes") {
      if (!all(c("MLDS_up", "MLDS_down") %in% names(module_list))) {
        stop("Module list must contain 'MLDS_up' and 'MLDS_down' for MLDS_topGenes.")
      }
      return(list(
        up = module_list[["MLDS_up"]],
        down = module_list[["MLDS_down"]]
      ))
    } else if (module_type == "MLDS_imprints_T21FL") {
      if (!all(c("T21_MLDS_up", "T21_MLDS_down") %in% names(module_list))) {
        stop("Module list must contain 'T21_MLDS_up' and 'T21_MLDS_down' for MLDS_imprints_T21FL")
      }
      return(list(
        up = module_list[["T21_MLDS_up"]],
        down = module_list[["T21_MLDS_down"]]
      ))
    } else if ("all" %in% names(module_list) &&
               all(c("all_up", "all_down") %in% names(module_list[["all"]]))) {
      return(list(
        up = module_list[["all"]][["all_up"]],
        down = module_list[["all"]][["all_down"]]
      ))
    } else {
      stop(paste0(
        "Unrecognized module_type '", module_type, 
        "' and default module structure ('all$all_up' / 'all$all_down') not found."
      ))
    }
  }
  
  ## ----- Helper function: check module validity -----
  check_module_validity <- function(module_list, module_type, min_genes) {
    gene_sets <- get_module_gene_sets(module_list, module_type)
    for (gs_name in names(gene_sets)) {
      genes <- gene_sets[[gs_name]]
      if (length(genes) == 0) {
        stop(sprintf("Gene set '%s' is empty.", gs_name))
      }
      if (length(genes) < min_genes) {
        stop(sprintf("Gene set '%s' has fewer than %d genes (actual: %d).", gs_name, min_genes, length(genes)))
      }
    }
    TRUE
  }
  
  module_type <- match.arg(module_type)

  ## ----- 1. Import bulk counts -----
  bulkSamples <- import_bulk_data(
    gns = gns,
    rm_henning_healthy = rm_henning_healthy,
    oxfordBulk = oxfordBulk,
    filter_lowExpr_gene = filter_lowExpr_gene,
    tpm_fp = tpm_fp,
    oxford_meta_fp = oxford_meta_fp,
    oxford_tpm_fp = oxford_tpm_fp,
    plot_filter = plot_filter
  )
  
  tpm_count <- bulkSamples[['tpm_count']]
  bulk_mdat <- bulkSamples[['mdat']]
  
  
  ## ----- 2. Prepare metadata -----
  bulk_mdat$os = bulk_mdat$OS
  bulk_mdat$os_group = '-'
  bulk_mdat$os_group[!is.na(bulk_mdat$os) & bulk_mdat$os != '-'] = ifelse(as.numeric(bulk_mdat$os[bulk_mdat$os != '-'])<2,'low',
                                                                          ifelse(as.numeric(bulk_mdat$os[bulk_mdat$os != '-'])<5,'mid',
                                                                                 ifelse(as.numeric(bulk_mdat$os[bulk_mdat$os != '-'])<10,'mid2','high')))
  bulk_mdat$category = ifelse(bulk_mdat$Subgroup %in% c('CMP','GMP','MEP','HSC','MPP','Bcell','ELP','LMPP','PreProB','ProB'),'Normal',
                              ifelse(bulk_mdat$Subgroup %in% c('TMD','MLDS'),'TAM / MLDS','Other leukaemia'))
  
  bulk_mdat$category = factor(bulk_mdat$category,c('Normal','TAM / MLDS','Other leukaemia'))
  bulk_mdat$sampleGroup = as.character(bulk_mdat$Subgroup)
  bulk_mdat$sampleGroup[grepl('^t',bulk_mdat$sampleGroup)] = 'MLL_rearrangement'
  bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD'] = 'TAM'
  bulk_mdat$sampleGroup = factor(bulk_mdat$sampleGroup,c('TAM','MLDS','AMKL','MLL','MLL_rearrangement',
                                                         'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))
  
  
  ## ----- 3. Filter and prepare matrix -----
  # Check sample match
  samp <- intersect(colnames(tpm_count)[-c(1:2)], bulk_mdat$Sample)
  if (length(samp) == 0) stop("No TPM samples match metadata Sample column.")
  mtx <- tpm_count[,colnames(tpm_count) %in% bulk_mdat$Sample[!grepl('^PDX',bulk_mdat$Sample)]]
  
  ## ----- 4. Rank genes -----
  ranked <- singscore::rankGenes(mtx)
  
  ## ----- 5. Check module validity -----
  check_module_validity(moduleList, module_type, min_genes)
  
  ## ----- 6. Get gene sets for scoring -----
  gene_sets <- get_module_gene_sets(moduleList, module_type)
  
  ## ----- 7. Calculate module scores -----
  allScore <- data.frame()
  
  for (i in seq_along(moduleList)) {
    module_name <- names(moduleList)[i]
    
    # For the 'all' module, use up/down gene sets per module_type
    if (module_name == 'all') {
      dims <- singscore::simpleScore(bulk_ranked, upSet = gene_sets$up, downSet = gene_sets$down)
      dims <- dims[, c('TotalScore', 'TotalDispersion'), drop = FALSE]
    } else {
      # For other modules, simpleScore with only upSet (no downSet)
      # Validate module genes count here optionally
      if (length(moduleList[[i]]) == 0) {
        warning(sprintf("Module '%s' is empty. Skipping.", module_name))
        next
      }
      if (length(moduleList[[i]]) < min_genes) {
        warning(sprintf("Module '%s' has fewer than %d genes. Skipping.", module_name, min_genes))
        next
      }
      dims <- singscore::simpleScore(bulk_ranked, upSet = moduleList[[i]])
    }
    
    # Merge scores with metadata
    scoredf <- merge(bulk_mdat, dims, by.x = 'Sample', by.y = 0)
    scoredf$moduleType <- paste0(module_type, '_', module_name)
    
    # Append to allScore
    allScore <- rbind(allScore, scoredf)
  }
  
  print(table(allScore$moduleType))
  
  return(allScore)
}







figS6_gata1s_moduleScore_inBulkSamples = function(plotDir,prefix,title=NULL,allScore,module_type = c("GATA1s_topGenes", "MLDS_topGenes","MLDS_imprints_T21FL")){
  library(ggbeeswarm)
  
  if(is.null(title)){
    if(module_type == 'GATA1s_topGenes'){
      title = 'GATA1s module'
    }else if(module_type == 'MLDS_imprints_T21FL'){
      title = 'MLDS imprints in T21 FL'
    }else if(module_type == 'MLDS_topGenes'){
      title = 'MLDS module'
    }
  }
  
  #allScore = bulkRNA_moduleScore_singScore(moduleList=moduleList,module_type=module_type,gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T)
  
  
  # # Group by Events
  # bulk_mdat$sampleGroup = as.character(bulk_mdat$Subgroup)
  # bulk_mdat$sampleGroup[grepl('^t',bulk_mdat$sampleGroup)] = 'MLL_rearrangement'
  # # bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD' & bulk_mdat$Event == 1] = 'TMD:1'
  # # bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD' & bulk_mdat$Event == 0] = 'TMD:0'
  # # bulk_mdat$sampleGroup = factor(bulk_mdat$sampleGroup,c('TMD','TMD:0','TMD:1','MLDS','AMKL','MLL','MLL_rearrangement',
  # #                                                        'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))
  # # bulk_mdat$Subgroup = factor(bulk_mdat$Subgroup,c('TMD','MLDS','AMKL','MLL',
  # #                                                  'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell',
  # #                                                  't(8;21)','t(10;11)','t(6;11)','t(9;11)'))
  #
  
  
  # create a dataframe with the data required: scores and sample group
  allScore$group_facet_hor = allScore$category
  allScore$group_facet_ver = allScore$moduleType
  allScore$sampleGroup[allScore$sampleGroup == 'MLL_rearrangement'] = 'MLL'
  
  if(module_type=='MLDS_topGenes'){
    allScore$sampleGroup = as.character(allScore$Category)
    allScore$sampleGroup[allScore$sampleGroup %in% c('TAM_MLDS','TAM_progressive')] = 'Progressive TAM'
    allScore$sampleGroup[allScore$sampleGroup == 'TAM_conventional'] = 'Conventional TAM'
    
    #allScore$sampleGroup[allScore$sampleGroup == 'TAM' & allScore$Event == 1] = 'TAM:1'
    #allScore$sampleGroup[allScore$sampleGroup == 'TAM' & allScore$Event == 0] = 'TAM:0'
    allScore$sampleGroup = factor(allScore$sampleGroup,c('Conventional TAM','Progressive TAM','TAM_earlyDeath','TAM_unknown','TAM_MLDS','unclassified','MLDS',
                                                         'AMKL','MLL','MLL_rearrangement',
                                                         'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))
    
  }
  
  allScore$group_fill = allScore$sampleGroup
  
  plotFun_GATA1s_moduleScore_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){
    # allScore$moduleType = factor(allScore$moduleType,c('GATA1s_Prog_up','GATA1s_MK_up','GATA1s_Mast_up','GATA1s_Ery_up',
    #                                                    'GATA1s_Myeloid_up','GATA1s_Lymphoid_up',
    #                                                    'GATA1s_Mast/MK/Prog_down','GATA1s_EE/MK/Prog_down','GATA1s_Ery_down',
    #                                                    'GATA1s_Myeloid_down','GATA1s_Bcells_down','GATA1s_NK_T_down',
    #                                                    'GATA1s_all_up','GATA1s_all_down','GATA1s_all'))
    
    p1 = ggplot(allScore, aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.25,lty=2)+
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.25,col='black') +
      geom_quasirandom(aes(col=Event),size=0.2,width = 0.18,alpha=1,col='black')+
      #scale_fill_manual(values = c(rep(col25[4],2),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      scale_fill_manual(values = c(rep(geno_cols[['T21']],2),rep(grey(0.6),7))) +
      #scale_color_manual(values = c('-'=grey(0),'0'='#2D4372','1'=col25[5])) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic()+
      ggtitle(title)+xlab('')+ylab('Module score')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.2),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            strip.text = element_text(size = 6,colour = 'black'),
            axis.ticks = element_line(colour = 'black',linewidth = 0.2),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text = element_text(colour = 'black'))
    
    p2 = ggplot(allScore, aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.3,lty=2)+
      geom_quasirandom(aes(color=group_fill),size=0.5,width = 0.2,alpha=1)+
      stat_summary(
        aes(group = group_fill), fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", color = "black", width = 0.6, lwd = 0.2,
        
        # add this bit here to your stat_summary function
        position=position_dodge(width=0.75)
      )+
      scale_color_manual(values = c(rep(geno_cols[['T21']],2),rep(grey(0.6),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      xlab('')+ylab('Module score')+ggtitle(title)+
      theme_classic(base_size = 12)+
      theme(axis.text.x = element_text(size=10,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
            strip.text = element_text(size = 11,colour = 'black'),
            axis.text.y = element_text(size=10,colour = 'black'),
            panel.border = element_rect(fill=F),axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.4),
            strip.background=element_rect(linewidth=0),
            legend.position = 'none',legend.text = element_text(size=9,colour = 'black'))
    
    print(p1)
    
  }
  
  
  
  
  plotFun_MLDS.degs_moduleScore_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){
    dd = allScore[!allScore$Category %in% c('TAM_earlyDeath','TAM_unknown','unclassified'),]
    
    # allScore$moduleType = factor(allScore$moduleType,c('GATA1s_Prog_up','GATA1s_MK_up','GATA1s_Mast_up','GATA1s_Ery_up',
    #                                                    'GATA1s_Myeloid_up','GATA1s_Lymphoid_up',
    #                                                    'GATA1s_Mast/MK/Prog_down','GATA1s_EE/MK/Prog_down','GATA1s_Ery_down',
    #                                                    'GATA1s_Myeloid_down','GATA1s_Bcells_down','GATA1s_NK_T_down',
    #                                                    'GATA1s_all_up','GATA1s_all_down','GATA1s_all'))
    
    p1 = ggplot(dd[dd$group_fill != 'TAM',], aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.2,lty=2)+
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.2,color='black') +
      geom_quasirandom(aes(col=Event),size=0.25,width = 0.15,alpha=1,col=grey(0))+
      scale_fill_manual(values = c('white',col25[2],rep('white',15))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic()+
      #scale_y_continuous(labels = c(0,'',0.1,'',0.2))+
      ggtitle(title)+xlab('')+ylab('Module score')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.2),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            strip.text = element_text(size = 11,colour = 'black'),
            axis.ticks = element_line(colour = 'black',linewidth = 0.2),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text = element_text(colour = 'black'))
    
    
    
    p2 = ggplot(dd[dd$group_fill != 'TAM',], aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.3,lty=2)+
      geom_quasirandom(aes(color=group_fill),size=0.5,width = 0.2,alpha=1)+
      stat_summary(
        aes(group = group_fill), fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", color = "black", width = 0.6, lwd = 0.2,
        
        # add this bit here to your stat_summary function
        position=position_dodge(width=0.75)
      )+
      scale_color_manual(values = c(rep(colAlpha(geno_cols[['T21']],1),1),'#c7065a',colAlpha(geno_cols[['T21']],0.7),rep(rep(grey(0.7),7)))) +
      #scale_color_manual(values = c(rep(geno_cols[['T21']],2),rep(grey(0.4),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      xlab('')+ylab('Module score')+ggtitle(title)+
      theme_classic(base_size = 12)+
      theme(axis.text.x = element_text(size=10,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
            strip.text = element_text(size = 11,colour = 'black'),
            axis.text.y = element_text(size=10,colour = 'black'),
            panel.border = element_rect(fill=F),axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.4),
            strip.background=element_rect(linewidth=0),
            legend.position = 'none',legend.text = element_text(size=9,colour = 'black'))
    
    
    print(p1)
    
  }
  
  if(module_type == 'GATA1s_topGenes'){
    saveFig(file.path(plotDir,prefix),plotFun_GATA1s_moduleScore_inBulkSamples,rawData=allScore,width = 4.3,height = 7.1,res = 500)
  }else if(module_type == 'MLDS_imprints_T21FL'){
    saveFig(file.path(plotDir,prefix),plotFun_GATA1s_moduleScore_inBulkSamples,rawData=allScore,width = 4.3,height = 5,res = 500)
  }else if(module_type == 'MLDS_topGenes'){
    saveFig(file.path(plotDir,prefix),plotFun_MLDS.degs_moduleScore_inBulkSamples,rawData=dd,width = 4.8,height = 4,res = 500)
    
  }
  
}
