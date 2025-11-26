# Deep-dive into L076 [relapse ML-DS]

# Libraries --------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(zeallot)
source("~/ML-DS/utils/misc.R")
source("~/ML-DS/utils/sc_utils.R")

outDir <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L076"
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = T)
}

plotDir <- "~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/"

# Global params ----------------------------------------------------------------
mlds_srat_fp <- "~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns.RDS"
mlds_mdat_fp <- "~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns_mdat_2508.csv"

generate_l076_srat_obj <- function(mlds_srat_fp, mlds_mdat_fp, outDir) {
  if (!file.exists(file.path(outDir, "L076_sratObj.RDS"))) {
    mlds <- readRDS(mlds_srat_fp)
    mlds_mdat <- read.csv(mlds_mdat_fp, row.names = 1)
    mlds_mdat$annot <- mlds_mdat$annot_aug24_new
    mlds_mdat$finalAnn_broad <- mlds_mdat$annot_aug24_new
    mlds_mdat$broadLineage[mlds_mdat$finalAnn_broad == "Tumour" & mlds_mdat$broadLineage != "Tumour"] <- "Tumour"
    mlds_mdat$broadLineage[mlds_mdat$finalAnn_broad %in% c("Tum_MK?", "Tumour_WT")] <- "Tumour_unsure"
    checkmate::assert_true(!any(is.na(match(rownames(mlds@meta.data), rownames(mlds_mdat)))))
    mlds@meta.data <- mlds_mdat[match(rownames(mlds@meta.data), rownames(mlds_mdat)), ]

    l076 <- subset(mlds, donorID == "L076")
    l076 <- standard_clustering(l076)
    saveRDS(l076, file.path(outDir, "L076_sratObj.RDS"))

    # L076 - bone marrow samples only
    l076_bm <- subset(l076, subset = tissue == "BM")
    l076_bm <- standard_clustering(l076_bm, clusteringRes = 0.2)
    saveRDS(l076_bm, file.path(outDir, "L076_BM.samples.only_sratObj.RDS"))
  } else {
    l076 <- readRDS(file.path(outDir, "L076_sratObj.RDS"))
    l076_bm <- readRDS(file.path(outDir, "L076_BM.samples.only_sratObj.RDS"))
  }

  list(l076 = l076, l076_bm = l076_bm)
}

# Sub-clustering L076 ----------------------------------------------------------
c(l076, l076_bm) %<-% generate_l076_srat_obj(
  mlds_srat_fp = mlds_srat_fp,
  mlds_mdat_fp = mlds_mdat_fp,
  outDir = outDir
)
fig4B <- function() {
  l076_bm$GATA1s_status <- l076_bm$GATA1_status
  l076_bm$GATA1s_status[l076_bm$GATA1s_status == "no_GATA1_expr"] <- "No GATA1 expression"
  l076_bm$GATA1s_status[l076_bm$GATA1s_status == "GATA1s_WT"] <- "GATA1 wild type"
  l076_bm$GATA1s_status[l076_bm$GATA1s_status == "GATA1s_mutant"] <- "GATA1s mutation"
  l076_bm$GATA1s_status[l076_bm$GATA1s_status %in% c("unsure", "noCov", "uninformative")] <- "Uninformative"

  dd <- cbind(
    l076_bm@meta.data[, c(
      "cellID", "orig.ident", "donorID", "tissue",
      "disease", "timePoint", "clinicalOutcome", "annot",
      "finalAnn_broad", "GATA1s_status"
    )],
    l076_bm@reductions$umap@cell.embeddings
  )

  plotFun_GATA1status <- function(noFrame = FALSE, noPlot = FALSE) {
    par(mar = c(0.1, 0.1, 1, 0.1))
    ccs <- c(
      "No GATA1 expression" = grey(0.9),
      "Uninformative" = grey(0.55),
      "GATA1s mutation" = "#A92821",
      "GATA1 wild type" = "#005579"
    )
    p <- ggplot(dd, aes(UMAP_1, UMAP_2)) +
      geom_point(size = 0.001, aes(col = GATA1s_status), alpha = 0.4) +
      scale_color_manual(values = ccs) +
      theme_classic(base_size = 8.3) +
      theme(
        panel.border = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 3)
      ) +
      xlab("") +
      ylab("")

    print(p)
  }
  saveFig(file.path(plotDir, "Figure4", "Fig4C_L076_GATA1s_UMAP"),
    plotFun_GATA1status,
    rawData = dd, width = 3.65, height = 2.65, res = 500
  )

  #  L076 - BM Tumour cells only
  l076.tum <- subset(l076_bm, subset = cellID %in%
    l076_bm$cellID[l076_bm$annot == "Tumour"])
  l076.tum <- standard_clustering(l076.tum, clusteringRes = 0.2)

  dd <- cbind(
    l076.tum@meta.data[, c(
      "cellID", "orig.ident", "donorID", "tissue",
      "disease", "timePoint", "clinicalOutcome", "annot",
      "finalAnn_broad", "GATA1s_status"
    )],
    l076.tum@reductions$umap@cell.embeddings
  )

  plotFun_timePoint <- function(noFrame = FALSE, noPlot = FALSE) {
    par(mar = c(0.1, 0.1, 1, 0.1))
    ccs <- c(
      "D.Relapse" = col25[5],
      "D.Relapse2" = col25[3],
      "Diagnostic" = pal34H[34],
      "Normal" = grey(0.9)
    )
    dd$group <- ifelse(dd$annot == "Tumour", dd$timePoint, "Normal")
    p <- ggplot(dd, aes(UMAP_1, UMAP_2)) +
      geom_point(size = 0.0001, aes(col = group)) +
      scale_color_manual(values = ccs) +
      theme_classic(base_size = 8.3) +
      theme(
        panel.border = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 3)
      ) +
      xlab("") +
      ylab("")

    print(p)
  }

  saveFig(file.path(plotDir, "Figure4", "Fig4C_L076.BM.tum_timePoint_UMAP"),
    plotFun_timePoint,
    rawData = dd, width = 4, height = 3.65, res = 500
  )
}


# Run alleleIntegrator on L076 -------------------------------------------------
