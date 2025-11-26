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


# Run alleleIntegrator on L076 -------------------------------------------------
# Import results
l076_aiRes = read.csv(file.path('~/ML-DS/Results/06_MLDS_refractory_relapse/L076/L076_alleleIntegrator/','L076_AI_res.csv'),row.names = 1)
l076@meta.data = l076_aiRes[match(rownames(l076@meta.data),rownames(l076_aiRes)),] %>% 
  dplyr::mutate(GATA1s_status = dplyr::case_when(GATA1_status == 'no_GATA1_expr' ~ 'No GATA1 expression',
                                                 GATA1_status == 'GATA1s_WT' ~ 'GATA1 wild type',
                                                 GATA1_status == 'GATA1s_mutant' ~ 'GATA1s mutation',
                                                 GATA1_status %in% c('unsure','noCov','uninformative','GATA1s_unsure') ~ 'Uninformative',
                                                 .default = 'others'),
                AIres = dplyr::case_when(AI_output %in% c('?','Uncalled') | broadLineage != 'Tumour' ~ 'Uninformative',
                                         AI_output %in% c('normFrac') ~ 'without CNA',
                                         AI_output %in% c('abbFrac') ~ 'with CNA',
                                         .default = 'others'),
                AIres = factor(AIres,c('Uninformative','with CNA','without CNA')),
                AI_output_pp = ifelse(broadLineage != 'Tumour',NA,AI_output_pp))

checkmate::assert_true(all(!is.na(l076@meta.data$AIres)))
checkmate::assert_true((!'others' %in% l076$AIres) && (!'others' %in% l076$GATA1s_status))

l076_bm$GATA1s_status = l076$GATA1s_status[match(l076_bm$cellID,l076$cellID)]
l076_bm$AIres = l076$AIres[match(l076_bm$cellID,l076$cellID)]
l076_bm$AI_output_pp = l076$AI_output_pp[match(l076_bm$cellID,l076$cellID)]


fig4B <- function() {
  # l076_bm$GATA1s_status <- l076_bm$GATA1_status
  # l076_bm$GATA1s_status[l076_bm$GATA1s_status == "no_GATA1_expr"] <- "No GATA1 expression"
  # l076_bm$GATA1s_status[l076_bm$GATA1s_status == "GATA1s_WT"] <- "GATA1 wild type"
  # l076_bm$GATA1s_status[l076_bm$GATA1s_status == "GATA1s_mutant"] <- "GATA1s mutation"
  # l076_bm$GATA1s_status[l076_bm$GATA1s_status %in% c("unsure", "noCov", "uninformative","GATA1s_unsure")] <- "Uninformative"
  
  dd <- cbind(
    l076_bm@meta.data[, c(
      "cellID", "orig.ident", "donorID", "tissue",
      "disease", "timePoint", "clinicalOutcome", "annot",
      "finalAnn_broad", "GATA1s_status", "AIres", "AI_output_pp"
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
  saveFig(file.path(plotDir, "Figure4", "Fig4B_L076_GATA1s_UMAP"),
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
      "finalAnn_broad", "GATA1s_status", "AIres", "AI_output_pp"
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
  
  saveFig(file.path(plotDir, "Figure4", "Fig4B_L076.BM.tum_timePoint_UMAP"),
          plotFun_timePoint,
          rawData = dd, width = 4, height = 3.65, res = 500
  )
  
  # AI results -----------------------------------------------------------------
  
  plotFun_AIresult = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    p <- ggplot(dd,aes(UMAP_1,UMAP_2,col=AI_output_pp))+
      theme_classic()+
      scale_color_gradient2(low = 'black',high='red')+
      theme(panel.border = element_rect(fill=NA,colour = 'black'),
            strip.background = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.text = element_text(size=7))
    
    if(!noPlot){
      p <- p + geom_point(size = 0.001,alpha = 0.4)
    }
    print(p)
  }
  
  saveFig(file.path(plotDir,"Figure4",'Fig4B_L076.BM.tum_AIres_UMAP'),
          plotFun_AIresult,
          rawData=dd,width = 5.4,height = 2.3,res = 500)
  
  
  # Plot AI BAF ----------------------------------------------------------------
  ## plot BAF in each celltype
  l076$group = ifelse(l076$broadLineage == 'Tumour',paste0(l076$timePoint,':',l076$tissue),
                      ifelse(!grepl('Tumour',l076$broadLineage),'Normal','others'))
  df = plot_BAF_byCellClusters(mDat=l076@meta.data, cellID_column='cellID',
                               group = 'group', 
                               phCnts_fp = file.path(outDir,'L076_alleleIntegrator','PD62331a_phCnts.RDS'),
                               normalGroups = c('Normal'),
                               outDir=plotDir,patientID='PD62331a',PDID='PD62331a',
                               tgtChrs=tgtChrs,minRead = 3)
  data = df %>% dplyr::filter(totCount > 10 & 
                                chr %in% c("chr5", "chr8", "chr13", "chr21") & 
                                clusterID != 'others (n=17)') %>% 
    dplyr::mutate(clusterID = factor(clusterID,c("Diagnostic:Blood (n=1040)",
                                                 "Diagnostic:BM (n=2008)",
                                                 "D.Relapse:BM (n=695)",
                                                 "D.Relapse2:BM (n=12256)",
                                                 "Normal (n=9469)")))
  
  p_base <- ggplot(data, aes(pos / 1e6, altFreq)) +
    geom_point(data = data[data$phasingAssign=='uninformative',],aes(col = phasingAssign), size = 0.05, alpha = 0.5) +
    geom_point(data = data[data$phasingAssign!='uninformative',],aes(col = phasingAssign), size = 0.05, alpha = 0.5) +
    geom_hline(yintercept = 0.5, lty = 1, lwd = 0.3, col = "black") +
    scale_color_manual(values = c('minor_allele' = "black", 'major_allele' = "#FF4D00", "uninformative"=grey(0.7))) +
    scale_y_continuous(breaks = c(0.0, 0.5, 1.0), labels = c(0, 0.5, 1)) +
    facet_grid(clusterID ~ chr, scales = "free_x") +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_blank(),
      axis.text = element_text(size = 8, colour = "black"),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      legend.position = 'none'
    ) +
    xlab("Genomic position") + ylab("Aggregated Alt allele frequency")
    
    plotFun_rawBAF_L076 = function(noFrame=FALSE,noPlot=FALSE){
      par(mar=c(0.1,0.1,1,0.1))
      if (!noPlot && !noFrame) {
        p <- p_base +
          theme(panel.border = element_rect(fill = NA, colour = "black"))
        print(p)
      } else if (!noPlot && noFrame) {
        print(p_base)
      } 
    }
    
    saveFig(file.path(plotDir,"FigureS9",'FigS9B_L076_rawBAF'),
            plotFun_rawBAF_L076,
            rawData=data,width = 5.5,height = 4.5,res = 500)
}


# L076 relapse markers ---------------------------------------------------------

l076$group = ifelse(l076$annot_aug24 == 'Tumour' & l076$tissue=='Blood','D_Blood',
                    ifelse(l076$annot_aug24 == 'Tumour' & l076$tissue!='Blood',l076$timePoint,'others'))
l076$group[l076$group == 'others'] = 'Normal'
DimPlot(l076,group.by = 'group',cols = col25,label = T,repel = T,label.box = T) + NoAxes() + NoLegend()

Idents(l076) = l076$group

## How is relapse 1 different from diagnostic
markers_rd1_BM = FindMarkers(l076,ident.1 = 'D.Relapse',ident.2 = 'Diagnostic')
markers_rd1_BM$geneSym = rownames(markers_rd1_BM)
markers_rd1_BM$pct_diff = markers_rd1_BM$pct.1 - markers_rd1_BM$pct.2
markers_rd1_BM$comp = 'RD_vs_D'
markers_rd1_BM.sub = markers_rd1_BM[abs(markers_rd1_BM$pct_diff) > 0.5,]
markers_rd1_BM.sub = markers_rd1_BM.sub[order(markers_rd1_BM.sub$pct_diff,decreasing = T),]
rownames(markers_rd1_BM.sub) = geneMap$ensID[match(markers_rd1_BM.sub$geneSym,geneMap$geneSym)]
markers_rd1_BM.sub = annotateGenes(markers_rd1_BM.sub,geneMap = geneMap)
DotPlot(l076,group.by = 'group',features = markers_rd1_BM.sub$geneSym[markers_rd1_BM.sub$avg_log2FC < 0]) + RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))

FeaturePlot(l076,'EPS8')

## How is relapse 2 different from diagnostic
markers_rd2_BM = FindMarkers(l076,ident.1 = 'D.Relapse2',ident.2 = 'Diagnostic')
markers_rd2_BM$geneSym = rownames(markers_rd2_BM)
markers_rd2_BM$pct_diff = markers_rd2_BM$pct.1 - markers_rd2_BM$pct.2
markers_rd2_BM$comp = 'R2D_vs_D'
markers_rd2_BM.sub = markers_rd2_BM[abs(markers_rd2_BM$pct_diff) > 0.1,]
markers_rd2_BM.sub = markers_rd2_BM.sub[order(markers_rd2_BM.sub$pct_diff,decreasing = T),]
rownames(markers_rd2_BM.sub) = geneMap$ensID[match(markers_rd2_BM.sub$geneSym,geneMap$geneSym)]
markers_rd2_BM.sub = annotateGenes(markers_rd2_BM.sub,geneMap = geneMap)

markers = rbind(markers_rd1_BM,markers_rd2_BM)
write.csv(markers,file.path(outDir,'L076_relapse_FindMarkers.csv'))
