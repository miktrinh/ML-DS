# Deep-dive into L038 [refractory ML-DS]

# Libraries --------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(zeallot)
source("~/ML-DS/utils/misc.R")
source("~/ML-DS/utils/sc_utils.R")

outDir <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L038/"
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = T)
}

plotDir <- "~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/"

# Global params ----------------------------------------------------------------
mlds_srat_fp <- "~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns.RDS"
mlds_mdat_fp <- "~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns_mdat_2508.csv"

generate_l038_srat_obj <- function(mlds_srat_fp, mlds_mdat_fp, outDir) {
  if (!file.exists(file.path(outDir, "L038_sratObj.RDS"))) {
    mlds <- readRDS(mlds_srat_fp)
    mlds_mdat <- read.csv(mlds_mdat_fp, row.names = 1)
    mlds_mdat$annot <- mlds_mdat$annot_aug24_new
    mlds_mdat$finalAnn_broad <- mlds_mdat$annot_aug24_new
    mlds_mdat$broadLineage[mlds_mdat$finalAnn_broad == "Tumour" & mlds_mdat$broadLineage != "Tumour"] <- "Tumour"
    mlds_mdat$broadLineage[mlds_mdat$finalAnn_broad %in% c("Tum_MK?", "Tumour_WT")] <- "Tumour_unsure"
    checkmate::assert_true(!any(is.na(match(rownames(mlds@meta.data), rownames(mlds_mdat)))))
    mlds@meta.data <- mlds_mdat[match(rownames(mlds@meta.data), rownames(mlds_mdat)), !colnames(mlds_mdat) %in% c('UMAP_1.1','UMAP_2.1')]
    
    l038 <- subset(mlds, donorID == "L038")
    l038 <- standard_clustering(l038)
    saveRDS(l038, file.path(outDir, "L038_sratObj.RDS"))
    
    # L038 - Tumour + Ery cells only
    l038_tumour = subset(l038,subset = cellID %in% l038$cellID[grepl('Tum|MEP|EE|ME|LE',l038$annot) &
                                                                 !grepl('unsure',l038$annot)])
    l038_tumour = standard_clustering(l038_tumour,clusteringRes = 1.3)
    # DimPlot(l038_tumour,label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()
    # DimPlot(l038_tumour,group.by = 'group',label = T,repel = T,label.box = T,cols = c(col25,pal34H)) + NoLegend()

    df = l038_tumour@meta.data[,colnames(l038_tumour@meta.data) != 'group']
    df$umap_1 = l038_tumour@reductions$umap@cell.embeddings[,1]
    df$umap_2 = l038_tumour@reductions$umap@cell.embeddings[,2]

    write.csv(df,file.path(outDir,'L038_TumourEry_subClustering_mdat.csv'))

    
  } else {
    l038 <- readRDS(file.path(outDir, "L038_sratObj.RDS"))
  }
  
  l038
}

# Sub-clustering L038 ----------------------------------------------------------

l038 <- generate_l038_srat_obj(
  mlds_srat_fp = mlds_srat_fp,
  mlds_mdat_fp = mlds_mdat_fp,
  outDir = outDir
)



# Run alleleIntegrator on L038 -------------------------------------------------
# Import results
l038_aiRes = read.csv(file.path('~/ML-DS/Results/06_MLDS_refractory_relapse/L038/L038_alleleIntegrator/','L038_AI_res.csv'),row.names = 1)
l038@meta.data = l038_aiRes[match(rownames(l038@meta.data),rownames(l038_aiRes)),] %>% 
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

checkmate::assert_true(all(!is.na(l038@meta.data$AIres)))
checkmate::assert_true((!'others' %in% l038$AIres) && (!'others' %in% l038$GATA1s_status))



# sub-cluster by time point
l038_by_timepoint <- function(l038_srat, outDir){
  l038_d_mdat_fp = file.path(outDir,'L038_Diagnostic_AIresult_mdat.csv')
  if(!file.exists(l038_d_mdat_fp)){
    l038_d = subset(l038_srat,subset = timePoint == 'Diagnostic')
    l038_d = standard_clustering(l038_d)
    l038_d = FindClusters(l038_d,resolution = 1.5)
    l038_d$group = ifelse(l038_d$seurat_clusters %in% c(32,13,17,20,0,21),'clone2_D',
                          ifelse(l038_d$annot == 'Tumour','clone1_D','normal_D'))
    DimPlot(l038_d,group.by = 'group',label = T,repel = T,label.box = T) + NoLegend()
    l038_d_mdat = cbind(l038_d@meta.data[,!colnames(l038_d@meta.data) %in% c('UMAP_1','UMAP_2')],
                        l038_d@reductions$umap@cell.embeddings)
    write.csv(l038_d_mdat,l038_d_mdat_fp)
  }else{
    l038_d_mdat = read.csv(l038_d_mdat_fp,row.names = 1)
  }
  
  l038_tp1_mdat_fp = file.path(outDir,'L038_TP1_AIresult_mdat.csv')
  if(!file.exists(l038_tp1_mdat_fp)){
    l038_tp1 = subset(l038_srat,subset = timePoint == 'TP1')
    l038_tp1 = standard_clustering(l038_tp1)
    l038_tp1$group = ifelse(l038_tp1$annot == 'Tumour','clone2_TP1','normal_TP1')
    
    DimPlot(l038_tp1,group.by = 'group',label = T,repel = T,label.box = T) + NoLegend()
    DimPlot(l038_tp1,cells.highlight = l038_tp1$cellID[l038_tp1$AI_output == 'normFrac' & l038_tp1$annot == 'Tumour'])
    DimPlot(l038_tp1,cells.highlight = l038_tp1$cellID[!l038_tp1$seurat_clusters %in% c(1,11,18,13,1,23) & l038_tp1$annot == 'Tumour'])
    l038_tp1_mdat = cbind(l038_tp1@meta.data[,!colnames(l038_tp1@meta.data) %in% c('UMAP_1','UMAP_2')],
                          l038_tp1@reductions$umap@cell.embeddings)
    write.csv(l038_tp1_mdat,l038_tp1_mdat_fp)
  }else{
    l038_tp1_mdat = read.csv(l038_tp1_mdat_fp,row.names = 1)  
  }
  
  list(l038_d_mdat,l038_tp1_mdat)
}

c(l038_d_mdat,l038_tp1_mdat) %<-% l038_by_timepoint(l038,outDir=outDir)


# Plot Figure 4E ---------------------------------------------------------------

fig4E <- function() {
  shared_columns = intersect(colnames(l038_d_mdat),colnames(l038_tp1_mdat))
  dd=rbind(l038_d_mdat[,shared_columns],l038_tp1_mdat[,shared_columns])
  dd <- dd %>% 
    dplyr::mutate(AI_output_pp = ifelse(broadLineage != 'Tumour',NA,AI_output_pp)) %>% 
    dplyr::select(c(
    "cellID", "orig.ident", "donorID", "tissue",
    "disease", "timePoint", "clinicalOutcome", "annot",
    "finalAnn_broad", "GATA1s_status", "AIres", "AI_output_pp", 'UMAP_1',"UMAP_2"
  ))
  
  
  # Plot GATA1 status ----------------------------------------------------------
  
  plotFun_GATA1status <- function(noFrame = FALSE, noPlot = FALSE) {
    par(mar = c(0.1, 0.1, 1, 0.1))
    ccs <- c(
      "No GATA1 expression" = grey(0.9),
      "Uninformative" = grey(0.55),
      "GATA1s mutation" = "#A92821",
      "GATA1 wild type" = "#005579"
    )
    p <- ggplot(dd,aes(UMAP_1,UMAP_2,col=GATA1s_status))+
      scale_color_manual(values = ccs)+
      facet_wrap(vars(timePoint))+
      theme_classic()+
      theme(panel.border = element_rect(fill=NA,colour = 'black'),
            strip.background = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            legend.text = element_text(size=9))
    
    if(!noPlot){
      p <- p + geom_point(size = 0.001, aes(col = GATA1s_status), alpha = 0.4)
    }
    print(p)
  }
  
  saveFig(file.path(plotDir, "Figure4", "Fig4E_L038_GATA1s_UMAP"),
          plotFun_GATA1status,
          rawData = dd, width = 5.9, height = 2.3, res = 500
  )
  
  # AI results -----------------------------------------------------------------
  
  plotFun_AIresult = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    group_cols = c('Uninformative'=grey(0.8),
                   'without CNA' = grey(0.3),
                   'with CNA' = '#A92821')
    
    
    p <- ggplot(dd,aes(UMAP_1,UMAP_2,col=AIres))+
      scale_color_manual(values = group_cols)+
      facet_wrap(vars(timePoint))+
      theme_classic()+
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
  
  saveFig(file.path(plotDir,"Figure4",'Fig4E_L038_AIresult_UMAP'),
          plotFun_AIresult,
          rawData=dd,width = 5.4,height = 2.3,res = 500)
  
  
  # Plot AI BAF ----------------------------------------------------------------
  ## plot BAF in each celltype
  dd=rbind(l038_d_mdat[,shared_columns],l038_tp1_mdat[,shared_columns])
  l038$group = ifelse(l038$cellID %in% dd$cellID,dd$group[match(l038$cellID,dd$cellID)],'?')
  df = plot_BAF_byCellClusters(mDat=l038@meta.data, cellID_column='cellID',
                               group = 'group', 
                               phCnts_fp = file.path(outDir,'L038_alleleIntegrator','PD60301a.PD61846a_phCnts.RDS'),
                               normalGroups = c('normal_D','normal_TP1'),
                               outDir=plotDir,patientID='PD60301a',PDID='PD60301a',
                               tgtChrs=tgtChrs,minRead = 3)
  
  ## Config
  chr_cfg <- tibble::tibble(
    chrom = c("chr5", "chr17"),
    x_max = c(max(df$pos[df$chr == "chr5"], na.rm = TRUE) / 1e6, 23),
    highlight = list(NULL, data.frame(xmin = 7661779/1e6, xmax = 7687550/1e6))
  )
  
  base_plot <- function(data, size, alpha, highlight = NULL) {
    p <- ggplot(data, aes(pos / 1e6, altFreq)) +
      geom_rect(aes(xmin = 0, xmax = x_max, ymin = -0.02, ymax = 1.02),
                fill = grey(0.9), alpha = 0.2) +
      geom_point(data = data[data$phasingAssign=='uninformative',],aes(col = phasingAssign), size = size, alpha = alpha) +
      geom_point(data = data[data$phasingAssign!='uninformative',],aes(col = phasingAssign), size = size, alpha = alpha) +
      geom_hline(yintercept = 0.5, lty = 1, lwd = 0.3, col = "black") +
      scale_color_manual(values = c('minor_allele' = "black", 'major_allele' = "#FF4D00", "uninformative"=grey(0.7))) +
      scale_y_continuous(breaks = c(0.0, 0.5, 1.0), labels = c(0, 0.5, 1)) +
      facet_grid(chr ~ clusterID, scales = "free_x") +
      theme_classic(base_size = 13) +
      theme(
        strip.background = element_blank(),
        axis.text = element_text(size = 8, colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        legend.position = 'none'
      ) +
      xlab("Genomic position") + ylab("Aggregated Alt allele frequency")
    
    if (!is.null(highlight)) {
      p <- p + geom_rect(data = highlight,
                         aes(xmin = xmin, xmax = xmax, ymin = -0.02, ymax = 1.02),
                         col = "darkblue", inherit.aes = FALSE)
    }
    p
  }
  
  for (i in seq_len(nrow(chr_cfg))) {
    chrom <- chr_cfg$chrom[i]
    x_max <- chr_cfg$x_max[i]
    highlight <- chr_cfg$highlight[[i]]

    dd = df %>% 
      filter(totCount > 10 & 
               chr %in% chrom & 
               grepl('clone1_D|clone2_D|clone2_TP1',clusterID)) %>% 
      mutate(clusterID = case_when(grepl('clone1_D \\(',clusterID) ~ 'Clone 1 (D)',
                                   grepl('clone2_D \\(',clusterID) ~ 'Clone 2 (D)',
                                   grepl('clone2_TP1 \\(',clusterID) ~ 'Clone 2 (TP1)',
                                   .default = clusterID
                                   ),
             clusterID = factor(clusterID,c('Clone 1 (D)','Clone 2 (D)','Clone 2 (TP1)'))
             )
    
    
    plotFun_rawBAF_L038Tum = function(noFrame=FALSE,noPlot=FALSE){
      par(mar=c(0.1,0.1,1,0.1))
      if (!noPlot && !noFrame) {
        p <- base_plot(dd, size = 0.05, alpha = 0.5, highlight=highlight)+
          theme(panel.border = element_rect(fill = NA, colour = "black"))
        print(p)
      } else if (!noPlot && noFrame) {
        print(base_plot(dd, size = 0.05, alpha = 0.5, highlight=highlight))
      } else if (noPlot && !noFrame) {
        df2 <- df %>%
          group_by(clusterID) %>%
          slice_head(n = 2) %>%
          ungroup()
        if (chrom == "chr17") {
          df2 <- bind_rows(df2,
                           df %>% filter(chr == chrom, pos == max(pos[chr == chrom])) %>% slice(1))
        }
        print(base_plot(df2, size = 0.05, alpha = 0.5, highlight=highlight))
      }
    }
    saveFig(file.path(plotDir,"Figure4",paste0('Fig4E_L038Tum_rawBAF_',chrom)),
            plotFun_rawBAF_L038Tum,
            rawData=dd,width = 7.3,height = 2,res = 500)
  }
}
    

# L038 Refractory markers ------------------------------------------------------
l038$group = ifelse(l038$annot_aug24 == 'Tumour' & l038$cellID %in% l038_d_mdat$cellID[l038_d_mdat$group == 'clone2_D'], 'clone2_D',
                    ifelse(l038$annot_aug24 == 'Tumour' & l038$timePoint == 'Diagnostic','clone1',
                           ifelse(l038$annot_aug24 == 'Tumour' & l038$timePoint == 'TP1' & l038$cellID %in% l038_tp1_mdat$cellID[l038_tp1_mdat$group == 'clone2_TP1'] & l038$seurat_clusters != 7,'clone2',
                                  ifelse(l038$annot_aug24 == 'Tumour' & l038$timePoint == 'TP1','clone2_ery','Normal'))))

DimPlot(l038,group.by = 'group',label = T)
DimPlot(l038,group.by = 'seurat_clusters',label = T)
Idents(l038) = l038$group

# Refractory clone vs diagnostic clone
l038_tp1.vs.clone1 = FindMarkers(l038,ident.1 = 'clone2',ident.2 = 'clone1')
l038_tp1.vs.clone1$geneSym = rownames(l038_tp1.vs.clone1)
l038_tp1.vs.clone1$pct_diff = l038_tp1.vs.clone1$pct.1 - l038_tp1.vs.clone1$pct.2
l038_tp1.vs.clone1$comp = 'TP1_vs_Dclone1'
rownames(l038_tp1.vs.clone1) = geneMap$ensID[match(l038_tp1.vs.clone1$geneSym,geneMap$geneSym)]
l038_tp1.vs.clone1 = annotateGenes(l038_tp1.vs.clone1,geneMap = geneMap)
l038_tp1.vs.clone1$direction = ifelse(l038_tp1.vs.clone1$avg_log2FC > 0,'clone2_up','clone2_down')

write.csv(l038_tp1.vs.clone1,file.path(outDir,'L038_tp1.vs.clone1.D_markers.csv'))

FeaturePlot(l038,'CD82')
    
    
    
      
      
    






















