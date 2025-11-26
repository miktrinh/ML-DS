##----   Deriving TAM-GATA1s module, using good TAM only    -----##

outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/oct24'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)


##---------------##
#### Libraries ####
##---------------##
library(Seurat)
library(GenomicFeatures)
library(tidyverse)
library(UpSetR)
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")



##----------------------------##
##   Set Global parameters  ####
##----------------------------##
max_cellfrac = 10/100
min_l2FC = 0.5

#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf.gz'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)
## Import Gene Map
geneMap = read.delim('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_45842_SB_Leuk13104278_GRCh38-2020-A/filtered_feature_bc_matrix/features.tsv.gz',header = F)
colnames(geneMap) = c('ensID','geneSym','GEX')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))



##------------------------------------------------##
#### 1. Import relevant scRNA datasets: MLDS    ####
##-----------------------------------------------##
## Import MLDS object
mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS')
mlds$annot = as.character(mlds$annot_aug24)




##-------------------------------------------------------##
##   2. DEGs between good TAM vs gosh_T21 MEMP_MEP     ####
##-------------------------------------------------------##
## I will use the results from FM goodTAM.vs.allMEMP, because:
#  1. There is only 1 good TAM sample --> no replicate --> should not use DESeq2 or EdgeR
#  2. for consistancy with MLDS.vs.MEMP comparison

##------- Using Find Markers -----------
mlds$group_tmp = ifelse(mlds$disease == 'TAM' & mlds$annot == 'Tumour' & mlds$donorID == 'L075','TAM',
                        ifelse(mlds$annot == 'MEP' & mlds$donorID %in% c('L019','L038','L039','L040') ,'MEP','others'))
Idents(mlds) = mlds$group_tmp
tam.vs.goshMEP.markers = FindMarkers(mlds,ident.1 = 'TAM',ident.2 = 'MEP')
tam.vs.goshMEP.markers$geneSym = rownames(tam.vs.goshMEP.markers)
tam.vs.goshMEP.markers$ensID = geneMap$ensID[match(tam.vs.goshMEP.markers$geneSym,geneMap$geneSym)]
rownames(tam.vs.goshMEP.markers) = tam.vs.goshMEP.markers$ensID
tam.vs.goshMEP.markers = annotateGenes(tam.vs.goshMEP.markers,geneMap = geneMap)

tam.vs.goshMEP.markers$pct.diff = tam.vs.goshMEP.markers$pct.1 - tam.vs.goshMEP.markers$pct.2
tam.vs.goshMEP.markers$direction = ifelse(tam.vs.goshMEP.markers$avg_log2FC > 0 ,'TAM_up','TAM_down')
write.csv(tam.vs.goshMEP.markers,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/FM_goodTAM.vs.goshMEP_allGene.csv')

tam.vs.goshMEP.markers = read.csv('FM_goodTAM.vs.goshMEP_allGene.csv')

# Subset for Significant Markers only
tam.vs.goshMEP.markers.deg = tam.vs.goshMEP.markers[tam.vs.goshMEP.markers$p_val_adj < 0.01 & abs(tam.vs.goshMEP.markers$pct.diff) >= 0.1,]
tam.vs.goshMEP.markers.deg = tam.vs.goshMEP.markers.deg[order(abs(tam.vs.goshMEP.markers.deg$avg_log2FC),decreasing = T),]
table(tam.vs.goshMEP.markers.deg$direction)
View(tam.vs.goshMEP.markers.deg[tam.vs.goshMEP.markers.deg$avg_log2FC < 0,])


## Perform enrichR on these top DEGs
upDEG_enriched <- enrichr(unique(tam.vs.goshMEP.markers.deg$geneSym[tam.vs.goshMEP.markers.deg$direction == 'TAM_up']), dbs)
downDEG_enriched <- enrichr(unique(tam.vs.goshMEP.markers.deg$geneSym[tam.vs.goshMEP.markers.deg$direction == 'TAM_down']), dbs)

View(upDEG_enriched[[6]])
downDEG_enriched[[6]][1,]
plotEnrich(upDEG_enriched[[5]], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")
plotEnrich(downDEG_enriched[[5]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")



mlds$group_tmp = ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID !='L041',paste0(mlds$disease,':',mlds$donorID),mlds$annot)
Idents(mlds) = mlds$group_tmp

DotPlot(mlds,idents = unique(mlds$group_tmp[grepl('MLDS:|TAM:|MEP',mlds$group_tmp)]),
        features = tam.vs.goshMEP.markers.deg$geneSym[tam.vs.goshMEP.markers.deg$avg_log2FC < 0][1:30]) + RotatedAxis()



##------- Using pseudobulk -----------

## Version 1: only canonical / good TAM
cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in% 
                                                c(mlds$cellID[mlds$donorID %in% c('CC1','CC2','L075','L182','L114','CC6','CC8','CC7') & mlds$annot == 'Tumour'],
                                                  mlds$cellID[mlds$donorID %in% c('L019','L039','L040','L178') & mlds$timePoint %in% c('TP1','TP2','TP4') & mlds$annot == 'MEP'])]]

# ## Version 2: using All TAM cases
# cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in% 
#                                                 c(mlds$cellID[mlds$disease == 'TAM' & mlds$annot == 'Tumour'],
#                                                   mlds$cellID[mlds$donorID %in% c('L019','L038','L039','L040') & mlds$annot == 'MEP'])]]


rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex')]

mDat$donorID = as.character(mDat$donorID)
mDat = mDat[match(colnames(cnt_mtx),mDat$cellID),]

tam.vs.gMEP_res = compareCell_simplified(toc = cnt_mtx,
                                         mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                                         coords = gns[rownames(cnt_mtx)],
                                         cellTypeName='TAM',
                                         formula='~ %s',tgtChrs = paste0('chr',c(1:22)),
                                         donorID='donorID',groupID='annot')

saveRDS(tam.vs.gMEP_res,'DESeq2_goodTAM.vs.goshMEP_res.RDS')
tam.vs.gMEP_res = readRDS('DESeq2_goodTAM.vs.goshMEP_res.RDS')

#saveRDS(tam.vs.gMEP_res,'DESeq2_goodTAM.vs.goshMEP_res.RDS')

tam.vs.gMEP_dds = tam.vs.gMEP_res[['dds']]
tam.vs.gMEP_deg = tam.vs.gMEP_res[['mainDE']]
tam.vs.gMEP_deg$direction = ifelse(tam.vs.gMEP_deg$log2FoldChange > 0, 'TAM_up','TAM_down')
tam.vs.gMEP_deg$cellFrac.diff = tam.vs.gMEP_deg$cellFrac_g1 - tam.vs.gMEP_deg$cellFrac_g2
tam.vs.gMEP_deg = tam.vs.gMEP_deg[tam.vs.gMEP_deg$padj < 0.05 & abs(tam.vs.gMEP_deg$log2FoldChange) >= min_l2FC & 
                                    abs(tam.vs.gMEP_deg$cellFrac.diff) >= max_cellfrac,]
table(tam.vs.gMEP_deg$direction)
dim(tam.vs.gMEP_deg)
tam.vs.gMEP_deg = tam.vs.gMEP_deg[order(abs(tam.vs.gMEP_deg$log2FoldChange),decreasing = T),]

write.csv(tam.vs.gMEP_deg,'DESeq2_goodTAM.vs.goshMEP_topDEGs.csv')
tam.vs.gMEP_deg = read.csv('DESeq2_goodTAM.vs.goshMEP_topDEGs.csv')





##---------   Using edgeR   -----------####



toc = cnt_mtx[rownames(cnt_mtx) %in% geneMap$ensID,]
mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),]
coords = gns[rownames(cnt_mtx)]
cellTypeName='TAM'
formula='~ %s'
tgtChrs = paste0('chr',c(1:22))
donorID='donorID'
group='annot'

## Drop irrelevant genes
#Uninteresting genes
w = which(as.character(seqnames(coords)) %in% tgtChrs)
coords = coords[w]
seqlevels(coords) = tgtChrs
toc = toc[coords$gene_id,]
#Non-expressed genes
w = which(rowSums(toc>0)>=3)
coords = coords[w]
toc = toc[w,]

##=========================##
# Get genomic coordinates
#Order aa and bb by genomic coords
#If it's positive stranded or unknown use start, else end
coords$TSS = ifelse(strand(coords)=='-',end(coords),start(coords))
o = order(seqnames(coords),coords$TSS)
coords = coords[o]
toc = toc[o,]

##=========================##
# Create pseudobulk 
pb = do.call(cbind,lapply(split(colnames(toc),mDat[,donorID]),function(e) rowSums(toc[,e,drop=FALSE])))
colDat = mDat[match(colnames(pb),mDat[,donorID]),]
colDat[[group]] = factor(colDat[[group]],c('MEP','Tumour'))
rownames(colDat) = colDat[,donorID]


##=========================##
# Fit EDGER model
pdf(file.path(outDir,paste0(tgtCell,'_',ageGroup,'_pbEdgeR_plots.pdf')))
tam.vs.gMEP_res = fit_model(pb,colDat,formula = '~ %s',geneMap,groupID=group,donorID=donorID,MDS_groups = c('annot','donorID','sex'),coef=2)
dev.off()



fit_model <- function(pb,colDat,formula,geneMap,groupID='group',donorID='donorID',MDS_groups = c('Genotype','donorID','sex'),coef=NULL){
  require(edgeR)
  if(!grepl('%s',formula,fixed=TRUE) || grepl('%s',sub('%s','',formula,fixed=TRUE),fixed=TRUE))
    stop("Invalid formula, must contain one instance of %s")
  #Convert it now in case it is still invalid in more subtle ways
  formula = as.formula(sprintf(formula,groupID))
  
  
  # create an edgeR object with counts and grouping factor
  y <- DGEList(pb, group = colDat[[groupID]])
  y$samples = cbind(y$samples,colDat[match(rownames(y$samples),colDat[[donorID]]),!colnames(colDat) %in% colnames(y$samples)])
  
  ## Plot library sizes
  df = y$samples
  df$donorID = rownames(df)
  df$group = df[[groupID]]
  p = ggplot(df,aes(reorder(donorID,`lib.size`),fill=group,y=`lib.size`))+
    geom_col()+
    scale_fill_manual(values = col25)+
    theme_classic(base_size = 12)+xlab('')+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          axis.line = element_blank(),axis.text = element_text(color='black'),
          panel.border = element_rect(fill=F,colour = 'black'))
  
  
  
  print(p)
  
  
  
  # filter out genes with low counts
  print("Dimensions before subsetting:")
  print(dim(y))
  print("")
  all(y$samples$lib.size>1e5)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  print("Dimensions after subsetting:")
  print(dim(y))
  print("")
  
  
  
  ## Normalization 
  y <- calcNormFactors(y)
  # MDS plot
  cluster<-as.factor(colDat[[groupID]][match(rownames(y$samples),colDat[[donorID]])])
  df = plotMDS(y, pch=16,col=col25[cluster], main="MDS") 
  df = data.frame(x = df$x,y=df$y,donorID = names(df$x))
  df = cbind(df,colDat[match(df$donorID,colDat$donorID),!colnames(colDat) %in% colnames(df)])
  for(f in MDS_groups){
    df$group = df[[f]]
    p = ggplot(df,aes(x,y,col=group))  +
      geom_point()+
      theme_classic(base_size = 10)+
      scale_color_manual(values = col25)+
      theme(axis.line = element_blank(),axis.text = element_text(color='black'),
            panel.border = element_rect(fill=F,colour = 'black'))
    print(p)
  }
  
  
  # create a design matrix: here we have multiple donors so also consider that in the design matrix
  design <- model.matrix(formula,data=y$samples)
  
  # estimate dispersion
  y <- estimateDisp(y, design = design,robust = T)
  # fit the model
  fit <- glmQLFit(y, design)
  plotQLDisp(fit)
  
  # Extract coefficients
  if(is.null(coef)){
    contrast = makeContrasts(groupT21-groupdiploid, levels=colnames(design))  
    qlf<-glmQLFTest(fit, contrast=contrast)
  }else{
    qlf<-glmQLFTest(fit, coef = coef)
  }
  
  
  # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
  tt <- topTags(qlf, n = Inf,p.value = 0.05,)
  tt <- tt$table
  tt = annotateGenes(tt,geneMap = geneMap)
  
  ## Calculate logCPM
  y$logCPM <- cpm(y, log=TRUE, prior.count = 1)
  
  return(list("fit"=fit, "design"=design, "y"=y, 'tt'=tt))
}






##---------------------------------------------------##
##   3. DEGs between MLDS vs gosh_T21 MEMP_MEP     ####
##---------------------------------------------------##
##------- Using Find Markers -----------
#DimPlot(big.srat,cells.highlight = big.srat$cellID[big.srat$finalAnn_broad == 'MEP' & big.srat$disease == 'BALL'])

##----    FM_1: comparing all_MEMP vs all_MLDS in 1 comparison
mlds$group_tmp = ifelse(mlds$disease == 'MLDS' & mlds$timePoint == 'Diagnostic' & mlds$annot =='Tumour' & mlds$donorID != 'L041','MLDS',
                        ifelse(mlds$annot == 'MEP' & mlds$donorID %in% c('L019','L038','L039','L040'),'MEP','others'))
table(mlds$group_tmp,mlds$donorID)

Idents(mlds) = mlds$group_tmp
mlds.vs.goshMEP.markers = FindMarkers(mlds,ident.1 = 'MLDS',ident.2 = 'MEP',logfc.threshold = 0.1,min.pct = 0.01)
mlds.vs.goshMEP.markers$geneSym = rownames(mlds.vs.goshMEP.markers)
mlds.vs.goshMEP.markers$ensID = geneMap$ensID[match(mlds.vs.goshMEP.markers$geneSym,geneMap$geneSym)]
rownames(mlds.vs.goshMEP.markers) = mlds.vs.goshMEP.markers$ensID
mlds.vs.goshMEP.markers = annotateGenes(mlds.vs.goshMEP.markers,geneMap = geneMap)

mlds.vs.goshMEP.markers$pct.diff = mlds.vs.goshMEP.markers$pct.1 - mlds.vs.goshMEP.markers$pct.2
mlds.vs.goshMEP.markers$direction = ifelse(mlds.vs.goshMEP.markers$avg_log2FC > 0 ,'MLDS_up','MLDS_down')
#mlds.vs.goshMEP.markers = mlds.vs.goshMEP.markers[,!colnames(mlds.vs.goshMEP.markers) %in% c('X.1','X')]
rownames(mlds.vs.goshMEP.markers) = mlds.vs.goshMEP.markers$ensID
write.csv(mlds.vs.goshMEP.markers,'FM_allMLDS.vs.goshMEP_allGenev2.csv',row.names = T)



##------- Using pseudobulk -----------
cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in% 
                                                c(mlds$cellID[mlds$disease == 'MLDS' & mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & !mlds$donorID %in% c('L041','CC3') & mlds$tissue == 'BM'],
                                                  mlds$cellID[mlds$donorID %in% c('L019','L039','L040','L178') & mlds$timePoint %in% c('TP1','TP2','TP4') & mlds$annot == 'MEP'])]]
rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex')]
mDat$annot = mDat$annot
mDat$donorID = as.character(mDat$donorID)
mDat$donorID[mDat$annot == 'MEP'] = paste0(mDat$donorID[mDat$annot == 'MEP'],':MEP')
mDat = mDat[match(colnames(cnt_mtx),mDat$cellID),]
table(mDat$donorID,mDat$annot)

mlds.vs.gMEP_res = compareCell_simplified(toc = cnt_mtx,
                                         mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                                         coords = gns[rownames(cnt_mtx)],
                                         cellTypeName='MLDS',
                                         formula='~ %s',
                                         tgtChrs = paste0('chr',c(1:22)),
                                         donorID='donorID',groupID='annot')

saveRDS(mlds.vs.gMEP_res,'DESeq2_dMLDS.vs.goshMEP_res.RDS')
mlds.vs.gMEP_res = readRDS('DESeq2_dMLDS.vs.goshMEP_res.RDS')


mlds.vs.gMEP_dds = mlds.vs.gMEP_res[['dds']]
mlds.vs.gMEP_deg = mlds.vs.gMEP_res[['mainDE']]
mlds.vs.gMEP_deg$direction = ifelse(mlds.vs.gMEP_deg$log2FoldChange > 0, 'MLDS_up','MLDS_down')
mlds.vs.gMEP_deg$cellFrac.diff = mlds.vs.gMEP_deg$cellFrac_g1 - mlds.vs.gMEP_deg$cellFrac_g2
mlds.vs.gMEP_deg$cellFrac_max = pmax(mlds.vs.gMEP_deg$cellFrac_g1, mlds.vs.gMEP_deg$cellFrac_g2)

# max_cellfrac = 20
# min_l2FC = 0.5

mlds.vs.gMEP_deg = mlds.vs.gMEP_deg[mlds.vs.gMEP_deg$padj < 0.05 & 
                                      abs(mlds.vs.gMEP_deg$log2FoldChange) >= min_l2FC &
                                      abs(mlds.vs.gMEP_deg$cellFrac_max) >= max_cellfrac,]

# mlds.vs.gMEP_deg = mlds.vs.gMEP_deg[mlds.vs.gMEP_deg$padj < 0.05 & 
#                                       abs(mlds.vs.gMEP_deg$log2FoldChange) >= min_l2FC &
#                                       abs(mlds.vs.gMEP_deg$pct.diff) >= 10,]

table(mlds.vs.gMEP_deg$direction)
dim(mlds.vs.gMEP_deg)
mlds.vs.gMEP_deg = mlds.vs.gMEP_deg[order(abs(mlds.vs.gMEP_deg$log2FoldChange),decreasing = T),]

write.csv(mlds.vs.gMEP_deg,'DESeq2_dMLDS.vs.goshMEP_topDEGs.csv')
mlds.vs.gMEP_deg = read.csv('DESeq2_dMLDS.vs.goshMEP_topDEGs.csv')



##---------   Using edgeR   -----------####

toc = cnt_mtx[rownames(cnt_mtx) %in% geneMap$ensID,]
mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),]
coords = gns[rownames(cnt_mtx)]
cellTypeName='MLDS'
formula='~ %s'
tgtChrs = paste0('chr',c(1:22))
donorID='donorID'
group='annot'

## Drop irrelevant genes
#Uninteresting genes
w = which(as.character(seqnames(coords)) %in% tgtChrs)
coords = coords[w]
seqlevels(coords) = tgtChrs
toc = toc[coords$gene_id,]
#Non-expressed genes
w = which(rowSums(toc>0)>=3)
coords = coords[w]
toc = toc[w,]

##=========================##
# Get genomic coordinates
#Order aa and bb by genomic coords
#If it's positive stranded or unknown use start, else end
coords$TSS = ifelse(strand(coords)=='-',end(coords),start(coords))
o = order(seqnames(coords),coords$TSS)
coords = coords[o]
toc = toc[o,]

##=========================##
# Create pseudobulk 
pb = do.call(cbind,lapply(split(colnames(toc),mDat[,donorID]),function(e) rowSums(toc[,e,drop=FALSE])))
colDat = mDat[match(colnames(pb),mDat[,donorID]),]
colDat[[group]] = factor(colDat[[group]],c('MEP','Tumour'))
rownames(colDat) = colDat[,donorID]


##=========================##
# Fit EDGER model
pdf(file.path(outDir,paste0(tgtCell,'_',ageGroup,'_pbEdgeR_plots.pdf')))
mlds.vs.gMEP_res_edgeR = fit_model(pb,colDat,formula = '~ %s',geneMap,groupID=group,donorID=donorID,MDS_groups = c('annot','donorID','sex'),coef=2)
dev.off()

## compare to DESeq2 results
library(UpSetR)
upset(fromList(list('edgeR_up' = mlds.vs.gMEP_res_edgeR[['tt']]$geneSym[mlds.vs.gMEP_res_edgeR[['tt']]$logFC > 0],
                    'edgeR_down' = mlds.vs.gMEP_res_edgeR[['tt']]$geneSym[mlds.vs.gMEP_res_edgeR[['tt']]$logFC < 0],
                    'deseq2_up' = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_up'],
                    'deseq2_down' = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_down'])))



##----  Gene set enrichment analysis with limma::camera() ####
## Obtain the gene sets
C2_genesets <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.c2.all.v7.1.entrez.rds"))
HM_genesets <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds")) 

library(org.Hs.eg.db)
entrez_geneid <- mapIds(org.Hs.eg.db, keys=rownames(fit), 
                        column="ENTREZID", 
                        keytype="ENSEMBL")
idx_cam <- ids2indices(c(C2_genesets,HM_genesets),id=entrez_geneid)

## Import the DEG object
out = readRDS(file.path(outDir,'pb_edgeR_byCT.RDS'))
fit = mlds.vs.gMEP_res_edgeR[['fit']]
fit<-glmQLFTest(fit, coef = 2)

summary(decideTests(fit))


## Run camera
cam_MLDS = limma::camera(mlds.vs.gMEP_res_edgeR[['y']]$logCPM, index=idx_cam, design=mlds.vs.gMEP_res_edgeR[['y']]$design, contrast=2) 
cam_MLDS$group = 'MLDS'
# cam_fTFC2 <- limma::camera(out[[2]][['y']]$logCPM, index=idx_cam, design=out[[2]][['y']]$design, contrast=contrast) 
# cam_fTFC2$group = 'fTFC2'

library(UpSetR)
upset(fromList(list('fTFC1_up' = rownames(cam_MLDS)[cam_MLDS$FDR < 0.05 & cam_MLDS$Direction == 'Up'],
                    'fTFC1_down' = rownames(cam_MLDS)[cam_MLDS$PValue < 0.05 & cam_MLDS$Direction == 'Down'])))








##-----------------------------##
##  4. GATA1s module in MLDS ####
##-----------------------------##
mlds.vs.gMEP_deg = read.csv('DESeq2_dMLDS.vs.goshMEP_topDEGs.csv')

## How much overlap is there in the list of TAM.vs.T21 and MLDS.vs.T21?
upset(fromList(list(tam_up = tam.vs.gMEP_deg$geneSym[tam.vs.gMEP_deg$direction == 'TAM_up'],
                    tam_down = tam.vs.gMEP_deg$geneSym[tam.vs.gMEP_deg$direction == 'TAM_down'],
                    mlds_up = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_up'],
                    mlds_down = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_down'])),
      nsets = 6,text.scale = 1.5)



## Of the MLDS-DEGs, which ones were also DEGs in TAM?

#-- Import raw gTAM vs gMEP result
tam.vs.gMEP_res = readRDS('DESeq2_goodTAM.vs.goshMEP_res.RDS')
tam.vs.gMEP_results = tam.vs.gMEP_res[['res']]
tam.vs.gMEP_results = annotateGenes(tam.vs.gMEP_results,geneMap = geneMap)

cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in% 
                                                c(mlds$cellID[mlds$donorID %in% c('CC1','CC2','L075','L182','L114','CC6','CC8','CC7') & mlds$annot == 'Tumour'],
                                                  mlds$cellID[mlds$donorID %in% c('L019','L039','L040','L178') & mlds$timePoint %in% c('TP1','TP2','TP4') & mlds$annot == 'MEP'])]]
rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]
mDat = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex')]
mDat$donorID = as.character(mDat$donorID)
mDat = mDat[match(colnames(cnt_mtx),mDat$cellID),]

tam.vs.gMEP_results$cellFrac_g1 = NA
tam.vs.gMEP_results$cellFrac_g2 = NA
for(i in split(seq(1:nrow(tam.vs.gMEP_results)),seq(1:nrow(tam.vs.gMEP_results)) %/% 1000)){
  print(max(i))
  tam.vs.gMEP_results$cellFrac_g1[i] = apply(cnt_mtx[match(tam.vs.gMEP_results$ensID[i],rownames(cnt_mtx)),colnames(cnt_mtx) %in% mDat$cellID[mDat$annot == 'Tumour']],1,function(x){sum(x>0)})
  tam.vs.gMEP_results$cellFrac_g2[i] = apply(cnt_mtx[match(tam.vs.gMEP_results$ensID[i],rownames(cnt_mtx)),colnames(cnt_mtx) %in% mDat$cellID[mDat$annot == 'MEP']],1,function(x){sum(x>0)})
}

tam.vs.gMEP_results$cellFrac_g1 = tam.vs.gMEP_results$cellFrac_g1 / ncol(cnt_mtx[,colnames(cnt_mtx) %in% mDat$cellID[mDat$annot == 'Tumour']])
tam.vs.gMEP_results$cellFrac_g2 = tam.vs.gMEP_results$cellFrac_g2 / ncol(cnt_mtx[,colnames(cnt_mtx) %in% mDat$cellID[mDat$annot == 'MEP']])  

tam.vs.gMEP_results$cellFrac.diff = tam.vs.gMEP_results$cellFrac_g1 - tam.vs.gMEP_results$cellFrac_g2
tam.vs.gMEP_results$cellFrac_max = pmax(tam.vs.gMEP_results$cellFrac_g1, tam.vs.gMEP_results$cellFrac_g2)


write.csv(tam.vs.gMEP_results,'DESeq2_goodTAM.vs.goshMEP_result_allGene.csv')
tam.vs.gMEP_results = read.csv('DESeq2_goodTAM.vs.goshMEP_result_allGene.csv')


#-- Subset to keep only DEGs from previous comparison 
# gata1s_module = rbind(tam.vs.gMEP_results[abs(tam.vs.gMEP_results$cellFrac_max) >= max_cellfrac & tam.vs.gMEP_results$log2FoldChange > 0 & tam.vs.gMEP_results$geneSym %in% mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_up'] & !is.na(tam.vs.gMEP_results$padj),],
#                       tam.vs.gMEP_results[abs(tam.vs.gMEP_results$cellFrac_max) >= max_cellfrac & tam.vs.gMEP_results$log2FoldChange < 0 & tam.vs.gMEP_results$geneSym %in% mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_down'] & !is.na(tam.vs.gMEP_results$padj),])
gata1s_module = rbind(tam.vs.gMEP_results[tam.vs.gMEP_results$log2FoldChange > 0 & tam.vs.gMEP_results$geneSym %in% mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_up'] & !is.na(tam.vs.gMEP_results$padj),],
                      tam.vs.gMEP_results[tam.vs.gMEP_results$log2FoldChange < 0 & tam.vs.gMEP_results$geneSym %in% mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_down'] & !is.na(tam.vs.gMEP_results$padj),])

#-- re-calculate padj
gata1s_module$padj = p.adjust(gata1s_module$pvalue,method = 'BH')
table(gata1s_module$padj < 0.05)
gata1s_module$direction = ifelse(gata1s_module$log2FoldChange > 0 & gata1s_module$padj < 0.05,'TAM_up',
                                       ifelse(gata1s_module$log2FoldChange < 0 & gata1s_module$padj < 0.05,'TAM_down','nonDE'))

table(gata1s_module$chr,gata1s_module$log2FoldChange > 0)
table(gata1s_module$log2FoldChange > 0,gata1s_module$direction)
table(gata1s_module$padj < 0.05)

gata1s_module = gata1s_module[order(abs(gata1s_module$log2FoldChange),decreasing = T),]
write.csv(gata1s_module,'DESeq2_goodTAM.vs.goshMEP_MLDSgenes.csv',row.names = T)
gata1s_module = read.csv('DESeq2_goodTAM.vs.goshMEP_MLDSgenes.csv',row.names = 1)




# # table(abs(gata1s_module$log2FoldChange) >= 0.5)
# # gata1s_module$direction[abs(gata1s_module$log2FoldChange) < 0.5] = 'nonDE'
# o1 = gata1s_module
# o2 = gata1s_module
# 
# dim(o2)
# 
# gata1s_module = tam.vs.gMEP_results[rownames(tam.vs.gMEP_results) %in% mlds.vs.gMEP_deg$ensID,]
# gata1s_module$padj = p.adjust(gata1s_module$pvalue,method = 'BH')
# table(gata1s_module$padj < 0.05)
# gata1s_module$direction = ifelse(gata1s_module$log2FoldChange > 0 & gata1s_module$padj < 0.05,'TAM_up',
#                                        ifelse(gata1s_module$log2FoldChange < 0 & gata1s_module$padj < 0.05,'TAM_down','nonDE'))
# 
# table(abs(gata1s_module$log2FoldChange) >= 0.5)
# gata1s_module$direction[abs(gata1s_module$log2FoldChange) < 0.5] = 'nonDE'
# 
# 
# 
# write.csv(tam.vs.gMEP_results,'DESeq2_goodTAM.vs.goshMEP_MLDSgenes.csv',row.names = T)
# tam.vs.gMEP_results = read.csv('DESeq2_goodTAM.vs.goshMEP_MLDSgenes.csv',row.names = 1)






#---- DotPlot expression of these GATA1s genes  ----####
genes_toPlot = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 0.3 & !grepl('^LINC\\d+|^AC\\d+',gata1s_module$geneSym),]
genes_toPlot = mlds.vs.gMEP_deg[abs(mlds.vs.gMEP_deg$cellFrac.diff) >= 0.2 & 
                                  abs(mlds.vs.gMEP_deg$cellFrac.diff) < 0.3 &
                                  mlds.vs.gMEP_deg$tam_vs_mempT21_group !='notDE',]
genes_toPlot$group = genes_toPlot$tam_vs_mempT21_group

genes_toPlot = mlds_vs_goodTAM

genes_toPlot = gata1s_module


## Plot expression of top GATA1s modules - DotPlot
mlds$group = ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L038',
                    paste0(mlds$annot,':',mlds$disease,':',mlds$donorID,':',mlds$timePoint),
                    ifelse(mlds$annot %in% c('MEP','Tumour'),paste0(mlds$annot,':',mlds$disease,':',mlds$donorID),
                    mlds$annot))
Idents(mlds) = mlds$group
mlds.vs.gMEP_deg = mlds.vs.gMEP_deg[order(abs(mlds.vs.gMEP_deg$cellFrac.diff),decreasing = T),]

DotPlot(mlds,#idents = unique(mlds$group[grepl('^Tumour:TAM|^Tumour:MLDS|^MEP:MLDS:L019|^MEP:MLDS:L038|^MEP:MLDS:L039|^MEP:MLDS:L040',mlds$group) & mlds$donorID != 'L041']),
        idents = unique(mlds$group[!grepl('NA|unsure_|\\?',mlds$group)]),
        features = c(genes_toPlot$geneSym[genes_toPlot$group == 'TAM.MLDS.up'][1:30],
                     genes_toPlot$geneSym[genes_toPlot$group == 'TAM.MLDS.down'][1:30]),
        #features = c(gata1s_module$geneSym[gata1s_module$GATA2_target == '1']),
        #features = c(gata1_network$target[gata1_network$mor == -1]),
        #features = c(gata1s_module$geneSym[gata1s_module$low_FLexpr == T & gata1s_module$isCSM==T & gata1s_module$direction == 'MLDS_up']),
        
        group.by = 'group',scale = T)+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 


##----- In otherLeukaemia seurat object
big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/combinedLeuk_2408.RDS')
DotPlot(big.srat,#idents = unique(big.srat$group_4_dotPlot[grepl('^Tumour.*RD|Tumour.*D',big.srat$group_4_dotPlot)]),
        idents = unique(big.srat$group_4_dotPlot[!grepl('other',big.srat$group_4_dotPlot)]),
        # features = c(genes_toPlot$geneSym[genes_toPlot$group == 'TAM.MLDS.up'][1:100],
        #              genes_toPlot$geneSym[genes_toPlot$group == 'TAM.MLDS.down'][1:32]))+
        features = c(genes_toPlot$geneSym[genes_toPlot$direction == 'MLDS_down'],
                     genes_toPlot$geneSym[genes_toPlot$direction == 'MLDS_up']))+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 











#--- Plot fraction of MLDS-DEGs being DEGs in TAM
mlds.vs.gMEP_deg = mlds.vs.gMEP_deg[order(abs(mlds.vs.gMEP_deg$log2FoldChange),decreasing = T),]
mlds.vs.gMEP_deg$tam_vs_mempT21 = gata1s_module$direction[match(mlds.vs.gMEP_deg$geneSym,gata1s_module$geneSym)]
mlds.vs.gMEP_deg$tam_vs_mempT21_group = 'notDE'
mlds.vs.gMEP_deg$tam_vs_mempT21_group[mlds.vs.gMEP_deg$direction == 'MLDS_down' & mlds.vs.gMEP_deg$tam_vs_mempT21 == 'TAM_down'] = 'TAM.MLDS.down'
mlds.vs.gMEP_deg$tam_vs_mempT21_group[mlds.vs.gMEP_deg$direction == 'MLDS_up' & mlds.vs.gMEP_deg$tam_vs_mempT21 == 'TAM_up'] = 'TAM.MLDS.up'
table(mlds.vs.gMEP_deg$tam_vs_mempT21_group)

write.csv(mlds.vs.gMEP_deg,'DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv')
mlds.vs.gMEP_deg = read.csv('DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv')

# write.csv(mlds.vs.gMEP_deg,'DESeq2_goshMEP.allTAM.MLDS_topGeneModule.csv')
# mlds.vs.gMEP_deg = read.csv('DESeq2_goshMEP.allTAM.MLDS_topGeneModule.csv')



ggplot(mlds.vs.gMEP_deg,aes(direction,fill = tam_vs_mempT21_group))+
  geom_bar(position = 'fill',width = 0.7)+
  scale_fill_manual(values = c(grey(0.6),col25[1],col25[2]),name='TAM vs GOSH_T21_MEP')+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme_classic(base_size = 15)+xlab('')+ylab('Fraction of DEGs')+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank(),#legend.position = 'bottom',
        legend.title = element_text(size=10),legend.text = element_text(size=8),legend.key.size = unit(0.5,'cm'),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

table(mlds.vs.gMEP_deg$tam_vs_mempT21_group,mlds.vs.gMEP_deg$direction)

## Make a strange plot
df = as.data.frame(table(mlds.vs.gMEP_deg$direction,mlds.vs.gMEP_deg$tam_vs_mempT21_group))
colnames(df) = c('MLDS_deg','TAM_deg','nGene')
df2 = df %>% group_by(MLDS_deg) %>% summarise(nGene = sum(nGene))
df2$group = 'MLDS vs T21_MEMP'
colnames(df2)[1] = 'direction'
df = df[df$TAM_deg == 'TAM.MLDS.down' & df$MLDS_deg == 'MLDS_down' |
          df$TAM_deg == 'TAM.MLDS.up' & df$MLDS_deg == 'MLDS_up',]
df = df[,colnames(df) != 'MLDS_deg']
df$group = 'TAM vs T21_MEMP'
colnames(df)[1] = 'direction'
df = rbind(df,df2)
df$nGene[grepl('down',df$direction)] = -df$nGene[grepl('down',df$direction)]
df$nGene_log10 = log10(abs(df$nGene))
df$nGene_log10[df$direction %in% c('TAM.MLDS.down','MLDS_down')] = -df$nGene_log10[df$direction %in% c('TAM.MLDS.down','MLDS_down')]


plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots/noCC3'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
plotFun_MLDS.TAM_DEGs = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(df,aes(group,nGene_log10,fill=direction))+
    geom_col(width = 0.6)+
    #scale_fill_manual(values = c(col25[1],col25[2],colAlpha(col25[1:2],alphas = 0.6)))+
    scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
    geom_hline(yintercept = 0,col='black',lwd=0.5)+
    geom_hline(yintercept = df$nGene_log10,col='black',lty=2)+
    scale_y_continuous(limits = c(-log10(300),log10(650)),
                       breaks = c(-log10(c(200,10)),0,log10(c(10,200,400,600))),
                       labels = as.character(c(-200,10,0,10,200,400,600)))+
    theme_classic(base_size = 11)+xlab('')+ylab('# DEGs')+
    theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
          axis.line = element_blank(),#legend.position = 'bottom',
          axis.ticks = element_line(colour = 'black'),
          legend.title = element_text(size=10,colour = 'black'),
          legend.text = element_text(size=8,colour = 'black'),legend.key.size = unit(0.5,'cm'),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = 'black'),
          axis.text = element_text(color='black'))
  print(p1)
}

#saveFig(file.path(plotDir,paste0('Fig3x_MLDS.goodTAM.DEGs_barplot')),plotFun_MLDS.TAM_DEGs,rawData=df,width = 3.5,height = 4,res = 500,useDingbats = T)
saveFig(file.path(plotDir,paste0('Fig3B_MLDS.goodTAM.DEGs_barplot')),plotFun_MLDS.TAM_DEGs,rawData=df,width = 3.5,height = 4,res = 500)





df = df[df$group == 'TAM vs T21_MEMP',]
plotFun_MLDS.TAM_DEGs_barplot_onlyTAM = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(df,aes(group,nGene,fill=direction))+
    geom_col(width = 0.58)+
    #scale_fill_manual(values = c(col25[1],col25[2],colAlpha(col25[1:2],alphas = 0.6)))+
    scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
    geom_hline(yintercept = 0,col='black',lwd=0.5)+
    # scale_y_continuous(limits = c(-log10(250),log10(600)),
    #                    breaks = c(-log10(c(200,10)),0,log10(c(10,200,400,600))),
    #                    labels = as.character(c(-200,10,0,10,200,400,600)))+
    scale_y_continuous(limits = c(-250,600),
                       breaks = c(-200,0,200,400,600),
                       labels = c(-200,0,200,400,600))+
    #geom_hline(yintercept = df$nGene_log10,col='black',lty=2)+
    theme_classic(base_size = 11)+xlab('')+ylab('# DEGs')+
    theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
          axis.line = element_blank(),#legend.position = 'bottom',
          axis.ticks = element_line(colour = 'black'),
          legend.title = element_text(size=10,colour = 'black'),
          legend.text = element_text(size=8,colour = 'black'),legend.key.size = unit(0.5,'cm'),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = 'black'),
          axis.text = element_text(color='black'))
  print(p1)
}

saveFig(file.path(plotDir,paste0('Fig3B_MLDS.goodTAM.DEGs_barplot_simplified')),plotFun_MLDS.TAM_DEGs_barplot_onlyTAM,rawData=df,width = 2.8,height = 5,res = 500)









##---------------------------------------##
##  5. UCell Scoring in other scLeuk   ####
##---------------------------------------##
source('~/lustre_mt22/Aneuploidy/scripts/finalScripts/xx01_moduleScore_helperFunctions.R')

## Define gene module
#gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/oct24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)

table(gata1s_module$tam_vs_mempT21_group)
gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)

## Select top genes only
gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 30/100,]
geneList = split(gata1s_module$geneSym,gata1s_module$group)
module_type = 'GATA1s'


## Score the module - using xx01_moduleScore_terminal
## Define gene module
module_type = 'GATA1s_topGenes2' 
data = import_UCell_result(outDir= '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/UCell_moduleScore',
                           module_type=module_type)


module_name = colnames(data)[grepl('_UCell$',colnames(data))]#paste0(names(geneList)[1],'_UCell')
if(module_type %in% c('GATA1s','GATA1s_topGenes','GATA1s_topGenes2')){
  data$module = data[[module_name[4]]]
}else if(module_type %in% c('badTAM','badTAM_topGenes','MLDS_T21imprint')){
  data$module = data[[module_name[2]]]
}else{
  data$module = data[[module_name[1]]]
}

data = data[data$dataset_ID %in% c('fBM','fBM.T21','fLiver','fLiverMuz','MLDS','otherLeuk'),]

data$dataset_ID_og = data$dataset_ID
data$dataset_ID[data$dataset_ID == 'fLiver' & data$Genotype == 'diploid' & !grepl('^Hsb',data$donorID)] = 'fLiver_ext'
data$dataset_ID[data$dataset_ID == 'fLiverMuz'] = 'fLiver_ext'
data$dataset_ID[data$dataset_ID == 'fBM.T21'] = 'fBM'
data$dataset_ID[data$dataset_ID == 'fBM.T21'] = 'fBM'
data$dataset_ID[data$dataset_ID == 'MLDS' & data$donorID %in% c('CC1','CC2','L075','L114','L182','L156','CC6','CC7','CC8')] = 'TAM'
#data$dataset_ID[data$dataset_ID == 'infantALL' & !is.na(data$donorID) & data$donorID == 'P9_iAML/Cancer'] = 'AMKL'
data$dataset_ID[data$dataset_ID == 'otherLeuk'] = data$disease[data$dataset_ID == 'otherLeuk']
# data$dataset_ID[data$dataset_ID == 'otherLeuk' & data$donorID == 'P9'] = 'AMKL'
# data$dataset_ID[data$orig.ident %in% c('NB14406184','NB14406185')] = 'L062_LPD'
# data$dataset_ID[data$orig.ident %in% c('NB14406184','NB14406185')] = 'L062_LPD'
#data$donorID = gsub('_iA.*\\/Cancer$','',data$donorID)
data$finalAnn_broad[data$finalAnn_broad %in% c('Tumour','TumourCNA')] = 'Cancer'
data$finalAnn_broad[data$finalAnn_broad %in% c('MEP','MEMP')] = 'MEMP_MEP'
data$dataset_ID = factor(data$dataset_ID,c('fLiver','fLiver_ext','fBM','TAM','MLDS','MDS','AMKL','infantALL','pAML','pBALL','LPD'))

# Define baseline = median score in fLiver_MEMP_inhouse2n + fLiver_MEMP_T21
baseLine = data[data$dataset_ID %in% c('fLiver') & data$Genotype %in% c('diploid','T21'),] %>% group_by(Genotype) %>% summarise(medScore = median(module))
# Only keep things we want to plots
data$finalAnn_broad[data$cellID %in% mlds$cellID] = mlds$annot[match(data$cellID[data$cellID %in% mlds$cellID],mlds$cellID)]
data$finalAnn_broad[data$finalAnn_broad == 'Tumour' & !is.na(data$finalAnn_broad)]  = 'Cancer'
df = data[data$finalAnn_broad %in% c('Cancer','MEMP_MEP'),]
df = df[!(df$donorID == 'L041' & df$finalAnn_broad == 'Cancer'),]
df = df[!(df$finalAnn_broad == 'MEMP_MEP' & df$dataset_ID %in% c('pBALL','pAML','MLDS','TAM','MDS','AMKL','infantALL','LPD')),]

df$group = ifelse(df$finalAnn_broad != 'Cancer',as.character(df$finalAnn_broad),
                  ifelse(df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Cancer','TAM / MLDS',
                         ifelse(df$finalAnn_broad == 'Cancer','Other leukaemia','others')))
df$group = factor(df$group,c('MEMP_MEP','TAM / MLDS','Other leukaemia'))
df$Genotype[df$Genotype == 'complete_trisomy'] = 'Triploid'
df$Genotype[df$Genotype == 'diploid'] = 'Diploid'
df$Genotype = factor(df$Genotype,c('Diploid','T21','T18','T22','MX','Triploid'))

# Remove infantALL datast
df = df[df$dataset_ID != 'infantALL',]

df$groupID = ifelse(!df$group %in% c('Other leukaemia','TAM / MLDS'),paste0(df$group,'_',df$dataset_ID,':',df$Genotype),paste0(df$group,'_',df$dataset_ID,':',df$Genotype,':',df$donorID))
d = df
d = d %>% group_by(groupID) %>% mutate(nCell_perGroupID = n())
table(d$groupID[d$nCell_perGroupID <= 30])
## Remove groups with <= 30 cells
d = d[d$nCell_perGroupID >= 30,]

## Calculate boxplot quantiles
quant_df = data.frame()
for(i in 1:n_distinct(d$groupID)){
  dd = d[d$groupID == unique(d$groupID)[i],]
  quant = quantile(dd$module) %>% t() %>% as.data.frame()
  quant$groupID = unique(d$groupID)[i]
  quant_df = rbind(quant_df,quant)
}
quant_df$group = d$group[match(quant_df$groupID,d$groupID)]
quant_df$dataset_ID = d$dataset_ID[match(quant_df$groupID,d$groupID)]
quant_df$Genotype = d$Genotype[match(quant_df$groupID,d$groupID)]




##--------------------------------------------------------------------##
##  Statistical test for the module score across different groups   ####
##--------------------------------------------------------------------##
# I think the appropriate statistical test for the single-cell module score is a nested ANOVA or mixed-model ANOVA
# Response var Y = module score for individual cell
# Predictor X = Disease_group, with dataset --> Patient_ID being nested within the each of the disease group
library(lme4)
dd = df[,c('cellID',"donorID",'group','dataset_ID','module')]
lmm1 = nlme::lme(module ~ group, random = ~ 1|dataset_ID:donorID, data = dd, method="REML")
summary(lmm1)
AIC(lmm1)
anova(lmm1)

## compare the median of each group
t.test(quant_df$`50%`[quant_df$group == 'TAM / MLDS'],quant_df$`50%`[quant_df$group == 'Other leukaemia'],alternative = 'greater')
t.test(quant_df$`50%`[quant_df$group == 'TAM / MLDS'],quant_df$`50%`[quant_df$group == 'MEMP_MEP'],alternative = 'greater')
t.test(quant_df$`50%`[quant_df$dataset_ID == 'TAM'],quant_df$`50%`[quant_df$dataset_ID == 'MLDS'])





##----------------------------------------------------##
##  6. Functional assessment of the GATA1s module   ####
##----------------------------------------------------##
gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/oct24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
table(gata1s_module$tam_vs_mempT21_group)
gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)

##----    EnrichR   -----## --> totally NOT informative
upDEG_enriched <- enrichr(gata1s_module$geneSym[gata1s_module$group == 'TAM.MLDS.up'], dbs)
View(upDEG_enriched[[6]][upDEG_enriched[[6]]$Adjusted.P.value < 0.05,])
plotEnrich(upDEG_enriched[[3]][upDEG_enriched[[3]]$Adjusted.P.value < 0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")

downDEG_enriched <- enrichr(unique(gata1s_module$geneSym[gata1s_module$group == 'TAM.MLDS.down']), dbs)
View(downDEG_enriched[[6]][downDEG_enriched[[6]]$Adjusted.P.value < 0.05,])
plotEnrich(downDEG_enriched[[5]][downDEG_enriched[[5]]$Adjusted.P.value < 0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")


## UP regulated module
allTerms_up = extract_enrichR.res(upDEG_enriched,module_direction = 'up',pVal_cutoff = 0.1,min_overlapFrac = 0.05,
                                  db_toExtract = names(upDEG_enriched)[c(3,5,6)])

## Plot up-MSigDB_Hallmark_2020
plotFun_enrichR.up_pPTC_DEGs = function(noFrame=FALSE,noPlot=FALSE){
  p3=ggplot(allTerms_up,aes(Combined.Score,Term))+
    geom_segment(aes(x=0,xend=(Combined.Score),y=Term,yend=Term),col=grey(0.7))+
    geom_point(aes(size = nGene,col=overlap2))+
    facet_grid(db~.,scales = 'free_y',space = 'free_y')+
    xlab('Combined Score')+ylab('') + 
    ggtitle('GATA1 up-regulated DEGs in TAM / ML-DS')+
    scale_color_gradient(low='#eba9a9',high = '#a10505')+
    theme_classic(base_size = 11)+
    theme(panel.border = element_rect(fill=F),axis.line = element_blank())
  print(p3)
}

saveFig(file.path(outDir,'gata1s_module_enrichR.up_DEGs'),plotFun_enrichR.up_pPTC_DEGs,rawData=allTerms_up,width = 8.5,height = 9.5,res = 500)


## UP regulated module
allTerms_down = extract_enrichR.res(downDEG_enriched,module_direction = 'down',pVal_cutoff = 0.1,min_overlapFrac = 0.05,
                                  db_toExtract = names(downDEG_enriched)[c(3,5,6)])

plotFun_enrichR.down_pPTC_DEGs = function(noFrame=FALSE,noPlot=FALSE){
  p3=ggplot(allTerms_down,aes(Combined.Score,Term))+
    geom_segment(aes(x=0,xend=(Combined.Score),y=Term,yend=Term),col=grey(0.7))+
    geom_point(aes(size = nGene,col=overlap2))+
    facet_grid(db~.,scales = 'free_y',space = 'free_y')+
    xlab('Combined Score')+ylab('') + 
    ggtitle('GATA1 down-regulated DEGs in TAM / ML-DS')+
    scale_color_gradient(low='#eba9a9',high = '#a10505')+
    theme_classic(base_size = 11)+
    theme(panel.border = element_rect(fill=F),axis.line = element_blank())
  print(p3)
}

saveFig(file.path(outDir,'gata1s_module_enrichR.down_DEGs'),plotFun_enrichR.up_pPTC_DEGs,rawData=allTerms_up,width = 8.5,height = 6.5,res = 500)



##----    Chr21 / GATA1 downstream targets   -----##
#gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
gata1s_module = read.csv('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
table(gata1s_module$tam_vs_mempT21_group)
gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
# Subset to keep just GATA1s gene module
gata1s_module = gata1s_module[gata1s_module$tam_vs_mempT21_group != 'notDE',]
table(gata1s_module$group)
table(gata1s_module$chr,gata1s_module$direction)
View(gata1s_module[gata1s_module$chr == 'chr21',])

## Which ones are GATA1 targets 
library(decoupleR)
net = get_collectri(organism='human', split_complexes=FALSE)
head(net)

gata1_network = net[net$source == 'GATA1',]
gata1s_module$GATA1_target = gata1_network$mor[match(gata1s_module$geneSym,gata1_network$target)]
gata1s_module$GATA1_target[is.na(gata1s_module$GATA1_target)] = '-'


gata2_network = net[net$source == 'GATA2',]
gata1s_module$GATA2_target = gata2_network$mor[match(gata1s_module$geneSym,gata2_network$target)]
gata1s_module$GATA2_target[is.na(gata1s_module$GATA2_target)] = '-'

##  Which ones are not expressed in fLiver  ##
gata1s_module$low_FLexpr = !(gata1s_module$geneSym %in% rownames(mtx))

##----    Which cell types utilises these genes   -----##
## Plot expression of top GATA1s modules - Average expression plot

##---- In MLDS
mlds$timePoint[mlds$orig.ident == 'MY.200531.14635833'] = 'D.Relapse'
mlds$group_4avgExpr = ifelse(mlds$annot == 'MEP' & mlds$donorID %in% c('L019','L038','L039','L040'),'MEMP_MEP',
                             ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID !='L041' & mlds$tissue == 'BM',paste0(mlds$annot,'.',mlds$disease,'.',mlds$donorID),
                                    ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID !='L041' & mlds$disease == 'TAM',paste0(mlds$annot,'.',mlds$disease,'.',mlds$donorID),
                                           ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L038' & mlds$timePoint == 'TP1','Tumour.MLDS.L038_TP1',
                                                  ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L076' & mlds$timePoint == 'D.Relapse','Tumour.MLDS.L076_RD',
                                                         ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L076' & mlds$tissue == 'Blood','Tumour.MLDS.L076_Blood',
                                                                mlds$broadLineage))))))

mlds$group_4avgExpr = factor(mlds$group_4avgExpr,levels = c('MEMP_MEP',
                                                            unique(mlds$group_4avgExpr[grepl('TAM',mlds$group_4avgExpr)]),
                                                            unique(mlds$group_4avgExpr[grepl('MLDS',mlds$group_4avgExpr)]),
                                                            unique(mlds$group_4avgExpr[!grepl('TAM|MLDS|MEMP_MEP',mlds$group_4avgExpr)])))




mlds$disease_group = ifelse(mlds$group_4avgExpr == 'MEMP_MEP','MEMP_MEP',
                            ifelse(grepl('TAM',mlds$group_4avgExpr),'TAM',
                                   ifelse(grepl('MLDS',mlds$group_4avgExpr),'MLDS','others')))
avgExpr = AverageExpression(mlds,group.by = 'group_4avgExpr')
avgExpr = avgExpr[['RNA']]




##---- In fLiver
fLiver = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS')

genes_toPlot = gata1s_module[order(gata1s_module$fLiver_geneGroup),]
genes_toPlot = gata1s_module[gata1s_module$chr == 'chr21',]
genes_toPlot = gata1s_module[gata1s_module$GATA1_target != '-',]

fLiver$group_4_dotPlot = as.character(fLiver$annot)
fLiver$group_4_dotPlot = factor(fLiver$group_4_dotPlot,c('HSC_MPP','LMPP_ELP','MEMP_MEP','CMP_GMP',
                                                         'earlyMK','MK','Mast.cell',
                                                         'EE','ME','LE',
                                                         'proMono','Monocyte','promyelocyte','myelocyte',
                                                         'DC1','DC2','Macrophage','Kupffer.cell','pDC',
                                                         'ILC.precursor','NK_T','T.cell',
                                                         'pro.B.cell','pre.B.cell','B.cell',
                                                         'Endo','Fibroblast','Hepatocyte','Mesenchyme','NPC','doublets'))
Idents(fLiver) = fLiver$Genotype
DotPlot(fLiver,#idents = unique(fLiver$annot[grepl('B|MEP',fLiver$annot)]),
        idents = 'T21', group.by = 'group_4_dotPlot',
        #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) >= 30]
        #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) <= 30 & abs(mlds.vs.gMEP_deg$pct.diff) >= 20]
        # features = c(genes_toPlot$geneSym[genes_toPlot$group == 'TAM.MLDS.up'][1:30],
        #              genes_toPlot$geneSym[genes_toPlot$group == 'TAM.MLDS.down'][1:30]),
        #features = genes_toPlot$geneSym#[!genes_toPlot$fLiver_geneGroup %in% c('others','-')],
        features = c(gata1s_module$geneSym[gata1s_module$low_FLexpr == T & gata1s_module$isCSM==T & gata1s_module$direction == 'MLDS_up']),
        #features = c(rownames(mtx2),rownames(mtx)[!rownames(mtx) %in% rownames(mtx2)])
        #features = gata1s_module$geneSym[gata1s_module$fLiver_geneGroup == 'MK' & abs(gata1s_module$pct.diff) >= 30][1:50]
        )+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 



plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_0424'



figS3_GATA1s_moduleExpr_inFLiver = function(){
  
  #gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  #gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.allTAM.MLDS_topGeneModule.csv',row.names = 1)
  gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/jul24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  table(gata1s_module$tam_vs_mempT21_group)
  gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
  
  
  ##---- In fLiver
  akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'
  akLiv_mdat_fp = gsub('_0824.RDS','_0824_mdat.csv',akLiv_srat_fp)
  
  fLiver = readRDS(akLiv_srat_fp)
  fLiver_mdat = read.csv(akLiv_mdat_fp,row.names = 1)
  fLiver$broadLineage = fLiver_mdat$broadLineage[match(fLiver$cellID,fLiver_mdat$cellID)]
  fLiver$finalAnn_broad = fLiver$annot
  fLiver$Sex[fLiver$Sex %in% c('F','Female','female','XX')] = 'F'
  fLiver$Sex[fLiver$Sex %in% c('M','Male','male','XY')] = 'M'
  fLiver$sex = fLiver$Sex
  
  
  
  
  fLiver$group_4_avgExpr = as.character(fLiver$finalAnn_broad)
  fLiver$group_4_avgExpr[fLiver$group_4_avgExpr == 'promyelocyte'] = 'myelocyte'
  fLiver$group_4_avgExpr[fLiver$group_4_avgExpr %in% c('proMono','Monocyte','Kupffer.cell','Macrophage')] = 'Mono.Mac'
  fLiver$group_4_avgExpr[fLiver$broadLineage == 'Dendritic cells'] = 'Dendritic.cells'
  fLiver$group_4_avgExpr[fLiver$group_4_avgExpr == 'T.cell'] = 'NK_T'
  fLiver$group_4_avgExpr[fLiver$group_4_avgExpr %in% c('pro.B.cell','pre.B.cell')] = 'early.B.cell'
  
  fLiver$group_4_avgExpr = paste0(fLiver$group_4_avgExpr,'.',fLiver$Genotype)
  avgExpr = AverageExpression(fLiver,group.by = 'group_4_avgExpr')
  avgExpr = avgExpr[['RNA']]
  
  
  ## Subset to genes of interest
  topGenes = gata1s_module[gata1s_module$group %in% c('TAM.MLDS.down','TAM.MLDS.up'),]
  table(topGenes$group)
  
  
  
  ## Make the heatmap
  mtx = avgExpr[rownames(avgExpr) %in% topGenes$geneSym[topGenes$group %in% c('TAM.MLDS.down','TAM.MLDS.up')],grepl('diploid',colnames(avgExpr)) & !grepl('NPC|Endo|Fibroblast|Hepa|Mesen|doublets|Trophoblast|Cholangiocyte|Mesothelial|Neuron',colnames(avgExpr))]
  #mtx = avgExpr[rownames(avgExpr) %in% topGenes$geneSym[topGenes$group %in% c('TAM.MLDS.down','TAM.MLDS.up')],grepl('diploid',colnames(avgExpr)) & !grepl('NPC|doublets',colnames(avgExpr))]
  mtx=mtx[rowSums(mtx) > 0,]
  
  # Remove lowly expressed genes?
  #mtx = mtx[apply(mtx,1,max) > 0.5,]
  
  library(ComplexHeatmap)
  library(circlize)
  col_fun = colorRamp2(c(-3, 0, 3), c('#1a4a87','white','#a4282c'))
  
  avgExpr_cols  = c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
  bCols_avgExpr = circlize::colorRamp2(seq(0.05,5,length.out=length(avgExpr_cols)),avgExpr_cols)
  maxAvgExpr = apply(mtx,1,max)
  maxAvgExpr = as.numeric(maxAvgExpr>0.5)
  rowAnno = rowAnnotation(maxAvgExpr = maxAvgExpr,
                          col=list(maxAvgExpr=c('1'='red','0'='white')))
  hm = Heatmap(t(scale(t(mtx))),name='zScaled Avg. Expression',na_col = 'grey',
               show_row_dend = F,show_column_dend = F,
               show_row_names = T,show_column_names = T,
               cluster_rows = T,cluster_columns = T,
               col = col_fun,
               right_annotation = rowAnno,
               row_names_gp = gpar(fontsize=6),column_names_gp = gpar(fontsize=6),column_title_gp = gpar(fontsize=10),column_title_rot = 90,km=8,
               #split = topGenes$group[match(rownames(mtx),topGenes$geneSym)],
               row_title_rot = 0,
               column_split = fLiver$broadLineage[match(colnames(mtx),fLiver$group_4_avgExpr)])
  
  # hm = Heatmap(t(scale(t(mtx))),name='zScaled Avg. Expression',na_col = 'grey',
  #              show_row_dend = F,show_column_dend = F,
  #              show_row_names = T,show_column_names = T,
  #              cluster_rows = T,cluster_columns = T,
  #              col = col_fun,
  #              right_annotation = rowAnno,
  #              row_names_gp = gpar(fontsize=6),column_names_gp = gpar(fontsize=6),column_title_gp = gpar(fontsize=10),column_title_rot = 90,km=8,
  #              #split = topGenes$group[match(rownames(mtx),topGenes$geneSym)],
  #              column_split = fLiver$broadLineage[match(colnames(mtx),fLiver$group_4_avgExpr)])
  # 
  ht = draw(hm)
  ## Add gene group
  gata1s_module$fLiver_geneGroup = '-'
  
  # gata1s_module$fLiver_geneGroup_withStromal = ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['1']],row_order(ht)[['2']])],'Ery',
  #                                         ifelse(gata1s_module$geneSym %in% rownames(mtx)[row_order(ht)[['3']]],'MK',
  #                                                ifelse(gata1s_module$geneSym %in% rownames(mtx)[row_order(ht)[['4']]],'Prog',
  #                                                       ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[[c('5')]])],'NK.T',
  #                                                              ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['9']],row_order(ht)[['8']],row_order(ht)[['10']])],'Stromal',
  #                                                                     ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['6']])],'Myeloid',
  #                                                                            ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['7']])],'B','others')))))))
  # 
  
  gata1s_module$fLiver_geneGroup = ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['1']])],'Ery',
                                          ifelse(gata1s_module$geneSym %in% rownames(mtx)[row_order(ht)[['4']]],'MK',
                                                 ifelse(gata1s_module$geneSym %in% rownames(mtx)[row_order(ht)[['2']]],'Prog',
                                                        ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[[c('7')]])],'NK.T',
                                                               ifelse(gata1s_module$geneSym %in% rownames(mtx)[row_order(ht)[['3']]],'Mast',
                                                                      ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['5']],row_order(ht)[[c('6')]])],'Myeloid',
                                                                             ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['8']])],'B',
                                                                                    ifelse(!gata1s_module$geneSym %in% rownames(mtx),'-','others'))))))))
  
  table(gata1s_module$fLiver_geneGroup,gata1s_module$group)
  gata1s_module$fLiver_geneGroup = factor(gata1s_module$fLiver_geneGroup,c('Prog','MK','Mast','Ery','Myeloid','NK.T','B','others','-'))
  
  fLiver$broadLineage = factor(fLiver$broadLineage,levels = c( "HSC & prog.","Megakaryocytes","Mast.cell","Erythroblasts","Monocyte/Macrophage", "Dendritic cells","Myelocytes","T/NK lineage","B lineage","Stromal",'others'))
  Heatmap(t((t(mtx[rownames(mtx) %in% gata1s_module$geneSym[gata1s_module$fLiver_geneGroup == 'MK'],]))),name='zScaled Avg. Expression',
          show_row_dend = F,show_column_dend = F,
          show_row_names = T,show_column_names = T,
          cluster_rows = T,cluster_columns = T,
          #col = col_fun,
          row_names_gp = gpar(fontsize=6),column_names_gp = gpar(fontsize=6),column_title_gp = gpar(fontsize=10),column_title_rot = 90,km=6,
          #split = topGenes$group[match(rownames(mtx),topGenes$geneSym)],
          column_split = fLiver$broadLineage[match(colnames(mtx),fLiver$group_4_avgExpr)])
  
  
  
  
  
  ##------ Pre-determine the order of genes to plot -------##
  ## Make the heatmap
  mtx = avgExpr[rownames(avgExpr) %in% topGenes$geneSym[topGenes$group %in% c('TAM.MLDS.down','TAM.MLDS.up')],grepl('diploid',colnames(avgExpr)) & !grepl('NPC|Endo|Fibroblast|Hepa|Mesen|doublets',colnames(avgExpr))]
  mtx=mtx[rowSums(mtx) > 0,]
  # Remove lowly expressed genes?
  #mtx = mtx[apply(mtx,1,max) > 0.5,]
  
  mtx2 = mtx
  colnames(mtx2) = gsub('\\.diploid$','',colnames(mtx2))
  mtx2 = mtx2[,c("HSC_MPP","MEMP_MEP","CMP_GMP","LMPP_ELP","earlyMK","MK","Mast.cell",'EE','ME','LE','Mono.Mac','Dendritic.cells','myelocyte',"ILC.precursor",'NK_T',"early.B.cell",'B.cell')]
  mtx = mtx[,paste0(colnames(mtx2),'.diploid')]
  
  colnames(mtx2) = gsub('_',' / ',colnames(mtx2))
  colnames(mtx2) = gsub('\\.',' ',colnames(mtx2))
  colnames(mtx2) = gsub('cells$','cell',colnames(mtx2))
  colnames(mtx2)[colnames(mtx2) == 'EE'] = 'early Ery'
  colnames(mtx2)[colnames(mtx2) == 'ME'] = 'mid Ery'
  colnames(mtx2)[colnames(mtx2) == 'LE'] = 'late Ery'
  colnames(mtx2)[colnames(mtx2) == 'earlyMK'] = 'early MK'
  colnames(mtx2)[colnames(mtx2) == 'myelocyte'] = 'Myelocyte'
  
  
  
  genes_toMark = gata1s_module$geneSym[gata1s_module$isTF==T | gata1s_module$isCSM==T | gata1s_module$isCosmic==T | gata1s_module$chr == 'chr21']
  rownames(mtx2) = ifelse(rownames(mtx2) %in% genes_toMark,rownames(mtx2),'')
  
  
  gene_order = gata1s_module[gata1s_module$geneSym %in% rownames(mtx),]
  gene_order$direction = factor(gene_order$direction,c('MLDS_up','MLDS_down'))
  gene_order$FL_group = gene_order$fLiver_geneGroup
  gene_order$gene_split_order = paste0(gene_order$group,':',gene_order$FL_group)
  gene_order$gene_split_order = factor(gene_order$gene_split_order,c('TAM.MLDS.up:Prog',
                                                                     'TAM.MLDS.up:MK',
                                                                     'TAM.MLDS.up:Mast',
                                                                     'TAM.MLDS.up:Ery',
                                                                     'TAM.MLDS.up:Myeloid',
                                                                     'TAM.MLDS.up:NK.T',
                                                                     'TAM.MLDS.up:B',
                                                                     
                                                                     'TAM.MLDS.down:Prog',
                                                                     'TAM.MLDS.down:MK',
                                                                     'TAM.MLDS.down:Mast',
                                                                     'TAM.MLDS.down:Ery',
                                                                     'TAM.MLDS.down:Myeloid',
                                                                     'TAM.MLDS.down:NK.T',
                                                                     'TAM.MLDS.down:B'
  ))
  
  
  
  
  ## Heatmap colours
  pR2Cols  = c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
  bCols = circlize::colorRamp2(seq(log1p(0),log1p(3),length.out=length(pR2Cols)),pR2Cols)
  l2FC_cols = colorRamp2(c(-2, -1, 0,1,2), c('#1a4a87','#7592b7','white','#bf686b','#a4282c'))
  
  
  
  data_forPlot = mtx %>% as.data.frame()
  data_forPlot = merge(gene_order,data_forPlot,by.x='geneSym',by.y=0)
  
  data_forPlot = read.delim('~/lustre_mt22/Aneuploidy/manuscriptDraft_0424/Plots/Supp.FigSxx_GATA1s_moduleExpr_in2nFLiver_heatmap_vertical_rawData.tsv',sep = '\t')
  
  
  
  
  plotFun_GATA1s_moduleExpr_inFLiver_vertical = function(noFrame=FALSE,noPlot=FALSE){
    
    rowAnno = rowAnnotation(df=data.frame(MEMP_expr = log1p(mtx[,'MEMP_MEP.diploid']),
                                          avg_l2FC = gene_order$log2FoldChange[match(rownames(mtx),gene_order$geneSym)]),
                            col=list(MEMP_expr=bCols,
                                     avg_l2FC=l2FC_cols
                            ))
    # rowAnno = rowAnnotation(df=data.frame(MEMP_expr = log1p(mtx[,'MEMP_MEP.diploid']),
    #                                       avg_l2FC = gata1s_module$log2FoldChange[match(rownames(mtx2),gata1s_module$geneSym)]),
    #                         col=list(MEMP_expr=bCols,
    #                                  avg_l2FC=l2FC_cols
    #                         ),annotation_legend_param = list(MEMP_expr = list(direction = "horizontal"),
    #                                                          avg_l2FC = list(direction = "horizontal")))
    
    rownames(mtx2) = ifelse(rownames(mtx2) %in% gene_order$geneSym[gene_order$isTF==T | gene_order$chr == 'chr21'],rownames(mtx2),'')
    
    
    
    hm_toPrint = Heatmap(t(scale(t(mtx2))),name='zScaled Avg. Expression',
                         show_row_dend = F,show_column_dend = F,
                         show_row_names = T,show_column_names = T,
                         cluster_rows = T,cluster_columns = F,cluster_column_slices = F,cluster_row_slices = F,
                         col = col_fun,
                         row_names_gp = gpar(fontsize=6,col = ifelse(gene_order$chr[match(rownames(mtx2),gene_order$geneSym)] == 'chr21','purple','black'),
                                             face =  ifelse(gene_order$isTF[match(rownames(mtx2),gene_order$geneSym)] == T,'bold','plain')),
                         column_names_gp = gpar(fontsize=7.5),
                         column_title_gp = gpar(fontsize=10),column_title_rot = 90,
                         row_title_gp = gpar(fontsize=0),
                         left_annotation = rowAnno,
                         row_gap = unit(0.12,'cm'),column_gap = unit(0.12,'cm'),
                         row_names_side = "right", 
                         #right_annotation = ha,
                         #split = gene_order$group[match(rownames(mtx),gene_order$geneSym)],
                         split = gene_order$gene_split_order[match(rownames(mtx),gene_order$geneSym)],
                         column_split = fLiver$broadLineage[match(colnames(mtx),fLiver$group_4_avgExpr)])
    
    
    draw(hm_toPrint)
    
  }
  
  saveFig(file.path(plotDir,paste0('FigSuppXX_GATA1s_moduleExpr_in2nFLiver_heatmap_vertical')),plotFun_GATA1s_moduleExpr_inFLiver_vertical,rawData=data_forPlot,width = 6,height = 11,res = 500)  
  
  
  
  
  
  plotFun_GATA1s_moduleExpr_inFLiver_horizontal = function(noFrame=FALSE,noPlot=FALSE){
    
    botAnno = HeatmapAnnotation(df=data.frame(MEMP_expr = log1p(mtx[,'MEMP_MEP.diploid']),
                                              avg_l2FC = gata1s_module$log2FoldChange[match(rownames(mtx),gata1s_module$geneSym)]),
                                col=list(MEMP_expr=bCols,
                                         avg_l2FC=l2FC_cols),
                                annotation_legend_param = list(MEMP_expr = list(direction = "horizontal"),
                                                               avg_l2FC = list(direction = "horizontal")),annotation_name_side = 'left')
    
    # ha = rowAnnotation(geneInterest = anno_mark(at = which(rownames(mtx2) %in% genes_toMark), 
    #                                    labels = rownames(mtx2)[rownames(mtx2) %in% genes_toMark],which = 'column'))
    
    hm_toPrint_hor = Heatmap((scale(t(mtx2))),name='zScaled Avg. Expression',
                             show_row_dend = F,show_column_dend = F,
                             show_row_names = T,show_column_names = T,
                             cluster_rows = F,cluster_columns = T,cluster_column_slices = F,cluster_row_slices = F,
                             col = col_fun,
                             column_names_gp = gpar(fontsize=6, face = ifelse(gata1s_module$chr[match(rownames(mtx2),gata1s_module$geneSym)] == 'chr21','bold','plain'),
                                                    col =  ifelse(gata1s_module$isTF[match(rownames(mtx2),gata1s_module$geneSym)] == T,'purple',
                                                                     ifelse(gata1s_module$isCosmic[match(rownames(mtx2),gata1s_module$geneSym)] == T,'darkgreen',
                                                                            ifelse(gata1s_module$isCSM[match(rownames(mtx2),gata1s_module$geneSym)] == T,'darkblue','black')))),
                             row_names_gp = gpar(fontsize=6),
                             row_title_gp = gpar(fontsize=10),row_title_rot = 0,
                             column_title_gp = gpar(fontsize=0),column_title_rot = 90,
                             top_annotation = botAnno,
                             row_gap = unit(0.12,'cm'),column_gap = unit(0.15,'cm'),
                             row_names_side = "right", 
                             #bottom_annotation = ha,
                             column_split = gene_order$gene_split_order[match(rownames(mtx),gene_order$geneSym)],
                             
                             split = fLiver$broadLineage[match(colnames(mtx),fLiver$group_4_avgExpr)],
                             heatmap_legend_param = list(legend_direction = "horizontal"))
    ht = draw(hm_toPrint_hor,heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  }
  
  saveFig(file.path(plotDir,paste0('FigSupXX_GATA1s_moduleExpr_in2nFLiver_heatmap_horizontal_allTAM')),plotFun_GATA1s_moduleExpr_inFLiver_horizontal,rawData=data_forPlot,width = 12.4,height = 5,res = 500,useDingbats = F)  
}





##----    Overlap with published results in Marta Byrska-Bishop et al 2015 JCI  -----##
##  Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4362246/
##  Published on my birthday!


bishop_geneSig = read.delim('~/lustre_mt22/Aneuploidy/bishop_etal2015_GATA1s_geneSig.txt',sep = ' ',skip = 1)
bishop_geneSig$direction = c(rep('GATA1s_down',which(bishop_geneSig$Gene == 'JHDM1D')),
                             rep('GATA1s_up',(nrow(bishop_geneSig) - which(bishop_geneSig$Gene == 'JHDM1D'))))

bishop_geneSig$in_GATA1s_module = '?'
bishop_geneSig$in_GATA1s_module[bishop_geneSig$direction == 'GATA1s_down'] = (bishop_geneSig$Gene[bishop_geneSig$direction == 'GATA1s_down'] %in% gata1s_module$geneSym[gata1s_module$group == 'TAM.MLDS.down'])
bishop_geneSig$in_GATA1s_module[bishop_geneSig$direction == 'GATA1s_up'] = (bishop_geneSig$Gene[bishop_geneSig$direction == 'GATA1s_up'] %in% gata1s_module$geneSym[gata1s_module$group == 'TAM.MLDS.up'])
table(bishop_geneSig$in_GATA1s_module)





