## Markers of T_CD4 B vs BM
Idents(l076) = paste0(l076$annot,':',l076$tissue)
markers_T = FindMarkers(l076,ident.1 ='T_CD4:BM',ident.2 = 'T_CD4:Blood')
markers_T$geneSym = rownames(markers_T)
markers_T = markers_T[markers_T$p_val_adj < 0.05,]
markers_T = markers_T[order(abs(markers_T$avg_log2FC),decreasing = T),]

DotPlot(l076,idents = c(#'T_CD4:BM','T_CD4:Blood',
  'Tumour:BM','Tumour:Blood'),
  features = markers_T$geneSym[markers_T$avg_log2FC >0][1:50]) + RotatedAxis()



markers_Tumour = FindMarkers(l076,ident.1 ='Tumour:BM',ident.2 = 'Tumour:Blood')
markers_Tumour$geneSym = rownames(markers_Tumour)
markers_Tumour = markers_Tumour[markers_Tumour$p_val_adj < 0.05,]
markers_Tumour = markers_Tumour[order(abs(markers_Tumour$avg_log2FC),decreasing = T),]
markers_Tumour$also_Tcell_markers = ifelse(markers_Tumour$avg_log2FC >0 & markers_Tumour$geneSym %in% markers_T$geneSym[markers_T$avg_log2FC >0],T,
                                           ifelse(markers_Tumour$avg_log2FC < 0 & markers_Tumour$geneSym %in% markers_T$geneSym[markers_T$avg_log2FC < 0],T,F))
DotPlot(l076,idents = c('T_CD4:BM','T_CD4:Blood',
                        'Tumour:BM','Tumour:Blood'),
        features = markers_Tumour$geneSym[markers_Tumour$avg_log2FC >0 & markers_Tumour$also_Tcell_markers == F][1:50]) + RotatedAxis()









## Implement some sort of rolling mean aggregated BAF across all SNPs along chromosomes
binSize = 1e7
tgtChrs = paste0('chr',c(1:22,'X'))
genomicBins = GRanges(rep(tgtChrs,1e9/binSize),
                      IRanges(sapply(seq(1,1e9,binSize),function(i){rep(i,length(tgtChrs))}),
                              sapply(seq(binSize,1e9,binSize),function(i){rep(i,length(tgtChrs))})))
genomicBins$binID = paste0(seqnames(genomicBins),':',seq(1:length(genomicBins)))
# genomicBins$chr = as.character(seqnames(genomicBins))
# genomicBins$pos = start(genomicBins)
# View(as.data.frame(mcols(genomicBins[seqnames(genomicBins) == 'chr21',])))
clusterCnts = clusterCnts.segs
clusterCnts_gr = GRanges(clusterCnts$chr,IRanges(clusterCnts$pos,clusterCnts$pos),
                         clusterID = clusterCnts$cellID,altCount = clusterCnts$altCount,refCount = clusterCnts$refCount,
                         regionID = clusterCnts$regionID,cov=clusterCnts$cov)

clusterCnts_gr$varID = seq(1:length(clusterCnts_gr))
tmp = mergeByOverlaps(clusterCnts_gr,genomicBins)
clusterCnts_gr$binID = tmp$binID[match(clusterCnts_gr$varID,tmp$varID)]
clusterCnts_gr = clusterCnts_gr[clusterCnts_gr$cov != '<5']

gCnts.segs.byChrBins = aggregateByLists(gCnts = clusterCnts_gr, assays = c("altCount", "refCount"), cellList = clusterCnts_gr$clusterID, regionList = clusterCnts_gr$binID)
gCnts.segs.byChrBins$chr = as.character(seqnames(genomicBins))[match(gCnts.segs.byChrBins$regionID,genomicBins$binID)]
gCnts.segs.byChrBins$pos = as.numeric(start(genomicBins))[match(gCnts.segs.byChrBins$regionID,genomicBins$binID)]

gCnts.segs.byChrBins$totCount = gCnts.segs.byChrBins$altCount + gCnts.segs.byChrBins$refCount
gCnts.segs.byChrBins$majorFreq = ifelse(gCnts.segs.byChrBins$altCount > gCnts.segs.byChrBins$refCount,
                                        gCnts.segs.byChrBins$altCount/gCnts.segs.byChrBins$totCount,
                                        gCnts.segs.byChrBins$refCount/gCnts.segs.byChrBins$totCount)
gCnts.segs.byChrBins$altFreq = gCnts.segs.byChrBins$altCount / (gCnts.segs.byChrBins$altCount + gCnts.segs.byChrBins$refCount)

gCnts.segs.byChrBins$cov = ifelse(gCnts.segs.byChrBins$totCount < 5,'<5',
                                  ifelse(gCnts.segs.byChrBins$totCount < 10,'<10',
                                         ifelse(gCnts.segs.byChrBins$totCount < 20,'<20','>=20')))
gCnts.segs.byChrBins$cov = factor(gCnts.segs.byChrBins$cov,c('<5','<10','<20','>=20'))
gCnts.segs.byChrBins$tissue = gsub('^.*:','',gCnts.segs.byChrBins$cellID)
gCnts.segs.byChrBins$celltype = gsub(':.*$','',gCnts.segs.byChrBins$cellID)


ggplot(gCnts.segs.byChrBins[gCnts.segs.byChrBins$totCount > 5 &
                              gCnts.segs.byChrBins$celltype %in%  ct_toPlot &
                              gCnts.segs.byChrBins$chr %in% chr_toPlot,],aes(pos,majorFreq))+
  geom_point(aes(col=cov),size=0.8,alpha=1)+
  scale_color_manual(values = c(grey(0.7),col25[c(3,4)]))+
  geom_hline(yintercept = 0.5,lwd=0.4)+
  facet_grid(celltype + tissue ~ chr,scales = 'free_x')+
  theme_classic()+theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        panel.border = element_rect(fill=F),
                        axis.line = element_blank())+ylim(0,1)


ggplot(clusterCnts.segs[clusterCnts.segs$totCount > 5 &
                          clusterCnts.segs$celltype %in%  ct_toPlot &
                          clusterCnts.segs$chr %in% chr_toPlot,],aes(pos,altFreq))+
  geom_point(aes(col=cov),size=0.2,alpha=1)+
  geom_point(data = gCnts.segs.byChrBins[gCnts.segs.byChrBins$totCount > 5 &
                                           gCnts.segs.byChrBins$celltype %in%  ct_toPlot &
                                           gCnts.segs.byChrBins$chr %in% chr_toPlot,],
             aes(pos,majorFreq,col=cov),size=0.8,alpha=1,col='red')+
  scale_color_manual(values = c(grey(0.7),col25[c(3,4)]))+
  geom_hline(yintercept = 0.5,lwd=0.4)+
  geom_hline(yintercept = c(1/3,2/3),lwd=0.2,col='red')+
  facet_grid(celltype + tissue ~ chr,scales = 'free_x')+
  theme_classic()+theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        panel.border = element_rect(fill=F),
                        axis.line = element_blank())

