# Annotate ML-DS scRNAseq dataset, using:
outDir = '~/ML-DS/Results/03_MLDS_annotation'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}


## Import ML-DS srat object ------------------------------------------------------------------
mlds=readRDS('~/ML-DS/Results/02_MLDS_scQC/MLDS/MLDS_clean_noMTCells.RDS')

## Import existing annotation (August 2024)
mlds_mdat = read.csv('~/ML-DS/Results/MLDS_annotation_250527.csv')
mlds_mdat = mlds_mdat[,!colnames(mlds_mdat) %in% c('X','X.1')]

# Check that cellIDs are the same
table(mlds$cellID %in% mlds_mdat$cellID)
table(mlds_mdat$cellID %in% mlds$cellID)
table(mlds_mdat$annot[!mlds_mdat$cellID %in% mlds$cellID])

table(mlds_mdat$cellID %in% mlds@misc$preQC$cellID)

# Subset from preQC object
mlds_new = subset(mlds@misc$preQC,subset = cellID %in% mlds_mdat$cellID)
colnames(mlds_new@meta.data)[!colnames(mlds_new@meta.data) %in% colnames(mlds_mdat)]
table(mlds_mdat$cellID %in% mlds_new$cellID)
colnames(mlds_mdat)[!colnames(mlds_mdat) %in% colnames(mlds_new@meta.data)]

mlds_new@meta.data = cbind(mlds_new@meta.data,mlds_mdat[match(mlds_new$cellID,mlds_mdat$cellID),!colnames(mlds_mdat) %in% colnames(mlds_new@meta.data)])
mlds_new@reductions$umap@cell.embeddings[,'UMAP_1'] = mlds_mdat$UMAP_1[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new@reductions$umap@cell.embeddings[,'UMAP_2'] = mlds_mdat$UMAP_2[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$reasonForFail = mlds_mdat$reasonForFail[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$Phase = mlds_mdat$Phase[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$S.Score = mlds_mdat$S.Score[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$G2M.Score = mlds_mdat$G2M.Score[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$PASS_doublet = mlds_mdat$PASS_doublet[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$PASS_soupX = mlds_mdat$PASS_soupX[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$PASS_cluster = mlds_mdat$PASS_cluster[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$PASS = mlds_mdat$PASS[match(mlds_new$cellID,mlds_mdat$cellID)]
mlds_new$PASS_withMT = mlds_mdat$PASS_withMT[match(mlds_new$cellID,mlds_mdat$cellID)]



DimPlot(mlds_new,group.by = 'annot',label = T,repel = T,label.box = T)+NoLegend()

table(mlds_new$reasonForFail,mlds_new$annot)
DimPlot(mlds_new,cells.highlight = mlds_new$cellID[mlds_new$reasonForFail != 'Pass'])



saveRDS(mlds_new,'~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns.RDS')



