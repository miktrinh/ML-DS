##----   Deriving MLDS imprint in T21 foetal MEMP, using good TAM only    -----##


outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/jul24'
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
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')


##----------------------------##
##   Set Global parameters  ####
##----------------------------##
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





#Define genomic coordinates
# gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf'
# txdb = makeTxDbFromGFF(gtf)
# gns = genes(txdb)
# 
# geneMap = read.delim('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data/Donor15680/liver/fLiver_MY_200531_10043298/filtered_feature_bc_matrix/features.tsv.gz',header = F)
# colnames(geneMap) = c('ensID','geneSym','GEX')
# geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
# geneMap$geneSym = gsub('_','-',geneMap$geneSym)
# geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))


##---------------------------------------------------------##
#### 1. Import relevant scRNA datasets: MLDS + fLiver    ####
##---------------------------------------------------------##
## Import MLDS object
#mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_noMTCells.RDS')
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2406_mdat.csv')
mlds$annot_jun24 = mdat$annot_jun24[match(mlds$cellID,mdat$cellID)]
# Add Henning's samples clinical outcome
mlds$clinicalOutcome[mlds$donorID %in% c('CC4','CC5','CC1','CC2')] = 'Remission'
mlds$clinicalOutcome[mlds$donorID %in% c('CC3')] = 'Relapse'
DimPlot(mlds,group.by = 'annot_jun24',cols = c(col25,pal34H),label.size = 3,label = T,repel = T,label.box = T) + NoLegend()


## Import fLiver object
akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0724.RDS'
akLiv_mdat_fp = gsub('_0724.RDS','_0724_mdat.csv',akLiv_srat_fp)

fLiver = readRDS(akLiv_srat_fp)
fLiver_mdat = read.csv(akLiv_mdat_fp,row.names = 1)
fLiver$broadLineage = fLiver_mdat$broadLineage[match(fLiver$cellID,fLiver_mdat$cellID)]
all(fLiver$finalAnn_broad == fLiver$annot_jun24)
if(sum(is.na(fLiver$Sex) | fLiver$Sex == '?') > 0){
  ## Determine sex
  avg.expr.sexMarkers = AverageExpression(fLiver,group.by = 'donorID',features = c('XIST','RPS4Y1'))
  avg.expr.sexMarkers = as.data.frame(t(avg.expr.sexMarkers[[1]]))
  avg.expr.sexMarkers$sex = fLiver@meta.data$sex[match(rownames(avg.expr.sexMarkers),fLiver@meta.data$donorID)]
  avg.expr.sexMarkers$pred.sex = ifelse(avg.expr.sexMarkers$XIST > avg.expr.sexMarkers$RPS4Y1,'F',
                                        ifelse(avg.expr.sexMarkers$XIST < avg.expr.sexMarkers$RPS4Y1,'M','??'))
  
  avg.expr.sexMarkers[avg.expr.sexMarkers$sex != avg.expr.sexMarkers$pred.sex,]
  
  fLiver$Sex[fLiver$donorID %in% c('Hsb36','Hsb37')] = 'XX'
  
}


fLiver$Sex[fLiver$Sex %in% c('F','Female','female','XX')] = 'F'
fLiver$Sex[fLiver$Sex %in% c('M','Male','male','XY')] = 'M'
fLiver$sex = fLiver$Sex



##-------------------------------------------------##
##   2. DEGs between T21_fLiver vs 2n_fLiver     ####
##-------------------------------------------------##
## NOTE: formula for comparison
##       1. I have compared the results from this "compareCell_simplified" function to the normal "compareCell" function when perforing perCT_perGeno DESeq2 comparison - and the results are similar to one another
##       2. The question is - should I add sex as a co-variate in the model?
# if yes --> results give more DEGs (~ 40 AK_up and ~ 4 AK_down, in addition to 42/50 shared DEGs). 
# on one hand, I do not trust this very much given that we only have 1/3 diploid sample being male. although we do have 3 male samples in T21 cohort, 2 of these are the lower quality samples when non-sorted and have fewer cells...
# additionally, add sex now makes the results in consistent with that done with the fLiver analysis.
#   --> on this point, I think we can say that to make it more comparable with other genotypes, given that the other genotypes don't have enough sex representation to do this? 
#       that being said - for fLiver analysis, should I do 1 cell type, all Geno at once, and specify the controlGenes to be non-trisomy chromosomes? 
#       --> this does not work well for Triploid case --> easier to just do them seperately
#       Can justify including sex here to make the comparison more specific/sensitive?
# 
# OR we can just not add the sex here?

##------- Using pseudobulk -----------
cnt_mtx = fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$cellID %in% c(fLiver$cellID[fLiver$Genotype %in% c('diploid','T21') & fLiver$annot_jun24 == 'MEMP_MEP'])]]
rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = fLiver@meta.data[fLiver$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot_jun24','Sex','Genotype','assay')]
mDat$sex = mDat$Sex
mDat$geno = mDat$Genotype

table(mDat$sex,mDat$donorID,mDat$Genotype)


# t21_vs_2n_res = compareCell_simplified(toc = cnt_mtx,
#                                        mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
#                                        coords = gns[rownames(cnt_mtx)],
#                                        cellTypeName='MEMP_MEP',
#                                        formula='~ %s + sex',tgtChrs = paste0(c(1:22)),
#                                        donorID='donorID',groupID='geno')

## To keep it consistent with the other AK analysis - use compareCell function instead. The only real difference is the use of "controlGenes"
## this requires mDat to have the "group" column
mDat$group = mDat$Genotype
controlGenes = gns[!seqnames(gns) %in% c('X','Y','21')]$gene_id

t21_vs_2n_res = compareCell(toc = cnt_mtx,
                            mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                            coords = gns[rownames(cnt_mtx)],
                            cellTypeName='MEMP_MEP',
                            formula='~ %s + assay',
                            tgtChrs = paste0('chr',c(1:22)),
                            donorID='donorID',groupID='geno',
                            geneMap = geneMap,
                            controlGenes = controlGenes)

saveRDS(t21_vs_2n_res,'DESeq2_T21.vs.diploid.fLiver.MEMP_geno.assay_res.RDS')
t21_vs_2n_res = readRDS('DESeq2_T21.vs.diploid.fLiver.MEMP_geno.assay_res.RDS')
t21_vs_2n_results = t21_vs_2n_res[['res']] %>% as.data.frame()
t21_vs_2n_results = annotateGenes(t21_vs_2n_results,geneMap = geneMap)

# mDat$group = mDat$Genotype

# t21_vs_2n_res4 = compareCell(toc = cnt_mtx,
#                              mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
#                              coords = gns[rownames(cnt_mtx)],
#                              cellTypeName='MEMP_MEP',
#                              formula='~ %s',
#                              tgtChrs = paste0(c(1:22)),
#                              donorID='donorID',groupID='geno',
#                              geneMap = geneMap,
#                              controlGenes = controlGenes)
# 
# 
# t21_vs_2n_deg = t21_vs_2n_res[['mainDE']]
# t21_vs_2n_deg$direction = ifelse(t21_vs_2n_deg$log2FoldChange > 0 ,'AK_up','AK_down')
# t21_vs_2n_deg = t21_vs_2n_deg2
# t21_vs_2n_deg = t21_vs_2n_deg[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',t21_vs_2n_deg$geneSym),]
# t21_vs_2n_deg$max_pct = pmax(t21_vs_2n_deg$cellFrac_g1,t21_vs_2n_deg$cellFrac_g2)
# t21_vs_2n_deg = t21_vs_2n_deg[t21_vs_2n_deg$max_pct > 10,]
# ## compared to AK fLiver DESeq results
# allDEGs.sub.filtered = read.csv('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/apr24/without_published2n/geno/liver_allDEGs_filtered.csv')
# mempT21 = allDEGs.sub.filtered[allDEGs.sub.filtered$geno == 'T21' & allDEGs.sub.filtered$ct == 'MEMP_MEP',]
# t21_vs_2n_deg$inBoth = ifelse(t21_vs_2n_deg$log2FoldChange > 0 & t21_vs_2n_deg$geneSym %in% mempT21$geneSym[mempT21$log2FoldChange >0],T,
#                               ifelse(t21_vs_2n_deg$log2FoldChange < 0 & t21_vs_2n_deg$geneSym %in% mempT21$geneSym[mempT21$log2FoldChange < 0],T,F))
# t21_vs_2n_deg$in_noSex = ifelse(t21_vs_2n_deg$log2FoldChange > 0 & t21_vs_2n_deg$geneSym %in% t21_vs_2n_deg2$geneSym[t21_vs_2n_deg2$log2FoldChange >0],T,
#                               ifelse(t21_vs_2n_deg$log2FoldChange < 0 & t21_vs_2n_deg$geneSym %in% t21_vs_2n_deg2$geneSym[t21_vs_2n_deg2$log2FoldChange < 0],T,F))
# table(t21_vs_2n_deg$inBoth,t21_vs_2n_deg$direction)
# 
# t21_vs_2n_deg = t21_vs_2n_deg[order(abs(t21_vs_2n_deg$log2FoldChange),decreasing = T),]
# s$tmp = paste0(s$Genotype,':',s$Sex)
# DotPlot(s,group.by = 'tmp',
#         cols = c(colAlpha(grey(0.95),0.8),'black'),
#         # features = c(mempDEGs$geneSym[mempDEGs$geno == 'T21' & mempDEGs$direction == 'AK_up'][1:30],
#         #              mempDEGs$geneSym[mempDEGs$geno == 'T21' & mempDEGs$direction == 'AK_down'])
#         features = c(t21_vs_2n_deg$geneSym[t21_vs_2n_deg$direction == 'AK_up' & t21_vs_2n_deg$inBoth==F & t21_vs_2n_deg$in_noSex ==T & t21_vs_2n_deg$cellFrac_overall >= 6]
#                      #t21_vs_2n_deg$geneSym[t21_vs_2n_deg$direction == 'AK_down' & t21_vs_2n_deg$inBoth==F]
#                      )
# ) + RotatedAxis()+
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7),
#         axis.text.y = element_text(size=11),
#         legend.position = 'top',legend.text = element_text(size=9),legend.title = element_text(size=8.5),legend.key.size = unit(0.6,'cm')) + xlab('') + ylab('') 
# 
# out = readRDS('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/apr24/without_published2n/geno/liver/T21_liver_pseudoBulkDESeq2out.RDS')
# out = out[['MEMP_MEP_unmerged']]
# out = out[['res']]
# out = annotateGenes(out,geneMap = geneMap)
# out = as.data.frame(out)
# View(out[out$geneSym %in% t21_vs_2n_deg$geneSym[t21_vs_2n_deg$direction == 'AK_up' & t21_vs_2n_deg$inBoth==F & t21_vs_2n_deg$cellFrac_overall >= 6],])
# 
# res3 = as.data.frame(t21_vs_2n_res3[['res']])
# res3 = annotateGenes(res3,geneMap = geneMap)
# View(res3[res3$geneSym %in% t21_vs_2n_deg$geneSym[t21_vs_2n_deg$direction == 'AK_up' & t21_vs_2n_deg$inBoth==F & t21_vs_2n_deg$cellFrac_overall >= 6],])
# res3.sex = results(t21_vs_2n_res3[['dds']],name = 'sex_XY_vs_XX')
# res3.sex = annotateGenes(as.data.frame(res3.sex),geneMap = geneMap)




##------- Using pseudobulk to calculate assay differences  -----------
cnt_mtx = fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$cellID %in% c(fLiver$cellID[fLiver$Genotype %in% c('diploid','T21') & fLiver$annot_jun24 == 'MEMP_MEP'])]]
rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = fLiver@meta.data[fLiver$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot_jun24','Sex','Genotype','assay')]
mDat$sex = mDat$Sex
mDat$geno = mDat$Genotype

table(mDat$sex,mDat$donorID,mDat$Genotype)


# t21_vs_2n_res = compareCell_simplified(toc = cnt_mtx,
#                                        mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
#                                        coords = gns[rownames(cnt_mtx)],
#                                        cellTypeName='MEMP_MEP',
#                                        formula='~ %s + sex',tgtChrs = paste0(c(1:22)),
#                                        donorID='donorID',groupID='geno')

## To keep it consistent with the other AK analysis - use compareCell function instead. The only real difference is the use of "controlGenes"
## this requires mDat to have the "group" column
controlGenes = gns[!seqnames(gns) %in% c('X','Y','21')]$gene_id

t21_vs_2n_res = compareCell(toc = cnt_mtx,
                            mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                            coords = gns[rownames(cnt_mtx)],
                            cellTypeName='MEMP_MEP',
                            formula='~ geno + assay',
                            tgtChrs = paste0('chr',c(1:22)),
                            donorID='donorID',groupID='assay',
                            geneMap = geneMap,
                            controlGenes = controlGenes)



##----------------------------------------------------##
##   3. DEGs between MLDS vs diploid 2n fLiver      ####
##----------------------------------------------------##
## Cluster MLDS / TAM blasts + MEP + fLiver cells
fLiver = subset(fLiver,subset = cellID %in% fLiver$cellID[fLiver$annot_jun24 %in% c('MEMP_MEP','HSC_MPP','Mast.cell','MK') &
                                                            fLiver$Genotype %in% c('diploid','T21')])

mlds.sub = subset(mlds,subset = cellID %in% mlds$cellID[mlds$annot_jun24 %in% c('HSC_MPP','EE','MEP','MK','Tumour','unsure_Tum_MK?')])
srat = merge_seurat_objects(fLiver,mlds.sub,keepAllGenes = F,genomeVersions = c('v38','v38'))
srat$tissue[srat$cellID %in% fLiver$cellID] = 'fLiver'
srat = standard_clustering(srat,runHarmony = T,harmonyVar = 'tissue')
DimPlot(srat,group.by = 'annot_jun24',cols = col25,label = T,repel = T,label.size = 3)
FeaturePlot(srat,'GATA1')
Idents(srat) = paste0(srat$annot_jun24,'_',srat$Genotype,'_',srat$tissue,'_',srat$donorID)
p = DotPlot(srat,
            # idents = unique(Idents(srat)[grepl('MEMP_MEP_T21_fLiver|MEMP_MEP_diploid_fLiver|Tumour_T21_Blood|Tumour_T21_BM',Idents(srat)) &
            #                                     !grepl('L041|CC3',Idents(srat))]),
            idents = unique(Idents(srat)[grepl('MEMP_MEP_diploid_fLiver|Tumour_T21_BM',Idents(srat)) &
                                           !grepl('L041|CC3',Idents(srat))]),
            group.by = 'donorID',
            #idents = c('MEMP_MEP_T21_fLiver','MEMP_MEP_diploid_fLiver','Tumour_T21_Blood','Tumour_T21_BM'),
            #idents = unique(mlds$broadLineage[mlds$broadLineage != 'others']),
            #idents = c('Tumour','MEP','HSC_MPP','MK'),
            #cols = c("#EBFFE5", "#244b05"),
            #cols = c(colAlpha('#F1F5FA',1),'#425580'),
            #cols = c(colAlpha(grey(0.95),0.8),'black'),
            #cols = c(grey(0.99), grey(0.2)),
            #group.by = 'seurat_clusters',
            #idents = unique(srat$finalAnn[srat$finalAnn != 'others']),
            #features = genes
            #features = unique(c('RUNX1','MIR99AHG','SON',
            #                    'APP','CYYR1',mlds_vs_2n_deg$geneSym[mlds_vs_2n_deg$direction == 'MLDS_up' & mlds_vs_2n_deg$chr == 'chr21']))
            # features = geneMap$geneSym[geneMap$chr == 'chr21' & geneMap$geneSym %in% de$geneSym[abs(de$log2FoldChange) >= 1] &
            #                              !geneMap$geneSym %in% de_assay_fLiver2n$geneSym][1:100]
            #features = geneMap$geneSym[geneMap$chr == 'chr21' & geneMap$geneSym %in% de_assay_fLiver2nT21$geneSym[!de_assay_fLiver2nT21$geneSym %in% c(de_assay_fLiver2n$geneSym,de_assay_fLiverT21$geneSym)]]
            #features = mlds_vs_2n_deg$geneSym[mlds_vs_2n_deg$chr=='chr21']
            features = geneMap$geneSym[geneMap$chr == 'chr21' & !geneMap$geneSym %in% mlds_vs_2n_deg$geneSym]
            #features = c()
)+RotatedAxis() +
  theme(axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'top') + xlab('') + ylab('')
p

DimPlot(srat,cells.highlight = srat$cellID[srat$annot_jun24 == 'MEMP_MEP'])

##------- Using pseudobulk -----------
cnt_mtx.FL = fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$cellID %in% c(fLiver$cellID[fLiver$Genotype %in% c('diploid') & fLiver$annot_jun24 == 'MEMP_MEP'])]]
cnt_mtx.mlds = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in%
                                                     c(mlds$cellID[mlds$disease == 'MLDS' & mlds$annot_jun24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$tissue == 'BM' & !mlds$donorID %in% c('L041','CC3')])]]
genes = intersect(rownames(cnt_mtx.mlds),rownames(cnt_mtx.FL))
cnt_mtx=cbind(cnt_mtx.FL[genes,],cnt_mtx.mlds[genes,])

rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat.mlds = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot_jun24','sex','assay')]
colnames(mDat.mlds)[colnames(mDat.mlds) == 'annot_jun24'] = 'annot'
table(is.na(mDat.mlds$sex))
mDat.mlds$sex[mDat.mlds$donorID == 'L178'] = 'F'
mDat.fl = fLiver@meta.data[fLiver$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot_jun24','Sex','assay')]
colnames(mDat.fl)[colnames(mDat.fl) == 'Sex'] = 'sex'
colnames(mDat.fl)[colnames(mDat.fl) == 'annot_jun24'] = 'annot'
mDat = rbind(mDat.mlds,mDat.fl)
mDat$dataset = ifelse(mDat$cellID %in% mDat.mlds$cellID,'MLDS','fLiver')
mDat$sex[mDat$sex == 'XX'] = 'F'
mDat$sex[mDat$sex == 'XY'] = 'M'

# mlds_vs_2n_res = compareCell_simplified(toc = cnt_mtx,
#                                          mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
#                                          coords = gns[rownames(cnt_mtx)],
#                                          cellTypeName='MLDS',
#                                          formula='~ %s + sex',tgtChrs = paste0(c(1:22)),
#                                          donorID='donorID',groupID='dataset')

mDat$group = mDat$dataset
controlGenes = gns[!seqnames(gns) %in% c('X','Y','21')]$gene_id

mlds_vs_2n_res = compareCell(toc = cnt_mtx,
                             mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                             coords = gns[rownames(cnt_mtx)],
                             cellTypeName='MLDS',
                             formula='~ %s',
                             tgtChrs = paste0('chr',c(1:22)),
                             donorID='donorID',groupID='dataset',
                             geneMap = geneMap,
                             controlGenes = controlGenes)


saveRDS(mlds_vs_2n_res,'DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_res.RDS')
mlds_vs_2n_res = readRDS('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_res.RDS')

mlds_vs_2n_dds = mlds_vs_2n_res[['dds']]
mlds_vs_2n_results = mlds_vs_2n_res[['res']] %>% as.data.frame()
mlds_vs_2n_results = annotateGenes(mlds_vs_2n_results,geneMap = geneMap)

mlds_vs_2n_deg = mlds_vs_2n_res[['de']]
mlds_vs_2n_deg$direction = ifelse(mlds_vs_2n_deg$log2FoldChange > 0,'MLDS_up','MLDS_down')
mlds_vs_2n_deg$pct.diff = mlds_vs_2n_deg$cellFrac_g1 - mlds_vs_2n_deg$cellFrac_g2
mlds_vs_2n_deg$cellFrac_max = pmax(mlds_vs_2n_deg$cellFrac_g1, mlds_vs_2n_deg$cellFrac_g2)

max_cellfrac = 10
min_l2FC = 0.5


#mlds_vs_2n_deg = mlds_vs_2n_deg[abs(mlds_vs_2n_deg$log2FoldChange) >= min_l2FC & abs(mlds_vs_2n_deg$pct.diff) >= 10,]
mlds_vs_2n_deg = mlds_vs_2n_deg[abs(mlds_vs_2n_deg$log2FoldChange) >= min_l2FC & abs(mlds_vs_2n_deg$cellFrac_max) >= max_cellfrac,]
#mlds_vs_2n_deg = mlds_vs_2n_deg[abs(mlds_vs_2n_deg$log2FoldChange) >= 1 & abs(mlds_vs_2n_deg$cellFrac_max) >= 20,]

table(mlds_vs_2n_deg$direction)
dim(mlds_vs_2n_deg)
mlds_vs_2n_deg = mlds_vs_2n_deg[order(abs(mlds_vs_2n_deg$pct.diff),decreasing = T),]

## Local vs Global impact
table(mlds_vs_2n_deg$direction,mlds_vs_2n_deg$chr)
write.csv(mlds_vs_2n_deg,'DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_topDEGs.csv')
mlds_vs_2n_deg = read.csv('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_topDEGs.csv')















##--------------------------------------------------------------------------##
##   3.2 DEGs between MLDS vs diploid T21 / 2n fLiver in 1 comparison     ####
##--------------------------------------------------------------------------##
##------- Using pseudobulk -----------
cnt_mtx.FL = fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$cellID %in% c(fLiver$cellID[fLiver$Genotype %in% c('diploid','T21') & fLiver$annot_jun24 == 'MEMP_MEP'])]]
cnt_mtx.mlds = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in%
                                                     c(mlds$cellID[mlds$disease == 'MLDS' & mlds$annot_jun24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$tissue == 'BM' & !mlds$donorID %in% c('L041','CC3')])]]
genes = intersect(rownames(cnt_mtx.mlds),rownames(cnt_mtx.FL))
cnt_mtx=cbind(cnt_mtx.FL[genes,],cnt_mtx.mlds[genes,])

rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat.mlds = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot_jun24','sex','assay','Genotype')]
colnames(mDat.mlds)[colnames(mDat.mlds) == 'annot_jun24'] = 'annot'
table(is.na(mDat.mlds$sex))
mDat.mlds$sex[mDat.mlds$donorID == 'L178'] = 'F'
mDat.fl = fLiver@meta.data[fLiver$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot_jun24','Sex','assay','Genotype')]
colnames(mDat.fl)[colnames(mDat.fl) == 'Sex'] = 'sex'
colnames(mDat.fl)[colnames(mDat.fl) == 'annot_jun24'] = 'annot'
mDat = rbind(mDat.mlds,mDat.fl)
mDat$dataset = ifelse(mDat$cellID %in% mDat.mlds$cellID,'MLDS','fLiver')
mDat$sex[mDat$sex == 'XX'] = 'F'
mDat$sex[mDat$sex == 'XY'] = 'M'
mDat$group = ifelse(mDat$dataset == 'fLiver',paste0(mDat$dataset,'_',mDat$Genotype),mDat$dataset)

controlGenes = gns[!seqnames(gns) %in% c('X','Y','21')]$gene_id

mlds_vs_2n_res = compareCell(toc = cnt_mtx,
                             mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                             coords = gns[rownames(cnt_mtx)],
                             cellTypeName='MLDS',
                             formula='~ %s + assay',
                             tgtChrs = paste0('chr',c(1:22)),
                             donorID='donorID',groupID='group',
                             geneMap = geneMap,
                             controlGenes = controlGenes)


mDat$cellID = rownames(mDat)
if(!all(rownames(mDat)==colnames(toc)))
  stop("Rownames of mDat must equal colnames of toc")
if(!all(names(coords)==rownames(toc)))
  stop("Names of coords must equal rownames of toc")
if(!donorID %in% colnames(mDat))
  stop("Donor ID not found in mDat")
if(!groupID %in% colnames(mDat))
  stop("Group ID not found in mDat")
if(!is.factor(mDat[,groupID])){
  warning("groupID column in mDat not a factor, contrast will be based on alphabetical order")
  mDat[,groupID] = factor(mDat[,groupID])
}
#Formula should contain one and only one instance of %s
if(!grepl('%s',formula,fixed=TRUE) || grepl('%s',sub('%s','',formula,fixed=TRUE),fixed=TRUE))
  stop("Invalid formula, must contain one instance of %s")
#Convert it now in case it is still invalid in more subtle ways
formula = as.formula(sprintf(formula,groupID))
#Drop bullshit
#Uninteresting genes
w = which(as.character(seqnames(coords)) %in% tgtChrs)
coords = coords[w]
seqlevels(coords) = tgtChrs
toc = toc[w,]
#Non-expressed genes
w = which(rowSums(toc)>0)
coords = coords[w]
toc = toc[w,]
##########################
# Get genomic coordinates
#Order aa and bb by genomic coords
#If it's positive stranded or unknown use start, else end
coords$TSS = ifelse(strand(coords)=='-',end(coords),start(coords))
o = order(seqnames(coords),coords$TSS)
coords = coords[o]
toc = toc[o,]
#Genomic position linearised
xpos = coords$TSS
chrLens = sapply(split(xpos,as.character(seqnames(coords))),max)
chrLens = chrLens[seqlevels(coords)]
offset = c(0,cumsum(as.numeric(chrLens))[-length(chrLens)])
names(offset) = names(chrLens)
xpos = xpos + offset[as.character(seqnames(coords))]+1e6
mergeMap=NULL

#####################
# Create pseudobulk 
pb = do.call(cbind,lapply(split(colnames(toc),mDat[,donorID]),function(e) rowSums(toc[,e,drop=FALSE])))
colDat = mDat[match(colnames(pb),mDat[,donorID]),]
rownames(colDat) = colDat[,donorID]
dds = DESeqDataSetFromMatrix(round(pb),
                             colData=colDat,
                             design = formula)
#####################
# General comparison
rld = rlog(dds,blind=TRUE)
rldMat = assay(rld)
rldCor = cor(rldMat)

if(!is.null(controlGenes)){
  ctrlGenes_id = which(rownames(dds) %in% controlGenes)  
  dds <- estimateSizeFactors(dds,controlGenes=ctrlGenes_id)
}

dds = DESeq(dds)
if(doPlot){
  #PCA
  pca_dat = plotPCA_inhouse(object = rld,intgroup=groupID,ntop = 5000,returnData = TRUE)  
  plot(pca_dat[[1]] + ggtitle(cellTypeName))
  
  pca_data = as.data.frame(pca_dat[[2]])
  percent_var = pca_dat[[3]]
  for(col_factor in colnames(colDat)){
    if(col_factor == 'nCells'){
      next
    }
    p = ggplot(pca_data,aes(x=PC1,y=PC2,col=.data[[col_factor]]))+
      geom_point() + 
      scale_color_manual(values = c(brewer.pal(8,'Dark2'),brewer.pal(12,'Paired')[-11],col25,col22,col25,col22))+
      theme_bw() + 
      xlab(paste0('PC1: ',round(percent_var[1],digits = 3)*100,'% variance'))+
      ylab(paste0('PC2: ',round(percent_var[2],digits = 3)*100,'% variance'))+
      labs(colour=col_factor)+
      ggtitle(cellTypeName)
    
    print(p)
  }
  
  
  #Heatmap
  hm = Heatmap(rldCor,
               name='Correlation',
               top_annotation = columnAnnotation(group=colDat[,groupID]),
               right_annotation = rowAnnotation(group=colDat[,groupID])
  )
  draw(hm)
  plotDispEsts(dds)
  
  # Plot distribution of estimated dispersion 
  est_disp = data.frame(disp = dispersions(dds))
  est_disp$ensID = rownames(dds)
  est_disp$chr = geneMap$chr[match(est_disp$ensID,geneMap$ensID)]
  p1 = ggplot(est_disp,aes(log10(disp)))+
    geom_density() + ggtitle('Distribution of estimated dispersion') +
    facet_wrap(vars(chr))+
    theme_classic() + theme(panel.border = element_rect(fill=F),axis.line = element_blank())
  print(p1)
  
}

##########################
# Differential expression
#Syntax is the factor to look at, numerator of logFC, denominator of logFC
resultsNames(dds)
contrast1 = sprintf('%s_%s_vs_%s',groupID,levels(colDat[,groupID])[2],levels(colDat[,groupID])[1])
contrast2 = sprintf('%s_%s_vs_%s',groupID,levels(colDat[,groupID])[3],levels(colDat[,groupID])[1])
contrast3 = 'assay_scRNA.3..v3.1_vs_GEX5p'
de_assay = de
de_assay_fLiver2n = de
de_assay_fLiverT21 = de
de_assay_fLiver2nT21 = de
table(de$chr)
table(de_assay$chr)

lfcCut=0
pCut=0.05
kMerge = NULL
tfMap=NULL
membrane_protein=NULL
cosmicGenes=NULL
tsg=NULL
doPlot=T
for(contrast in c(contrast1,contrast2)){
  groupID2 = gsub('group_|_vs_fLiver_diploid','',contrast)
  res = results(dds,name = contrast, alpha = pCut,lfcThreshold=lfcCut)  
  #Shrink the LFC estimates
  res = tryCatch({lfcShrink(dds,coef = contrast,res=res)
  }, error = function(e){
    print(e)
    message('\nTrying ashr instead...')
    res = lfcShrink(dds,coef = contrast,res=res,type = 'ashr')
    return(res)
  }, finally = {message('Successfully performed lfcShrinkage!')})
  
  #res = lfcShrink(dds,coef = contrast,res=res)
  res = res[order(res$pvalue),]
  de = as.data.frame(res[which(res$padj<pCut),])
  if(nrow(de) == 0){
    de = NULL
  }else{
    de$contrast = contrast
    if(is.null(kMerge)){
      de = annotateGenes(de,geneMap=geneMap,tfMap=tfMap,membrane_protein=membrane_protein,cosmicGenes=cosmicGenes,tsg=tsg)  
      
      # Add cell fraction
      if(nrow(de) >1){
        de$cellFrac_overall = apply(toc[match(de$ensID,rownames(toc)),],1,function(x){sum(x>0.55)})  
        de$cellFrac_g1 = apply(toc[match(de$ensID,rownames(toc)),colnames(toc) %in% mDat$cellID[mDat[,groupID] %in% levels(colDat[,groupID])[1]]],1,function(x){sum(x>0.55)})
        de$cellFrac_g2 = apply(toc[match(de$ensID,rownames(toc)),colnames(toc) %in% mDat$cellID[mDat[,groupID] %in% groupID2]],1,function(x){sum(x>0.55)})
      }else{
        de$cellFrac_overall = sum(toc[match(de$ensID,rownames(toc)),] > 0.55)
        de$cellFrac_g1 = sum(toc[match(de$ensID,rownames(toc)),colnames(toc) %in% mDat$cellID[mDat[,groupID] %in% levels(colDat[,groupID])[1]]] > 0.55)
        de$cellFrac_g2 = sum(toc[match(de$ensID,rownames(toc)),colnames(toc) %in% mDat$cellID[mDat[,groupID] %in% groupID2]] > 0.55)
      }
      
      de$cellFrac_overall = 100*de$cellFrac_overall / ncol(toc)
      
      de$cellFrac_g1 = 100*de$cellFrac_g1 / ncol(toc[,colnames(toc) %in% mDat$cellID[mDat[,groupID] %in% levels(colDat[,groupID])[1]]])
      
      de$cellFrac_g2 = 100*de$cellFrac_g2 / ncol(toc[,colnames(toc) %in% mDat$cellID[mDat[,groupID] %in% groupID2]])
      
    }
  }
  
  
  normCnts = counts(dds,normalized=TRUE)
  if(doPlot & !is.null(de)){
    #Scatterplot of top 20
    d = melt(normCnts,varnames=c('genes','donors'))
    d$group = colDat$group[match(d$donors,colDat$donor)]
    gg = ggplot(d[d$genes %in% head(rownames(de),n=20),],aes(genes,value,colour=group)) +
      geom_point(position=position_jitter(w=0.1,h=0)) +
      scale_y_log10() + 
      xlab('Genes') +
      ylab('Normalised counts') +
      ggtitle(cellTypeName,subtitle = as.character(formula))
    plot(gg)
    #Heatmap, with z-transform on rows
    if(length(rownames(de))>1){
      hm = Heatmap(t(scale(t(normCnts[rownames(de),]))),
                   name='z-scaled expression',
                   top_annotation = columnAnnotation(group=colDat[,groupID]),
                   show_row_names=FALSE)
      draw(hm)
    }else{
      message(sprintf('only %d DE genes detected for %s with %s',nrow(de),cellTypeName,unique(colDat[,groupID])[unique(colDat[,groupID]) != 'diploid']))
    }
    
    #Volcano!
    plot(res$log2FoldChange,-log10(res$padj),
         col = ifelse(res$padj<pCut & abs(res$log2FoldChange)>lfcCut,'red','black'),
         cex=0.5,
         xlab='logFC',
         ylab='-log10(adjusted p-value)',
         main=cellTypeName
    )
  }
  
  ######################
  # Coordinate analysis
  #Visualise which chromosome things are on
  x = table(factor(de$log2FoldChange>0,levels=c(TRUE,FALSE)),
            factor(as.character(seqnames(coords[rownames(de)])),levels=tgtChrs)
  )
  class(x) = 'matrix'
  #Normalise by number of expressed genes
  nGenes = table(as.character(seqnames(coords)))
  fracDE = t(t(x)/as.numeric(nGenes[colnames(x)]))
  dd = data.frame(chr = tgtChrs,
                  nExpressed = as.numeric(nGenes[tgtChrs]),
                  nDE = colSums(x)[tgtChrs],
                  nUp = x[1,tgtChrs],
                  nDown = x[2,tgtChrs])
  dd$fracDE = dd$nDE/dd$nExpressed
  dd$fracUp = dd$nUp/dd$nDE
  dd$fracDown = dd$nDown/dd$nDE
  deTab = dd
  if(doPlot){
    #Fraction DE (and in which direction) by chromosome
    plot(seq_along(tgtChrs),rep(0,length(tgtChrs)),
         type='n',
         ylim=c(-max(fracDE,na.rm = T),max(fracDE,na.rm = T)),
         xaxt='n',
         bty='n',
         ylab='Frac genes DE',
         xlab='Chromosome',
         main=cellTypeName,sub = as.character(formula)
    )
    axis(1,at=seq_along(tgtChrs),labels=tgtChrs)
    abline(h=0,col='black')
    for(i in seq_along(tgtChrs)){
      if(!tgtChrs[i] %in% colnames(fracDE))
        next
      m = match(tgtChrs[i],colnames(fracDE))
      lines(x=rep(i,2),
            y=c(0,fracDE[1,m]),
            col='red',
            lwd=2
      )
      lines(x=rep(i,2),
            y=c(0,-fracDE[2,m]),
            col='blue',
            lwd=2
      )
    }
    y = res[rownames(toc),'log2FoldChange']
    #Boxplots of logFC
    #We will discard the outliers as they are the DE genes and already covered by the above
    log2FC_perChr = split(y,as.character(seqnames(coords)))[tgtChrs]
    boxplot(log2FC_perChr,
            outline=FALSE,
            xlab='Chromosome',
            ylab='logFC'
    )
    abline(h=0,col='red')
    #As close to raw as we can get
    yMax=1
    plot(xpos,y,
         type='n',
         xaxt='n',
         bty='n',
         xlab='Genomic Position',
         ylab='logFC',
         #ylim=c(-yMax,yMax),
         main=cellTypeName)
    abline(h=0,col='red')
    a = sapply(split(xpos,as.character(seqnames(coords))),max)[tgtChrs]
    b = sapply(split(xpos,as.character(seqnames(coords))),min)[tgtChrs]
    bounds = a + c((b[-1]-a[-length(a)])/2,0)
    abline(v=c(0,bounds))
    axis(1,at=sapply(split(xpos,as.character(seqnames(coords))),mean)[tgtChrs],labels=tgtChrs)
    points(xpos,y,
           col=ifelse(rownames(toc) %in% rownames(de),'red','black'),
           pch=19,
           cex=0.1)
    #arrows(x0=xpos,
    #       x1=xpos, 
    #       y0=y-res[rownames(toc),'lfcSE'], 
    #       y1=y+res[rownames(toc),'lfcSE'], 
    #       col = ifelse(rownames(toc) %in% rownames(de),'red','black'),
    #       lwd=0.1,
    #       code=3, 
    #       angle=90, 
    #       length=0.1)
    #Rolling mean, but confined to chromosomes
    for(tgtChr in tgtChrs){
      w = which(as.character(seqnames(coords))==tgtChr)
      if(length(w)<=kSmooth){
        yy = rep(mean(y[w]),length(w))
      }else{
        yy = rollmean(y[w],k=kSmooth,fill='extend')
      }
      lines(xpos[w],yy,col='orange')
    }
  }
  return(list(res=res,log2FC_perChr=log2FC_perChr,de=de,deTable=deTab,mergeMap=mergeMap,dds=dds,rldCor = rldCor))
}




saveRDS(mlds_vs_2n_res,'DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_res.RDS')
mlds_vs_2n_res = readRDS('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_res.RDS')

mlds_vs_2n_dds = mlds_vs_2n_res[['dds']]
mlds_vs_2n_results = mlds_vs_2n_res[['res']] %>% as.data.frame()
mlds_vs_2n_results = annotateGenes(mlds_vs_2n_results,geneMap = geneMap)

mlds_vs_2n_deg = mlds_vs_2n_res[['de']]
mlds_vs_2n_deg$direction = ifelse(mlds_vs_2n_deg$log2FoldChange > 0,'MLDS_up','MLDS_down')
mlds_vs_2n_deg$pct.diff = mlds_vs_2n_deg$cellFrac_g1 - mlds_vs_2n_deg$cellFrac_g2
mlds_vs_2n_deg$cellFrac_max = pmax(mlds_vs_2n_deg$cellFrac_g1, mlds_vs_2n_deg$cellFrac_g2)

max_cellfrac = 10
min_l2FC = 0.1


#mlds_vs_2n_deg = mlds_vs_2n_deg[abs(mlds_vs_2n_deg$log2FoldChange) >= min_l2FC & abs(mlds_vs_2n_deg$pct.diff) >= 10,]
mlds_vs_2n_deg = mlds_vs_2n_deg[abs(mlds_vs_2n_deg$log2FoldChange) >= min_l2FC & abs(mlds_vs_2n_deg$cellFrac_max) >= max_cellfrac,]
#mlds_vs_2n_deg = mlds_vs_2n_deg[abs(mlds_vs_2n_deg$log2FoldChange) >= 1 & abs(mlds_vs_2n_deg$cellFrac_max) >= 20,]

table(mlds_vs_2n_deg$direction)
dim(mlds_vs_2n_deg)
mlds_vs_2n_deg = mlds_vs_2n_deg[order(abs(mlds_vs_2n_deg$pct.diff),decreasing = T),]

## Local vs Global impact
table(mlds_vs_2n_deg$direction,mlds_vs_2n_deg$chr)
write.csv(mlds_vs_2n_deg,'DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_topDEGs.csv')
mlds_vs_2n_deg = read.csv('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_topDEGs.csv')









##-----------------------------------------------------##
##    4. Extract MLDS imprints in T21  fLiver MEMP   ####
##----------------------------------------------------##

### 1. This is how one go from absolute normal diploid fetal liver to full-blown MLDS (alter the expression of these genes in the direction suggested)
# The universe of DE genes = Tum-vs-Dip
mlds_vs_2n_deg = read.csv('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset.assay_topDEGs.csv',row.names = 1)

# For a gene to be considered as DE, it must be expressed in at least 20% cells from 1 group
mlds_vs_2n_deg = mlds_vs_2n_deg %>% filter(!is.na(padj)) %>%  filter(padj < 0.05)
mlds_vs_2n_deg$DE = ifelse((mlds_vs_2n_deg$cellFrac_g1 >= max_cellfrac | mlds_vs_2n_deg$cellFrac_g2 >= max_cellfrac),T,F)
#mlds_vs_2n_deg$DE = T
table(mlds_vs_2n_deg$DE)
sum(mlds_vs_2n_deg$log2FoldChange[mlds_vs_2n_deg$DE == T]>0)
sum(mlds_vs_2n_deg$log2FoldChange[mlds_vs_2n_deg$DE == T]<0)


# 2. Now, which of these genes were altered due to T21
#t21_vs_2n_res = read.csv('DESeq2_dip.vs.T21_MEMP_allGeneRes.csv',row.names = 1)

t21_vs_2n_res = readRDS('DESeq2_T21.vs.diploid.fLiver.MEMP_geno.assay_res.RDS')
t21_vs_2n_result = as.data.frame(t21_vs_2n_res[['res']])
t21_vs_2n_result = annotateGenes(t21_vs_2n_result,geneMap = geneMap)
t21_vs_2n_result$cellFrac_g1 = apply(cnt_mtx[match(t21_vs_2n_result$ensID,rownames(cnt_mtx)),colnames(cnt_mtx) %in% mDat$cellID[mDat$Genotype == 'diploid']],1,function(x){sum(x>0.55)})
t21_vs_2n_result$cellFrac_g1 = 100*t21_vs_2n_result$cellFrac_g1 / ncol(cnt_mtx[,colnames(cnt_mtx) %in% mDat$cellID[mDat$Genotype == 'diploid']])

t21_vs_2n_result$cellFrac_g2 = apply(cnt_mtx[match(t21_vs_2n_result$ensID,rownames(cnt_mtx)),colnames(cnt_mtx) %in% mDat$cellID[mDat$Genotype == 'T21']],1,function(x){sum(x>0.55)})
t21_vs_2n_result$cellFrac_g2 = 100*t21_vs_2n_result$cellFrac_g2 / ncol(cnt_mtx[,colnames(cnt_mtx) %in% mDat$cellID[mDat$Genotype == 'T21']])  

t21_vs_2n_result$pct.diff = t21_vs_2n_result$cellFrac_g1 - t21_vs_2n_result$cellFrac_g2
t21_vs_2n_result$cellFrac_max = pmax(t21_vs_2n_result$cellFrac_g1, t21_vs_2n_result$cellFrac_g2)



## New criteria to select gene module:
## 1. DE in MLDS_vs_Diploid comparison
## 2. also changed in the same direction in T21_vs_diploid comparison
mlds_degs = mlds_vs_2n_deg[mlds_vs_2n_deg$DE == T,]
#mlds_upGenes_T21up = t21_vs_2n_result[t21_vs_2n_result$log2FoldChange > 0.25 & t21_vs_2n_result$geneSym %in% mlds_upGenes$geneSym,]
mlds_imprints_in_T21 = rbind(t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_up'] & t21_vs_2n_result$log2FoldChange > 0 & !is.na(t21_vs_2n_result$padj),],
                             t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_down'] & t21_vs_2n_result$log2FoldChange < 0 & !is.na(t21_vs_2n_result$padj),])

mlds_imprints_in_T21 = rbind(t21_vs_2n_result[t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_up'] & t21_vs_2n_result$log2FoldChange > 0 & !is.na(t21_vs_2n_result$padj),],
                             t21_vs_2n_result[t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_down'] & t21_vs_2n_result$log2FoldChange < 0 & !is.na(t21_vs_2n_result$padj),])


mlds_imprints_in_T21$padj = p.adjust(mlds_imprints_in_T21$pvalue,method = 'BH')
table(mlds_imprints_in_T21$padj < 0.05,mlds_imprints_in_T21$log2FoldChange > 0)
mlds_imprints_in_T21 = mlds_imprints_in_T21[mlds_imprints_in_T21$padj < 0.05,]
mlds_imprints_in_T21 = mlds_imprints_in_T21[order(abs(mlds_imprints_in_T21$log2FoldChange),decreasing = T),]
mlds_imprints_in_T21 = annotateGenes(mlds_imprints_in_T21,geneMap = geneMap)
table(mlds_imprints_in_T21$chr,mlds_imprints_in_T21$log2FoldChange > 0)
table(mlds_imprints_in_T21$log2FoldChange > 0)

mlds_imprints_in_T21$direction = ifelse(mlds_imprints_in_T21$log2FoldChange > 0 , 'T21_MLDS_up','T21_MLDS_down')
write.csv(mlds_imprints_in_T21,'MLDS_imprint_in_T21.fLiver.MEMP_v2.csv',row.names = T)
mlds_imprints_in_T21 = read.csv('MLDS_imprint_in_T21.fLiver.MEMP.csv',row.names = 1)

## Having done the comparison, I have ordered these in order of most relaxed --> most stringent.
## And the more stringent options do not identify new genes, just remove a few genes from the option before that
## I think I'm gonna go with the cellFrac20 option - not too strict, not too relaxed?
genes_cellFrac10 = mlds_imprints_in_T21
genes_cellFrac20 = mlds_imprints_in_T21
genes_pctDiff10 = mlds_imprints_in_T21
genes_pctDiff20 = mlds_imprints_in_T21


Idents(fLiver) = fLiver$annot_jun24
DotPlot(fLiver,idents = 'MEMP_MEP',group.by = 'Genotype',
        features = c(mlds_imprints_in_T21$geneSym[mlds_imprints_in_T21$log2FoldChange < 0],
                     mlds_imprints_in_T21$geneSym[mlds_imprints_in_T21$log2FoldChange > 0]))+
  #features = genes_cellFrac10$geneSym[genes_cellFrac10$log2FoldChange > 0 & genes_cellFrac10$geneSym %in% c(genes_cellFrac20$geneSym)])+
  RotatedAxis()






##--------------------------------------------------##
##    Deriving a module - excluding chr21 genes   ####
##--------------------------------------------------##
mlds_degs = mlds_vs_2n_deg[mlds_vs_2n_deg$DE == T & mlds_vs_2n_deg$chr != '21',]
#mlds_upGenes_T21up = t21_vs_2n_result[t21_vs_2n_result$log2FoldChange > 0.25 & t21_vs_2n_result$geneSym %in% mlds_upGenes$geneSym,]
mlds_imprints_in_T21 = rbind(t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_up'] & t21_vs_2n_result$log2FoldChange > 0 & !is.na(t21_vs_2n_result$padj),],
                             t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_down'] & t21_vs_2n_result$log2FoldChange < 0 & !is.na(t21_vs_2n_result$padj),])

mlds_imprints_in_T21$padj = p.adjust(mlds_imprints_in_T21$pvalue,method = 'BH')
table(mlds_imprints_in_T21$padj < 0.05,mlds_imprints_in_T21$log2FoldChange > 0)
mlds_imprints_in_T21 = mlds_imprints_in_T21[mlds_imprints_in_T21$padj < 0.05,]
mlds_imprints_in_T21 = mlds_imprints_in_T21[order(abs(mlds_imprints_in_T21$log2FoldChange),decreasing = T),]
mlds_imprints_in_T21 = annotateGenes(mlds_imprints_in_T21,geneMap = geneMap)
table(mlds_imprints_in_T21$chr,mlds_imprints_in_T21$log2FoldChange > 0)
table(mlds_imprints_in_T21$log2FoldChange > 0)

## Check that by removing Chr21 genes before re-calculating p-value, the result yield the same set of non-chr21 genes as the gene list above (correcting p-value before removal of chr21 genes)
table(mlds_imprints_in_T21$geneSym %in% genes_cellFrac20$geneSym)


mlds_imprints_in_T21$direction = ifelse(mlds_imprints_in_T21$log2FoldChange > 0 , 'T21_MLDS_up','T21_MLDS_down')
write.csv(mlds_imprints_in_T21,'MLDS_imprint_in_T21.fLiver.MEMP_noChr21.csv',row.names = T)
mlds_imprints_in_T21 = read.csv('MLDS_imprint_in_T21.fLiver.MEMP_noChr21.csv',row.names = 1)



##-------------------------------------##
##   Plot expression of gene module  ####
##-------------------------------------##
genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/MLDS_imprint_in_T21.fLiver.MEMP_v2.csv',row.names = 1)
genes = mlds_imprints_in_T21
genes = rbind(genes[genes$direction == 'T21_MLDS_up',],
              genes[genes$direction == 'T21_MLDS_down',])


##---- In MLDS
mlds_mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
mlds$broadLineage = mlds_mdat$broadLineage[match(mlds$cellID,mlds_mdat$cellID)]

mlds$group_4avgExpr = ifelse(mlds$annot_mar24 == 'MEP' & mlds$donorID %in% c('L019','L038','L039','L040'),'MEMP_MEP',
                             ifelse(mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID !='L041' & mlds$tissue == 'BM',paste0(mlds$annot_mar24,'.',mlds$disease,'.',mlds$donorID),
                                    ifelse(mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID !='L041' & mlds$disease == 'TAM',paste0(mlds$annot_mar24,'.',mlds$disease,'.',mlds$donorID),
                                           ifelse(mlds$annot_mar24 == 'Tumour' & mlds$donorID == 'L038' & mlds$timePoint == 'TP1','Tumour.MLDS.L038_TP1',
                                                  ifelse(mlds$annot_mar24 == 'Tumour' & mlds$donorID == 'L076' & mlds$timePoint == 'D.Relapse','Tumour.MLDS.L076_RD',
                                                         ifelse(mlds$annot_mar24 == 'Tumour' & mlds$donorID == 'L076' & mlds$tissue == 'Blood','Tumour.MLDS.L076_Blood',
                                                                mlds$broadLineage))))))

mlds$group_4avgExpr = factor(mlds$group_4avgExpr,levels = c('MEMP_MEP',
                                                            unique(mlds$group_4avgExpr[grepl('TAM',mlds$group_4avgExpr)]),
                                                            unique(mlds$group_4avgExpr[grepl('MLDS',mlds$group_4avgExpr)]),
                                                            unique(mlds$group_4avgExpr[!grepl('TAM|MLDS|MEMP_MEP',mlds$group_4avgExpr)])))

DotPlot(mlds,group.by = 'group_4avgExpr',
        #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) >= 30]
        #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) <= 30 & abs(mlds.vs.gMEP_deg$pct.diff) >= 20]
        features = genes$geneSym,
        #features = gata1s_module$geneSym[gata1s_module$fLiver_geneGroup == 'MK' & abs(gata1s_module$pct.diff) >= 30][1:50]
)+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                   face=ifelse(genes$chr == '21','bold','plain'),
                                   colour =ifelse(genes$isTF == T,'purple',
                                                  ifelse(genes$isCSM == T,'blue','black'))),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 






##---- In fLiver
akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_annotated_0424.RDS'
akLiv_mdat_fp = gsub('_0424.RDS','_0424_mdat.csv',akLiv_srat_fp)

fLiver = readRDS(akLiv_srat_fp)
fLiver_mdat = read.csv(akLiv_mdat_fp,row.names = 1)
fLiver$broadLineage = fLiver_mdat$broadLineage[match(fLiver$cellID,fLiver_mdat$cellID)]
fLiver$Genotype = fLiver_mdat$Genotype[match(fLiver$cellID,fLiver_mdat$cellID)]
fLiver$Genotype[fLiver$Genotype == 'complete_trisomy'] = 'Triploid'
fLiver$Genotype[fLiver$Genotype == 'diploid'] = 'Diploid'
fLiver$finalAnn_broad = fLiver$annot_jun24
fLiver$Sex[fLiver$donorID %in% c('Hsb36','Hsb37')] = 'XX'

fLiver$Genotype = factor(fLiver$Genotype,c('Diploid','T21','MX','T18','T22','Triploid'))
Idents(fLiver) = paste0(fLiver$Genotype,':',fLiver$annot_jun24)
Idents(fLiver) = fLiver$annot_jun24

plotFun_Sup.Fig.3x_MLDSimprintsT21_FLexpr_dotplot=function(noFrame=FALSE,noPlot=FALSE){
  p = DotPlot(fLiver,idents = 'MEMP_MEP',group.by = 'Genotype',
              cols = c(colAlpha(grey(0.98),0.85),'black'),
              #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) >= 30]
              #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) <= 30 & abs(mlds.vs.gMEP_deg$pct.diff) >= 20]
              features = c(genes$geneSym[genes$direction == 'T21_MLDS_down'],
                           genes$geneSym[genes$direction == 'T21_MLDS_up']),
              #features = gata1s_module$geneSym[gata1s_module$fLiver_geneGroup == 'MK' & abs(gata1s_module$pct.diff) >= 30][1:50]
  )+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                     face = ifelse(genes$chr == '21','bold','plain'),
                                     colour=ifelse(genes$isTF[match(genes$geneSym,genes$geneSym)],'red',
                                                   ifelse(genes$isCSM[match(genes$geneSym,genes$geneSym)],'blue','black'))),
          axis.text.y = element_text(size=11),
          legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 
  print(p)  
  
}


saveFig(file.path(plotDir,'Sup.Fig3x_MLDS.imprints.inT21_FLexpr_dotplot'),plotFun_Sup.Fig.3x_MLDSimprintsT21_FLexpr_dotplot,rawData=df,width = 11,height = 4,res = 500,useDingbats = T)






##---- In old fLiver with published fLiver
fLiver_wPublishedData = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver_2n/jan24/liver_liverREFmerged_clean_processed_annotated_noUnknowns_0124.RDS')
fLiver_wPublishedData$finalAnn_broad = fLiver_wPublishedData$annot_jan24

fLiver$finalAnn_broad = fLiver$annot_jun24
fLiver$Sex[fLiver$donorID %in% c('Hsb36','Hsb37')] = 'XX'

Idents(fLiver) = fLiver$annot_jun24
DotPlot(fLiver,idents = unique(fLiver$annot_jun24[grepl('MEP',fLiver$annot_jun24)]),group.by = 'Genotype',
        #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) >= 30]
        #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) <= 30 & abs(mlds.vs.gMEP_deg$pct.diff) >= 20]
        features = genes$geneSym
        #features = gata1s_module$geneSym[gata1s_module$fLiver_geneGroup == 'MK' & abs(gata1s_module$pct.diff) >= 30][1:50]
)+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 




##------- In big.srat
big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/combinedLeuk_2404.RDS')

## Dot plot of MLDS-sepcific genes in other leukaemia
big.srat$timePoint[big.srat$timePoint == 'Diagnostic'] = 'D0'
big.srat$broadLineage[big.srat$finalAnn_broad == 'Tumour' & !is.na(big.srat$finalAnn_broad)] = 'Tumour'
big.srat$group_4_dotPlot = ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID != 'L041' & big.srat$timePoint %in% c('D0','Diagnostic','RelapseD0') & !big.srat$disease %in% c('TAM','MLDS'), 
                                  paste0('Tumour.',big.srat$disease,'.',big.srat$timePoint,'.',big.srat$tissue),
                                  ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID != 'L041' & big.srat$timePoint %in% c('D0','Diagnostic','RelapseD0') & big.srat$disease %in% c('TAM','MLDS'),
                                         paste0('Tumour.',big.srat$disease,'.',big.srat$timePoint,'.',big.srat$tissue,'.',big.srat$donorID),
                                         ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID == 'L038' & big.srat$timePoint %in% c('TP1') & big.srat$disease %in% c('MLDS'),
                                                'Tumour.MLDS.TP1.BM.L038',
                                                ifelse(big.srat$broadLineage == 'Tumour' & big.srat$disease == 'Lymphoma','Tumour.Lymphoma.D0.L062','others'))))
big.srat$group_4_dotPlot[grepl('^Tumour\\.BALL\\.D0',big.srat$group_4_dotPlot) & big.srat$Genotype == 'T21'] = gsub('Tumour.BALL.D0','Tumour.BALL.T21.D0',big.srat$group_4_dotPlot[grepl('^Tumour\\.BALL\\.D0',big.srat$group_4_dotPlot) & big.srat$Genotype == 'T21'])
big.srat$group_4_dotPlot[is.na(big.srat$group_4_dotPlot)] = 'others'
big.srat$group_4_dotPlot = factor(big.srat$group_4_dotPlot,levels = c(unique(big.srat$group_4_dotPlot[grepl('Tumour.MLDS',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.TAM',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.MDS',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.AMKL',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.BALL\\.T21',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.BALL\\.D0',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.BALL.Relapse',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.iALL',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.Lymphoma',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.pAML',big.srat$group_4_dotPlot)]),
                                                                      'others'
))

big.srat$group_4_dotPlot = paste0(big.srat$group_4_dotPlot,':',big.srat$donorID)
Idents(big.srat) = big.srat$group_4_dotPlot

DotPlot(big.srat,idents = unique(big.srat$group_4_dotPlot[grepl('^Tumour.*Diag|Tumour.*D0',big.srat$group_4_dotPlot)]),
        features = genes$geneSym)+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                   face = ifelse(genes$chr == '21','bold','plain'),
                                   colour=ifelse(genes$isTF[match(genes$geneSym,genes$geneSym)],'red',
                                                 ifelse(genes$isCSM[match(genes$geneSym,genes$geneSym)],'blue','black'))),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 






##--------    Plots  ----------##
## Make a strange plot
mlds_degs$diploid_vs_t21 = ifelse(mlds_degs$direction == 'MLDS_up' & mlds_degs$geneSym %in% mlds_imprints_in_T21$geneSym[mlds_imprints_in_T21$direction == 'T21_MLDS_up'],'T21_up',
                                  ifelse(mlds_degs$direction == 'MLDS_down' & mlds_degs$geneSym %in% mlds_imprints_in_T21$geneSym[mlds_imprints_in_T21$direction == 'T21_MLDS_down'],'T21_down',F))

df = as.data.frame(table(mlds_degs$direction,mlds_degs$diploid_vs_t21))
colnames(df) = c('MLDS_deg','T21_deg','nGene')
df2 = df %>% group_by(MLDS_deg) %>% filter(MLDS_deg != 'nonDE') %>% summarise(nGene = sum(nGene))
df2$group = 'MLDS vs 2n_MEMP_fLiver'
colnames(df2)[1] = 'direction'
df = df[df$T21_deg == 'T21_down' & df$MLDS_deg == 'MLDS_down' |
          df$T21_deg == 'T21_up' & df$MLDS_deg == 'MLDS_up',]
df = df[,colnames(df) != 'MLDS_deg']
df$group = 'T21 vs 2n_MEMP_fLiver'
colnames(df)[1] = 'direction'
df = rbind(df,df2)
df$nGene[grepl('down',df$direction)] = -df$nGene[grepl('down',df$direction)]
# df$nGene2 = log10(abs(df$nGene))
# df$nGene2[df$nGene < 0] = -df$nGene2[df$nGene < 0]
plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_0424/Plots/'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
plotFun_MLDS.T21_DEGs_barplot = function(noFrame=FALSE,noPlot=FALSE){
  library(ggbreak) 
  
  p1 = ggplot(df,aes(group,nGene,fill=direction))+
    geom_col(width = 0.58)+
    #scale_fill_manual(values = c(col25[1],col25[2],colAlpha(col25[1:2],alphas = 0.6)))+
    scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
    geom_hline(yintercept = 0,col='black',lwd=0.5)+
    #scale_y_log10()+
    scale_y_break(c(40,1560),scales = 0.8,expand = T) +
    scale_y_break(c(-40,-1870),scales = 4) +
    #scale_y_break(c(-40,-1880,40,1550),ticklabels = c(1550,1580,40,20,0,-20,-40,-1893)) +
    #geom_hline(yintercept = df$nGene,col='black',lty=2,lwd=0.3)+
    theme_classic(base_size = 11)+xlab('')+ylab('# DEGs')+
    theme(#panel.border = element_rect(fill=F),
      axis.line = element_blank(),#legend.position = 'bottom',
      legend.title = element_text(size=10),legend.text = element_text(size=8),legend.key.size = unit(0.5,'cm'),
      axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
      axis.ticks = element_line(colour='black'),
      axis.text = element_text(color='black'))
  print(p1)
}

saveFig(file.path(plotDir,paste0('Fig3x_MLDS.T21.DEGs_barplot')),plotFun_MLDS.T21_DEGs_barplot,rawData=df,width = 3.2,height = 4,res = 500,useDingbats = T)



df = df[df$group == 'T21 vs 2n_MEMP_fLiver',]
plotFun_MLDS.T21_DEGs_barplot_onlyT21imprints = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(df,aes(group,nGene,fill=direction))+
    geom_col(width = 0.58)+
    scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
    geom_hline(yintercept = 0,col='black',lwd=0.5)+
    
    # scale_y_break(c(40,1560),scales = 0.8,expand = T) +
    # scale_y_break(c(-40,-1870),scales = 4) +
    #scale_y_break(c(-40,-1880,40,1550),ticklabels = c(1550,1580,40,20,0,-20,-40,-1893)) +
    #scale_y_continuous(limits = c(-20,50))+
    ylim(-200,300)+
    theme_classic(base_size = 11)+xlab('')+ylab('# DEGs')+
    theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
          axis.line = element_blank(),#legend.position = 'bottom',
          legend.title = element_text(size=10),legend.text = element_text(size=8),legend.key.size = unit(0.5,'cm'),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
          axis.text = element_text(color='black'))
  print(p1)
}

#saveFig(file.path(plotDir,paste0('Fig3x_MLDS.T21.DEGs_barplot_simplified')),plotFun_MLDS.T21_DEGs_barplot_onlyT21imprints,rawData=df,width = 2.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,paste0('Fig3x_MLDS.T21.DEGs_barplot_simplified_sameScaleGATA1smodule')),plotFun_MLDS.T21_DEGs_barplot_onlyT21imprints,rawData=df,width = 2.45,height = 5.3,res = 500,useDingbats = T)





