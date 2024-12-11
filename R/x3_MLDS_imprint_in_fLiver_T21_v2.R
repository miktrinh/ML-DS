##----   Deriving MLDS imprint in T21 foetal MEMP, using good TAM only    -----##


outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/3_MLDSimprint_in_fLiverT21/jul24'
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




##---------------------------------------------------------##
#### 1. Import relevant scRNA datasets: MLDS + fLiver    ####
##---------------------------------------------------------##

## Import MLDS object
mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS')
# mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_mdat.csv',row.names = 1)
# table(mlds$cellID %in% mdat$cellID)
# mlds@meta.data = mdat[match(mlds$cellID,mdat$cellID),]
# rownames(mlds@meta.data) = mlds$cellID
mlds$annot = as.character(mlds$annot_aug24)
DimPlot(mlds,group.by = 'donorID',cols = c(col25,pal34H),label.size = 3,label = T,repel = T,label.box = T) + NoLegend() + NoAxes()
DimPlot(mlds,group.by = 'annot',cols = c(col25,pal34H),label.size = 3,label = T,repel = T,label.box = T) + NoLegend() + NoAxes()
ee_cells = WhichCells(mlds,cells = mlds$cellID[mlds$annot_aug24 == 'MEP'],expression =(IL1B == 0 & APOC1 > 2 & ITGA2B < 0.5))
DimPlot(mlds,cells.highlight = ee_cells)
mlds$group = ifelse(mlds$cellID %in% ee_cells,'EE_1',mlds$annot_aug24)
Idents(mlds)=mlds$group
mlds$annot[mlds$annot == 'MEP' & mlds$cellID %in% ee_cells] = 'EE'


## Import fLiver object
akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'
akLiv_mdat_fp = gsub('_0824.RDS','_0824_mdat.csv',akLiv_srat_fp)

fLiver = readRDS(akLiv_srat_fp)
fLiver$annot = fLiver$annot_aug24
# fLiver_mdat = read.csv(akLiv_mdat_fp,row.names = 1)
# fLiver$broadLineage = fLiver_mdat$broadLineage[match(fLiver$cellID,fLiver_mdat$cellID)]
all(fLiver$finalAnn_broad == fLiver$annot)
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
cnt_mtx = fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$cellID %in% c(fLiver$cellID[fLiver$Genotype %in% c('diploid','T21') & fLiver$annot == 'MEMP_MEP'])]]
rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = fLiver@meta.data[fLiver$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','Sex','Genotype','assay')]
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
t21_vs_2n_result = as.data.frame(t21_vs_2n_res[['res']])
t21_vs_2n_result = annotateGenes(t21_vs_2n_result,geneMap = geneMap)
t21_vs_2n_result$cellFrac_g1 = apply(cnt_mtx[match(t21_vs_2n_result$ensID,rownames(cnt_mtx)),colnames(cnt_mtx) %in% mDat$cellID[mDat$Genotype == 'diploid']],1,function(x){sum(x>0)})
t21_vs_2n_result$cellFrac_g1 = t21_vs_2n_result$cellFrac_g1 / ncol(cnt_mtx[,colnames(cnt_mtx) %in% mDat$cellID[mDat$Genotype == 'diploid']])

t21_vs_2n_result$cellFrac_g2 = apply(cnt_mtx[match(t21_vs_2n_result$ensID,rownames(cnt_mtx)),colnames(cnt_mtx) %in% mDat$cellID[mDat$Genotype == 'T21']],1,function(x){sum(x>0)})
t21_vs_2n_result$cellFrac_g2 = t21_vs_2n_result$cellFrac_g2 / ncol(cnt_mtx[,colnames(cnt_mtx) %in% mDat$cellID[mDat$Genotype == 'T21']])  

t21_vs_2n_result$pct.diff = t21_vs_2n_result$cellFrac_g1 - t21_vs_2n_result$cellFrac_g2
t21_vs_2n_result$cellFrac_max = pmax(t21_vs_2n_result$cellFrac_g1, t21_vs_2n_result$cellFrac_g2)

write.csv(t21_vs_2n_result,'DESeq2_T21.vs.diploid.fLiver.MEMP_geno.assay_results_allGenes.csv')







##----------------------------------------------------##
##   3. DEGs between MLDS vs diploid 2n fLiver      ####
##----------------------------------------------------##

##---- Cluster MLDS / TAM blasts + MEP + fLiver cells   -----------
# NOTE: try Harmony over "assay" --> created a mess...
fLiver = subset(fLiver,subset = cellID %in% fLiver$cellID[fLiver$annot %in% c('MEMP_MEP','HSC_MPP','Mast.cell','MK') &
                                                            fLiver$Genotype %in% c('diploid','T21')])

mlds.sub = subset(mlds,subset = cellID %in% mlds$cellID[mlds$annot %in% c('HSC_MPP','EE','MEP','MK','Tumour','unsure_Tum_MK?')])
srat = merge_seurat_objects(fLiver,mlds.sub,keepAllGenes = F,genomeVersions = c('v38','v38'))
srat$tissue[srat$cellID %in% fLiver$cellID] = 'fLiver'
srat = standard_clustering(srat,runHarmony = T,harmonyVar = 'tissue')
DimPlot(srat,group.by = 'annot',cols = col25,label = T,repel = T,label.size = 3)
FeaturePlot(srat,'GATA1')
Idents(srat) = paste0(srat$annot,'_',srat$Genotype,'_',srat$tissue,'_',srat$donorID)
p = DotPlot(srat,
            idents = unique(Idents(srat)[grepl('MEMP_MEP_diploid_fLiver|Tumour_T21_BM',Idents(srat)) &
                                           !grepl('L041|CC3',Idents(srat))]),
            group.by = 'donorID',
            features = geneMap$geneSym[geneMap$chr == 'chr21' & !geneMap$geneSym %in% mlds_vs_2n_deg$geneSym]
)+RotatedAxis() +
  theme(axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = 'top') + xlab('') + ylab('')
p

DimPlot(srat,cells.highlight = srat$cellID[srat$annot == 'MEMP_MEP'])






##------- Using pseudobulk -----------
#cnt_mtx.FL = fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$cellID %in% c(fLiver$cellID[fLiver$Genotype %in% c('diploid') & fLiver$annot == 'MEMP_MEP'])]]
cnt_mtx.FL = fLiver@assays$RNA@counts[,fLiver$cellID[fLiver$cellID %in% c(fLiver$cellID[fLiver$Genotype %in% c('diploid') & fLiver$annot == 'MEMP_MEP' & fLiver$assay == 'GEX5p'])]]
cnt_mtx.mlds = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in%
                                                     c(mlds$cellID[mlds$disease == 'MLDS' & mlds$annot == 'Tumour' & 
                                                                     #mlds$timePoint %in% c('Diagnostic','D.Relapse','D.Relapse2') &
                                                                     mlds$timePoint %in% c('Diagnostic') & 
                                                                     mlds$tissue == 'BM' & !mlds$donorID %in% c('L041')])]]
genes = intersect(rownames(cnt_mtx.mlds),rownames(cnt_mtx.FL))
cnt_mtx=cbind(cnt_mtx.FL[genes,],cnt_mtx.mlds[genes,])

rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat.mlds = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex','assay')]
table(is.na(mDat.mlds$sex))
mDat.fl = fLiver@meta.data[fLiver$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex','assay')]
mDat = rbind(mDat.mlds,mDat.fl)
mDat$dataset = ifelse(mDat$cellID %in% mDat.mlds$cellID,'MLDS','fLiver')
mDat$sex[mDat$sex == 'XX'] = 'F'
mDat$sex[mDat$sex == 'XY'] = 'M'
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


#saveRDS(mlds_vs_2n_res,'DESeq2_diag.relapse.MLDS.vs.2nfLiverMEMP_dataset_res.RDS')
#mlds_vs_2n_res = readRDS('DESeq2_diag.relapse.MLDS.vs.2nfLiverMEMP_dataset_res.RDS')
saveRDS(mlds_vs_2n_res,'DESeq2_dMLDS.vs.2nfLiverMEMP_dataset_res.RDS')
mlds_vs_2n_res = readRDS('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset_res.RDS')

# mlds_vs_2n_dds = mlds_vs_2n_res[['dds']]
# mlds_vs_2n_results = mlds_vs_2n_res[['res']] %>% as.data.frame()
# mlds_vs_2n_results = annotateGenes(mlds_vs_2n_results,geneMap = geneMap)

mlds_vs_2n_deg = mlds_vs_2n_res[['de']]
mlds_vs_2n_deg$direction = ifelse(mlds_vs_2n_deg$log2FoldChange > 0,'MLDS_up','MLDS_down')
mlds_vs_2n_deg$pct.diff = mlds_vs_2n_deg$cellFrac_g1 - mlds_vs_2n_deg$cellFrac_g2
mlds_vs_2n_deg$cellFrac_max = pmax(mlds_vs_2n_deg$cellFrac_g1, mlds_vs_2n_deg$cellFrac_g2)

max_cellfrac = 10/100
min_l2FC = 0.5

mlds_vs_2n_deg = mlds_vs_2n_deg[abs(mlds_vs_2n_deg$log2FoldChange) >= min_l2FC & abs(mlds_vs_2n_deg$cellFrac_max) >= max_cellfrac,]

table(mlds_vs_2n_deg$direction)
dim(mlds_vs_2n_deg)
mlds_vs_2n_deg = mlds_vs_2n_deg[order(abs(mlds_vs_2n_deg$pct.diff),decreasing = T),]

## Local vs Global impact
table(mlds_vs_2n_deg$direction,mlds_vs_2n_deg$chr)
# write.csv(mlds_vs_2n_deg,'DESeq2_diag.relapse.MLDS.vs.2nfLiverMEMP_dataset_topDEGs.csv')
# mlds_vs_2n_deg = read.csv('DESeq2_diag.relapse.MLDS.vs.2nfLiverMEMP_dataset_topDEGs.csv')
write.csv(mlds_vs_2n_deg,'DESeq2_dMLDS.vs.2nfLiverMEMP_dataset_topDEGs.csv')
mlds_vs_2n_deg = read.csv('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset_topDEGs.csv')

##----------------------------------------------##
##    Gene Expression stats in MLDS blasts    ####
##----------------------------------------------##
group = paste0(mlds$annot_aug24,':',mlds$donorID,':',mlds$timePoint)
names(group) = mlds$cellID
mldsLeuk_mtx_list = lapply(split(mlds$cellID[mlds$annot_aug24 == 'Tumour'],group[match(mlds$cellID[mlds$annot_aug24 == 'Tumour'],names(group))]),function(e){
  mtx = mlds@assays$RNA@counts[,e]
  g = unique(group[match(colnames(mtx),names(group))])
  tmp_nCell = as.data.frame(rowSums(mtx>0))
  colnames(tmp_nCell) = g
  tmp_percCell = as.data.frame(100*rowSums(mtx>0) / ncol(mtx))
  colnames(tmp_percCell) = g
  return(list('nCell'=tmp_nCell,'percCell'=tmp_percCell))
  })

# Extract percentage of cells expressing each genes
mlds_geneStats_percCells = do.call(cbind,sapply(mldsLeuk_mtx_list,'[[',2)) %>%  as.data.frame()
rownames(mlds_geneStats_percCells) = rownames(mldsLeuk_mtx_list[[1]][[2]])
write.csv(mlds_geneStats_percCells,'geneStats_percentageCellsExpressed.csv')
mlds_geneStats_nCells = do.call(cbind,sapply(mldsLeuk_mtx_list,'[[',1)) %>%  as.data.frame()
rownames(mlds_geneStats_nCells) = rownames(mldsLeuk_mtx_list[[1]][[1]])
write.csv(mlds_geneStats_nCells,'geneStats_fractionCellsExpressed.csv')

## How many genes expressed in >= 10% cells per group
table(rowSums(mlds_geneStats_percCells[,grepl('Diagnostic',colnames(mlds_geneStats_percCells)) & !grepl('CC3|L041|CC1|CC2|CC6|CC7|CC8|L075|L114|L182|L156',colnames(mlds_geneStats_percCells))] > 10))
genes_highlyExpressed = rowSums(mlds_geneStats_percCells[,grepl('Diagnostic',colnames(mlds_geneStats_percCells)) & !grepl('CC3|L041|CC1|CC2|CC6|CC7|CC8|L075|L114|L182|L156',colnames(mlds_geneStats_percCells))] > 10)
## There are 4053 genes expressed in >= 10% blasts across all 10 MLDS diagnostic samples
table(mlds_vs_2n_deg$geneSym %in% names(genes_highlyExpressed[genes_highlyExpressed >= 3]),mlds_vs_2n_deg$direction)
## Of 4053 genes expressed in >= 10% blasts across all 10 MLDS diagnostic samples, only 93 genes were marked as DEGs (90 up and 3 down)
## Of X, only 175 up and 16 down DEGs





##-----------------------------------------------------##
##    4. Extract MLDS imprints in T21  fLiver MEMP   ####
##----------------------------------------------------##

### 1. This is how one go from absolute normal diploid fetal liver to full-blown MLDS (alter the expression of these genes in the direction suggested)
# The universe of DE genes = Tum-vs-Dip
mlds_vs_2n_deg = read.csv('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset_topDEGs.csv',row.names = 1)

# For a gene to be considered as DE, it must be expressed in at least 20% cells from 1 group
mlds_vs_2n_deg = mlds_vs_2n_deg %>% filter(!is.na(padj)) %>%  filter(padj < 0.05)
mlds_vs_2n_deg$DE = ifelse((mlds_vs_2n_deg$cellFrac_g1 >= max_cellfrac | mlds_vs_2n_deg$cellFrac_g2 >= max_cellfrac),T,F)
#mlds_vs_2n_deg$DE = T
table(mlds_vs_2n_deg$DE)
sum(mlds_vs_2n_deg$log2FoldChange[mlds_vs_2n_deg$DE == T]>0)
sum(mlds_vs_2n_deg$log2FoldChange[mlds_vs_2n_deg$DE == T]<0)


# 2. Now, which of these genes were altered due to T21
#t21_vs_2n_res = read.csv('DESeq2_dip.vs.T21_MEMP_allGeneRes.csv',row.names = 1)

t21_vs_2n_result = read.csv('DESeq2_T21.vs.diploid.fLiver.MEMP_geno.assay_results_allGenes.csv',row.names = 1)


## New criteria to select gene module:
## 1. DE in MLDS_vs_Diploid comparison
## 2. also changed in the same direction in T21_vs_diploid comparison
## 3. expressed in at least max_cellfrac % of 1 group T21 or 2n fLiver
mlds_degs = mlds_vs_2n_deg[mlds_vs_2n_deg$DE == T,]

# mlds_imprints_in_T21 = rbind(t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_up'] & t21_vs_2n_result$log2FoldChange > 0 & !is.na(t21_vs_2n_result$pvalue),],
#                              t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_down'] & t21_vs_2n_result$log2FoldChange < 0 & !is.na(t21_vs_2n_result$pvalue),])
mlds_imprints_in_T21 = rbind(t21_vs_2n_result[t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_up'] & t21_vs_2n_result$log2FoldChange > 0 & !is.na(t21_vs_2n_result$pvalue),],
                             t21_vs_2n_result[t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_down'] & t21_vs_2n_result$log2FoldChange < 0 & !is.na(t21_vs_2n_result$pvalue),])


mlds_imprints_in_T21$padj = p.adjust(mlds_imprints_in_T21$pvalue,method = 'BH')
table(mlds_imprints_in_T21$padj < 0.05,mlds_imprints_in_T21$log2FoldChange > 0)
mlds_imprints_in_T21 = mlds_imprints_in_T21[mlds_imprints_in_T21$padj < 0.05,]
mlds_imprints_in_T21 = mlds_imprints_in_T21[order(abs(mlds_imprints_in_T21$log2FoldChange),decreasing = T),]
mlds_imprints_in_T21 = annotateGenes(mlds_imprints_in_T21,geneMap = geneMap)
table(mlds_imprints_in_T21$chr,mlds_imprints_in_T21$log2FoldChange > 0)
table(mlds_imprints_in_T21$log2FoldChange > 0)

mlds_imprints_in_T21$direction = ifelse(mlds_imprints_in_T21$log2FoldChange > 0 , 'T21_MLDS_up','T21_MLDS_down')
write.csv(mlds_imprints_in_T21,'MLDS_imprint_in_T21.fLiver.MEMP_2411.csv',row.names = T)
mlds_imprints_in_T21 = read.csv('MLDS_imprint_in_T21.fLiver.MEMP_2411.csv',row.names = 1)

## Having done the comparison, I have ordered these in order of most relaxed --> most stringent.
## And the more stringent options do not identify new genes, just remove a few genes from the option before that
## I think I'm gonna go with the cellFrac20 option - not too strict, not too relaxed?


Idents(fLiver) = fLiver$annot
DotPlot(fLiver,idents = 'MEMP_MEP',group.by = 'Genotype',
        features = c(mlds_imprints_in_T21$geneSym[mlds_imprints_in_T21$log2FoldChange < 0],
                     mlds_imprints_in_T21$geneSym[mlds_imprints_in_T21$log2FoldChange > 0]))+
  #features = genes_cellFrac10$geneSym[genes_cellFrac10$log2FoldChange > 0 & genes_cellFrac10$geneSym %in% c(genes_cellFrac20$geneSym)])+
  RotatedAxis()










##--------------------------------------------------##
##    Deriving a module - excluding chr21 genes   ####
##--------------------------------------------------##
mlds_degs = mlds_vs_2n_deg[mlds_vs_2n_deg$DE == T & mlds_vs_2n_deg$chr != 'chr21',]
mlds_imprints_in_T21 = rbind(t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_up'] & t21_vs_2n_result$log2FoldChange > 0 & !is.na(t21_vs_2n_result$padj),],
                             t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_down'] & t21_vs_2n_result$log2FoldChange < 0 & !is.na(t21_vs_2n_result$padj),])

mlds_imprints_in_T21$padj = p.adjust(mlds_imprints_in_T21$pvalue,method = 'BH')
table(mlds_imprints_in_T21$padj < 0.05,mlds_imprints_in_T21$log2FoldChange > 0)
mlds_imprints_in_T21 = mlds_imprints_in_T21[mlds_imprints_in_T21$padj < 0.05,]
mlds_imprints_in_T21 = mlds_imprints_in_T21[order(abs(mlds_imprints_in_T21$log2FoldChange),decreasing = T),]
mlds_imprints_in_T21 = annotateGenes(mlds_imprints_in_T21,geneMap = geneMap)
table(mlds_imprints_in_T21$chr,mlds_imprints_in_T21$log2FoldChange > 0)
table(mlds_imprints_in_T21$log2FoldChange > 0)


mlds_imprints_in_T21$direction = ifelse(mlds_imprints_in_T21$log2FoldChange > 0 , 'T21_MLDS_up','T21_MLDS_down')
write.csv(mlds_imprints_in_T21,'MLDS_imprint_in_T21.fLiver.MEMP_noChr21.csv',row.names = T)
mlds_imprints_in_T21 = read.csv('MLDS_imprint_in_T21.fLiver.MEMP_noChr21.csv',row.names = 1)



##-------------------------------------##
##   Plot expression of gene module  ####
##-------------------------------------##
genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/MLDS_imprint_in_T21.fLiver.MEMP.csv',row.names = 1)
genes = mlds_imprints_in_T21
genes = rbind(genes[genes$direction == 'T21_MLDS_up',],
              genes[genes$direction == 'T21_MLDS_down',])


##---- In MLDS
# mlds_mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
# mlds$broadLineage = mlds_mdat$broadLineage[match(mlds$cellID,mlds_mdat$cellID)]

mlds$group_4avgExpr = ifelse(mlds$annot == 'MEP' & mlds$donorID %in% c('L019','L038','L039','L040'),'MEMP_MEP',
                             ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID !='L041' & mlds$tissue == 'BM',paste0(mlds$annot,'.',mlds$disease,'.',mlds$donorID),
                                    ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID !='L041' & mlds$disease == 'TAM',paste0(mlds$annot,'.',mlds$disease,'.',mlds$donorID),
                                           ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L038' & mlds$timePoint == 'TP1','Tumour.MLDS.L038_TP1',
                                                  ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L156' & mlds$timePoint == 'D.Relapse','Tumour.TAM.L156_RD',
                                                         ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L076' & mlds$timePoint == 'D.Relapse','Tumour.MLDS.L076_RD',
                                                                ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L076' & mlds$timePoint == 'D.Relapse2','Tumour.MLDS.L076_RD2',
                                                                       ifelse(mlds$annot == 'Tumour' & mlds$donorID == 'L076' & mlds$tissue == 'Blood','Tumour.MLDS.L076_Blood',
                                                                              mlds$broadLineage))))))))


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
fLiver$Genotype[fLiver$Genotype == 'complete_trisomy'] = 'Triploid'
fLiver$Genotype[fLiver$Genotype == 'diploid'] = 'Diploid'
fLiver$finalAnn_broad = fLiver$annot

fLiver$Genotype = factor(fLiver$Genotype,c('Diploid','T21','MX','T18','T22','Triploid'))
Idents(fLiver) = fLiver$annot

genes = rbind(genes[genes$direction == 'T21_MLDS_down',],
              genes[genes$direction == 'T21_MLDS_up',])
plotFun_Sup.Fig.XX_MLDSimprintsT21_FLexpr_dotplot=function(noFrame=FALSE,noPlot=FALSE){
  p = DotPlot(fLiver,idents = 'MEMP_MEP',group.by = 'Genotype',
              cols = c(colAlpha(grey(0.98),0.85),'black'),
              #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) >= 30]
              #features = mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$tam_vs_mempT21_group == 'TAM.MLDS.up' & abs(mlds.vs.gMEP_deg$pct.diff) <= 30 & abs(mlds.vs.gMEP_deg$pct.diff) >= 20]
              features = genes$geneSym,
              #features = gata1s_module$geneSym[gata1s_module$fLiver_geneGroup == 'MK' & abs(gata1s_module$pct.diff) >= 30][1:50]
  )+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                     face = ifelse(genes$chr == 'chr21','bold','plain'),
                                     colour=ifelse(genes$isTF,'red',
                                                   ifelse(genes$isCSM,'blue','black'))),
          axis.text.y = element_text(size=11,colour = 'black'),
          axis.line = element_blank(),
          panel.border = element_rect(fill=F,colour = 'black'),
          axis.ticks = element_line(color='black'),
          legend.position = 'right',legend.text = element_text(size=9),
          legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 
  print(p)  
  
}


saveFig(file.path(plotDir,'Sup.FigXX_MLDS.imprints.inT21_FLexpr_dotplot'),plotFun_Sup.Fig.XX_MLDSimprintsT21_FLexpr_dotplot,rawData=genes,width = 7,height = 4,res = 500,useDingbats = T)






##---- In old fLiver with published fLiver
fLiver_wPublishedData = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver_2n/jan24/liver_liverREFmerged_clean_processed_annotated_noUnknowns_0124.RDS')
fLiver_wPublishedData$finalAnn_broad = fLiver_wPublishedData$annot_jan24

fLiver$finalAnn_broad = fLiver$annot
fLiver$Sex[fLiver$donorID %in% c('Hsb36','Hsb37')] = 'XX'

Idents(fLiver) = fLiver$annot
DotPlot(fLiver,idents = unique(fLiver$annot[grepl('MEP',fLiver$annot)]),group.by = 'Genotype',
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
df$nGene_log10 = log10(abs(df$nGene))
df$nGene_log10[df$direction %in% c('T21_down','MLDS_down')] = -df$nGene_log10[df$direction %in% c('T21_down','MLDS_down')]

plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/Plots'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')

plotFun_MLDS.T21_DEGs_barplot = function(noFrame=FALSE,noPlot=FALSE){
  library(ggbreak) 
  
  p1 = ggplot(df,aes(group,nGene_log10,fill=direction))+
    geom_col(width = 0.58)+
    #scale_fill_manual(values = c(col25[1],col25[2],colAlpha(col25[1:2],alphas = 0.6)))+
    scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
    geom_hline(yintercept = 0,col='black',lwd=0.5)+
    #scale_y_log10()+
    # scale_y_break(c(40,1560),scales = 0.8,expand = T) +
    # scale_y_break(c(-40,-1870),scales = 4) +
    scale_y_continuous(limits = c(-log10(250),log10(600)),
                       breaks = c(-log10(c(200,10)),0,log10(c(10,200,400,600))),
                       labels = as.character(c(-200,10,0,10,200,400,600)))+
    #scale_y_break(c(-40,-1880,40,1550),ticklabels = c(1550,1580,40,20,0,-20,-40,-1893)) +
    geom_hline(yintercept = df$nGene_log10,col='black',lty=2,lwd=0.3)+
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

saveFig(file.path(plotDir,paste0('Fig3A_MLDS.T21.DEGs_barplot')),plotFun_MLDS.T21_DEGs_barplot,rawData=df,width = 3.2,height = 4,res = 500,useDingbats = T)



df = df[df$group == 'T21 vs 2n_MEMP_fLiver',]

df = read.delim(file.path(plotDir,'Fig3A_MLDS.T21.DEGs_barplot_simplified_sameScaleGATA1smodule_rawData.tsv'),sep = '\t')
plotFun_MLDS.T21_DEGs_barplot_onlyT21imprints = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(df,aes(group,nGene,fill=direction))+
    geom_col(width = 0.58)+
    scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
    geom_hline(yintercept = 0,col='black',lwd=0.5)+
    
    # scale_y_break(c(40,1560),scales = 0.8,expand = T) +
    # scale_y_break(c(-40,-1870),scales = 4) +
    #scale_y_break(c(-40,-1880,40,1550),ticklabels = c(1550,1580,40,20,0,-20,-40,-1893)) +
    #scale_y_continuous(limits = c(-20,50))+
    #ylim(-log10(250),log10(600))+
    # scale_y_continuous(limits = c(-log10(250),log10(600)),
    #                    breaks = c(-log10(c(200,10)),0,log10(c(10,200,400,600))),
    #                    labels = as.character(c(-200,10,0,10,200,400,600)))+
    scale_y_continuous(limits = c(-250,600),
                       breaks = c(-200,0,200,400,600),
                       labels = c(-200,0,200,400,600))+
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

#saveFig(file.path(plotDir,paste0('Fig3x_MLDS.T21.DEGs_barplot_simplified')),plotFun_MLDS.T21_DEGs_barplot_onlyT21imprints,rawData=df,width = 2.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,paste0('Fig3A_MLDS.T21.DEGs_barplot_simplified_sameScaleGATA1smodule')),plotFun_MLDS.T21_DEGs_barplot_onlyT21imprints,rawData=df,width = 2.45,height = 5.3,res = 500,useDingbats = T)



