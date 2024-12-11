##----   Deriving TAM vs MLDS differences module, using good TAM only    -----##


outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24'
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
# mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv',row.names = 1)
# all(mlds$annot_aug24 == mdat$annot[match(mlds$cellID,mdat$cellID)])
# # Convert to h5ad for cell2drug (in python)
# devtools::install_github("cellgeni/sceasy")
# library(sceasy)
# use_condaenv('base',required = F,conda = '/opt/conda/etc/profile.d/conda.sh')
# # Converting Raw counts to adata
# convertFormat(mlds, from="seurat", to="anndata",assay = "RNA", main_layer = "counts",
#              outFile='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.h5ad')

##---------------------------------------------------------------##
##   2a. DEGs between good TAM vs all ML-DS BM diagnostic      ####
##---------------------------------------------------------------##
##------- Using pseudobulk -----------
cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in%
                                                c(mlds$cellID[mlds$donorID %in% c('CC1','CC2','L075','L182','L114','CC6','CC8') & mlds$annot == 'Tumour'],
                                                  mlds$cellID[mlds$disease == 'MLDS' & mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & !mlds$donorID %in% c('L041','CC3') & mlds$tissue == 'BM'])]]
rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex','disease','tissue')]
mDat$group = factor(mDat$disease,levels = c('TAM','MLDS'))

tam.vs.all.mlds_res = compareCell_simplified(toc = cnt_mtx,
                                         mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                                         coords = gns[rownames(cnt_mtx)],
                                         cellTypeName='TAM',
                                         formula='~ %s',tgtChrs = paste0('chr',c(1:22)),
                                         donorID='donorID',groupID='group')

saveRDS(tam.vs.all.mlds_res,'DESeq2_goodTAM.vs.all.dMLDS_res_2410.RDS')
tam.vs.all.mlds_res = readRDS('DESeq2_goodTAM.vs.all.dMLDS_res_2410.RDS')

tam.vs.all.mlds_deg = tam.vs.all.mlds_res[['mainDE']]
tam.vs.all.mlds_deg$direction = ifelse(tam.vs.all.mlds_deg$log2FoldChange > 0, 'MLDS_up','MLDS_down')
tam.vs.all.mlds_deg$cellFrac.diff = tam.vs.all.mlds_deg$cellFrac_g1 - tam.vs.all.mlds_deg$cellFrac_g2
tam.vs.all.mlds_deg$cellFrac_max = pmax(tam.vs.all.mlds_deg$cellFrac_g1, tam.vs.all.mlds_deg$cellFrac_g2)


tam.vs.all.mlds_deg = tam.vs.all.mlds_deg[tam.vs.all.mlds_deg$padj < 0.05 & 
                                      abs(tam.vs.all.mlds_deg$log2FoldChange) >= min_l2FC &
                                      abs(tam.vs.all.mlds_deg$cellFrac_max) >= max_cellfrac,]


table(tam.vs.all.mlds_deg$direction)
dim(tam.vs.all.mlds_deg)
tam.vs.all.mlds_deg = tam.vs.all.mlds_deg[order(abs(tam.vs.all.mlds_deg$log2FoldChange),decreasing = T),]


write.csv(tam.vs.all.mlds_deg,'DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
tam.vs.all.mlds_deg = read.csv('DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')


##----  Bar plot of number of DEGs -------
plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots/noCC3'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')

df = tam.vs.all.mlds_deg %>% group_by(direction) %>% summarise(nGene = n())
df$nGene[df$direction == 'MLDS_down'] = -df$nGene[df$direction == 'MLDS_down']
df$nGene_log10 = log10(abs(df$nGene))
df$nGene_log10[df$direction %in% c('MLDS_down')] = -df$nGene_log10[df$direction %in% c('MLDS_down')]

df$group = 'MLDS_vs_canonicalTAM'

df = read.delim(file.path(plotDir,'Fig3C_gTAM.allMLDS_barplot_rawData.tsv'),sep = '\t')
plotFun_TAM.vs.allMLDS_DEGs_barplot = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(df,aes(group,nGene,fill=direction))+
    geom_col(width = 0.58)+
    #scale_fill_manual(values = c(col25[1],col25[2],colAlpha(col25[1:2],alphas = 0.6)))+
    scale_fill_manual(values = c('#1a4a87','#a4282c','#5878a1','#b36d6e'))+
    geom_hline(yintercept = 0,col='black',lwd=0.5)+
    #geom_hline(yintercept = df$nGene,col='black',lty=2)+
    #scale_y_break(c(40,200),scales = 0.8) +
    # scale_y_continuous(limits = c(-log10(300),log10(650)),
    #                    breaks = c(-log10(c(200,10)),0,log10(c(10,200,400,600))),
    #                    labels = as.character(c(-200,10,0,10,200,400,600)))+
    scale_y_continuous(limits = c(-250,600),
                       breaks = c(-200,0,200,400,600),
                       labels = c(-200,0,200,400,600))+
    theme_classic(base_size = 11)+xlab('')+ylab('# DEGs')+
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

saveFig(file.path(plotDir,paste0('Fig3C_gTAM.allMLDS_barplot')),plotFun_TAM.vs.allMLDS_DEGs_barplot,rawData=df,width = 2.7,height = 5.9,res = 500)





##----------------------------------------------------------------------##
##   2b. DEGs between good TAM vs individual ML-DS BM diagnostic      ####
##----------------------------------------------------------------------##

##----    FM: comparing L075 vs each indv_MLDS
gTAM.vs.dMLDS.markers = data.frame()
for(donor in unique(mlds$donorID[mlds$disease == 'MLDS' & mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic'])){
  print(donor)
  mlds$group_tmp = ifelse(mlds$disease == 'MLDS' & mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$tissue == 'BM' & mlds$donorID == donor,'MLDS',
                          ifelse(mlds$donorID %in% c('L075','CC1','CC2') & mlds$annot == 'Tumour','TAM','others'))
  Idents(mlds) = mlds$group_tmp
  
  markers = FindMarkers(mlds,ident.1 = 'TAM',ident.2 = 'MLDS')
  markers$geneSym = rownames(markers)
  markers$ensID = geneMap$ensID[match(markers$geneSym,geneMap$geneSym)]
  rownames(markers) = markers$ensID
  markers = annotateGenes(markers,geneMap = geneMap)
  markers$donorID = donor
  # Remove rubbish genes
  markers = markers[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',markers$geneSym),]
  gTAM.vs.dMLDS.markers = rbind(gTAM.vs.dMLDS.markers,markers)
  
  
  if(donor == 'L076'){
    mlds$group_tmp = ifelse(mlds$disease == 'MLDS' & mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$tissue != 'BM' & mlds$donorID == donor,'MLDS',
                            ifelse(mlds$donorID %in% c('L075','CC1','CC2') & mlds$annot == 'Tumour','TAM','others'))
    Idents(mlds) = mlds$group_tmp
    
    markers = FindMarkers(mlds,ident.1 = 'TAM',ident.2 = 'MLDS')
    markers$geneSym = rownames(markers)
    markers$ensID = geneMap$ensID[match(markers$geneSym,geneMap$geneSym)]
    rownames(markers) = markers$ensID
    markers = annotateGenes(markers,geneMap = geneMap)
    markers$donorID = 'L076_Blood'
    # Remove rubbish genes
    markers = markers[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',markers$geneSym),]
    gTAM.vs.dMLDS.markers = rbind(gTAM.vs.dMLDS.markers,markers)
  }
  
}


gTAM.vs.dMLDS.markers$pct.diff = gTAM.vs.dMLDS.markers$pct.1 - gTAM.vs.dMLDS.markers$pct.2
gTAM.vs.dMLDS.markers$direction = ifelse(gTAM.vs.dMLDS.markers$avg_log2FC > 0 ,'TAM_up','TAM_down')

write.csv(gTAM.vs.dMLDS.markers,'FM_gTAM.vs.indiv_dMLDS_allGene.csv')
# gTAM.vs.dMLDS.markers = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/18_TAM_analysis/MLDS_vs_goodTAM/FM_L075.vs.indiv_MLDS_allGene.csv',row.names = 1)
# markers_summary = gTAM.vs.dMLDS.markers %>% group_by(direction,geneSym) %>% summarise(nDonor_tested=n_distinct(donorID))
# ## Filter the output to keep top genes
# gTAM.vs.dMLDS.markers_deg = gTAM.vs.dMLDS.markers[gTAM.vs.dMLDS.markers$p_val_adj < 0.01 & abs(gTAM.vs.dMLDS.markers$avg_log2FC) >= 0.3 &
#                                                   abs(gTAM.vs.dMLDS.markers$pct.diff) >= 0.2,]
# gTAM.vs.dMLDS.markers_deg = gTAM.vs.dMLDS.markers_deg[order(abs(gTAM.vs.dMLDS.markers_deg$avg_log2FC),decreasing = T),]
# table(gTAM.vs.dMLDS.markers_deg$direction)
# 
# #View(gTAM.vs.dMLDS.markers_deg[gTAM.vs.dMLDS.markers_deg$avg_log2FC < 0,])
# gTAM.vs.dMLDS.markers_summary = gTAM.vs.dMLDS.markers_deg %>% group_by(direction,geneSym) %>% summarise(nDonor = n_distinct(donorID))
# gTAM.vs.dMLDS.markers_summary = merge(gTAM.vs.dMLDS.markers_summary,markers_summary,by=c('direction','geneSym'),all.x=T)
# gTAM.vs.dMLDS.markers_summary$frac = gTAM.vs.dMLDS.markers_summary$nDonor/gTAM.vs.dMLDS.markers_summary$nDonor_tested
# table(gTAM.vs.dMLDS.markers_summary$direction,gTAM.vs.dMLDS.markers_summary$nDonor)
# 
# gTAM.vs.dMLDS.markers_topDEG = rbind(gTAM.vs.dMLDS.markers_deg[gTAM.vs.dMLDS.markers_deg$avg_log2FC > 0 & gTAM.vs.dMLDS.markers_deg$geneSym %in% gTAM.vs.dMLDS.markers_summary$geneSym[gTAM.vs.dMLDS.markers_summary$direction == 'MLDS_up' & gTAM.vs.dMLDS.markers_summary$nDonor >= 5],],
#                                     gTAM.vs.dMLDS.markers_deg[gTAM.vs.dMLDS.markers_deg$avg_log2FC < 0 & gTAM.vs.dMLDS.markers_deg$geneSym %in% gTAM.vs.dMLDS.markers_summary$geneSym[gTAM.vs.dMLDS.markers_summary$direction == 'MLDS_down' & gTAM.vs.dMLDS.markers_summary$nDonor >= 5],])
# 
# write.csv(gTAM.vs.dMLDS.markers_summary[gTAM.vs.dMLDS.markers_summary$nDonor >= 6,],'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/18_TAM_analysis/MLDS_vs_goodTAM/FM_L075.vs.indiv_MLDS_FM_MLDS.vs.goshT21_topMarkers.csv')
# gTAM.vs.dMLDS.markers_summary = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/18_TAM_analysis/MLDS_vs_goodTAM/FM_L075.vs.indiv_MLDS_FM_MLDS.vs.goshT21_topMarkers.csv',row.names = 1)
# 
# 
# ## Plot expression of these genes
# mlds$group = ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID != 'L041',paste0(mlds$disease,':',mlds$donorID),mlds$annot)
# mlds$group[mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076'] = paste0('MLDS:L076:',mlds$tissue[mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076'])
# mlds$group = factor(mlds$group,c(unique(mlds$group[!grepl('TAM:|MLDS:',mlds$group)]),
#                                  unique(mlds$group[grepl('MLDS:',mlds$group) & !grepl('L076',mlds$group)]),
#                                  unique(mlds$group[grepl('MLDS:L076:BM',mlds$group)]),
#                                  unique(mlds$group[grepl('MLDS:L076:Blood',mlds$group)]),
#                                  unique(mlds$group[grepl('TAM:L075',mlds$group)]),
#                                  unique(mlds$group[grepl('TAM:CC2',mlds$group)]),
#                                  unique(mlds$group[grepl('TAM:CC1',mlds$group)]),
#                                  unique(mlds$group[grepl('TAM:L156',mlds$group)])
# ))
# Idents(mlds) = mlds$group
# DotPlot(mlds,idents = unique(Idents(mlds)[grepl('MLDS:|TAM:|HSC_MPP|Mono_CD14|MEP|MK|EE|CMP',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT',Idents(mlds))]),
#         #group.by = 'donorID',
#         scale = T,
#         #features = gTAM.vs.dMLDS.markers_summary$geneSym[gTAM.vs.dMLDS.markers_summary$nDonor == 6 & gTAM.vs.dMLDS.markers_summary$direction == 'L075_down'][1:50]
#         #features = c('SEC31A','NFKBIA','RALGAPA2','KIAA1211')
#         #features = c('IFNGR1','IFNGR2')
#         #features = c('JUN','MAPK8','PIK3CA','FOS','SOS1','NFKB1','TRAF3','ITPR2'),
#         #features = geneOrder$geneSym[geneOrder$group == '3']
#         #features = good.vs.bad.TAM_topDEGs$geneSym[good.vs.bad.TAM_topDEGs$direction == 'good_down']
#         # features = c(gTam.vs.mlds_deg$geneSym[gTam.vs.mlds_deg$is_topDEG == T & gTam.vs.mlds_deg$direction == 'TAM_up'],
#         #              gTam.vs.mlds_deg$geneSym[gTam.vs.mlds_deg$is_topDEG == T & gTam.vs.mlds_deg$direction == 'TAM_down'])
#         features = c(good.vs.bad.TAM_deg$geneSym[good.vs.bad.TAM_deg$is_topDEG == T & good.vs.bad.TAM_deg$direction == 'good_up'],
#                      good.vs.bad.TAM_deg$geneSym[good.vs.bad.TAM_deg$is_topDEG == T & good.vs.bad.TAM_deg$direction == 'good_down']
#         )
#         #features = good.vs.bad.TAM_MLDS_deg$geneSym
#         #features = df$geneSym[df$module == 'MLDS.TAM.up'][1:60]
#         
# )+RotatedAxis()+
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
#         axis.text.y = element_text(size=11)) + xlab('') + ylab('')
# 
# 
# 
# 
##------- Using pseudobulk -----------
cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in%
                                                c(mlds$cellID[mlds$donorID %in% c('CC1','CC2','L075','L182','L114') & mlds$annot == 'Tumour'],
                                                  mlds$cellID[mlds$disease == 'MLDS' & mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & !mlds$donorID %in% c('L041','CC3') & mlds$tissue == 'BM'])]]

rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot','sex','disease','tissue')]
mDat$donorID = as.character(mDat$donorID)
mDat = mDat[match(colnames(cnt_mtx),mDat$cellID),]
table(mDat$disease,mDat$donorID)
table(mDat$disease,mDat$tissue)
mDat$group = ifelse(mDat$disease == 'MLDS',mDat$donorID,mDat$disease)
mDat$group = factor(mDat$group,levels = c('TAM',unique(mDat$group[mDat$group !='TAM'])))

tam.vs.mlds_res = compareCell_simplified(toc = cnt_mtx,
                                         mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                                         coords = gns[rownames(cnt_mtx)],
                                         cellTypeName='TAM',
                                         formula='~ %s',tgtChrs = paste0('chr',c(1:22)),
                                         donorID='donorID',groupID='group')

saveRDS(tam.vs.mlds_res,'DESeq2_goodTAM.vs.indiv_dMLDS_res.RDS')
tam.vs.mlds_res = readRDS('DESeq2_goodTAM.vs.indiv_dMLDS_res.RDS')

# tam.vs.mlds_dds = tam.vs.mlds_res[['dds']]
# res_name = results(tam.vs.mlds_dds,name = 'group_L019_vs_TAM',alpha = 0.05,lfcThreshold=0) %>% as.data.frame()
# res_contrast = results(tam.vs.mlds_dds,contrast = c('group','L019','TAM'),alpha = 0.05,lfcThreshold=0) %>% as.data.frame()
# res_contrast = res_contrast[match(rownames(res_name),rownames(res_contrast)),]
# 
# View(res_name[res_contrast$log2FoldChange != res_name$log2FoldChange,])
# res = results(tam.vs.mlds_dds,contrast = list(c(resultsNames(tam.vs.mlds_dds)[resultsNames(tam.vs.mlds_dds) != 'Intercept'])),
#               alpha = 0.05,lfcThreshold=0)
# res = annotateGenes(res,geneMap = geneMap)
# View(as.data.frame(res))
# tam.vs.mlds_deg = tam.vs.mlds_deg[tam.vs.mlds_deg$padj < 0.05 & abs(tam.vs.mlds_deg$log2FoldChange) >= 0.5 &
#                                     abs(tam.vs.mlds_deg$pct.diff) >= 10,]
# table(tam.vs.mlds_deg$direction)
# dim(tam.vs.mlds_deg)
# tam.vs.mlds_deg = tam.vs.mlds_deg[order(abs(tam.vs.mlds_deg$log2FoldChange),decreasing = T),]
# 
# write.csv(tam.vs.mlds_deg,'DESeq2_goodTAM.vs.dMLDS_topDEGs.csv')
# tam.vs.mlds_deg = read.csv('DESeq2_goodTAM.vs.dMLDS_topDEGs.csv')

res = lfcShrink(dds,coef = contrast,res=res,type = 'ashr')
res_contrast <- lfcShrink(tam.vs.mlds_dds, coef='group_L019_vs_TAM', res=res_contrast)
plotMA(res_contrast)


## old dds
dds = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
dds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.dMLDS_res.RDS')
resultsNames(dds[['dds']])
res_old = results(dds[['dds']],name = 'group_L019_vs_TAM',alpha = 0.05,lfcThreshold = 0) %>% as.data.frame()
res_old = annotateGenes(res_old,geneMap = geneMap)



tam.vs.mlds_res = readRDS('DESeq2_goodTAM.vs.indiv_dMLDS_res.RDS')
## Extracting DEGs when comparing TAM to individual MLDS
tam.vs.mlds_allDEG = tam.vs.mlds_res[['all_de']]
tam.vs.mlds_allDEG$direction = ifelse(tam.vs.mlds_allDEG$log2FoldChange > 0, 'TAM_down','TAM_up')
tam.vs.mlds_allDEG$cellFrac.diff = tam.vs.mlds_allDEG$cellFrac_g1 - tam.vs.mlds_allDEG$cellFrac_g2
tam.vs.mlds_allDEG$cellFrac_max = pmax(tam.vs.mlds_allDEG$cellFrac_g1, tam.vs.mlds_allDEG$cellFrac_g2)
# ## Remove DEGs overlapping with sex
# sex_deg = tam.vs.mlds_allDEG[tam.vs.mlds_allDEG$contrast == 'sex_M_vs_F',]
# #tam.vs.mlds_allDEG = tam.vs.mlds_allDEG[tam.vs.mlds_allDEG$contrast != 'sex_M_vs_F' & !tam.vs.mlds_allDEG$geneSym %in% sex_deg$geneSym,]
# tam.vs.mlds_allDEG = tam.vs.mlds_allDEG[tam.vs.mlds_allDEG$contrast != 'sex_M_vs_F',]

## Only keep genes with abs(log2FC) >= 0.5 & abs(pct.diff) >= 10
tam.vs.mlds_allDEG = tam.vs.mlds_allDEG[abs(tam.vs.mlds_allDEG$log2FoldChange) >= min_l2FC & 
                                          abs(tam.vs.mlds_allDEG$cellFrac_max) >= max_cellfrac,]
                                          #abs(tam.vs.mlds_allDEG$pct.diff) >= 30,]
tam.vs.mlds_allDEG = tam.vs.mlds_allDEG[order(abs(tam.vs.mlds_allDEG$cellFrac.diff),decreasing = T),]

## Tally up number of times each gene was identified as DEG
tam.vs.mlds_allDEG = tam.vs.mlds_allDEG %>% group_by(geneSym,direction) %>% mutate(nContrast = n_distinct(contrast))
deg_summary = tam.vs.mlds_allDEG %>% group_by(geneSym,direction) %>% summarise(nContrast = n_distinct(contrast))



old_fig3a_MLDS_DEGs = function(){
  
  # tam.vs.mlds_allDEG.sub = tam.vs.mlds_allDEG %>% filter(abs(pct.diff) >= 30 & 
  #                                                          abs(log2FoldChange) >= 1) %>%  
  #   group_by(geneSym,direction) %>% mutate(nContrast = n_distinct(contrast))
  
  tam.vs.mlds_allDEG.sub = tam.vs.mlds_allDEG %>%  
    group_by(geneSym,direction) %>% mutate(nContrast = n_distinct(contrast))
  degFiltered_summary = tam.vs.mlds_allDEG.sub %>% group_by(geneSym,direction) %>% summarise(nContrast = n_distinct(contrast))
  table(degFiltered_summary$nContrast,degFiltered_summary$direction)
  tam.vs.mlds_allDEG.sub = tam.vs.mlds_allDEG.sub[order(tam.vs.mlds_allDEG.sub$nContrast,decreasing = T),]
  
  
  
  # Plot frequency of genes being DEGs
  df = as.data.frame(table(degFiltered_summary$nContrast,degFiltered_summary$direction))
  colnames(df) = c('nContrast','direction','freq')
  df$direction = as.character(df$direction)
  df$direction[df$direction == 'TAM_up'] = 'MLDS_down'
  df$direction[df$direction == 'TAM_down'] = 'MLDS_up'
  df$freq[df$direction == 'MLDS_down'] = -df$freq[df$direction == 'MLDS_down']
  df$log10_nGene = ifelse(df$freq > 0,log10(df$freq),-log10(-df$freq))
  library(ggbreak) 
  
  plotFun_Fig.3a_MLDS.DEGs_barplot = function(noFrame=FALSE,noPlot=FALSE){
    p1 = ggplot(df,aes(nContrast,(freq),fill=direction))+
      geom_col(width = 0.85) + 
      geom_hline(yintercept = 0,col='black',lwd=0.1)+
      scale_fill_manual(values = c('#1a4a87','#a4282c'))+
      theme_classic(base_size = 7) + 
      scale_y_break(c(200,600,-300,-1000))+
      #scale_y_break(c(), scales = 3,expand = F)+
      #ylim(-920,550)+
      ylim(-1100,800)+
      xlab('Number of significant comparison')+ylab('# DEGs')+
      theme(axis.line = element_blank(),axis.text = element_text(color='black'))
    
    
    p2 = ggplot(df[is.finite(df$log10_nGene),],aes(nContrast,log10_nGene,fill=direction))+
      geom_col(width = 0.85) + 
      geom_hline(yintercept = 0,col='black',lwd=0.5)+
      scale_fill_manual(values = c('#1a4a87','#a4282c'))+
      theme_classic(base_size = 12) + 
      #scale_y_break(c(200,300,-300,-900))+
      #scale_y_break(c(), scales = 3,expand = F)+
      #ylim(-920,550)+
      #ylim(-1100,800)+
      xlab('Number of significant comparison')+ylab('log10 # DEGs')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),#legend.position = 'bottom',
             legend.title = element_text(size=10),legend.text = element_text(size=8),legend.key.size = unit(0.5,'cm'),
             axis.ticks = element_line(colour='black'),
             axis.text = element_text(color='black'))
    
    print(p1)    
  }
  saveFig(file.path(plotDir,'Fig3a_MLDS.DEG_barplot_v2'),plotFun_Fig.3a_MLDS.DEGs_barplot,rawData=df,width = 4.5,height = 4,res = 500)
  saveFig(file.path(plotDir,'Fig3a_MLDS.DEG_barplot_logScale_v2'),plotFun_Fig.3a_MLDS.DEGs_barplot,rawData=df,width = 4,height = 4,res = 500)
  
  
  ## WITHOUT CC3
  
  ## Version 1: stringent filtering criteria
  # tam.vs.mlds_allDEG.sub = tam.vs.mlds_allDEG %>% filter(contrast != 'group_CC3_vs_TAM' &
  #                                                          abs(pct.diff) >= 30 & 
  #                                                          abs(log2FoldChange) >= 1) %>%  
  #   group_by(geneSym,direction) %>% mutate(nContrast = n_distinct(contrast))
  # 
  
  ## Version 2: less stringent, consistent with the other comparisons
  tam.vs.mlds_allDEG.sub = tam.vs.mlds_allDEG %>% filter(contrast != 'group_CC3_vs_TAM') %>%  
    group_by(geneSym,direction) %>% mutate(nContrast = n_distinct(contrast))
  
  degFiltered_summary = tam.vs.mlds_allDEG.sub %>% group_by(geneSym,direction) %>% summarise(nContrast = n_distinct(contrast))
  table(degFiltered_summary$nContrast,degFiltered_summary$direction)
  tam.vs.mlds_allDEG.sub = tam.vs.mlds_allDEG.sub[order(tam.vs.mlds_allDEG.sub$nContrast,decreasing = T),]
  
  
  
  # Plot frequency of genes being DEGs
  df = as.data.frame(table(degFiltered_summary$nContrast,degFiltered_summary$direction))
  colnames(df) = c('nContrast','direction','freq')
  df$direction = as.character(df$direction)
  df$direction[df$direction == 'TAM_up'] = 'MLDS_down'
  df$direction[df$direction == 'TAM_down'] = 'MLDS_up'
  df$freq[df$direction == 'MLDS_down'] = -df$freq[df$direction == 'MLDS_down']
  library(ggbreak) 
  
  
  plotFun_Fig.3a_MLDS.DEGs_barplot_noCC3 = function(noFrame=FALSE,noPlot=FALSE){
    p1 = ggplot(df,aes(nContrast,(freq),fill=direction))+
      geom_col(width = 0.85) + 
      geom_hline(yintercept = 0,col='black',lwd=0.1)+
      scale_fill_manual(values = c('#1a4a87','#a4282c'))+
      theme_classic(base_size = 7) + 
      scale_y_break(c(100,300))+
      #scale_y_break(c(), scales = 3,expand = F)+
      #ylim(-20,350)+
      ylim(-400,700)+
      xlab('Number of significant comparison')+ylab('# DEGs')+
      theme(axis.line = element_blank(),axis.text = element_text(color='black'))
    
    print(p1)    
  }
  
  saveFig(file.path(plotDir,'Fig3a_MLDS.DEG_barplot_noCC3'),plotFun_Fig.3a_MLDS.DEGs_barplot_noCC3,rawData=df,width = 3.5,height = 3,res = 500)
  
}

## With the new type of plot, showing number of DEGs for each patient, we dont need to remove CC3 anymore  
fig3E_MLDS_DEGs_byDonor = function(){
  tam.vs.mlds_allDEG.sub = tam.vs.mlds_allDEG %>%  
    group_by(geneSym,direction) %>% mutate(nContrast = n_distinct(contrast))
  tam.vs.mlds_allDEG.sub$nComparison = ifelse(tam.vs.mlds_allDEG.sub$nContrast %in% c(1),1,
                                              ifelse(tam.vs.mlds_allDEG.sub$nContrast %in% c(2,3),'2-3','>=4'))
  
  donorLevel_summary = tam.vs.mlds_allDEG.sub %>% group_by(contrast,direction,nComparison) %>% summarise(nGene = n())
  donorLevel_summary$donorID = gsub('group_|_vs_.*$','',donorLevel_summary$contrast)
  donorLevel_summary = donorLevel_summary[,c('donorID','direction','nComparison','nGene')]
  donorLevel_summary$direction[donorLevel_summary$direction == 'TAM_up'] = 'MLDS_down'
  donorLevel_summary$direction[donorLevel_summary$direction == 'TAM_down'] = 'MLDS_up'
  donorLevel_summary$nGene[donorLevel_summary$direction == 'MLDS_down'] = -donorLevel_summary$nGene[donorLevel_summary$direction == 'MLDS_down']
  donorLevel_summary$donorID = factor(donorLevel_summary$donorID,c('CC3','L076','L178','L091','L019','L042','CC5','L038','L040','CC4','L039'))
  donorLevel_summary = donorLevel_summary[order(donorLevel_summary$donorID),]
  
  plotFun_Fig.3C_MLDS.DEGs_byDonor_barplot = function(noFrame=FALSE,noPlot=FALSE){
    donorLevel_summary$fillGroup = paste0(donorLevel_summary$direction,':',donorLevel_summary$nComparison)
    # donorLevel_summary$fillGroup = factor(donorLevel_summary$fillGroup,c('MLDS_up:>=3','MLDS_up:2','MLDS_up:1',
    #                                                                      'MLDS_down:>=3','MLDS_down:2','MLDS_down:1'))
    # 
    donorLevel_summary$fillGroup = factor(donorLevel_summary$fillGroup,c('MLDS_up:>=4','MLDS_up:2-3','MLDS_up:1',
                                                                         'MLDS_down:>=4','MLDS_down:2-3','MLDS_down:1'))
    # donorLevel_summary$log10_nGene = ifelse(donorLevel_summary$nGene > 0, log10(donorLevel_summary$nGene),-log10(-donorLevel_summary$nGene))  
    # 
    # p1 = ggplot(donorLevel_summary,aes(donorID,log10_nGene,fill=fillGroup))+
    #   geom_col(width = 0.85) + 
    #   geom_hline(yintercept = 0,col='black',lwd=0.5)+
    #   scale_fill_manual(values = c('#a4282c','#b36d6e','#d6b0b0','#1a4a87','#5878a1','#a9bad0'))+
    #   theme_classic(base_size = 12) + 
    #   scale_y_continuous(breaks = seq(-8,6))+
    #   #scale_y_break(c(), scales = 3,expand = F)+
    #   #ylim(-20,350)+
    #   #ylim(-1400,500)+
    #   xlab('')+ylab('log10 # DEGs')+
    # theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank(),#legend.position = 'bottom',
    #       legend.title = element_text(size=10),legend.text = element_text(size=8),legend.key.size = unit(0.5,'cm'),
    #       axis.ticks = element_line(colour='black'),
    #       axis.text = element_text(color='black'))
    # 
    # 
    p2 = ggplot(donorLevel_summary,aes(donorID,nGene,fill=fillGroup))+
      geom_col(width = 0.85) + 
      geom_hline(yintercept = 0,col='black',lwd=0.5)+
      scale_fill_manual(values = c('#a4282c','#b36d6e','#d6b0b0','#1a4a87','#5878a1','#a9bad0'))+
      theme_classic(base_size = 12) + 
      #scale_y_continuous(breaks = seq(-8,6))+
      scale_y_continuous(limits = c(-500,600),
                         breaks = c(-c(500,250),0,300,600),
                         labels = as.character(c(-500,-250,0,300,600)))+
      
      #ylim(-20,350)+
      #ylim(-1400,500)+
      xlab('')+ylab('# DEGs')+
      theme_classic(base_size = 11)+xlab('')+ylab('# DEGs')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 1),
            axis.line = element_blank(),#legend.position = 'bottom',
            axis.ticks = element_line(colour = 'black'),
            legend.title = element_text(size=10,colour = 'black'),
            legend.text = element_text(size=8,colour = 'black'),legend.key.size = unit(0.5,'cm'),
            axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text = element_text(color='black'))
    
    print(p2)    
  }
  
  saveFig(file.path(plotDir,'Fig3C_MLDS.DEG_byDonor_barplot'),plotFun_Fig.3C_MLDS.DEGs_byDonor_barplot,rawData=donorLevel_summary,width = 6,height = 6,res = 500)
  
}




# Plot frequency of genes being DEGs
df = as.data.frame(table(degFiltered_summary$nContrast,degFiltered_summary$direction))
colnames(df) = c('nContrast','direction','freq')




##------- Plot expression of these genes ------####

mlds$group = ifelse(mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID != 'L041',paste0(mlds$disease,':',mlds$donorID),mlds$annot)
mlds$group[mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076'] = paste0('MLDS:L076:',mlds$tissue[mlds$annot == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076'])
mlds$group[mlds$annot == 'Tumour' & mlds$timePoint == 'D.Relapse' & mlds$donorID == 'L076'] = 'MLDS:L076:DR'
mlds$group[mlds$annot == 'Tumour' & mlds$timePoint == 'D.Relapse2' & mlds$donorID == 'L076'] = 'MLDS:L076:DR2'
mlds$group[mlds$annot == 'Tumour' & mlds$timePoint == 'TP1' & mlds$donorID == 'L038'] = 'MLDS:L038:TP1'
mlds$group[mlds$annot == 'Tumour' & mlds$timePoint == 'D.Relapse' & mlds$donorID == 'L156'] = 'TAM:L156:DR'
mlds$group = gsub('MLDS:|TAM:','',mlds$group)
mlds$group = factor(mlds$group,c(unique(mlds$group[grepl('CC',mlds$group) & !grepl('L076|L038',mlds$group) & mlds$disease == 'MLDS']),
                                 unique(mlds$group[grepl('L\\d+',mlds$group) & !grepl('L076|L038',mlds$group) & mlds$disease == 'MLDS']),
                                 unique(mlds$group[grepl('L038',mlds$group)]),
                                 unique(mlds$group[grepl('L076:DR',mlds$group)]),
                                 unique(mlds$group[grepl('L076:BM',mlds$group)]),
                                 unique(mlds$group[grepl('L076:Blood',mlds$group)]),
                                 unique(mlds$group[grepl('L075',mlds$group)]),
                                 unique(mlds$group[grepl('CC2',mlds$group)]),
                                 unique(mlds$group[grepl('CC1',mlds$group)]),
                                 unique(mlds$group[grepl('CC6',mlds$group)]),
                                 unique(mlds$group[grepl('CC7',mlds$group)]),
                                 unique(mlds$group[grepl('CC8',mlds$group)]),
                                 unique(mlds$group[grepl('L114',mlds$group)]),
                                 unique(mlds$group[grepl('L182',mlds$group)]),
                                 unique(mlds$group[grepl('L156',mlds$group)]),
                                 unique(mlds$group[mlds$annot != 'Tumour']),
                                 'Tumour'
                                 
))

Idents(mlds) = mlds$group
DotPlot(mlds,
        #idents = unique(Idents(mlds)[grepl('MLDS:|TAM:|HSC_MPP|Mono_CD14|MEP|MK|EE|CMP',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT',Idents(mlds))]),
        idents = unique(Idents(mlds)[grepl('CC\\d+|L\\d+',Idents(mlds)) & 
                                       !grepl('CC3',Idents(mlds)) & 
                                       !grepl('unsure|Tum_MK|MK_WT|Tumour',Idents(mlds))]),
        
        scale = T,

        # features = c('CA1',unique(tam.vs.mlds_allDEG.sub$geneSym[tam.vs.mlds_allDEG.sub$direction == 'TAM_up' & tam.vs.mlds_allDEG.sub$nContrast >= 4])[1:30],
        #              unique(tam.vs.mlds_allDEG.sub$geneSym[tam.vs.mlds_allDEG.sub$direction == 'TAM_down' & tam.vs.mlds_allDEG.sub$nContrast >= 4])[1:50]
        # )
        features = c(tam.vs.all.mlds_deg$geneSym[tam.vs.all.mlds_deg$direction == 'MLDS_up'],
                     tam.vs.all.mlds_deg$geneSym[tam.vs.all.mlds_deg$direction == 'MLDS_down'])

)+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('')



library(decoupleR)
## 1. Import TF-target network from DOROTHEA 
#dorothea_regulon_human <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
net <- get_collectri(organism='human', split_complexes=FALSE)
head(net)

##---- Further select strongest top DEGs  -----####
tam.vs.mlds_allDEG.sub = tam.vs.all.mlds_deg %>% filter(abs(cellFrac.diff) >= 30/100 &
                                                         abs(log2FoldChange) >= 1) %>%
  group_by(geneSym,direction) %>% mutate(nContrast = n_distinct(contrast))
degFiltered_summary = tam.vs.mlds_allDEG.sub %>% group_by(geneSym,direction) %>% summarise(nContrast = n_distinct(contrast))


mlds_markers = as.data.frame(
  rbind(degFiltered_summary[degFiltered_summary$nContrast >= 3 & degFiltered_summary$direction == 'TAM_up',],
        degFiltered_summary[degFiltered_summary$nContrast >= 3 & degFiltered_summary$direction == 'TAM_down',]))
                     

# mlds.vs.L075.markers_summary[mlds.vs.L075.markers_summary$nDonor >= 7 & mlds.vs.L075.markers_summary$direction == 'L075_down' &
#                                                     mlds.vs.L075.markers_summary$geneSym %in% c('BST2','CA1','CDKN2C','FCGRT','IFI6','TBC1D32'),],
#                      mlds.vs.L075.markers_summary[mlds.vs.L075.markers_summary$nDonor >= 7 & mlds.vs.L075.markers_summary$direction == 'L075_down' &
#                                                     !mlds.vs.L075.markers_summary$geneSym %in% c('BST2','CA1','CDKN2C','FCGRT','IFI6','TBC1D32'),]
# )
mlds_markers$ensID = geneMap$ensID[match(mlds_markers$geneSym,geneMap$geneSym)]
rownames(mlds_markers) = mlds_markers$ensID
mlds_markers = annotateGenes(mlds_markers,geneMap = geneMap)
mlds_markers$isTF[mlds_markers$isTF == F] = (mlds_markers$geneSym[mlds_markers$isTF == F] %in% net$source)

#write.csv(mlds_markers,'DESeq2_goodTAM.vs.individual.dMLDS_topDEGs_noCC3.csv')
write.csv(mlds_markers,'DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')



# OLD- used to derive gene module from individual comparisons
# mlds_markers = read.csv('DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
# mlds_markers = mlds_markers[order(mlds_markers$nContrast,decreasing = T),]
# mlds_markers = rbind(mlds_markers[mlds_markers$direction == 'TAM_up',],
#                      mlds_markers[mlds_markers$direction == 'TAM_down',])



mlds_markers = read.csv('DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')
mlds_markers = mlds_markers[abs(mlds_markers$cellFrac.diff) >= 0.2,]
mlds_markers = mlds_markers[order(abs(mlds_markers$log2FoldChange),decreasing = T),]
mlds_markers = rbind(mlds_markers[mlds_markers$direction == 'MLDS_down',],
                     mlds_markers[mlds_markers$direction == 'MLDS_up',])


plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots/'
fig3a_MLDSmarkers_dotplot = function(){
  
  
  plotFun_MLDS.module_dotplot=function(noFrame=FALSE,noPlot=FALSE){
    p1 = DotPlot(mlds,idents = unique(Idents(mlds)[grepl('CC\\d+|L\\d+',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT|Tumour',Idents(mlds)) & !grepl('L041',Idents(mlds))]),
                 cols = c(colAlpha(grey(0.98),0.85),'black'),
                 features = c(mlds_markers$geneSym)
    )+
      RotatedAxis()+
      theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,face=ifelse(mlds_markers$chr == 'chr21','bold','plain'),
                                       colour =ifelse(mlds_markers$isTF==T,'purple',ifelse(mlds_markers$isCosmic==T,'darkgreen',
                                                                                           ifelse(mlds_markers$isCSM == T,'darkblue','black')))),
            axis.text.y = element_text(size=11),
            legend.position = 'top',legend.text = element_text(size=9),legend.title = element_text(size=8.5),legend.key.size = unit(0.6,'cm')) + xlab('') + ylab('') 
    
    print(p1)  
  }
  #saveFig(file.path(plotDir,'FigSupXX_MLDS.module_dotplot_noCC3'),plotFun_MLDS.module_dotplot,rawData=mlds_markers,width = 9,height = 6,res = 500)
  saveFig(file.path(plotDir,'FigSupXX_MLDS.module_dotplot_v2'),plotFun_MLDS.module_dotplot,rawData=mlds_markers,width = 12,height = 6,res = 500)
  
  
  
  plotFun_mainFig.3_MLDS.module_dotplot_vertical=function(noFrame=FALSE,noPlot=FALSE){
    p1 = DotPlot(mlds,idents = unique(Idents(mlds)[grepl('MLDS:|TAM:',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT',Idents(mlds))]),
                 cols = c(colAlpha(grey(0.98),0.85),'black'),
                 features = c(mlds_markers$geneSym[mlds_markers$direction == 'TAM_down'],
                              mlds_markers$geneSym[mlds_markers$direction == 'TAM_up'])
    )+
      RotatedAxis()+coord_flip()+
      theme(axis.text.y = element_text(size=9,face=ifelse(mlds_markers$chr == 'chr21','bold','plain'),
                                       colour =ifelse(mlds_markers$isTF==T,'purple',ifelse(mlds_markers$isCosmic==T,'darkgreen',
                                                                                           ifelse(mlds_markers$isCSM == T,'darkblue','black')))),
            axis.text.x = element_text(size=11,angle = 90,vjust = 0.5,hjust = 1),
            legend.position = 'right',legend.text = element_text(size=9),legend.title = element_text(size=8.5),legend.key.size = unit(0.6,'cm')) + xlab('') + ylab('') 
    
    
    print(p1)  
  }
  saveFig(file.path(plotDir,'Fig3a_MLDS.module_dotplot_vertical'),plotFun_mainFig.3_MLDS.module_dotplot_vertical,rawData=df,width = 6.5,height = 10,res = 500)
  
  
  
  
  
  plotFun_Sup.Fig.3_MLDS_module_dotplot_combinedLeuk=function(noFrame=FALSE,noPlot=FALSE){
    p1 = DotPlot(big.srat,
                 #idents = unique(big.srat$group_4_dotPlot[grepl('^Tumour.*Diag|Tumour.*D0|Tumour.*L038|Tumour.*L076',big.srat$group_4_dotPlot) & big.srat$donorID !='CC3']),
                 idents = unique(big.srat$group_4_dotPlot[!grepl('other',big.srat$group_4_dotPlot)]),
                 cols = c(colAlpha(grey(0.98),0.85),'black'),
                 features = mlds_markers$geneSym
                 
    )+
      
      RotatedAxis()+
      theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,face=ifelse(mlds_markers$chr == 'chr21','bold','plain'),
                                       colour =ifelse(mlds_markers$isTF==T,'purple',ifelse(mlds_markers$isCosmic==T,'darkgreen',
                                                                                           ifelse(mlds_markers$isCSM == T,'darkblue','black')))),
            axis.text.y = element_text(size=11),
            legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 
    
    
    print(p1)  
  }
  saveFig(file.path(plotDir,'FigSupXX_MLDSmodule_dotplot_combinedLeuk'),plotFun_Sup.Fig.3_MLDS_module_dotplot_combinedLeuk,rawData=df,width = 12,height = 7.7,res = 500)
  saveFig(file.path(plotDir,'FigSupXX_MLDSmodule_dotplot_combinedLeuk_v2'),plotFun_Sup.Fig.3_MLDS_module_dotplot_combinedLeuk,rawData=df,width = 14.7,height = 7.7,res = 500)
}





##------- In big.srat
tam.vs.all.mlds_deg = read.csv('DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')
tam.vs.all.mlds_deg = tam.vs.all.mlds_deg[abs(tam.vs.all.mlds_deg$cellFrac.diff) >= 0.2,]
mlds_markers = rbind(tam.vs.all.mlds_deg[tam.vs.all.mlds_deg$direction == 'MLDS_down',],
                     tam.vs.all.mlds_deg[tam.vs.all.mlds_deg$direction == 'MLDS_up',])

mlds_srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS'
otherLeuk_srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2408.RDS'
mlds = readRDS(mlds_srat_fp)
mlds$finalAnn_broad = mlds$annot_aug24
DimPlot(mlds,cells.highlight = mlds$cellID[mlds$finalAnn_broad != mlds$annot_aug24 & mlds$finalAnn_broad == 'LE'])
## Plot expression of these genes in ALL leukaemia datasets
big.srat = merge_seurat_objects(mlds,otherLeuk,keepAllGenes = T,genomeVersions = c('v38','v38'))
big.srat = standard_clustering(big.srat)
big.srat$finalAnn_broad = big.srat$annot_aug24
big.srat$broadLineage[!is.na(big.srat$GATA1_status) & big.srat$GATA1_status == 'GATA1s_mutant' & big.srat$broadLineage != 'Tumour'] = 'Tumour'

saveRDS(big.srat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/combinedLeuk_2408.RDS')
big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/combinedLeuk_2408.RDS')





#big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/combinedLeuk_2404.RDS')

## Dot plot of MLDS-sepcific genes in other leukaemia
big.srat$timePoint[big.srat$timePoint == 'Diagnostic'] = 'D'
big.srat$timePoint[big.srat$timePoint %in% c('D.Relapse','RelapseD0')] = 'RD'
big.srat$timePoint[big.srat$timePoint == 'D.Relapse2'] = 'R2D'
big.srat$disease[big.srat$disease == 'pAML'] = 'AML'
big.srat$disease[big.srat$disease == 'pBALL'] = 'B-ALL'
big.srat$tissue[is.na(big.srat$tissue) & big.srat$disease == 'B-ALL'] = 'BM'
big.srat$tissue[is.na(big.srat$tissue) & big.srat$disease == 'LPD'] = 'Brain'
big.srat$tissue[big.srat$tissue %in% c('Blood','PBMC')] = 'PB'
big.srat$broadLineage[is.na(big.srat$broadLineage) & big.srat$finalAnn_broad == 'Tumour'] = 'Tumour'
big.srat$broadLineage[is.na(big.srat$broadLineage) & big.srat$finalAnn_broad != 'Tumour'] = 'others'
big.srat$group_4_dotPlot = ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID != 'L041' & big.srat$timePoint %in% c('D0','Diagnostic','RelapseD0','D','RD','R2D') & !big.srat$disease %in% c('TAM','MLDS'), 
                                  paste0('Tumour.',big.srat$disease,'.',big.srat$timePoint,'.',big.srat$tissue),
                                  ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID != 'L041' & big.srat$timePoint %in% c('D0','Diagnostic','RelapseD0','D','RD','R2D') & big.srat$disease %in% c('TAM','MLDS'),
                                         paste0('Tumour.',big.srat$disease,'.',big.srat$tissue,'.',big.srat$donorID,'.',big.srat$timePoint),
                                         ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID == 'L038' & big.srat$timePoint %in% c('TP1') & big.srat$disease %in% c('MLDS'),
                                                'Tumour.MLDS.BM.L038.TP1',
                                                # ifelse(big.srat$broadLineage == 'Tumour' & big.srat$disease == 'Lymphoma','Tumour.Lymphoma.D0.L062',
                                                #        ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID == 'L076' & big.srat$timePoint %in% c('D.Relapse') & big.srat$disease %in% c('MLDS'),'Tumour.MLDS.DR.BM.L076',
                                                              'others')))#))
big.srat$group_4_dotPlot[grepl('^Tumour\\.B-ALL\\.D',big.srat$group_4_dotPlot) & big.srat$Genotype == 'T21'] = gsub('Tumour.B-ALL.D','Tumour.B-ALL.T21.D',big.srat$group_4_dotPlot[grepl('^Tumour\\.B-ALL\\.D',big.srat$group_4_dotPlot) & big.srat$Genotype == 'T21'])
big.srat$group_4_dotPlot[big.srat$group_4_dotPlot == 'Tumour.LPD.D.Brain'] = 'Tumour.LPD.T21.D.Brain'

big.srat$group_4_dotPlot = gsub('Tumour\\.|BM\\.','',big.srat$group_4_dotPlot)
big.srat$group_4_dotPlot[big.srat$finalAnn_broad %in% c('MK','MEP')] = big.srat$finalAnn_broad[big.srat$finalAnn_broad %in% c('MK','MEP')]

nCells = table(big.srat$group_4_dotPlot)
big.srat$group_4_dotPlot = paste0(big.srat$group_4_dotPlot,' (n=',nCells[match(big.srat$group_4_dotPlot,names(nCells))],')')
big.srat$group_4_dotPlot = gsub('TAM\\.|MLDS\\.|\\.D|\\.BM','',big.srat$group_4_dotPlot)
big.srat$group_4_dotPlot = gsub('PB\\.|\\.PB','*',big.srat$group_4_dotPlot)

big.srat$group_4_dotPlot = factor(big.srat$group_4_dotPlot,levels = c(unique(big.srat$group_4_dotPlot[grepl('.*L075',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*CC1',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*CC2',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*L114',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*L182',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*CC6',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*CC7',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*CC8',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*L156',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*L076 ',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*L076\\.RD',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*L076\\.R2D',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*L038 ',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('.*L038\\.TP1',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('',big.srat$group_4_dotPlot) & !big.srat$donorID %in% c('L076','L038') & big.srat$disease %in% c('MLDS') & !grepl('others|MEP|MK',big.srat$group_4_dotPlot)]),
                                                                      
                                                                      unique(big.srat$group_4_dotPlot[grepl('MDS',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('AMKL',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('AML',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('B-ALL\\.T21',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('B-ALL |B-ALL\\* ',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('B-ALL.RD',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('iALL',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('LPD',big.srat$group_4_dotPlot)]),
                                                                      
                                                                      unique(big.srat$group_4_dotPlot[grepl('MEP',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('^MK',big.srat$group_4_dotPlot)]),
                                                                      
                                                                      unique(big.srat$group_4_dotPlot[grepl('others',big.srat$group_4_dotPlot)])
))
table(is.na(big.srat$group_4_dotPlot))


Idents(big.srat) = big.srat$group_4_dotPlot



driver = read.csv('~/lustre_mt22/generalResources/MLDS_drivers.csv')
library(decoupleR)
net = get_collectri(organism='human', split_complexes=FALSE)
head(net)
myc=net[net$source == 'MYC',]

DotPlot(big.srat,#idents = unique(big.srat$group_4_dotPlot[grepl('^Tumour.*RD|Tumour.*D',big.srat$group_4_dotPlot)]),
        idents = unique(big.srat$group_4_dotPlot[!grepl('other',big.srat$group_4_dotPlot)]),
        #features = mlds_markers$geneSym)+
        #features = unique(c(driver$Gene[order(driver$Gene)],'CD81')))+
        features = myc$target[myc$mor == -1][1:100])+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,
                                   colour=ifelse(mlds_markers$isTF[match(mlds_markers$geneSym,mlds_markers$geneSym)],'red',
                                                 ifelse(mlds_markers$isCSM[match(mlds_markers$geneSym,mlds_markers$geneSym)],'blue','black'))),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 



mlds_markers = mlds_vs_goodTAM
mlds_markers = rbind(mlds_markers[mlds_markers$direction == 'MLDS_down',],
                     mlds_markers[mlds_markers$direction == 'MLDS_up',])
mlds_markers = mlds_markers[abs(mlds_markers$cellFrac.diff) > 0.1 & abs(mlds_markers$log2FoldChange) >= 1,]
table(mlds_markers$direction)
#mlds_markers = gata1s_module[gata1s_module$tam_vs_mempT21_group == 'TAM.MLDS.up' & gata1s_module$isTF == T,]
plotFun_Sup.Fig.3_MLDSmarkers_dotplot=function(noFrame=FALSE,noPlot=FALSE){
  p1 = DotPlot(big.srat,
               #idents = unique(big.srat$group_4_dotPlot[grepl('^Tumour.*Diag|Tumour.*D0|MEP',big.srat$group_4_dotPlot)]),
               idents = unique(big.srat$group_4_dotPlot[!grepl('others|CC3',big.srat$group_4_dotPlot)]),
               cols = c(colAlpha(grey(0.98),0.85),'black'),
               features = mlds_markers$geneSym
               
  )+
    
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,face=ifelse(mlds_markers$isTF,'bold','plain'),
                                     colour =ifelse(mlds_markers$chr == 'chr21','purple','black')),
          axis.text.y = element_text(size=11),
          legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 
  
  
  print(p1)  
}
saveFig(file.path(plotDir,'Sup.Fig3x_MLDSmarkers_dotplot'),plotFun_Sup.Fig.3_MLDSmarkers_dotplot,rawData=df,width = 12,height = 8,res = 500,useDingbats = F)



##-----------------------------------------------------------------##
##  3. Expression of MLDS-specific gene module in bRNAseq data   ####
##-----------------------------------------------------------------##
source('~/lustre_mt22/Aneuploidy/scripts/finalScripts/xx01_moduleScore_helperFunctions.R')
library(edgeR)
library(readxl)
#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

##----- 1. Import bulk counts -----##
bulkSamples = import_HenningsBulkSamples(gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T,
                                         tpm_fp = c('~/lustre_mt22/Down_Leukemia/Down_Leukemia_Klusmann_StarSalmon_tpm.csv'))

## extract CPM count
cpmCnt = bulkSamples[['cpm_count']]
tpmCnt = bulkSamples[['tpm_count']]
bulk_mdat = bulkSamples[['mdat']]
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


bulk_mdat$sampleGroup = as.character(bulk_mdat$sampleGroup)
bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TAM' & bulk_mdat$Event == 1] = 'TAM:1'
bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TAM' & bulk_mdat$Event == 0] = 'TAM:0'
bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'MLL_rearrangement'] = 'MLL'
bulk_mdat$sampleGroup = factor(bulk_mdat$sampleGroup,c('TAM','TAM:0','TAM:1','MLDS','AMKL','MLL','MLL_rearrangement',
                                                       'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))


##----- 2. boxplot of log10CPM count for gene of interests ----##
tam.vs.all.mlds_deg
gene_toPlot = mlds_markers
gene_toPlot = tam.vs.all.mlds_deg[tam.vs.all.mlds_deg$direction == 'MLDS_up',]
mtx = tpmCnt[rownames(tpmCnt) %in% gene_toPlot$ensID[gene_toPlot$direction == 'MLDS_up'],]
mtx = mtx[,colnames(mtx) != 'ensID']
rownames(mtx) = gene_toPlot$geneSym[match(rownames(mtx),gene_toPlot$ensID)]
mtx = mtx[match(gene_toPlot$geneSym,rownames(mtx)),]
mtx = as.data.frame(t(as.matrix(mtx)))

mtx$sampleID = rownames(mtx)
mtx$group = bulk_mdat$sampleGroup[match(mtx$sampleID,bulk_mdat$Sample)]
mtx$group[mtx$group == 'MLL_rearrangement'] = 'MLL'

mtx = pivot_longer(mtx,cols = c(1:(which(colnames(mtx) == 'sampleID')-1)),names_to = 'gene',values_to = 'CPM')
mtx = mtx[mtx$group !='TAM',]
mtx = mtx[!is.na(mtx$CPM),]

mtx2 = mtx[mtx$group %in% c('TAM:0','TAM:1'),]
mtx2 = mtx2 %>% group_by(group,gene) %>% summarise(med = median(CPM))
mtx2 = pivot_wider(mtx2,id_cols = 'gene',names_from = 'group',values_from = 'med')
mtx2$med_diff = mtx2$`TAM:1` - mtx2$`TAM:0`



library(ggbeeswarm)
mtx$gene = factor(mtx$gene,gene_toPlot$geneSym)
ggplot(mtx[mtx$gene %in% c('CLEC2B','STXBP6', 'ANGPT1', 'CA1', 'PRAME','CD81','CD44','CD74','CPA3','ARHGEF12','CDKN2C','FCGRT') & 
             mtx$group %in% c('TAM:0','TAM:1','MLDS'),],aes(group,log(CPM)))+
  geom_boxplot(aes(fill=group),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.3)+
  geom_point(size=0.5,alpha=0.7,col=grey(0))+
  scale_fill_manual(values = c(rep(col25[4],1),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),2),rep(grey(0.7),5))) +
  facet_wrap(vars(gene),nrow=2,scales = 'free_y')+
  theme_classic()+
  #ggtitle(title)+xlab('')+ylab('Module score')+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
        strip.background=element_rect(linewidth=0),
        axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))


mtx = mtx[mtx$gene %in% c('CDKN2C', 'CA1','CD81','PTK2','TMX4','MSRB3','ZHX1','ZNF141','PRAME','CD44','CD74','CPA3','ARHGEF12','MAGI1','PGM1','RRAS2','DSTN','LNCAROD','SLC5A3','BST2','ZNF345','LINS1'),]
mtx$gene = factor(mtx$gene,c('CDKN2C', 'CA1','CD81','PTK2','TMX4','MSRB3','ZHX1','ZNF141','PRAME','CD44','CD74','CPA3','ARHGEF12','MAGI1','PGM1','RRAS2','DSTN','LNCAROD','SLC5A3','BST2','ZNF345','LINS1'))
ggplot(mtx[mtx$gene %in% mtx2$gene[mtx2$med_diff > 0] & 
  mtx$group %in% c('TAM:0','TAM:1','MLDS'),],aes(group,log(CPM)))+
  geom_boxplot(aes(fill=group),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.3)+
  geom_point(size=0.5,alpha=0.7,col=grey(0))+
  scale_fill_manual(values = c('white','#c7065a',col25[4],rep(colAlpha(col25[1],0.4),2),rep(grey(0.7),5))) +
  facet_wrap(vars(gene),nrow=2,scales = 'free_y')+
  theme_classic()+
  xlab('')+
  #ggtitle(title)+xlab('')+ylab('Module score')+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
        strip.background=element_rect(linewidth=0),
        axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))+
  ggtitle('Bulk RNA-seq dataset')


##--- Heatmap of CPM in bulk data -----##
mtx = cpmCnt[rownames(cpmCnt) %in% gene_toPlot$ensID,]
rownames(mtx) = gene_toPlot$geneSym[match(rownames(mtx),gene_toPlot$ensID)]
mtx = mtx[match(gene_toPlot$geneSym,rownames(mtx)),]
mtx = mtx[,rownames(mtx) ]
botAnno = HeatmapAnnotation(group = bulk_mdat$sampleGroup[match(colnames(mtx),bulk_mdat$Sample)],
                            col=list(group = c('TAM'=grey(0.7),'TAM:0'=col25[4],'TAM:1'='#c7065a','MLDS'=col25[4],
                                                 'AMKL' = col25[1],'MLL'=col25[1],
                                                 'HSC'=grey(0.7),'MPP'=grey(0.7),'MEP'=grey(0.7),'CMP'=grey(0.7),'GMP'=grey(0.7))))

mtx = log10(mtx)
mtx[!is.finite(mtx)] = NA
mtx = mtx[,colnames(mtx) %in% bulk_mdat$Sample[bulk_mdat$sampleGroup %in% c('MLDS','TAM:0','TAM:1')]]
Heatmap(t(scale(t(mtx))),show_column_names = T,show_row_names = T,cluster_rows = T,cluster_columns = T,
        row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=8),show_row_dend = F,show_column_dend = F,
        bottom_annotation = botAnno,
        split=gene_toPlot$direction[match(rownames(mtx),gene_toPlot$geneSym)],
        column_split = bulk_mdat$sampleGroup[match(colnames(mtx),bulk_mdat$Sample)])


##----------------------------------------------------##
##  6. Functional assessment of the GATA1s module   ####
##----------------------------------------------------##
mlds_markers = read.csv('DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')

##----    EnrichR   -----## --> totally NOT informative
upDEG_enriched <- enrichr(tam.vs.all.mlds_deg$geneSym[tam.vs.all.mlds_deg$direction == 'MLDS_up'], dbs)
View(upDEG_enriched[[i]][upDEG_enriched[[i]]$Adjusted.P.value < 0.05,])
View(upDEG_enriched[[5]][upDEG_enriched[[5]]$Adjusted.P.value < 0.05,])
i=3
plotEnrich(upDEG_enriched[[i]][upDEG_enriched[[i]]$Adjusted.P.value < 0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")

downDEG_enriched <- enrichr(unique(mlds_markers$geneSym[mlds_markers$direction == 'TAM_up']), dbs)
View(downDEG_enriched[[6]][downDEG_enriched[[6]]$Adjusted.P.value < 0.05,])
plotEnrich(downDEG_enriched[[1]][downDEG_enriched[[1]]$Adjusted.P.value < 0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")






##---------------------------------------------------------------------##
##  7. heatmap of average expression of TF in TAM / ML-DS clusters   ####
##---------------------------------------------------------------------##
mlds$group = ifelse(mlds$annot == 'Tumour',paste0(mlds$disease,'_',mlds$donorID,'_',mlds$timePoint),mlds$annot)

avgExpr = AverageExpression(mlds,group.by = 'group')
avgExpr_df = as.data.frame(avgExpr[['RNA']])
avgExpr_df = avgExpr_df[rowSums(avgExpr_df) > 0,]
avgExpr_df$geneSym = rownames(avgExpr_df)
avgExpr_df$ensID = geneMap$ensID[match(avgExpr_df$geneSym,geneMap$geneSym)]
avgExpr_df = annotateGenes(avgExpr_df,geneMap = geneMap)

library(ComplexHeatmap)
mtx = avgExpr_df[avgExpr_df$isCSM==T,colnames(avgExpr_df) %in% colnames(avgExpr[['RNA']])]

hm = Heatmap(t(scale(t(mtx))),show_row_dend = F,show_column_dend = F,
        km=10,
        show_row_names = F,show_column_names = T,column_split = mlds$broadLineage[match(colnames(mtx),mlds$group)])
ht = draw(hm)
genes = rownames(mtx)[row_order(ht)[['8']]]
length(genes)
View(mtx[genes,])


DotPlot(mlds,group.by = 'group',features = tam.vs.all.mlds_deg$geneSym[tam.vs.all.mlds_deg$direction == 'MLDS_down']) + 
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 

DotPlot(mlds,group.by = 'group',features = c('RUNX1','GATA1','MORC3')) + 
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 

plotCounts(dds, gene=geneMap$ensID[geneMap$geneSym == 'MORC3'], intgroup="dataset")
Idents(fLiver) = fLiver$annot
fLiver$group_tmp = paste0(fLiver$Genotype,'_',fLiver$donorID)
DotPlot(fLiver,idents = 'MEMP_MEP',group.by = 'group_tmp',features = c('RUNX1','GATA1','MORC3','DACH1')) + 
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 

DimPlot(mlds,group.by = 'donorID',label = T,repel = T,label.box = T,cols = col25) + NoLegend()
