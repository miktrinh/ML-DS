##----   Deriving badTAM module, using good TAM only    -----##


outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/goodTAM_vs_badTAM'
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



##------------------------------------------------##
#### 1. Import relevant scRNA datasets: MLDS    ####
##-----------------------------------------------##
## Import MLDS object
mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
mlds$broadLineage = mdat$broadLineage[match(mlds$cellID,mdat$cellID)]
mlds$timePoint[mlds$orig.ident %in% c('MY.200531.14635833')] = 'D.Relapse'
mlds$finalAnn_broad = as.character(mlds$annot_mar24)






##----------------------------------------------##
##   2. DEGs between good TAM vs bad TAM      ####
##----------------------------------------------##
##------- Using pseudobulk -----------
cnt_mtx = mlds@assays$RNA@counts[,mlds$cellID[mlds$cellID %in%
                                                c(mlds$cellID[mlds$donorID %in% c('CC1','CC2','L075') & mlds$annot_mar24 == 'Tumour'],
                                                  mlds$cellID[mlds$donorID %in% c('L156') & mlds$annot_mar24 == 'Tumour'])]]
rownames(cnt_mtx) = geneMap$ensID[match(rownames(cnt_mtx),geneMap$geneSym)]

mDat = mlds@meta.data[mlds$cellID %in% colnames(cnt_mtx),c('cellID','donorID','annot_mar24','sex','disease','tissue')]
mDat$donorID = as.character(mDat$donorID)
mDat = mDat[match(colnames(cnt_mtx),mDat$cellID),]
table(mDat$sex,mDat$donorID)
table(mDat$disease,mDat$tissue)
mDat$group = ifelse(mDat$donorID == 'L156','bad','good')
mDat$group = factor(mDat$group,levels = c('good','bad'))

bad.vs.good_res = compareCell_simplified(toc = cnt_mtx,
                                         mDat = mDat[match(colnames(cnt_mtx),rownames(mDat)),],
                                         coords = gns[rownames(cnt_mtx)],
                                         cellTypeName='TAM',
                                         formula='~ %s + sex',tgtChrs = paste0('chr',c(1:22)),
                                         donorID='donorID',groupID='group')

saveRDS(bad.vs.good_res,'DESeq2_good.vs.bad.TAM_res.RDS')
bad.vs.good_res = readRDS('DESeq2_good.vs.bad.TAM_res.RDS')

bad.vs.good_dds = bad.vs.good_res[['dds']]
bad.vs.good_deg = bad.vs.good_res[['mainDE']]
bad.vs.good_deg$direction = ifelse(bad.vs.good_deg$log2FoldChange > 0,'bad_up','bad_down')
bad.vs.good_deg$pct.diff = bad.vs.good_deg$cellFrac_g1 - bad.vs.good_deg$cellFrac_g2
bad.vs.good_deg = bad.vs.good_deg[abs(bad.vs.good_deg$log2FoldChange) >= 1 & abs(bad.vs.good_deg$pct.diff) >= 20,]
table(bad.vs.good_deg$direction)
dim(bad.vs.good_deg)
bad.vs.good_deg = bad.vs.good_deg[order(abs(bad.vs.good_deg$pct.diff),decreasing = T),]

write.csv(bad.vs.good_deg,'DESeq2_good.vs.bad.TAM_topDEGs.csv')
bad.vs.good_deg = read.csv('DESeq2_good.vs.bad.TAM_topDEGs.csv')


##---------------------------------------##
##    Compared to gTAM_vs_MLDS DEGs    ####
##---------------------------------------##
mlds_markers = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
mlds_markers = deg_summary %>% as.data.frame()
mlds_markers$ensID = geneMap$ensID[match(mlds_markers$geneSym,geneMap$geneSym)]
mlds_markers_up = mlds_markers[mlds_markers$direction == 'TAM_down',]
rownames(mlds_markers_up) = mlds_markers_up$ensID
mlds_markers_up = annotateGenes(mlds_markers_up,geneMap = geneMap)

mlds_markers_down = mlds_markers[mlds_markers$direction == 'TAM_up',]
rownames(mlds_markers_down) = mlds_markers_down$ensID
mlds_markers_down = annotateGenes(mlds_markers_down,geneMap = geneMap)

mlds_markers = rbind(mlds_markers_up,mlds_markers_down)


mlds_markers$in_badTAM = ifelse(mlds_markers$direction == 'TAM_up' & mlds_markers$geneSym %in% bad.vs.good_deg$geneSym[bad.vs.good_deg$direction == 'bad_down'],T,
                                ifelse(mlds_markers$direction == 'TAM_down' & mlds_markers$geneSym %in% bad.vs.good_deg$geneSym[bad.vs.good_deg$direction == 'bad_up'],T,F))

table(mlds_markers$in_badTAM,mlds_markers$direction)
tam.vs.mlds_allDEG.sub2 = rbind(tam.vs.mlds_allDEG.sub[tam.vs.mlds_allDEG.sub$direction == 'TAM_down' & tam.vs.mlds_allDEG.sub$geneSym %in% bad.vs.good_deg$geneSym[bad.vs.good_deg$direction == 'bad_up'],],
                                tam.vs.mlds_allDEG.sub[tam.vs.mlds_allDEG.sub$direction == 'TAM_up' & tam.vs.mlds_allDEG.sub$geneSym %in% bad.vs.good_deg$geneSym[bad.vs.good_deg$direction == 'bad_down'],])

degFiltered_summary2 = tam.vs.mlds_allDEG.sub2 %>% group_by(geneSym,direction) %>% summarise(nContrast = n_distinct(contrast))
table(degFiltered_summary2$direction)



##------- Plot expression of these genes ------##
mlds$group = ifelse(mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID != 'L041',paste0(mlds$disease,':',mlds$donorID),mlds$annot_mar24)
mlds$group[mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076'] = paste0('MLDS:L076:',mlds$tissue[mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'Diagnostic' & mlds$donorID == 'L076'])
mlds$group[mlds$annot_mar24 == 'Tumour' & mlds$timePoint == 'D.Relapse' & mlds$donorID == 'L076'] = 'MLDS:L076:DR'
mlds$group = factor(mlds$group,c(unique(mlds$group[grepl('MLDS:CC',mlds$group) & !grepl('L076',mlds$group)]),
                                 unique(mlds$group[grepl('MLDS:L',mlds$group) & !grepl('L076',mlds$group)]),
                                 unique(mlds$group[grepl('MLDS:L076:DR',mlds$group)]),
                                 unique(mlds$group[grepl('MLDS:L076:BM',mlds$group)]),
                                 unique(mlds$group[grepl('MLDS:L076:Blood',mlds$group)]),
                                 unique(mlds$group[grepl('TAM:L075',mlds$group)]),
                                 unique(mlds$group[grepl('TAM:CC2',mlds$group)]),
                                 unique(mlds$group[grepl('TAM:CC1',mlds$group)]),
                                 unique(mlds$group[grepl('TAM:CC6',mlds$group)]),
                                 unique(mlds$group[grepl('TAM:CC7',mlds$group)]),
                                 unique(mlds$group[grepl('TAM:CC8',mlds$group)]),
                                 unique(mlds$group[grepl('TAM:L156',mlds$group)]),
                                 unique(mlds$group[!grepl('TAM:|MLDS:',mlds$group)])
                                 
))
Idents(mlds) = mlds$group
DotPlot(mlds,
        #idents = unique(Idents(mlds)[grepl('MLDS:|TAM:|HSC_MPP|Mono_CD14|MEP|MK|EE|CMP',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT',Idents(mlds))]),
        idents = unique(Idents(mlds)[grepl('MLDS:|TAM:',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT|CC6|CC7|CC8',Idents(mlds))]),
        scale = T,
        
        features = c(bad.vs.good_deg$geneSym[bad.vs.good_deg$direction == 'bad_up' & !grepl('AC\\d+|LINC\\d+',bad.vs.good_deg$geneSym)],
                     bad.vs.good_deg$geneSym[bad.vs.good_deg$direction == 'bad_down' & !grepl('AC\\d+|LINC\\d+',bad.vs.good_deg$geneSym)],
                     'HOXA1'
        )
        
)+RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('')











plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_0424/Plots/'
figxx_badTAM.markers_dotplot = function(){
  badTAM_markers = rbind(bad.vs.good_deg[bad.vs.good_deg$direction == 'bad_up' & !grepl('AC\\d+|LINC\\d+',bad.vs.good_deg$geneSym),],
                         bad.vs.good_deg[bad.vs.good_deg$direction == 'bad_down' & !grepl('AC\\d+|LINC\\d+',bad.vs.good_deg$geneSym),])
  
  
  plotFun_figxx_badTAM.markers_dotplot=function(noFrame=FALSE,noPlot=FALSE){
    p1 = DotPlot(mlds,idents = unique(Idents(mlds)[grepl('MLDS:|TAM:',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT',Idents(mlds))]),
                 cols = c(colAlpha(grey(0.98),0.85),'black'),
                 features = badTAM_markers$geneSym
    )+
      RotatedAxis()+
      theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9,face=ifelse(badTAM_markers$isTF,'bold',
                                                                                           ifelse(badTAM_markers$isCosmic,'italic','plain')),
                                       colour =ifelse(badTAM_markers$chr == 'chr21','purple','black')),
            axis.text.y = element_text(size=11),
            legend.position = 'top',legend.text = element_text(size=9),legend.title = element_text(size=8.5),legend.key.size = unit(0.6,'cm')) + xlab('') + ylab('') 
    
    
    print(p1)  
  }
  saveFig(file.path(plotDir,'Figxx_badTAM_degs_dotplot'),plotFun_figxx_badTAM.markers_dotplot,rawData=df,width = 7.3,height = 6,res = 500)
  
  
  plotFun_mainFig.3_MLDS.module_dotplot_vertical=function(noFrame=FALSE,noPlot=FALSE){
    p1 = DotPlot(mlds,idents = unique(Idents(mlds)[grepl('MLDS:|TAM:',Idents(mlds)) & !grepl('unsure|Tum_MK|MK_WT',Idents(mlds))]),
                 cols = c(colAlpha(grey(0.98),0.85),'black'),
                 features = c(mlds_markers$geneSym[mlds_markers$direction == 'TAM_down'],
                              mlds_markers$geneSym[mlds_markers$direction == 'TAM_up'])
    )+
      RotatedAxis()+coord_flip()+
      theme(axis.text.y = element_text(size=9,face=ifelse(mlds_markers$isTF,'bold',ifelse(mlds_markers$isCosmic,'italic','plain')),
                                       colour =ifelse(mlds_markers$chr == 'chr21','purple','black')),
            axis.text.x = element_text(size=11,angle = 90,vjust = 0.5,hjust = 1),
            legend.position = 'right',legend.text = element_text(size=9),legend.title = element_text(size=8.5),legend.key.size = unit(0.6,'cm')) + xlab('') + ylab('') 
    
    
    print(p1)  
  }
  saveFig(file.path(plotDir,'Fig3a_MLDS.module_dotplot_vertical'),plotFun_mainFig.3_MLDS.module_dotplot_vertical,rawData=df,width = 6.5,height = 10,res = 500)
  
}








##------- In big.srat
big.srat = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24.RDS')
## Dot plot of MLDS-sepcific genes in other leukaemia
big.srat$timePoint[big.srat$timePoint == 'Diagnostic'] = 'D0'
big.srat$group_4_dotPlot = ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID != 'L041' & big.srat$timePoint %in% c('D0','Diagnostic','RelapseD0') & !big.srat$disease %in% c('TAM','MLDS'), 
                                  paste0('Tumour.',big.srat$disease,'.',big.srat$timePoint,'.',big.srat$tissue),
                                  ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID != 'L041' & big.srat$timePoint %in% c('D0','Diagnostic','RelapseD0') & big.srat$disease %in% c('TAM','MLDS'),
                                         paste0('Tumour.',big.srat$disease,'.',big.srat$timePoint,'.',big.srat$tissue,'.',big.srat$donorID),
                                         ifelse(big.srat$broadLineage == 'Tumour' & big.srat$donorID == 'L038' & big.srat$timePoint %in% c('TP1') & big.srat$disease %in% c('MLDS'),
                                                'Tumour.MLDS.TP1.BM.L038',
                                                ifelse(big.srat$broadLineage == 'Tumour' & big.srat$disease == 'Lymphoma','Tumour.Lymphoma.D0.L062','others'))))
big.srat$group_4_dotPlot[grepl('^Tumour\\.BALL\\.D0',big.srat$group_4_dotPlot) & big.srat$Genotype == 'T21'] = gsub('Tumour.BALL.D0','Tumour.BALL.T21.D0',big.srat$group_4_dotPlot[grepl('^Tumour\\.BALL\\.D0',big.srat$group_4_dotPlot) & big.srat$Genotype == 'T21'])

big.srat$group_4_dotPlot = factor(big.srat$group_4_dotPlot,levels = c(unique(big.srat$group_4_dotPlot[grepl('Tumour.MLDS',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.TAM',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.MDS',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.AMKL',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.BALL\\.T21',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.BALL\\.D0',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.BALL.Relapse',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.iALL',big.srat$group_4_dotPlot)]),
                                                                      unique(big.srat$group_4_dotPlot[grepl('Tumour.Lymphoma',big.srat$group_4_dotPlot)]),
                                                                      'others'
))
Idents(big.srat) = big.srat$group_4_dotPlot

DotPlot(big.srat,idents = unique(big.srat$group_4_dotPlot[grepl('^Tumour.*Diag|Tumour.*D0',big.srat$group_4_dotPlot)]),
        features = c(bad.vs.good_deg$geneSym[bad.vs.good_deg$direction == 'bad_up' & !grepl('AC\\d+|LINC\\d+',bad.vs.good_deg$geneSym)],
                                bad.vs.good_deg$geneSym[bad.vs.good_deg$direction == 'bad_down' & !grepl('AC\\d+|LINC\\d+',bad.vs.good_deg$geneSym)])
        )+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=9),
        axis.text.y = element_text(size=11),
        legend.position = 'top',legend.text = element_text(size=9),legend.key.size = unit(0.5,'cm')) + xlab('') + ylab('') 





