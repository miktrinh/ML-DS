##--- Module scoring in bulk / sc / sn RNAseq dataset ---##

source('~/lustre_mt22/Aneuploidy/scripts/finalScripts/xx01_moduleScore_helperFunctions.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')


##----------------------------##
##    Define gene module    ####
##----------------------------##

#module_type = 'MLDS_topGenes_2410_v2' 
module_type = 'L076_relapse_module' 
## 24-11 version: using the same list of genes, but bi-directional UCell
# MLDS_topGenes_2410_cf2
# GATA1s_topGenes2v2_2410
# MLDS_T21imprint
# L076_relapse_module
# L076_L038_badMLDS_module


## 24-10 version
#MLDS_topGenes_2410_v2
#GATA1s_topGenes2v2_2410


## 24-07 version
#MLDS_T21imprint
#GATA1s_topGenes3
#MLDS_topGenes2

##----- bad MLDS L076-L038 gene Module  -----##

if(module_type == 'L076_L038_badMLDS_module'){
  bad.vs.good.mlds_deg = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/7_bad.vs.good_MLDS/DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')
  table(bad.vs.good.mlds_deg$direction)
  geneList = split(bad.vs.good.mlds_deg$geneSym,bad.vs.good.mlds_deg$direction)
  geneList[['badMLDS_composite']] = c(paste0(geneList[["MLDS.good_down"]],'+'),paste0(geneList[['MLDS.good_up']],'-'))  
  
  sc_tittle = 'L076.L038 bad MLDS transcriptional module'
}

##----- Relapse module in L076 --------##
if(module_type == 'L076_relapse_module'){
  l076_resistance_geneModule = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/7_bad.vs.good_MLDS/L076_resistance_geneModules.csv')
  table(l076_resistance_geneModule$direction,l076_resistance_geneModule$comp)
  l076_resistance_geneModule$group = paste0(l076_resistance_geneModule$comp,':',l076_resistance_geneModule$direction)
  geneList = split(l076_resistance_geneModule$geneSym,l076_resistance_geneModule$group)
  for(comparison in unique(l076_resistance_geneModule$comp)){
    up = names(geneList)[grepl(comparison,names(geneList)) & grepl('up',names(geneList))]
    down = names(geneList)[grepl(comparison,names(geneList)) & grepl('down',names(geneList))]
    geneList[[paste0(comparison,'_composite')]] = c(paste0(geneList[[up]],'+'),paste0(geneList[[down]],'-'))  
  }
  
  sc_tittle = 'L076 relapse transcriptional module'
}



##----- MLDS imprints in T21 --------##
if(module_type == 'MLDS_T21imprint'){
  #genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/DESeq2_T21.vs.2n.fLiver.MEMP_MLDSgenes.csv') # where did this file come from?
  #genes = genes[genes$direction != 'nonDE',]
  genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/3_MLDSimprint_in_fLiverT21/jul24/MLDS_imprint_in_T21.fLiver.MEMP_2411.csv',row.names = 1)
  table(genes$direction)
  geneList = split(genes$geneSym,genes$direction)
  geneList[['T21imprint_composite']] = c(paste0(geneList[['T21_MLDS_up']],'+'),paste0(geneList[['T21_MLDS_down']],'-'))
  module_type = 'MLDS_T21imprint' 
  sc_tittle = 'T21 - ML-DS transcriptional module'
}

# if(module_type == 'MLDS_T21imprint_v2'){
#   #genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/DESeq2_T21.vs.2n.fLiver.MEMP_MLDSgenes.csv') # where did this file come from?
#   #genes = genes[genes$direction != 'nonDE',]
#   genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/MLDS_imprint_in_T21.fLiver.MEMP_v2.csv',row.names = 1)
#   geneList = split(genes$geneSym,genes$direction)
#   
#   module_type = 'MLDS_T21imprint_v2' 
#   sc_tittle = 'MLDS imprints in T21 fLiver MEMP - v2'
# }
# 
# 
# if(module_type == 'MLDS_T21imprint_noChr21' ){
#   #genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/DESeq2_T21.vs.2n.fLiver.MEMP_MLDSgenes.csv') # where did this file come from?
#   #genes = genes[genes$direction != 'nonDE',]
#   genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/MLDS_imprint_in_T21.fLiver.MEMP_noChr21.csv',row.names = 1)
#   geneList = split(genes$geneSym,genes$direction)
#   
#   module_type = 'MLDS_T21imprint_noChr21' 
#   sc_tittle = 'MLDS imprints in T21 fLiver MEMP (no Chr21)'
# }
# 
# if(module_type == 'MLDS_T21imprint_noChr21_v2' ){
#   #genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/DESeq2_T21.vs.2n.fLiver.MEMP_MLDSgenes.csv') # where did this file come from?
#   #genes = genes[genes$direction != 'nonDE',]
#   genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/MLDS_imprint_in_T21.fLiver.MEMP_v2.csv',row.names = 1)
#   genes = genes[!genes$chr %in% c('chr21','21'),]
#   geneList = split(genes$geneSym,genes$direction)
#   
#   module_type = 'MLDS_T21imprint_noChr21_v2' 
#   sc_tittle = 'MLDS imprints in T21 fLiver MEMP (no Chr21) - v2'
# }
# 
# 
# 
# if(module_type == 'MLDS_T21imprint_onlyChr21'){
#   #genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/DESeq2_T21.vs.2n.fLiver.MEMP_MLDSgenes.csv') # where did this file come from?
#   #genes = genes[genes$direction != 'nonDE',]
#   genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/MLDS_imprint_in_T21.fLiver.MEMP.csv',row.names = 1)
#   genes = genes[genes$chr == 21,]
#   geneList = split(genes$geneSym,genes$direction)
#   
#   module_type = 'MLDS_T21imprint_onlyChr21' 
#   sc_tittle = 'MLDS imprints in T21 fLiver MEMP (only Chr21)'
# }
# 
# 
# if(module_type == 'MLDS_T21imprint_onlyChr21_v2'){
#   #genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/DESeq2_T21.vs.2n.fLiver.MEMP_MLDSgenes.csv') # where did this file come from?
#   #genes = genes[genes$direction != 'nonDE',]
#   genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/MLDSimprint_in_fLiverT21/MLDS_imprint_in_T21.fLiver.MEMP_v2.csv',row.names = 1)
#   genes = genes[genes$chr == 21,]
#   geneList = split(genes$geneSym,genes$direction)
#   
#   module_type = 'MLDS_T21imprint_onlyChr21_v2' 
#   sc_tittle = 'MLDS imprints in T21 fLiver MEMP (only Chr21) - v2'
# }






##---   GATA1s module ----##
if(module_type == 'GATA1s' ){
  #gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/jul24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  table(gata1s_module$tam_vs_mempT21_group)
  gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
  gata1s_module = gata1s_module[!grepl('^LINC\\d+|AC\\d+',gata1s_module$geneSym),]
  sc_tittle = 'GATA1s upregulated module'
  geneList = split(gata1s_module$geneSym,gata1s_module$group)
  module_type = 'GATA1s'
}

if(module_type == 'GATA1s_topGenes2'){
  #gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/jul24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  table(gata1s_module$tam_vs_mempT21_group)
  gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
  ## Select top genes only
  gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 20/100 & !grepl('^LINC\\d+|AC\\d+',gata1s_module$geneSym),]
  geneList = split(gata1s_module$geneSym,gata1s_module$group)
  sc_tittle = 'GATA1s top upregulated gene module'
  module_type = 'GATA1s_topGenes2'
}

if(module_type == 'GATA1s_topGenes3'){
  #gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/jul24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  table(gata1s_module$tam_vs_mempT21_group)
  gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
  ## Select top genes only
  gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 30/100 & !grepl('^LINC\\d+|AC\\d+',gata1s_module$geneSym),]
  geneList = split(gata1s_module$geneSym,gata1s_module$group)
  sc_tittle = 'GATA1s top upregulated gene module'
  module_type = 'GATA1s_topGenes3'
}



##----- GATA1s top genes 24-10 version
if(module_type == 'GATA1s_topGenes2v2_2410'){
  gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/oct24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  table(gata1s_module$tam_vs_mempT21_group)
  gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
  ## Select top genes only
  gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 20/100,]
  table(gata1s_module$group)
  geneList = split(gata1s_module$geneSym,gata1s_module$group)
  geneList[['GATA1s_composite']] = c(paste0(geneList[['TAM.MLDS.up']],'+'),paste0(geneList[['TAM.MLDS.down']],'-'))
  sc_tittle = 'GATA1s top upregulated gene module'
}


if(module_type == 'GATA1s_topGenes3_2410'){
  gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/oct24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  table(gata1s_module$tam_vs_mempT21_group)
  gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
  ## Select top genes only
  gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 30/100 & !grepl('^LINC\\d+|AC\\d+',gata1s_module$geneSym),]
  geneList = split(gata1s_module$geneSym,gata1s_module$group)
  sc_tittle = 'GATA1s top upregulated gene module'
  module_type = 'GATA1s_topGenes3_2410'
}


## Version 2: with all TAM, not just canonical TAM
if(module_type == 'GATA1s_topGenes_v2'){
  gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.allTAM.MLDS_topGeneModule.csv',row.names = 1)
  table(gata1s_module$tam_vs_mempT21_group)
  gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
  ## Select top genes only
  gata1s_module = gata1s_module[abs(gata1s_module$pct.diff) >= 30 & !grepl('^LINC\\d+|^AC\\d+',gata1s_module$geneSym),]
  geneList = split(gata1s_module$geneSym,gata1s_module$group)
  sc_tittle = 'GATA1s top upregulated gene module - allTAM'
  module_type = 'GATA1s_topGenes_v2'
}


##---   MLDS vs goodTAM module ----##
if(module_type == 'MLDS_topGenes'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')
  
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  sc_tittle = 'MLDS transcriptional module (up-regulated genes only)'
}

if(module_type == 'MLDS_topGenes2'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.2,]
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only)'
}

if(module_type == 'MLDS_topGenes3'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.2 & !grepl('^HLA-',mlds_vs_goodTAM$geneSym),]
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only)'
}

if(module_type == 'MLDS_topGenes_2410'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only)'
}


if(module_type == 'MLDS_topGenes_2410_v1'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
  mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.1 & !grepl('^HLA-',mlds_vs_goodTAM$geneSym),]
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only)'
}

if(module_type == 'MLDS_topGenes_2410_noFilter'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  geneList[['MLDS_composite']] = c(paste0(geneList[['MLDS_up']],'+'),paste0(geneList[['MLDS_down']],'-'))
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only) - noFilter'
}

if(module_type == 'MLDS_topGenes_2410_cf1'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
  mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.1,]
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  geneList[['MLDS_composite']] = c(paste0(geneList[['MLDS_up']],'+'),paste0(geneList[['MLDS_down']],'-'))
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only) - CF >= 0.1'
}

if(module_type == 'MLDS_topGenes_2410_cf2'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
  mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.2,]
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  geneList[['MLDS_composite']] = c(paste0(geneList[['MLDS_up']],'+'),paste0(geneList[['MLDS_down']],'-'))
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only) - CF >= 0.2'
}

if(module_type == 'MLDS_topGenes_2410_v2'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
  mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.2,]
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only)'
}

if(module_type == 'MLDS_topGenes_2410_l2FC_1'){
  #mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/DESeq2_goodTAM.vs.individual.dMLDS_topDEGs.csv')
  mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
  mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.2 & abs(mlds_vs_goodTAM$log2FoldChange) >= 0.75,]
  mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.2 & !grepl('^HLA-',mlds_vs_goodTAM$geneSym),]
  geneList = split(mlds_vs_goodTAM$geneSym,mlds_vs_goodTAM$direction)
  sc_tittle = 'MLDS transcriptional module (top up-regulated genes only)'
}


if(module_type == 'badTAM'){
  bad.vs.good_deg = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/goodTAM_vs_badTAM/DESeq2_good.vs.bad.TAM_topDEGs.csv')
  geneList = split(bad.vs.good_deg$geneSym,bad.vs.good_deg$direction)
  
  moduleList = split(bad.vs.good_deg$ensID,bad.vs.good_deg$direction)
  
  module_type = 'badTAM'
}


if(module_type == 'badTAM_topGenes'){
  bad.vs.good_deg = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/goodTAM_vs_badTAM/DESeq2_good.vs.bad.TAM_topDEGs.csv')
  badTAM_markers = rbind(bad.vs.good_deg[bad.vs.good_deg$direction == 'bad_up' & !grepl('AC\\d+|LINC\\d+|^HLA-',bad.vs.good_deg$geneSym),],
                         bad.vs.good_deg[bad.vs.good_deg$direction == 'bad_down' & !grepl('AC\\d+|LINC\\d+|^HLA-',bad.vs.good_deg$geneSym),])
  geneList = split(badTAM_markers$geneSym,badTAM_markers$direction)
  
  moduleList = split(badTAM_markers$ensID,badTAM_markers$direction)
  
  module_type = 'badTAM_topGenes'
}










##--------------------------##
##    Score the module    ####
##--------------------------##

scRNA_moduleScore_UCell(geneList = geneList,
                        outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/UCell_moduleScore/nov24',
                        module_type=module_type,datasets = c('fLiver','fLiver_Muz','MLDS','otherLeuk','fBM_2n','fBM_T21','infantALL'),ncores=10)
                        #module_type=module_type,datasets = c('otherLeuk'),ncores=10)
  
# scRNA_moduleScore_UCell(geneList = geneList,
#                         outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/UCell_moduleScore',
#                         module_type=module_type,datasets = c('fLiver','fLiver_Muz','MLDS', 'MDS','BALL','fBM_2n','fBM_T21','infantALL','pAML'),ncores=10)
#
# scRNA_moduleScore_UCell(geneList = geneList,
#                         outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/UCell_moduleScore',
#                         module_type=module_type,datasets = c('fLiver_Muz'),ncores=10)
# 









##----------------------------------##
##    Import the module scores    ####
##----------------------------------##

# Import latest MLDS annotation
mlds_mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv')
mlds_mdat$finalAnn_broad = mlds_mdat$annot_aug24
mlds_mdat$annot = mlds_mdat$annot_aug24
# #table(data$finalAnn_broad[data$cellID %in% mlds_mdat$cellID])
# #table(mlds_mdat$annot_mar24[mlds_mdat$cellID %in% data$cellID])
# Import latest big.srat (combined Leukaemia) annotation
otherLeuk_mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2408_mdat.csv')
otherLeuk_mdat$annot[otherLeuk_mdat$donorID == 'L010' & otherLeuk_mdat$annot == 'T_cells'] = 'unsure'
otherLeuk_mdat$broadLineage[otherLeuk_mdat$donorID == 'L010' & otherLeuk_mdat$annot == 'unsure'] = 'unsure'
# # Remove MLDS / TAM samples
# combinedLeuk_mdat = combinedLeuk_mdat[!combinedLeuk_mdat$cellID %in% mlds_mdat$cellID,]
# combinedLeuk_mdat$timePoint[combinedLeuk_mdat$donorID == 'L062' & combinedLeuk_mdat$disease == 'Lymphoma'] = 'D0'


## Import moduleScore results
data = import_UCell_result(outDir= '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/UCell_moduleScore/nov24',
                           module_type=module_type)

## Add latest annotation
data$finalAnn_broad[data$cellID %in% mlds_mdat$cellID] = mlds_mdat$annot[match(data$cellID[data$cellID %in% mlds_mdat$cellID],mlds_mdat$cellID)]
data$finalAnn_broad[data$cellID %in% otherLeuk_mdat$cellID] = otherLeuk_mdat$annot[match(data$cellID[data$cellID %in% otherLeuk_mdat$cellID],otherLeuk_mdat$cellID)]

## Add timepoint to data
#data$donorID[is.na(data$donorID) & data$dataset == 'infantALL'] = '?'
#data$finalAnn_broad[data$cellID %in% mlds_mdat$cellID] = mlds_mdat$annot_mar24[match(data$cellID[data$cellID %in% mlds_mdat$cellID],mlds_mdat$cellID)]
data$timePoint[data$cellID %in% mlds_mdat$cellID] = mlds_mdat$timePoint[match(data$cellID[data$cellID %in% mlds_mdat$cellID],mlds_mdat$cellID)]
data$timePoint[data$cellID %in% otherLeuk_mdat$cellID] = otherLeuk_mdat$timePoint[match(data$cellID[data$cellID %in% otherLeuk_mdat$cellID],otherLeuk_mdat$cellID)]
data$timePoint[is.na(data$timePoint)] = 'NA'
data$timePoint[data$timePoint %in% c('D.Relapse','RelapseD0')] = 'RD'
data$timePoint[data$timePoint %in% c('D.Relapse2')] = 'R2D'
#data$timePoint[data$timePoint %in% c('RelapseD29')] = 'RD29'
data$timePoint[data$timePoint == 'Diagnostic'] = 'D'


module_name = colnames(data)[grepl('_UCell$',colnames(data))]#paste0(names(geneList)[1],'_UCell')
if(module_type %in% c('GATA1s','GATA1s_topGenes','GATA1s_topGenes2','GATA1s_topGenes3','GATA1s_topGenes_v2','GATA1s_topGenes2_2410','GATA1s_topGenes3_2410','GATA1s_topGenes2v2_2410')){
  data$module = data[[module_name[5]]]
}else if(module_type %in% c('badTAM','badTAM_topGenes','MLDS_T21imprint','MLDS_T21imprint_v2','MLDS_T21imprint_noChr21','MLDS_T21imprint_noChr21_v2','MLDS_T21imprint_onlyChr21_v2')){
  #data$module = data[[module_name[2]]]
  data$module = data[[module_name[grepl('composite',module_name)]]]
}else if(module_type %in% c('MLDS_topGenes','MLDS_topGenes2','MLDS_topGenes3','MLDS_topGenes_2410','MLDS_topGenes_2410_v1','MLDS_topGenes_2410_v2','MLDS_topGenes_2410_noFilter','MLDS_topGenes_2410_cf1','MLDS_topGenes_2410_cf2')){
  #data$module = data[['MLDS_up_UCell']]
  data$module = data[[module_name[grepl('MLDS_composite',module_name)]]]
}else{
  print('Unrecognised module type')
}


## Reassign data to different dataset categories
data = data[data$dataset_ID %in% c('fBM','fBM.T21','fLiver','fLiverMuz','MLDS','otherLeuk'),]
data$dataset_ID_og = data$dataset_ID
table(data$dataset_ID_og)
data$dataset_ID[data$dataset_ID == 'fLiver' & data$Genotype == 'diploid' & !grepl('^Hsb',data$donorID)] = 'FL_ext'
data$dataset_ID[data$dataset_ID == 'fLiverMuz'] = 'FL_ext'
data$dataset_ID[data$dataset_ID == 'fBM.T21'] = 'FBM'
data$dataset_ID[data$dataset_ID == 'fBM'] = 'FBM'
data$dataset_ID[data$dataset_ID == 'fLiver'] = 'FL'
#data$dataset_ID[data$dataset_ID == 'bALL'] = 'B-ALL'
data$dataset_ID[data$dataset_ID == 'otherLeuk'] = data$disease[data$dataset_ID == 'otherLeuk']
data$dataset_ID[data$dataset_ID_og == 'otherLeuk' & data$disease == 'pAML' & data$donorID == 'L069'] = 'AEL'
data$dataset_ID[data$dataset_ID_og == 'otherLeuk' & data$donorID == 'L010'] = 'AMKL'
data$dataset_ID[data$dataset_ID == 'MLDS' & data$donorID %in% c('CC1','CC2','L075','L114','L182','L156','CC6','CC7','CC8')] = 'TAM'
# data$dataset_ID[data$dataset_ID == 'infantALL' & !is.na(data$donorID) & data$donorID == 'P9_iAML/Cancer'] = 'AMKL'
# data$donorID[data$donorID == 'L062' & data$orig.ident %in% c('NB14406184','NB14406185')] = 'L062_LPD'
# data$dataset_ID[data$donorID == 'L062_LPD'] = 'LPD'

data$donorID[data$donorID == 'L076' & data$timePoint == 'RD'] = 'L076_RD_BM'
data$donorID[data$donorID == 'L076' & data$timePoint == 'R2D'] = 'L076_R2D_BM'
data$donorID[data$donorID == 'L076' & data$timePoint == 'D' & data$tissue == 'BM'] = 'L076_BM'
data$donorID[data$donorID == 'L076' & data$timePoint == 'D' & data$tissue == 'Blood'] = 'L076_PB'

data$donorID[data$donorID == 'L038' & data$timePoint == 'TP1'] = 'L038_TP1'
data$donorID[data$donorID == 'L038' & data$timePoint == 'D'] = 'L038'

data$donorID[data$donorID == 'L156' & data$timePoint == 'RD'] = 'L156_RD'

data$donorID[data$donorID == 'L080' & data$timePoint == 'RD'] = 'L080_RD'
data$donorID[data$donorID == 'L001' & data$timePoint == 'RD'] = 'L001_RD'

data$donorID = gsub('_iA.*\\/Cancer$','',data$donorID)
data$finalAnn_broad[data$finalAnn_broad %in% c('Tumour','TumourCNA')] = 'Cancer'
data$finalAnn_broad[data$finalAnn_broad %in% c('MEP','MEMP')] = 'MEMP_MEP'
data$finalAnn_broad[data$finalAnn_broad %in% c('B.cell','naive.B','naive B cell')] = 'B.cell'

data$dataset_ID = as.character(data$dataset_ID)
data$dataset_ID[data$dataset_ID %in% c('MDS','AMKL','AEL')] = 'MDS_AMKL_AEL'
data$dataset_ID[data$dataset_ID == 'TAM' & data$donorID %in% c('L156','L156_RD')] = 'Recurrent TAM'
data$dataset_ID = factor(data$dataset_ID,c('FL','FL_ext','FBM','TAM','Recurrent TAM','MLDS','MDS','AMKL','AEL','MDS_AMKL_AEL','infantALL','pAML','pBALL','LPD'))

# Define baseline = median score in fLiver_MEMP_inhouse2n + fLiver_MEMP_T21
baseLine = data[data$dataset_ID %in% c('FL') & data$Genotype %in% c('diploid','T21'),] %>% group_by(Genotype) %>% summarise(medScore = median(module))

##----- Only keep things we want to plots  ---------
df = data[data$finalAnn_broad %in% c('Cancer','MEMP_MEP','B.cell','earlyMK','HSC_MPP','HSC/MPP','mature NK','pro B progenitor'),]
data$finalAnn_broad[data$finalAnn_broad == 'Tumour' & !is.na(data$finalAnn_broad)]  = 'Cancer'
df = data[data$finalAnn_broad %in% c('Cancer','MEMP_MEP'),]
#df = df[!(df$dataset_ID %in% c('fLiver') & df$Genotype %in% c('diploid','T21')) & !(df$donorID == 'L041' & df$finalAnn_broad == 'Cancer'),]
df = df[!(df$donorID %in% c('L041','CC3') & df$finalAnn_broad == 'Cancer'),]
df = df[!(df$finalAnn_broad == 'MEMP_MEP' & df$dataset_ID %in% c('pBALL','pAML','MLDS','TAM','MDS','AMKL','AEL','infantALL','LPD')),]

# remove blasts from TP1/TP4, except for L076 and L038
df$toKeep = T
df$toKeep[!grepl('L076|L038|L156',df$donorID) & df$timePoint != 'D' & df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Cancer'] = F
#df$toKeep[!grepl('L076|L038|L156',df$donorID) & df$timePoint != 'D0' & df$dataset_ID %in% c('pBALL','pAML','MDS','AMKL','AEL','infantALL','LPD') & df$finalAnn_broad == 'Cancer'] = F
df = df[df$toKeep==T,]
#df$toKeep[!df$donorID %in% c('L076','L038') & !is.na(df$timePoint) & df$timePoint != 'Diagnostic' & df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Tumour'] = F


df$group = ifelse(df$finalAnn_broad != 'Cancer',as.character(df$finalAnn_broad),
                  ifelse(df$dataset_ID %in% c('MLDS','TAM','Recurrent TAM') & df$finalAnn_broad == 'Cancer','TAM / MLDS',
                         ifelse(df$finalAnn_broad == 'Cancer','Other leukaemia','others')))
df$group = factor(df$group,c('MEMP_MEP','TAM / MLDS','Other leukaemia'))
df$Genotype[df$Genotype == 'complete_trisomy'] = 'Triploid'
df$Genotype[df$Genotype == 'diploid'] = 'Diploid'
df$Genotype = factor(df$Genotype,c('Diploid','T21','T18','T22','MX','Triploid'))

## Add AlleleIntegrator call for L038 ##
# aiRes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x_alleleIntegrator/L038/L038_TumourEry_subClustering_AIRes_240131.csv')
# df$AIres = '?'
# df$AIres[df$cellID %in% aiRes$cellID] = aiRes$AI_output[match(df$cellID[df$cellID %in% aiRes$cellID],aiRes$cellID)]
#df$timePoint = big.srat$timePoint[match(df$cellID,big.srat$cellID)]

# Remove infantALL dataset
df = df[df$dataset_ID != 'infantALL',]


## Annotate L038_CNA cancer cells
df$donorID2 = paste0(df$donorID,':',df$AIres)
df$donorID2 = gsub(':\\?$','',df$donorID2)



# # seperate L038 w/wo CNA cancer cells
# g = 'TAM / MLDS'
# d = df[df$group == g,]
# d$donorID=factor(d$donorID,unique(d$donorID[order(as.numeric(d$Genotype))]))
# d = d[d$timePoint == 'Diagnostic' |
#         (d$timePoint == 'TP1' & d$donorID == 'L038'),]
# p1 = ggplot(d,aes(donorID2,(low_mid_high_UCell)))+
#   geom_boxplot(aes(fill = timePoint),width=0.5,outlier.size = 0.001)+
#   scale_fill_manual(values = geno_cols)+
#   scale_fill_manual(values = c(col25[c(5:7)],grey(0.8)))+
#   facet_grid(.~dataset_ID,scales = 'free_x',space = 'free')+
#   geom_hline(yintercept = baseLine$medScore[baseLine$Genotype == 'diploid'],color='grey')+
#   geom_hline(yintercept = baseLine$medScore[baseLine$Genotype == 'T21'],color='#6A3D9A')+
#   ylim(min(df$low_mid_high_UCell),max(df$low_mid_high_UCell))+
#   xlab('')+ylab('MLDS module score')+
#   theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + ggtitle(g)
# p1
#
# ggplot(d[d$donorID == 'L038',],aes(timePoint,(low_mid_high_UCell)))+
#   geom_boxplot(aes(fill = AIres),width=0.5,outlier.size = 0.001)+
#   scale_fill_manual(values = c(col25[c(1,2,3)],grey(0.8)))+
#   facet_grid(.~dataset_ID,scales = 'free_x',space = 'free')+
#   geom_hline(yintercept = baseLine$medScore[baseLine$Genotype == 'diploid'],color='grey')+
#   geom_hline(yintercept = baseLine$medScore[baseLine$Genotype == 'T21'],color='#6A3D9A')+
#   ylim(min(df$low_mid_high_UCell),max(df$low_mid_high_UCell))+
#   xlab('')+ylab('MLDS module score')+
#   theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + ggtitle(g,subtitle = 'L038')


# df$xlab = paste0(df$dataset_ID,':',df$Genotype)
# df$xlab = factor(df$xlab,c('fLiver:MX','fLiver:T18','fLiver:T22','fLiver:complete_trisomy','fLiver_ext:diploid',
#                            'fBM:diploid','fBM:T21',
#                            'MLDS:T21','TAM:T21',
#                            'MDS:diploid','AMKL:diploid','pAML:diploid','infantALL:diploid','bALL:diploid','bALL:T21'))

geno_cols = c('Diploid' = grey(0.2),
              'T21' = '#93221E',
              'T18' = '#3d5dad',
              'T22' = '#679551',
              'T13' = '#526691',
              'MX' = '#b18db8',
              'Triploid' = '#e07d26')

plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots/noCC3/'
plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/Plots'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')

plotFun_UCell_moduleScore_v2 = function(noFrame=FALSE,noPlot=FALSE){
  library(ggbeeswarm)

  df$groupID = ifelse(!df$group %in% c('Other leukaemia','TAM / MLDS'),paste0(df$group,'_',df$dataset_ID,':',df$Genotype),
                      paste0(df$group,'_',df$dataset_ID,':',df$Genotype,':',df$donorID))
                      # ifelse(df$donorID =='L076',paste0(df$group,'_',df$dataset_ID,':',df$Genotype,':',df$donorID,'_',df$timePoint,'_',df$tissue),
                      #        ifelse(df$donorID =='L038',paste0(df$group,'_',df$dataset_ID,':',df$Genotype,':',df$donorID,'_',df$timePoint),
                      #               )))

  #d = df[df$group %in% c('TAM / MLDS', 'Other leukaemia'),]
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

  ## Define the x-axis order
  groupID_order.pre = c(paste0(sapply(levels(d$group),function(x){rep(x,length(levels(d$dataset_ID)))}) %>% as.vector(),'_',levels(d$dataset_ID)))
  groupID_order = paste0(sapply(groupID_order.pre[grepl('MEMP_MEP',groupID_order.pre)],function(x){rep(x,length(levels(d$Genotype)))}) %>% as.vector(),
                         ':',levels(d$Genotype))

  groupID_order_byMedian = quant_df$groupID[order(quant_df$`50%`,decreasing = F)]
  groupID_order = c(groupID_order,groupID_order_byMedian[grepl('_TAM:',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Recurrent TAM:T21:L156$',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Recurrent TAM:T21:L156_RD',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('MLDS_MLDS:',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Other.*L067$',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Other.*L010$',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Other.*P9$',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Other.*L069$',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Other.*_pAML:',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Other.*_pBALL:Diploid',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Other.*_pBALL:T21',groupID_order_byMedian)],
                    groupID_order_byMedian[grepl('Other.*_LPD:T21',groupID_order_byMedian)])

  table(!df$groupID %in% groupID_order)

  if(!all(unique(d$groupID) %in% groupID_order)){
    notFound = unique(d$groupID)[!unique(d$groupID) %in% groupID_order]
    stop(sprintf('These groups were not found %s',paste(notFound,collapse = ', ')))
  }
  #groupID_order[!(groupID_order %in% unique(d$groupID))]
  groupID_order = groupID_order[groupID_order %in% unique(d$groupID)]
  length(groupID_order)
  groupID_order[duplicated(groupID_order)]

  quant_df$groupID_level = NA
  for(g in unique(quant_df$group)){
    for(g.sub in unique(quant_df$dataset_ID[quant_df$group == g])){
      quant_df$groupID_level[quant_df$group == g &
                               quant_df$dataset_ID == g.sub] = as.numeric(factor(quant_df$groupID[quant_df$group == g &
                                                                                                    quant_df$dataset_ID == g.sub],
                                                                                 levels = groupID_order[groupID_order %in% quant_df$groupID[quant_df$group == g &
                                                                                                                                              quant_df$dataset_ID == g.sub]]))
    }
  }


  ## Downsample group with >=2000 cells
  d.large = d[d$nCell_perGroupID > 500,]
  set.seed(1234)
  d.large = d.large %>% group_by(groupID) %>% mutate(id = seq(1:n()),
                                                     selected = id %in% sample(1:n(),500))
  d = rbind(d[d$nCell_perGroupID <= 500,],
            d.large[d.large$selected == T,!colnames(d.large) %in% c('id','selected')])

  d$genoAlpha_group = ifelse(d$nCell_perGroupID > 300,'low',
                             ifelse(d$nCell_perGroupID < 100,'high','mid'))
  d$genoAlpha_group[grepl('MLDS_TAM:|MLDS_MLDS:|MLDS_Recurrent TAM:',d$groupID)] = 'mid'
  d$genoAlpha_group = factor(d$genoAlpha_group,c('low','mid','high'))
  d$groupID = factor(d$groupID,groupID_order)
  d$groupID2 = as.character(d$groupID)
  d$groupID2 = factor(gsub('^.*:','',as.character(d$groupID)),unique(gsub('^.*:','',groupID_order)))

  alphaValue = c(0.007,0.02,0.15)

  ## Fix max ylim
  max_ylim = max(0.2,max(quant_df$`75%`)) + 0.01
  d$module[d$module > max_ylim] = max_ylim

  #d$donorID = factor(d$donorID,gsub('^.*:','',levels(d$groupID)))
  if(!noPlot & !noFrame){
    p1 = ggplot(d,aes(groupID2,module))+
      geom_quasirandom(aes(alpha = genoAlpha_group),width = 0.3,size=0.1)+
      scale_alpha_manual(values = alphaValue)+
      geom_segment(data=quant_df,aes(x=groupID_level-0.3,y=`50%`,xend=groupID_level+0.3,yend=`50%`,col=Genotype),lwd=0.55)+
      geom_segment(data=quant_df,aes(x=groupID_level,y=`25%`,xend=groupID_level,yend=`75%`,col=Genotype),lwd=0.45)+
      scale_fill_manual(values = geno_cols)+
      scale_color_manual(values = geno_cols)+
      facet_grid(.~group+dataset_ID,scales = 'free_x',space = 'free')+
      ylim(min(min(d$module),min(quant_df$`25%`)),max_ylim)+
      xlab('')+ylab('Module score')+
      theme_classic(base_size = 13)+
      theme(axis.text.x = element_text(size=7,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text.y = element_text(size=9,colour = 'black'),
            strip.text = element_text(size = 8,colour = 'black'),
            panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.4,colour = 'black'),
            strip.background=element_rect(linewidth=0,colour = 'black'),
            legend.position = 'bottom',legend.text = element_text(size=9,colour = 'black')) +
      ggtitle(sc_tittle)
      #ggtitle('bad TAM upregulated module')
      #ggtitle('GATA1s upregulated module')
    #ggtitle('MLDS upregulated module compared to canonical TAM')
    #ggtitle('MLDS imprints in T21 fLiver MEMP')

    #ggtitle('T21 upregulated module')
    #ggtitle('agg. upregulated module')

    print(p1)

  }

  if(!noPlot & noFrame){
    p1 = ggplot(d,aes(groupID2,module))+
      geom_quasirandom(aes(alpha = genoAlpha_group),width = 0.3,size=0.1)+
      scale_alpha_manual(values = alphaValue)+
      geom_segment(data=quant_df,aes(x=groupID_level-0.3,y=`50%`,xend=groupID_level+0.3,yend=`50%`,col=Genotype),lwd=0.65,alpha=0)+
      #geom_segment(data=quant_df,aes(x=groupID_level,y=`25%`,xend=groupID_level,yend=`75%`,col=Genotype),lwd=0.55)+
      scale_fill_manual(values = geno_cols)+
      scale_color_manual(values = geno_cols)+
      facet_grid(.~group+dataset_ID,scales = 'free_x',space = 'free')+
      ylim(min(d$module),max_ylim)+
      xlab('')+ylab('Module score')+
      theme_classic(base_size = 13)+
      theme(axis.text.x = element_text(size=7,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text.y = element_text(size=9,colour = 'black'),
            strip.text = element_text(size = 8,colour = 'black'),
            #panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.4,colour = 'black'),
            strip.background=element_rect(linewidth=0,colour = 'black'),
            legend.position = 'bottom',legend.text = element_text(size=9,colour = 'black')) +
      ggtitle(sc_tittle)

    print(p1)

  }



  if(noPlot & !noFrame){
    ## Downsample
    set.seed(1234)
    d.ds = d %>% group_by(groupID) %>% mutate(id = seq(1:n()),selected = id %in% sample(1:n(),1))
    d2 = rbind(d.ds[d.ds$selected == T,!colnames(d.ds) %in% c('id','selected')])

    p1 = ggplot(d2,aes(groupID2,module))+
      geom_quasirandom(aes(alpha = genoAlpha_group),width = 0.3,size=0.1)+
      scale_alpha_manual(values = alphaValue)+
      geom_segment(data=quant_df,aes(x=groupID_level-0.3,y=`50%`,xend=groupID_level+0.3,yend=`50%`,col=Genotype),lwd=0.65)+
      geom_segment(data=quant_df,aes(x=groupID_level,y=`25%`,xend=groupID_level,yend=`75%`,col=Genotype),lwd=0.55)+
      scale_fill_manual(values = geno_cols)+
      scale_color_manual(values = geno_cols)+
      facet_grid(.~group+dataset_ID,scales = 'free_x',space = 'free')+
      ylim(min(d$module),max_ylim)+
      xlab('')+ylab('Module score')+
      theme_classic(base_size = 13)+
      theme(axis.text.x = element_text(size=7,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text.y = element_text(size=9,colour = 'black'),
            strip.text = element_text(size = 8,colour = 'black'),
            panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.4,colour = 'black'),
            strip.background=element_rect(linewidth=0,colour = 'black'),
            legend.position = 'bottom',legend.text = element_text(size=9,colour = 'black')) +
      ggtitle(sc_tittle)

    print(p1)

  }


}

# saveFig(file.path(plotDir,paste0('Fig3x_',module_type)),plotFun_MLDS_moduleScore_v2,rawData=df,width = 7,height = 4.75,res = 500,useDingbats = T)
# saveFig(file.path(plotDir,paste0('Fig3x_',module_type,'_Frame')),plotFun_MLDS_moduleScore_v2,rawData=df,width = 7,height = 3.865,res = 500,useDingbats = T)
saveFig(file.path(plotDir,paste0('Fig3x_',module_type)),plotFun_UCell_moduleScore_v2,rawData=df,width = 7.5,height = 4,res = 500)
saveFig(file.path(plotDir,paste0('Fig3x_composite_',module_type)),plotFun_UCell_moduleScore_v2,rawData=df,width = 7.21,height = 3.96,res = 500)
#saveFig(file.path(plotDir,paste0('Fig3x_',module_type,'_D0only')),plotFun_MLDS_moduleScore_v2,rawData=df,width = 7.5,height = 4.75,res = 500)
#saveFig(file.path(plotDir,paste0('Fig3x_',module_type,'_Frame')),plotFun_MLDS_moduleScore_v2,rawData=df,width = 7.3,height = 4.75,res = 500)
#saveFig(file.path(plotDir,paste0('Fig3x_',module_type,'_Frame')),plotFun_MLDS_moduleScore_v2,rawData=df,width = 7.3,height = 3.865,res = 500)






fig3G_moduleScore_groupMedian = function(){
  df$dataset_ID = as.character(df$dataset_ID)
  df$dataset_ID[df$donorID == 'L067'] = 'MDS'
  df$dataset_ID[df$donorID == 'L069'] = 'AEL'
  df$dataset_ID[df$donorID %in% c('P9','L010')] = 'AMKL'
  df$dataset_ID = factor(df$dataset_ID,c('FL','FL_ext','FBM','TAM','Recurrent TAM','MLDS','MDS','AMKL','AEL','infantALL','pAML','pBALL','LPD'))

  df$groupID = ifelse(!df$group %in% c('Other leukaemia','TAM / MLDS'),paste0(df$group,'_',df$dataset_ID,':',df$Genotype),paste0(df$group,'_',df$dataset_ID,':',df$Genotype,':',df$donorID))

  d = df %>% group_by(groupID) %>% mutate(nCell_perGroupID = n())
  table(d$groupID[d$nCell_perGroupID <= 30])
  ## Remove groups with <= 30 cells
  d = d[d$nCell_perGroupID >= 30,]

  ## Calcuate median score per individual
  d = d %>% group_by(groupID,dataset_ID,group,Genotype) %>% summarise(meanScore = mean(module),
                                                                      medianScore = median(module))

  d$dataset_ID[d$dataset_ID == 'FL_ext'] = 'FL'



  plotFun_moduleScore_medianPerGroup = function(noFrame=FALSE,noPlot=FALSE){

    if(!noPlot & !noFrame){
      p1 = ggplot(d,aes(dataset_ID,medianScore,col=Genotype))+
        stat_summary(
          aes(group = group), fun = median, fun.min = median, fun.max = median,
          geom = "crossbar", color = "black", width = 0.9, lwd = 0.2,

          # add this bit here to your stat_summary function
          position=position_dodge(width=0.75)
        )+
        geom_quasirandom(data=d[d$Genotype == 'diploid',],width = 0.5,size=1.3,alpha=1)+
        geom_quasirandom(data=d[d$Genotype != 'diploid',],width = 0.5,size=1.3,alpha=1)+

        scale_color_manual(values = geno_cols)+
        facet_grid(.~group,scales = 'free_x',space = 'free')+
        ylim(0,max(d$medianScore)+0.05)+
        xlab('')+ylab('Module score')+
        theme_classic(base_size = 12)+
        theme(axis.text.x = element_text(size=10,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
              strip.text = element_text(size = 11,colour = 'black'),
              axis.text.y = element_text(size=10,colour = 'black'),
              panel.border = element_rect(fill=F),axis.line = element_blank(),
              axis.ticks = element_line(linewidth=0.4,colour='black'),
              strip.background=element_rect(linewidth=0),
              legend.position = 'none',legend.text = element_text(size=9,colour = 'black'))

      print(p1)

    }
  }
  saveFig(file.path(plotDir,paste0('Fig3x_medianScorePerGroup_',module_type)),plotFun_moduleScore_medianPerGroup,rawData=d,width = 4.2,height = 3.3,res = 500)

}

# saveFig(file.path(plotDir,paste0('Fig3x_',module_type)),plotFun_MLDS_moduleScore_v2,rawData=df,width = 7,height = 4.75,res = 500,useDingbats = T)
# saveFig(file.path(plotDir,paste0('Fig3x_',module_type,'_Frame')),plotFun_MLDS_moduleScore_v2,rawData=df,width = 7,height = 3.865,res = 500,useDingbats = T)






##--------------------------------------------------------------------##
##  Statistical test for the module score across different groups   ####
##--------------------------------------------------------------------##
## Calculate group mean

mean_df = d %>% group_by(groupID) %>% summarise(mu = mean(module))
mean_df$group = d$group[match(mean_df$groupID,d$groupID)]
mean_df$dataset_ID = d$dataset_ID[match(mean_df$groupID,d$groupID)]
mean_df$Genotype = d$Genotype[match(mean_df$groupID,d$groupID)]

## compare the median of each group
wilcox.test(mean_df$mu[mean_df$group == 'TAM / MLDS' & mean_df$dataset_ID == 'MLDS'],mean_df$mu[mean_df$group == 'Other leukaemia'])
wilcox.test(mean_df$mu[mean_df$group == 'TAM / MLDS' & mean_df$dataset_ID == 'TAM'],mean_df$mu[mean_df$group == 'Other leukaemia'],alternative = 'less')
wilcox.test(mean_df$mu[mean_df$group == 'TAM / MLDS' & mean_df$dataset_ID == 'TAM'],mean_df$mu[mean_df$group == 'MEMP_MEP'],alternative = 'less')
wilcox.test(mean_df$mu[mean_df$group == 'TAM / MLDS' & mean_df$dataset_ID == 'TAM'],mean_df$mu[mean_df$group == 'MEMP_MEP' & mean_df$dataset_ID == 'FL' & mean_df$Genotype == 'T21'])
wilcox.test(mean_df$mu[mean_df$dataset_ID == 'TAM'],mean_df$mu[mean_df$dataset_ID == 'MLDS'])

wilcox.test(mean_df$mu[mean_df$group == 'TAM / MLDS'],mean_df$mu[mean_df$group == 'MEMP_MEP'],alternative = 'greater')
wilcox.test(mean_df$mu[mean_df$group == 'TAM / MLDS'],mean_df$mu[mean_df$group == 'Other leukaemia' & mean_df$dataset_ID != 'MDS'],alternative = 'greater')
wilcox.test(mean_df$mu[mean_df$group == 'TAM / MLDS'],mean_df$mu[mean_df$dataset_ID == 'MDS'],alternative = 'greater')
wilcox.test(mean_df$mu[mean_df$group == 'TAM / MLDS' & mean_df$dataset_ID == 'TAM'],mean_df$mu[mean_df$group == 'TAM / MLDS' & mean_df$dataset_ID == 'MLDS'])


##----------------------------------------------##
##    Score module in bulkRNAseq dataset      ####
##----------------------------------------------##
library(readxl)
library(SummarizedExperiment)
library(GenomicFeatures)
source('~/lustre_mt22/Aneuploidy/scripts/finalScripts/xx01_moduleScore_helperFunctions.R')


#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)




##----- 1. Import gene module -----##
## MLDS imprints
genes = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/3_MLDSimprint_in_fLiverT21/jul24/MLDS_imprint_in_T21.fLiver.MEMP.csv',row.names = 1)
moduleList = split(genes$ensID,genes$direction)
module_type = 'MLDS_imprints_T21FL'
moduleList[['all']] = genes$ensID[genes$direction %in% c("T21_MLDS_down", "T21_MLDS_up")]
names(moduleList)
moduleList = moduleList[unlist(lapply(moduleList,function(x){length(x)>=5}))]




## GATA1s module
gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/17_GATA1s_module/goodTAM/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
table(gata1s_module$tam_vs_mempT21_group)
gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
## Select top genes only
gata1s_module = gata1s_module[abs(gata1s_module$pct.diff) >= 30,]

gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/jul24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
table(gata1s_module$tam_vs_mempT21_group)
gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
## Select top genes only
gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 20/100 & !grepl('^LINC\\d+|AC\\d+',gata1s_module$geneSym),]


moduleList = split(gata1s_module$ensID,gata1s_module$group)
module_type = 'GATA1s_topGenes3'
moduleList[['all']] = gata1s_module$ensID[gata1s_module$group %in% c("TAM.MLDS.down", "TAM.MLDS.up")]
names(moduleList)
length(moduleList[[4]])
moduleList = moduleList[unlist(lapply(moduleList,function(x){length(x)>=5}))]

## MLDS - vs - good TAM
# # mlds_vs_goodTAM = read_csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/18_TAM_analysis/MLDS_vs_goodTAM/FM_L075.vs.indiv_MLDS_FM_MLDS.vs.goshT21_topMarkers.csv')
# # mlds_vs_goodTAM$ensID = geneMap$ensID[match(mlds_vs_goodTAM$geneSym,geneMap$geneSym)]
# # module_up = mlds_vs_goodTAM$ensID[mlds_vs_goodTAM$direction == 'L075_up']
# # module_down = tam.vs.mepT21.down$ensID
# mlds_vs_goodTAM = mlds.vs.L075.markers_summary[mlds.vs.L075.markers_summary$nDonor >= 7 & !mlds.vs.L075.markers_summary$geneSym %in% c('LMO4',
#                                                                                                                                        'NFKBIA','RALGAPA2','KIAA1211',
#                                                                                                                                        'RFX7','CDC42BPA','SLC25A21',
#                                                                                                                                        'SEC31A'),]
# table(mlds_vs_goodTAM$direction)
# mlds_vs_goodTAM$ensID = geneMap$ensID[match(mlds_vs_goodTAM$geneSym,geneMap$geneSym)]
# moduleList = split(mlds_vs_goodTAM$ensID,mlds_vs_goodTAM$direction)
# moduleList[['all']] = mlds_vs_goodTAM$ensID
# names(moduleList)
# length(moduleList[[1]])
#
# moduleList[['badTAM_up']] = geneMap$ensID[match(geneOrder$geneSym[geneOrder$group == '3'],geneMap$geneSym)]
#mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs.csv')
mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.2,]
moduleList = split(mlds_vs_goodTAM$ensID,mlds_vs_goodTAM$direction)
moduleList[['all']] = mlds_vs_goodTAM$ensID
module_type = 'MLDS_topGenes'
# moduleList[['all']] = gata1s_module$ensID[gata1s_module$group %in% c("TAM.MLDS.down", "TAM.MLDS.up")]
# names(moduleList)
# length(moduleList[[1]])
# moduleList = moduleList[unlist(lapply(moduleList,function(x){length(x)>=5}))]


## good vs bad TAM module
# good.vs.bad.TAM_deg$ensID[good.vs.bad.TAM_deg$is_topDEG==T]
# moduleList = split(good.vs.bad.TAM_deg$ensID[good.vs.bad.TAM_deg$is_topDEG==T],good.vs.bad.TAM_deg$direction[good.vs.bad.TAM_deg$is_topDEG==T])



moduleList = split(resistanceModule$ensID[resistanceModule$chr %in% c('chr5','chr17')],resistanceModule$direction[resistanceModule$chr %in% c('chr5','chr17')])
module_type = 'resistance_topGenes'








plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots/noCC3'
plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/Plots'
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')








id_col = 'Sample'



figS3_gata1s_moduleScore_inBulkSamples = function(){
  library(ggbeeswarm)
  title = 'GATA1s module'
  title = 'MLDS imprints in T21 FL'
  allScore = bulkRNA_moduleScore_singScore(moduleList=moduleList,module_type=module_type,gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T)


  # # Group by Events
  # bulk_mdat$sampleGroup = as.character(bulk_mdat$Subgroup)
  # bulk_mdat$sampleGroup[grepl('^t',bulk_mdat$sampleGroup)] = 'MLL_rearrangement'
  # # bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD' & bulk_mdat$Event == 1] = 'TMD:1'
  # # bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD' & bulk_mdat$Event == 0] = 'TMD:0'
  # # bulk_mdat$sampleGroup = factor(bulk_mdat$sampleGroup,c('TMD','TMD:0','TMD:1','MLDS','AMKL','MLL','MLL_rearrangement',
  # #                                                        'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))
  # # bulk_mdat$Subgroup = factor(bulk_mdat$Subgroup,c('TMD','MLDS','AMKL','MLL',
  # #                                                  'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell',
  # #                                                  't(8;21)','t(10;11)','t(6;11)','t(9;11)'))
  #


  # create a dataframe with the data required: scores and sample group
  allScore$group_facet_hor = allScore$category
  allScore$group_facet_ver = allScore$moduleType
  allScore$group_fill = allScore$sampleGroup
  allScore$group_fill[allScore$group_fill == 'MLL_rearrangement'] = 'MLL'


  plotFun_GATA1s_moduleScore_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){
    # allScore$moduleType = factor(allScore$moduleType,c('GATA1s_Prog_up','GATA1s_MK_up','GATA1s_Mast_up','GATA1s_Ery_up',
    #                                                    'GATA1s_Myeloid_up','GATA1s_Lymphoid_up',
    #                                                    'GATA1s_Mast/MK/Prog_down','GATA1s_EE/MK/Prog_down','GATA1s_Ery_down',
    #                                                    'GATA1s_Myeloid_down','GATA1s_Bcells_down','GATA1s_NK_T_down',
    #                                                    'GATA1s_all_up','GATA1s_all_down','GATA1s_all'))

    p1 = ggplot(allScore, aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.25,lty=2)+
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.25,col='black') +
      geom_quasirandom(aes(col=Event),size=0.2,width = 0.18,alpha=1,col='black')+
      #scale_fill_manual(values = c(rep(col25[4],2),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      scale_fill_manual(values = c(rep(geno_cols[['T21']],2),rep(grey(0.6),7))) +
      #scale_color_manual(values = c('-'=grey(0),'0'='#2D4372','1'=col25[5])) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic()+
      ggtitle(title)+xlab('')+ylab('Module score')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.2),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            strip.text = element_text(size = 11,colour = 'black'),
            axis.ticks = element_line(colour = 'black',linewidth = 0.2),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text = element_text(colour = 'black'))

    p2 = ggplot(allScore, aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.3,lty=2)+
      geom_quasirandom(aes(color=group_fill),size=0.5,width = 0.2,alpha=1)+
      stat_summary(
        aes(group = group_fill), fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", color = "black", width = 0.6, lwd = 0.2,

        # add this bit here to your stat_summary function
        position=position_dodge(width=0.75)
      )+
      scale_color_manual(values = c(rep(geno_cols[['T21']],2),rep(grey(0.6),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      xlab('')+ylab('Module score')+ggtitle(title)+
      theme_classic(base_size = 12)+
      theme(axis.text.x = element_text(size=10,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
            strip.text = element_text(size = 11,colour = 'black'),
            axis.text.y = element_text(size=10,colour = 'black'),
            panel.border = element_rect(fill=F),axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.4),
            strip.background=element_rect(linewidth=0),
            legend.position = 'none',legend.text = element_text(size=9,colour = 'black'))

    print(p1)

  }
  saveFig(file.path(plotDir,'Sup.FigS6x_GATA1s_topGenes_moduleScore_bulkSamples_newBulkMdat_1124'),plotFun_GATA1s_moduleScore_inBulkSamples,rawData=allScore,width = 4.3,height = 7.1,res = 500,useDingbats = F)
  saveFig(file.path(plotDir,'Sup.FigS6x_MLDSimprints.T21FL_topGenes_moduleScore_bulkSamples_newBulkMdat_1124'),plotFun_GATA1s_moduleScore_inBulkSamples,rawData=allScore,width = 4.3,height = 5,res = 500,useDingbats = F)
  #saveFig(file.path(plotDir,paste0('Sup.FigS3x_',module_type,'_moduleScore_bulkSamples')),plotFun_GATA1s_moduleScore_inBulkSamples,rawData=allScore,width = 4.8,height = 3.4,res = 500,useDingbats = F)





}


figxx_MLDS.degs_moduleScore_inBulkSamples = function(){
  library(ggbeeswarm)
  title = 'MLDS module'
  allScore = bulkRNA_moduleScore_singScore(moduleList=moduleList,module_type=module_type,gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T)


  # # Group by Events
  # bulk_mdat$sampleGroup = as.character(bulk_mdat$Subgroup)
  # bulk_mdat$sampleGroup[grepl('^t',bulk_mdat$sampleGroup)] = 'MLL_rearrangement'
  # # bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD' & bulk_mdat$Event == 1] = 'TMD:1'
  # # bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD' & bulk_mdat$Event == 0] = 'TMD:0'
  # # bulk_mdat$sampleGroup = factor(bulk_mdat$sampleGroup,c('TMD','TMD:0','TMD:1','MLDS','AMKL','MLL','MLL_rearrangement',
  # #                                                        'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))
  # # bulk_mdat$Subgroup = factor(bulk_mdat$Subgroup,c('TMD','MLDS','AMKL','MLL',
  # #                                                  'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell',
  # #                                                  't(8;21)','t(10;11)','t(6;11)','t(9;11)'))
  #
  allScore_og=allScore
  allScore$sampleGroup = as.character(allScore$Category)
  allScore$sampleGroup[allScore$sampleGroup == 'TAM_MLDS'] = 'TAM_progressive'

  #allScore$sampleGroup[allScore$sampleGroup == 'TAM' & allScore$Event == 1] = 'TAM:1'
  #allScore$sampleGroup[allScore$sampleGroup == 'TAM' & allScore$Event == 0] = 'TAM:0'
  allScore$sampleGroup = factor(allScore$sampleGroup,c('TAM_conventional','TAM_progressive','TAM_earlyDeath','TAM_unknown','TAM_MLDS','unclassified','MLDS',
                                                       'AMKL','MLL','MLL_rearrangement',
                                                       'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))

  # create a dataframe with the data required: scores and sample group
  allScore$group_facet_hor = allScore$category
  allScore$group_facet_ver = allScore$moduleType
  allScore$group_fill = allScore$sampleGroup
  allScore$group_fill[allScore$group_fill == 'MLL_rearrangement'] = 'MLL'

  dd = allScore[!allScore$Category %in% c('TAM_earlyDeath','TAM_unknown','unclassified'),]

  plotFun_MLDS.degs_moduleScore_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){
    # allScore$moduleType = factor(allScore$moduleType,c('GATA1s_Prog_up','GATA1s_MK_up','GATA1s_Mast_up','GATA1s_Ery_up',
    #                                                    'GATA1s_Myeloid_up','GATA1s_Lymphoid_up',
    #                                                    'GATA1s_Mast/MK/Prog_down','GATA1s_EE/MK/Prog_down','GATA1s_Ery_down',
    #                                                    'GATA1s_Myeloid_down','GATA1s_Bcells_down','GATA1s_NK_T_down',
    #                                                    'GATA1s_all_up','GATA1s_all_down','GATA1s_all'))

    p1 = ggplot(dd[dd$group_fill != 'TAM',], aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.2,lty=2)+
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.2,color='black') +
      geom_quasirandom(aes(col=Event),size=0.25,width = 0.15,alpha=1,col=grey(0))+
      scale_fill_manual(values = c('white',col25[2],rep('white',15))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic()+
      #scale_y_continuous(labels = c(0,'',0.1,'',0.2))+
      ggtitle(title)+xlab('')+ylab('Module score')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.2),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            strip.text = element_text(size = 11,colour = 'black'),
            axis.ticks = element_line(colour = 'black',linewidth = 0.2),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text = element_text(colour = 'black'))



    p2 = ggplot(dd[dd$group_fill != 'TAM',], aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.3,lty=2)+
      geom_quasirandom(aes(color=group_fill),size=0.5,width = 0.2,alpha=1)+
      stat_summary(
        aes(group = group_fill), fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", color = "black", width = 0.6, lwd = 0.2,

        # add this bit here to your stat_summary function
        position=position_dodge(width=0.75)
      )+
      scale_color_manual(values = c(rep(colAlpha(geno_cols[['T21']],1),1),'#c7065a',colAlpha(geno_cols[['T21']],0.7),rep(rep(grey(0.7),7)))) +
      #scale_color_manual(values = c(rep(geno_cols[['T21']],2),rep(grey(0.4),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      xlab('')+ylab('Module score')+ggtitle(title)+
      theme_classic(base_size = 12)+
      theme(axis.text.x = element_text(size=10,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
            strip.text = element_text(size = 11,colour = 'black'),
            axis.text.y = element_text(size=10,colour = 'black'),
            panel.border = element_rect(fill=F),axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.4),
            strip.background=element_rect(linewidth=0),
            legend.position = 'none',legend.text = element_text(size=9,colour = 'black'))


    print(p1)

  }
  saveFig(file.path(plotDir,'Sup.FigSxx_MLDS_topGenes_moduleScore_bulkSamples'),plotFun_MLDS.degs_moduleScore_inBulkSamples,rawData=dd,width = 4.8,height = 4,res = 500)
  saveFig(file.path(plotDir,'Fig3_MLDS_topGenes_moduleScore_bulkSamples_v2_newBulkMdat'),plotFun_MLDS.degs_moduleScore_inBulkSamples,rawData=dd,width = 5.3,height = 6,res = 500)


  ## Perform statistical test:
  dd = read.delim(file.path(plotDir,'Fig3_MLDS_topGenes_moduleScore_bulkSamples_v2_newBulkMdat_rawData.tsv'),sep = '\t')
  dd.sub = dd[dd$moduleType == 'MLDS_topGenes_all' & dd$sampleGroup %in% c('TAM_conventional','TAM_progressive'),]
  wilcox.test(dd.sub$TotalScore[dd.sub$sampleGroup == 'TAM_conventional'],dd.sub$TotalScore[dd.sub$sampleGroup == 'TAM_progressive'],alternative = 'less')


  plotFun_MLDS_moduleScore_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){

    p1 = ggplot(allScore[allScore$moduleType == 'GATA1s_L075_down' & allScore$group_facet_hor != 'Normal',], aes(group_fill, TotalScore)) +
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 0.7,width=0.5,linewidth=0.3) +
      geom_quasirandom(aes(col=Event),size=0.7,width = 0.2,alpha=0.9,col=grey(0))+
      scale_fill_manual(values = c(rep(col25[4],2),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5)))+
      # c(rep(colAlpha(col25[4],0.3),2),colAlpha('#c7065a',0.7),colAlpha(col25[4],0.3),
      #          rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5))) +
      facet_grid(.~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic(base_size = 15)+
      scale_y_continuous(breaks = c(0.1,0.2,0.3))+
      ggtitle(title)+xlab('')+ylab('ML-DS module score')+
      theme(panel.border = element_rect(fill=F),
            axis.text.x = element_text(size=15,angle=90,vjust = 0.5,hjust = 1),
            strip.text = element_text(size = 10),
            axis.line = element_blank(),axis.ticks = element_line(linewidth=0.4),
            strip.background=element_rect(linewidth=0),
            legend.position = 'none')


    print(p1)
  }
  saveFig(file.path(plotDir,'Fig3c_MLDS_moduleScore_bulkSamples'),plotFun_MLDS_moduleScore_inBulkSamples,rawData=allScore,width = 7,height = 5,res = 500,useDingbats = T)
  saveFig(file.path(plotDir,'Fig3c_MLDS_moduleScore_bulkSamples_sub'),plotFun_MLDS_moduleScore_inBulkSamples,rawData=allScore,width = 4.7,height = 4.5,res = 500,useDingbats = F)

}


figxx_resistance.degs_moduleScore_inBulkSamples = function(){
  library(ggbeeswarm)
  title = 'MLDS module'
  allScore = bulkRNA_moduleScore_singScore(moduleList=moduleList,module_type=module_type,gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T)


  # # Group by Events
  # bulk_mdat$sampleGroup = as.character(bulk_mdat$Subgroup)
  # bulk_mdat$sampleGroup[grepl('^t',bulk_mdat$sampleGroup)] = 'MLL_rearrangement'
  # # bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD' & bulk_mdat$Event == 1] = 'TMD:1'
  # # bulk_mdat$sampleGroup[bulk_mdat$sampleGroup == 'TMD' & bulk_mdat$Event == 0] = 'TMD:0'
  # # bulk_mdat$sampleGroup = factor(bulk_mdat$sampleGroup,c('TMD','TMD:0','TMD:1','MLDS','AMKL','MLL','MLL_rearrangement',
  # #                                                        'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))
  # # bulk_mdat$Subgroup = factor(bulk_mdat$Subgroup,c('TMD','MLDS','AMKL','MLL',
  # #                                                  'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell',
  # #                                                  't(8;21)','t(10;11)','t(6;11)','t(9;11)'))
  #
  allScore_og=allScore
  allScore$sampleGroup = as.character(allScore$sampleGroup)
  allScore$sampleGroup[allScore$sampleGroup == 'TAM' & allScore$Event == 1] = 'TAM:1'
  allScore$sampleGroup[allScore$sampleGroup == 'TAM' & allScore$Event == 0] = 'TAM:0'
  allScore$sampleGroup[allScore$sampleGroup == 'MLDS' & allScore$Event == 1] = 'MLDS:1'
  allScore$sampleGroup[allScore$sampleGroup == 'MLDS' & allScore$Event == 0] = 'MLDS:0'
  allScore$sampleGroup = factor(allScore$sampleGroup,c('TAM','TAM:0','TAM:1','MLDS:1','MLDS:0','MLDS','AMKL','MLL','MLL_rearrangement',
                                                       'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))

  # create a dataframe with the data required: scores and sample group
  allScore$group_facet_hor = allScore$category
  allScore$group_facet_ver = allScore$moduleType
  allScore$group_fill = allScore$sampleGroup
  allScore$group_fill[allScore$group_fill == 'MLL_rearrangement'] = 'MLL'


  plotFun_resistance.degs_moduleScore_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){
    # allScore$moduleType = factor(allScore$moduleType,c('GATA1s_Prog_up','GATA1s_MK_up','GATA1s_Mast_up','GATA1s_Ery_up',
    #                                                    'GATA1s_Myeloid_up','GATA1s_Lymphoid_up',
    #                                                    'GATA1s_Mast/MK/Prog_down','GATA1s_EE/MK/Prog_down','GATA1s_Ery_down',
    #                                                    'GATA1s_Myeloid_down','GATA1s_Bcells_down','GATA1s_NK_T_down',
    #                                                    'GATA1s_all_up','GATA1s_all_down','GATA1s_all'))

    p1 = ggplot(allScore[allScore$group_fill != 'TAM',], aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.3,lty=2)+
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 0.7,width=0.5,linewidth=0.3) +
      geom_quasirandom(aes(col=Event),size=0.4,width = 0.15,alpha=0.5,col=grey(0))+
      scale_fill_manual(values = c(rep(col25[4],1),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),2),rep(grey(0.7),5))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic()+
      ggtitle(title)+xlab('')+ylab('Module score')+
      theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))



    p1 = ggplot(allScore[allScore$group_fill != 'TAM',], aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.3,lty=2)+
      geom_quasirandom(aes(color=group_fill),size=0.5,width = 0.2,alpha=1)+
      stat_summary(
        aes(group = group_fill), fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", color = "black", width = 0.6, lwd = 0.2,

        # add this bit here to your stat_summary function
        position=position_dodge(width=0.75)
      )+
      scale_color_manual(values = c(rep(colAlpha(geno_cols[['T21']],1),3),'#c7065a',colAlpha(geno_cols[['T21']],0.7),rep(rep(grey(0.7),7)))) +
      #scale_color_manual(values = c(rep(geno_cols[['T21']],2),rep(grey(0.4),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      xlab('')+ylab('Module score')+ggtitle(title)+
      theme_classic(base_size = 12)+
      theme(axis.text.x = element_text(size=10,angle=90,vjust = 0.5,hjust = 1,colour = 'black'),
            strip.text = element_text(size = 11,colour = 'black'),
            axis.text.y = element_text(size=10,colour = 'black'),
            panel.border = element_rect(fill=F),axis.line = element_blank(),
            axis.ticks = element_line(linewidth=0.4),
            strip.background=element_rect(linewidth=0),
            legend.position = 'none',legend.text = element_text(size=9,colour = 'black'))


    print(p1)

  }
  saveFig(file.path(plotDir,'Sup.FigSxx_MLDS_topGenes_moduleScore_bulkSamples'),plotFun_MLDS.degs_moduleScore_inBulkSamples,rawData=allScore,width = 4.8,height = 4,res = 500)
  saveFig(file.path(plotDir,'Sup.FigSxx_MLDS_topGenes_moduleScore_bulkSamples_v2'),plotFun_MLDS.degs_moduleScore_inBulkSamples,rawData=allScore,width = 4,height = 3.3,res = 500)







  plotFun_MLDS_moduleScore_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){

    p1 = ggplot(allScore[allScore$moduleType == 'GATA1s_L075_down' & allScore$group_facet_hor != 'Normal',], aes(group_fill, TotalScore)) +
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 0.7,width=0.5,linewidth=0.3) +
      geom_quasirandom(aes(col=Event),size=0.7,width = 0.2,alpha=0.9,col=grey(0))+
      scale_fill_manual(values = c(rep(col25[4],2),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5)))+
      # c(rep(colAlpha(col25[4],0.3),2),colAlpha('#c7065a',0.7),colAlpha(col25[4],0.3),
      #          rep(colAlpha(col25[1],0.4),3),rep(grey(0.7),5))) +
      facet_grid(.~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic(base_size = 15)+
      scale_y_continuous(breaks = c(0.1,0.2,0.3))+
      ggtitle(title)+xlab('')+ylab('ML-DS module score')+
      theme(panel.border = element_rect(fill=F),
            axis.text.x = element_text(size=15,angle=90,vjust = 0.5,hjust = 1),
            strip.text = element_text(size = 10),
            axis.line = element_blank(),axis.ticks = element_line(linewidth=0.4),
            strip.background=element_rect(linewidth=0),
            legend.position = 'none')


    print(p1)
  }
  saveFig(file.path(plotDir,'Fig3c_MLDS_moduleScore_bulkSamples'),plotFun_MLDS_moduleScore_inBulkSamples,rawData=allScore,width = 7,height = 5,res = 500,useDingbats = T)
  saveFig(file.path(plotDir,'Fig3c_MLDS_moduleScore_bulkSamples_sub'),plotFun_MLDS_moduleScore_inBulkSamples,rawData=allScore,width = 4.7,height = 4.5,res = 500,useDingbats = T)

}

figxx_badTAM.degs_moduleScore_inBulkSamples = function(){
  library(ggbeeswarm)
  title = 'bad TAM'
  allScore = bulkRNA_moduleScore_singScore(moduleList=moduleList,module_type=module_type,gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T)


  allScore_og=allScore
  allScore$sampleGroup = as.character(allScore$sampleGroup)
  allScore$sampleGroup[allScore$sampleGroup == 'TAM' & allScore$Event == 1] = 'TAM:1'
  allScore$sampleGroup[allScore$sampleGroup == 'TAM' & allScore$Event == 0] = 'TAM:0'
  allScore$sampleGroup = factor(allScore$sampleGroup,c('TAM','TAM:0','TAM:1','MLDS','AMKL','MLL','MLL_rearrangement',
                                                       'HSC','MPP','MEP','CMP','GMP','LMPP','ELP','PreProB','ProB','Bcell'))

  # create a dataframe with the data required: scores and sample group
  allScore$group_facet_hor = allScore$category
  allScore$group_facet_ver = allScore$moduleType
  allScore$group_fill = allScore$sampleGroup
  allScore$group_fill[allScore$group_fill == 'MLL_rearrangement'] = 'MLL'


  plotFun_badTAM.degs_moduleScore_inBulkSamples = function(noFrame=FALSE,noPlot=FALSE){

    p1 = ggplot(allScore[allScore$group_fill != 'TAM',], aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.3,lty=2)+
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 0.7,width=0.5,linewidth=0.3) +
      geom_quasirandom(aes(col=Event),size=0.4,width = 0.15,alpha=0.5,col=grey(0))+
      scale_fill_manual(values = c(rep(col25[4],1),'#c7065a',col25[4],rep(colAlpha(col25[1],0.4),2),rep(grey(0.7),5))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic()+
      ggtitle(title)+xlab('')+ylab('Module score')+
      theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))

    print(p1)

  }
  saveFig(file.path(plotDir,'Sup.FigSxx_badTAM_topGenes_moduleScore_bulkSamples'),plotFun_badTAM.degs_moduleScore_inBulkSamples,rawData=allScore,width = 5,height = 4,res = 500)
}






## Heatmap expression of these mlds genes
bulk_mdat$Event = as.character(bulk_mdat$Event)
bulk_mdat$Event[is.na(bulk_mdat$Event)] = '-'
mtx = tpm_count[rownames(tpm_count) %in% mlds_markers$ensID,bulk_mdat$Sample[bulk_mdat$Subgroup %in% c('TMD')]]
mtx = tpm_count[rownames(tpm_count) %in% good.vs.bad.TAM_deg$ensID[good.vs.bad.TAM_deg$is_topDEG==T & good.vs.bad.TAM_deg$direction == 'good_down'],bulk_mdat$Sample[bulk_mdat$sampleGroup %in% c('TMD:0','TMD:1')]]
mtx = logcounts[rownames(logcounts) %in% good.vs.bad.TAM_deg$ensID[good.vs.bad.TAM_deg$is_topDEG==T & good.vs.bad.TAM_deg$direction == 'good_down'],bulk_mdat$Sample[bulk_mdat$sampleGroup %in% c('TMD:0','TMD:1')]]

mtx = tpm_count[rownames(tpm_count) %in% good.vs.bad.TAM_deg$ensID[good.vs.bad.TAM_deg$is_topDEG==T],bulk_mdat$Sample[bulk_mdat$Subgroup %in% c('TMD')]]
mtx = tpm_count[rownames(tpm_count) %in% df$ensID,bulk_mdat$Sample[bulk_mdat$Subgroup %in% c('TMD','MLDS') & !grepl('PDX',bulk_mdat$Sample)]]
mtx = logc[rownames(tpm_count) %in% mlds_markers$ensID[mlds_markers$direction == 'L075_down'],bulk_mdat$Sample[bulk_mdat$sampleGroup %in% c('TMD:0','TMD:1')]]

mtx = bulk_ranked[rownames(bulk_ranked) %in% mlds_markers$ensID[mlds_markers$direction == 'L075_down'],colnames(bulk_ranked) %in% bulk_mdat$Sample[bulk_mdat$sampleGroup %in% c('TMD:0','TMD:1','MLDS') & bulk_mdat$Sample != 'Patient_TMD_1701']]
#rownames(mtx) = mlds_markers$geneSym[match(rownames(mtx),mlds_markers$ensID)]
rownames(mtx) = geneMap$geneSym[match(rownames(mtx),geneMap$ensID)]
botAnno = HeatmapAnnotation(event = bulk_mdat$Event[match(colnames(mtx),bulk_mdat$Sample)],
                            col = list(event = c('-'=grey(0.8),'NA'=grey(0.8),'1'='red','0'='blue')))
library(circlize)
col_fun = colorRamp2(c(-3, 0, 2), c('#1a4a87','white','#a4282c'))
hm = Heatmap(t((t(mtx))),
             #col = col_fun,
             show_column_dend = T,show_row_dend = F,
             row_names_gp = gpar(fontsize=7),
             column_title_rot = 90,bottom_annotation = botAnno,
             cluster_rows = T,
             column_split = bulk_mdat$sampleGroup[match(colnames(mtx),bulk_mdat$Sample)],
             #split = mlds_markers$direction[match(rownames(mtx),mlds_markers$geneSym)],
             split = good.vs.bad.TAM_deg$direction[match(rownames(mtx),good.vs.bad.TAM_deg$geneSym)]
             #split = df$direction[match(rownames(mtx),df$geneSym)]
             #km=4

)

ht = draw(hm)
## Extract bad TAM-specific genes
badTAM_markers_bulk = rownames(mtx)[c(row_order(ht)[['2']],
                                      row_order(ht)[['4']])]
moduleList = list('badTAM' = geneMap$ensID[geneMap$geneSym %in% badTAM_markers_bulk])

moduleList = split(df$ensID,df$direction)

# plot_singScore(adultUP_scores,annot = bulk_samples,group_col = 'cancerType',group_col2 = 'source',id_col = 'sampleID',title = 'upDEGs_AdultSpecific')
# plot_singScore(thyUP_scores,annot = bulk_samples,group_col = 'cancerType',group_col2 = 'source',id_col = 'sampleID',title = 'upDEGs_ThySpecific')
# plot_singScore(mesUP_scores,annot = bulk_samples,group_col = 'cancerType',group_col2 = 'source',id_col = 'sampleID',title = 'upDEGs_MesSpecific')
# plot_singScore(foetalUPsub_scores,annot = bulk_samples,group_col = 'group.sub',group_col2 = 'group.big',id_col = 'sampleID',title = 'upDEGs_FoetalSpecific.sub')

# bulk_samples$group.big = bulk_samples$group



plot_singScore = function(singScore_result,annot,group_col,group_col2,id_col,relSize = 0.9,title='singScore'){
  #relative size of text in the figure
  #relSize = 0.9

  #create annotation
  #bulk_samples$cancerType
  p1 = plotDispersion(singScore_result, annot = annot[[group_col]][match(rownames(singScore_result),annot[[id_col]])], textSize = relSize,size = 0.6)
  print(p1)

  library(reshape2)

  # ggplot theme
  rl = 0.9
  current_theme = theme_minimal() +
    theme(
      panel.border = element_rect(colour = 'black', fill = NA),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = rel(rl) * 1.1),
      axis.text = element_text(size = rel(rl)),
      plot.title = element_text(size = rel(rl)),
      strip.background = element_rect(fill = NA, colour = 'black'),
      strip.text = element_text(size = rel(rl)),
      legend.text = element_text(size = rel(rl)),
      legend.title = element_text(size = rel(rl), face = 'italic'),
      legend.position = 'bottom',
      legend.direction = 'horizontal'
    )

  #create a dataframe with the data required: scores and sample group
  scoredf = annot
  scoredf$group = scoredf[[group_col]]
  scoredf$group2 = scoredf[[group_col2]]
  scoredf$topDEGs_singScore = singScore_result$TotalScore[match(scoredf[[id_col]],rownames(singScore_result))]
  scoredf$topDEGs_up_singScore = singScore_result$UpScore[match(scoredf[[id_col]],rownames(singScore_result))]
  scoredf$topDEGs_down_singScore = singScore_result$DownScore[match(scoredf[[id_col]],rownames(singScore_result))]
  scoredf$Dispersion = singScore_result$TotalDispersion[match(scoredf[[id_col]],rownames(singScore_result))]

  p1 = ggplot(scoredf, aes(group, topDEGs_singScore,fill=group2)) +
    geom_boxplot(position = 'dodge', alpha = 0.9) +
    scale_fill_manual(values = c(col25,pal34H)) +
    facet_wrap(vars(group2),scales = 'free_x',nrow=1)+
    current_theme +
    #geom_hline(yintercept = 0.2,col='red')+
    ggtitle(title)+
    theme(axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1))

  print(p1)

}




#
# ## Tradeseq plot - not really working....
# df = data[!grepl('L041',data$donorID),]
# df$mk2ee = df$MK_UCell/df$EE_UCell
# df$mast2ee = df$Mast_UCell/df$EE_UCell
# table(df$dataset_ID)
# df$group = ifelse(df$dataset_ID %in% c('fLiver','fBM'),paste0(df$dataset_ID,':',df$Genotype),paste0(df$dataset_ID,':',df$donorID))
# ggplot(df[df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Cancer' ,],aes(mast2ee,EE_UCell))+
#   geom_point(size=0.001)+
#   facet_grid(group~finalAnn_broad)+
#   theme_classic()+
#   theme(panel.border = element_rect(fill=F))
#
# ggplot(df[df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Cancer' ,],aes(group,EE_UCell,fill=dataset_ID))+
#   geom_boxplot(outlier.size = 0.001)+
#   scale_y_log10()+
#   #facet_grid(group~finalAnn_broad)+
#   theme_classic()+
#   theme(panel.border = element_rect(fill=F))
#
# ggplot(df[df$dataset_ID %in% c('fLiver') & df$finalAnn_broad %in% c('MEMP_MEP','MK','HSC_MPP'),],aes(group,HSC_UCell,fill=finalAnn_broad))+
#   geom_boxplot(outlier.size = 0.001)+
#   scale_y_log10()+
#   #facet_grid(group~finalAnn_broad)+
#   theme_classic()+
#   theme(panel.border = element_rect(fill=F))
#





















#
# ## HBB plot
# df = data[!grepl('L041',data$donorID),]
# table(df$dataset_ID)
# df$group = ifelse(df$dataset_ID %in% c('fLiver','fBM'),paste0(df$dataset_ID,':',df$Genotype),paste0(df$dataset_ID,':',df$donorID))
# ggplot(df[df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Cancer' ,],aes(mast2ee,EE_UCell))+
#   geom_point(size=0.001)+
#   facet_grid(group~finalAnn_broad)+
#   theme_classic()+
#   theme(panel.border = element_rect(fill=F))
#
# ggplot(df[df$dataset_ID %in% c('MLDS','TAM') & df$finalAnn_broad == 'Cancer' ,],aes(group,HBZ_UCell,fill=dataset_ID))+
#   geom_boxplot(outlier.size = 0.001)+
#   scale_y_log10()+
#   #facet_grid(group~finalAnn_broad)+
#   theme_classic()+
#   theme(panel.border = element_rect(fill=F))
#
# ggplot(df[df$dataset_ID %in% c('fLiver') & df$finalAnn_broad %in% c('MEMP_MEP','MK','HSC_MPP','EE'),],aes(group,HBZ_UCell,fill=finalAnn_broad))+
#   geom_boxplot(outlier.size = 0.001)+
#   scale_y_log10()+
#   #facet_grid(group~finalAnn_broad)+
#   theme_classic()+
#   theme(panel.border = element_rect(fill=F))
#
#
#






