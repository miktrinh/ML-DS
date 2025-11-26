## Generate figures for Revision (NatComms) July/August 2025  ##


##----    Set working directory  -----##
setwd('~/ML-DS/')

#------------------------#
##      Libraries     ####
#------------------------#
library(tidyverse)
library(plotrix)
library(Seurat)
library(readxl)
library(alleleIntegrator)
library(RColorBrewer)
library(circlize)
source('R/utils/misc.R')

plotDir='~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots'
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}



mlds_srat_fp = '~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns.RDS'
mlds_mdat_fp = '~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns_mdat_2508.csv'
# mdat = cbind(mlds@meta.data,mlds@reductions$umap@cell.embeddings)
# write.csv(mdat,'~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns_mdat_2508.csv')

fig1b_MLDScohort = function(){
  ## MLDS dataset ##
  mdat = read.csv(mlds_mdat_fp,row.names = 1) 
  mdat$annot = mdat$annot_aug24_new
  mdat$finalAnn_broad = mdat$annot_aug24_new
  mdat$broadLineage[mdat$finalAnn_broad == 'Tumour' & mdat$broadLineage != 'Tumour'] = 'Tumour'
  mdat$broadLineage[mdat$finalAnn_broad %in% c('Tum_MK?','Tumour_WT')] = 'Tumour_unsure'
  #mdat = mdat[!mdat$broadLineage %in% c('doublets','unsure_others','others','Tumour?','Tumour_unsure','lowQual'),]
  
  if(!keep_CC3){
    # Remove CC3 and L041 Diagnostic samples
    mdat = mdat[!(mdat$donorID %in% c('CC3','L041') & mdat$timePoint == 'Diagnostic'),]
  }
  
  
  mdat$disease = factor(mdat$disease,rev(c('TAM','MLDS')))
  mdat$clinicalOutcome[mdat$donorID %in% c('L038')] = 'Refractory'
  mdat$clinicalOutcome[mdat$donorID %in% c('L156','L076','CC3')] = 'Relapse'
  mdat$clinicalOutcome[mdat$donorID %in% c('L019','L039','L040','L042','L091','L041','CC4','CC5','CC1','CC2','L178',
                                           'L075','L114','L182','CC6','CC7','CC8')] = 'Remission'
  mdat$clinicalOutcome[mdat$clinicalOutcome == '?'] = 'unknown'
  table(mdat$clinicalOutcome,mdat$donorID)
  
  mdat$donorID = factor(mdat$donorID,c('CC3','CC4','CC5',
                                       'L019','L039','L040','L041','L042','L091','L178','L076','L038',
                                       'CC1','CC2','CC6','CC7','CC8','L075','L114','L182','L156'))
  #mdat$tissue[mdat$orig.ident == 'Leuk13697519'] = 'Blood'
  
  mdat$timePoint = as.character(mdat$timePoint)
  mdat$timePoint[mdat$donorID == 'L156' & mdat$orig.ident %in% c('MY.200531.14784986','MY.200531.14784987')] = 'Recurrent'
  mdat$timePoint[mdat$timePoint == 'Diagnostic' & mdat$donorID == 'L076' & mdat$tissue == 'Blood'] = 'D (PBMC)'
  mdat$timePoint[mdat$timePoint == 'Diagnostic' & mdat$disease == 'TAM'] = 'D (PBMC)'
  mdat$timePoint[mdat$timePoint == 'Diagnostic'] = 'D (BM)'
  mdat$timePoint = factor(mdat$timePoint,(c('D (BM)','D (PBMC)','TP1','TP2','TP4','D.Relapse','D.Relapse2','Recurrent')))
  
  
  ## Save the annotation in shared folder
  # annot = mdat[,c('cellID',colnames(mdat)[c(1:19,26:29,43,30:35,37,45,55,56,57:62,53,65,42,63,64)])]
  # write.csv(annot,'/lustre/scratch126/casm/team274sb/project_folders/GOSH_Leuk/sample_metadata/TAM_MLDS_metadata_mt22_2410.csv')
  
  mdat$group = ifelse(mdat$broadLineage == 'Tumour','Tumour',ifelse(mdat$broadLineage %in% c('Tumour_unsure'),'Tumour_unsure','normal'))
  nLeuk = as.data.frame(table(mdat$donorID,mdat$orig.ident,mdat$group))
  colnames(nLeuk) = c('donorID','channelID','group','Freq')
  nLeuk = nLeuk[nLeuk$Freq >0,]
  
  #mdat = mdat[mdat$group != 'lowQual',]
  mlds_dataset = as.data.frame(table(mdat$donorID,mdat$Genotype,mdat$sex,mdat$age_yrs,mdat$disease,mdat$orig.ident,mdat$timePoint,mdat$tissue,mdat$clinicalOutcome,mdat$blastPerc))
  colnames(mlds_dataset) = c('donorID','Genotype','Sex','Age (year)','Disease','channelID','timePoint','Tissue','Clinical outcome','blastPerc','nCell')
  mlds_dataset = mlds_dataset[order(mlds_dataset$donorID),]
  
  
  
  
  mlds_dataset2 = mlds_dataset[mlds_dataset$nCell > 0,]
  mlds_dataset2$Sex = as.character(mlds_dataset2$Sex)
  mlds_dataset2$Sex[mlds_dataset2$Sex == 'F'] = 'Female'
  mlds_dataset2$Sex[mlds_dataset2$Sex == 'M'] = 'Male'
  mlds_dataset2$timePoint = gsub(' (.*)','',mlds_dataset2$timePoint)
  mlds_dataset2$timePoint = gsub('TP','Timepoint ',mlds_dataset2$timePoint)
  mlds_dataset2$timePoint[mlds_dataset2$timePoint == 'D'] = 'Diagnostic'
  mlds_dataset2$timePoint[mlds_dataset2$timePoint == 'D.Relapse'] = 'Relapse 1 diagnostic'
  mlds_dataset2$timePoint[mlds_dataset2$timePoint == 'D.Relapse2'] = 'Relapse 2 diagnostic'
  mlds_dataset2$Tissue = as.character(mlds_dataset2$Tissue)
  mlds_dataset2$Tissue[mlds_dataset2$Tissue == 'BM'] = 'Bone marrow'
  mlds_dataset2$Tissue[mlds_dataset2$Tissue == 'Blood'] = 'Peripheral blood'
  
  leuk = nLeuk[nLeuk$group == 'Tumour',]
  norm = nLeuk[nLeuk$group == 'normal',]
  mlds_dataset2$nLeuk = leuk$Freq[match(mlds_dataset2$channelID,leuk$channelID)]
  mlds_dataset2$nLeuk[is.na(mlds_dataset2$nLeuk)] = 0
  mlds_dataset2$nNorm = norm$Freq[match(mlds_dataset2$channelID,norm$channelID)]
  mlds_dataset2$nNorm[is.na(mlds_dataset2$nNorm)] = 0
  all((mlds_dataset2$nNorm + mlds_dataset2$nLeuk) == mlds_dataset2$nCell)
  which((mlds_dataset2$nNorm + mlds_dataset2$nLeuk) != mlds_dataset2$nCell)
  table(mlds_dataset2$donorID)
  
  
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$timePoint,decreasing = F),]
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$donorID),]
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$Disease),]
  mlds_dataset2$blastPerc = as.character(mlds_dataset2$blastPerc)
  mlds_dataset2$blastPerc[mlds_dataset2$blastPerc != '?'] = as.numeric(as.character(mlds_dataset2$blastPerc[mlds_dataset2$blastPerc != '?']))
  sum(mlds_dataset2$nLeuk)
  sum(mlds_dataset2$nNorm)
  
  ## Change Genotype
  mlds_dataset2$Genotype = as.character(mlds_dataset2$Genotype)
  mlds_dataset2$Genotype[mlds_dataset2$Genotype == 'T21'] = 'Trisomy 21'
  mlds_dataset2$Disease = as.character(mlds_dataset2$Disease)
  mlds_dataset2$Disease[mlds_dataset2$Disease == 'MLDS'] = 'ML-DS'
  mlds_dataset2$Disease = factor(mlds_dataset2$Disease,c('TAM','ML-DS'))
  
  ## Add PDID
  projectMani = read_excel('~/projectManifest.xlsx','MLDS_GOSH')
  projectMani = projectMani[!is.na(projectMani$donorID) & !is.na(projectMani$PDID) & projectMani$Disease %in% c('MLDS','TAM') & projectMani$assay != 't-NanoSeq',]
  projectMani$Disease[projectMani$Disease == 'MLDS'] = 'ML-DS'
  projectMani$`Point in treatment`=gsub('postChemo ','',projectMani$`Point in treatment`)
  projectMani$`Point in treatment` = gsub('TP','Timepoint ',projectMani$`Point in treatment`)
  projectMani$`Point in treatment`[projectMani$`Point in treatment` == 'Relapse'] = 'Relapse 1 diagnostic'
  projectMani$`Point in treatment`[projectMani$`Point in treatment` == 'Relapse2_D0'] = 'Relapse 2 diagnostic'
  projectMani$`Point in treatment`[projectMani$`Point in treatment` == 'Relapse 1 diagnostic' & projectMani$donorID == 'L156'] = 'Recurrent'
  projectMani$sampleID2 = paste0(projectMani$donorID,'_',projectMani$Disease,'_',projectMani$`Point in treatment`)
  projectMani$sampleID2[projectMani$sampleID2 == 'L091_MLDS_NA'] = 'L091_MLDS_Diagnostic'
  projectMani$sampleID2[projectMani$sampleID2 == 'L076_MLDS_Diagnostic' & projectMani$Tissue == 'Blood'] = 'L076_MLDS_Diagnostic_blood'
  
  mlds_dataset2$sampleID2 = paste0(mlds_dataset2$donorID,'_',mlds_dataset2$Disease,'_',mlds_dataset2$timePoint)
  mlds_dataset2$sampleID2[mlds_dataset2$sampleID2 == 'L076_MLDS_Diagnostic' & mlds_dataset2$Tissue == 'Peripheral blood'] = 'L076_MLDS_Diagnostic_blood'
  table(mlds_dataset2$sampleID2 %in% projectMani$sampleID2)
  table(projectMani$sampleID2 %in% mlds_dataset2$sampleID2)
  projectMani$sampleID2[!projectMani$sampleID2 %in% mlds_dataset2$sampleID2]
  table(mlds_dataset2$sampleID2[!mlds_dataset2$sampleID2 %in% projectMani$sampleID2])
  
  mlds_dataset2$PDID = projectMani$PDID[match(mlds_dataset2$sampleID2,projectMani$sampleID2)]
  mlds_dataset2$PDID[mlds_dataset2$PDID == 'no materials for WGS - Checked with Conor' & !is.na(mlds_dataset2$PDID)] = '-'
  mlds_dataset2 = mlds_dataset2[,colnames(mlds_dataset2) != 'sampleID2']
  
  ## Reorder samples
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$donorID),]
  mlds_dataset2$`Clinical outcome` = factor(mlds_dataset2$`Clinical outcome`,c('Remission','Relapse','Refractory'))
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$`Clinical outcome`),]
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$Disease),]
  
  mlds_dataset_by_donor = mlds_dataset2[,c('donorID','Genotype','Sex','Age (year)','Disease','timePoint','Tissue','Clinical outcome','blastPerc')]
  mlds_dataset_by_donor = mlds_dataset_by_donor[!duplicated(mlds_dataset_by_donor),]
  
  write.csv(mlds_dataset2,file.path(plotDir,'..','TableS1_MLDS_dataset.csv'),row.names = F)
  write.csv(mlds_dataset_by_donor,file.path(plotDir,'..','TableS1_MLDS_dataset_by_donor.csv'),row.names = F)
  
}

## Figure X: GATA1 status for L038 ##
figX_MLDS_EryCluster_GATA1s_status = function(){
  if(file.exists(file.path(plotDir,'Fig2c_MLDS_GATA1s_UMAP_rawData.tsv')) & skipIfExists){
    df = read.delim(file.path(plotDir,'Fig2c_MLDS_GATA1s_UMAP_rawData.tsv'),sep = '\t',header = T)
  }else{
    
    # Import the seurat object
    mdat = read.csv(mlds_mdat_fp,row.names = 1)

    dd = mdat[,c("cellID","donorID",'disease','annot',"broadLineage",'timePoint','tissue','GATA1_status','UMAP_1','UMAP_2','Phase')]
    dd = dd[!dd$broadLineage %in% c('doublets','unsure_others','others','lowQual'),]
    
    dd$GATA1_status[dd$GATA1_status %in% c('WT','GATA1s_WT')] = 'GATA1 wild type'
    dd$GATA1_status[dd$GATA1_status %in% c('Mut','GATA1s_mutant')] = 'GATA1s mutation'
    dd$GATA1_status[dd$GATA1_status %in% c('unsure','GATA1s_unsure')] = 'Not informative'
    dd$GATA1_status[dd$GATA1_status %in% c('noCov','uninformative')] = 'Not informative'#'No coverage at mutation site'
    dd$GATA1_status[dd$GATA1_status %in% c('noGATA1expr','no_GATA1_expr')] = 'No GATA1 expression'
    
    # dd$GATA1_status2 = dd$GATA1_status
    # dd$GATA1s_status2[dd$GATA1s_status2 %in% c('Not informative','No coverage at mutation site','uninformative')] = 'Uninformative'
    
    d = dd[dd$broadLineage == 'Erythroblasts' | (dd$annot == 'Tumour' & dd$donorID == 'L038'),]
    d = d[d$UMAP_1 > 3 & d$UMAP_2 < 10 & d$UMAP_2 > 3,]
    
  }
  
  plotFun_MLDS_EryCluster_byPDID_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    ccs_donorID = c('L156' = col25[4],
                    'L076' = col25[21],
                    'L038' = col25[23],
                    'L178' = col25[7],
                    'L019' = col25[12],
                    'CC8' = col25[5],
                    'L091' = col25[8],
                    'L040' = col25[6],
                    'CC5' = col25[9],
                    'CC1' = col25[3],
                    'L041' = col25[24],
                    'L042' = col25[18],
                    'CC2' = col25[10],
                    'CC6' = col25[17],
                    'L182'=col25[22],
                    'L039' = col25[16],
                    'L114' = col25[15],
                    'L075' = col25[1],
                    'CC4' = col25[20],
                    'CC7' = col25[13])
    cols = sapply(1:nrow(d),function(i){
      ifelse(d$donorID[i] == 'L038' & d$annot[i] == 'Tumour',ccs_donorID[d$donorID[i]],
             ifelse(d$donorID[i] == 'L038' & d$annot[i] != 'Tumour',colAlpha(ccs_donorID[d$donorID[i]],0.8),
                    colAlpha(ccs_donorID[d$donorID[i]],alphas = 0.3)))})
    
    plot(d$UMAP_1,d$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','TAM - MLDS'),
         frame.plot=F)
    
    if(!noPlot){
      points(d$UMAP_1,d$UMAP_2,
             col = cols,
             pch = 19,
             cex=0.08)
    }
  }
  
  saveFig(file.path(plotDir,'SuppFigXX_MLDS_EryCluster_donorID_UMAP'),plotFun_MLDS_EryCluster_byPDID_UMAP,rawData=d,width = 2.3,height = 1.5,res = 500)
  
  plotFun_MLDS_EryCluster_L038_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    d$cols = colAlpha(grey(0.8),0.5)
    d$cols[d$donorID == 'L038' & d$annot == 'Tumour' & d$timePoint == 'TP1'] = "black"#col25[2]
    d$cols[d$donorID == 'L038' & d$annot != 'Tumour' & d$timePoint == 'TP1'] = '#C16A36'
    #d$cols[d$donorID == 'L038' & d$annot == 'Tumour' & d$timePoint != 'TP1'] = "#C16A36"
    #d$cols[d$donorID == 'L038' & d$annot != 'Tumour' & d$timePoint != 'TP1'] = "#b5559b"
    #d$cols[d$donorID != 'L038' ] = 
    
    
    
    plot(d$UMAP_1,d$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','TAM - MLDS'),
         frame.plot=F)
    
    if(!noPlot){
      points(d$UMAP_1,d$UMAP_2,
             col = d$cols,
             pch = 19,
             cex=0.01)
    }
  }
  
  saveFig(file.path(plotDir,'SuppFigXX_MLDS_EryCluster_L038_UMAP'),plotFun_MLDS_EryCluster_L038_UMAP,rawData=d,width = 3.3,height = 2,res = 500)
  
  
  
  plotFun_MLDS_EryCluster_GATA1_status_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    gata1_status_ccs = c("GATA1 wild type"=colAlpha('#0666b5',0.5),
                         "GATA1s mutation"='#A92821',
                         "No GATA1 expression"=colAlpha(grey(0.15),0.3),
                         "Not informative"=colAlpha(grey(0.8),0.3),
                         "No coverage at mutation site"=colAlpha('#a17f45',1))
            
    plot(d$UMAP_1,d$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','TAM - MLDS'),
         frame.plot=F)
    
    if(!noPlot){
      
      points(d$UMAP_1,d$UMAP_2,
             col = gata1_status_ccs[d$GATA1_status],
             pch = 19,
             cex=0.01)
      
      d_gata1_mut = d[d$GATA1_status == 'GATA1s mutation',]
      points(d_gata1_mut$UMAP_1,d_gata1_mut$UMAP_2,
             col = gata1_status_ccs[d_gata1_mut$GATA1_status],
             pch = 19,
             cex=0.01)
      
    }
  }
  
  saveFig(file.path(plotDir,'SuppFigXX_MLDS_EryCluster_GATA1_status_UMAP'),plotFun_MLDS_EryCluster_GATA1_status_UMAP,rawData=d,width = 3.3,height = 2,res = 500)
  
  plotFun_MLDS_EryCluster_Phase_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    phase_ccs = c("G1"=colAlpha(grey(0.8),0.8),
                         "G2M"=colAlpha(grey(0.5),0.5),
                         "S"=colAlpha(grey(0.1),0.5))
    
    plot(d$UMAP_1,d$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','TAM - MLDS'),
         frame.plot=F)
    
    if(!noPlot){
      
      points(d$UMAP_1,d$UMAP_2,
             col = phase_ccs[d$Phase],
             pch = 19,
             cex=0.01)
      
    }
  }
  
  saveFig(file.path(plotDir,'SuppFigXX_MLDS_EryCluster_Phase_UMAP'),plotFun_MLDS_EryCluster_Phase_UMAP,rawData=d,width = 3.3,height = 2,res = 500)
  
  
  ggplot(d,aes(UMAP_1,UMAP_2))+
    #geom_point(aes(col=donorID),size=0.01)+
    scale_color_manual(values = as.character(ccs_donorID))+
    geom_point(aes(col=Phase),size=1)+
    #geom_point(data = d[d$donorID == 'L038' & d$annot == 'Tumour' & d$timePoint != 'Diagnostic',],aes(col=GATA1_status),size=0.1)+
    #geom_point(data = d[d$donorID == 'L038' & d$annot != 'Tumour',],size=0.01,col=grey(0.8))+
    theme_classic()
  
  ## Look at distribution of Cell Cycle state aross donors 
  ## --> L038 and L076 indeed have highest fraction of cells in G1
  a = mdat %>% filter(annot=='Tumour') %>% group_by(Phase,disease,donorID,timePoint) %>% 
    summarise(nCell = n()) %>% group_by(disease,donorID,timePoint) %>% mutate(totalCell = sum(nCell),
                                                                              frac = nCell/totalCell)
  ggplot(a,aes(timePoint,frac,fill=Phase))+
    geom_col()+facet_grid(.~donorID + disease,scales = 'free_x',space = 'free_x')+theme_bw()+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
}





## Replot figure 3D to highlight AMKL rearranged cases
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
  
  
  
  ## Perform statistical test:
  #dd = read.delim(file.path(plotDir,'Fig3_MLDS_topGenes_moduleScore_bulkSamples_v2_newBulkMdat_rawData.tsv'),sep = '\t')
  dd.sub = dd[dd$moduleType == 'MLDS_topGenes_all' & dd$sampleGroup %in% c('TAM_conventional','TAM_progressive'),]
  wilcox.test(dd.sub$TotalScore[dd.sub$sampleGroup == 'TAM_conventional'],dd.sub$TotalScore[dd.sub$sampleGroup == 'TAM_progressive'],alternative = 'less')
  t.test(dd.sub$TotalScore[dd.sub$sampleGroup == 'TAM_conventional'],dd.sub$TotalScore[dd.sub$sampleGroup == 'TAM_progressive'],alternative = 'less')
  
  
  
  
  plotFun_MLDS.degs_moduleScore_inBulkSamples_AMKLfusion = function(noFrame=FALSE,noPlot=FALSE){
    # allScore$moduleType = factor(allScore$moduleType,c('GATA1s_Prog_up','GATA1s_MK_up','GATA1s_Mast_up','GATA1s_Ery_up',
    #                                                    'GATA1s_Myeloid_up','GATA1s_Lymphoid_up',
    #                                                    'GATA1s_Mast/MK/Prog_down','GATA1s_EE/MK/Prog_down','GATA1s_Ery_down',
    #                                                    'GATA1s_Myeloid_down','GATA1s_Bcells_down','GATA1s_NK_T_down',
    #                                                    'GATA1s_all_up','GATA1s_all_down','GATA1s_all'))
    dd = dd[dd$group_fill != 'TAM' & dd$group_facet_ver == 'MLDS_topGenes_all',]
    dd$amkl_group = ifelse(dd$Sample %in% c('Patient_AMKL_3806','Patient_AMKL_3865'),'amkl_fusion','-')
    
    p1 = ggplot(dd, aes(group_fill, TotalScore)) +
      geom_hline(yintercept = 0,lwd=0.2,lty=2)+
      geom_boxplot(aes(fill=group_fill),outlier.colour = 'white',position = 'dodge', alpha = 1,width=0.5,linewidth=0.2,color='black') +
      geom_quasirandom(aes(col=amkl_group,size=amkl_group,shape=amkl_group),width = 0.15,alpha=0.8)+
      scale_size_manual(values = c('amkl_fusion'=1.2,'-'=0.1))+
      scale_shape_manual(values = c('amkl_fusion'=17,'-'=19))+
      scale_colour_manual(values = c('amkl_fusion'='#0F6FFF','-'='black'))+
      scale_fill_manual(values = c('white',col25[2],rep('white',15))) +
      #scale_fill_manual(values = c(rep(col25[4],2),rep(grey(0.7),7))) +
      facet_grid(group_facet_ver~group_facet_hor,scales = 'free',space = 'free_x')+
      theme_classic()+
      #scale_y_continuous(labels = c(0,'',0.1,'',0.2))+
      ggtitle(title)+xlab('')+ylab('Module score')+
      theme(panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.2),axis.line = element_blank(),
            strip.background=element_rect(linewidth=0),
            strip.text = element_text(size = 9,colour = 'black'),
            axis.ticks = element_line(colour = 'black',linewidth = 0.2),
            axis.text.x = element_text(size = 8,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
            legend.position = 'none',
            axis.text = element_text(colour = 'black')
            )
    
    print(p1)
    
  }
  saveFig(file.path(plotDir,'SuppFigXX_MLDS_topGenes_moduleScore_bulkSamples_AMKLfusion'),plotFun_MLDS.degs_moduleScore_inBulkSamples_AMKLfusion,rawData=dd,width = 3.8,height = 2.9,res = 500)
  
}


##------- Module scoring in bulk RNA-seq dataset ----------
source('~/ML-DS/R/xx01_moduleScore_bulkRNA_helperFunctions.R')


##------- T21 - Leukaemia module
## Regenerate T21-leukaemia and GATA1s-leukaemia modules ##
mlds_imprints_in_T21 = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/3_MLDSimprint_in_fLiverT21/jul24/MLDS_imprint_in_T21.fLiver.MEMP_2411.csv',row.names = 1)
# Save it for Supplementary Table
write.csv(mlds_imprints_in_T21,'~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/T21_leuk_module_37genes.csv')

moduleList = split(mlds_imprints_in_T21$ensID,mlds_imprints_in_T21$direction)
module_type = 'MLDS_imprints_T21FL'
moduleList[['all']] = mlds_imprints_in_T21$ensID[mlds_imprints_in_T21$direction %in% c("T21_MLDS_down", "T21_MLDS_up")]
names(moduleList)
moduleList = moduleList[unlist(lapply(moduleList,function(x){length(x)>=5}))]

allScore = compute_bulk_module_scores(moduleList=moduleList,
                                      module_type=module_type,gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T,
                                      tpm_fp = NULL,
                                      oxford_meta_fp = NULL,
                                      oxford_tpm_fp = NULL,
                                      plot_filter = TRUE,
                                      min_genes=5)

figS6_gata1s_moduleScore_inBulkSamples(plotDir=plotDir,prefix='Sup.FigS6x_MLDSimprints.T21FL_topGenes_moduleScore_bulkSamples_newBulkMdat_1124',title=NULL,allScore,module_type = module_type)

#**NOTE:
# I found a record of this file (on computer) but I think this was an older version due to lower number of DEGs
#t21_leuk_module = read.csv('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/TableS4_MLDSimprints_in_T21.fLiver.MEMP.csv')
#table(t21_leuk_module$is_DEG..MLDS.vs.diploid.FL.MEP.,t21_leuk_module$is_DEG..T21.FL.MEP.vs.diploid.FL.MEP.)
  




##------- GATA1s - Leukaemia module

#**NOTE:
# I found a record of this file (on computer) but I think this was an older version due to lower number of DEGs
# gata1s_module = read.csv('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/TableS5_GATA1s_gene_module.csv',skip = 3)
# table(gata1s_module$is_DEG..TAM.vs.T21.MEP.,gata1s_module$direction..MLDS.vs.T21.MEP.)


gata1s_module = read.csv('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
table(gata1s_module$tam_vs_mempT21,gata1s_module$direction)
# Save it for Supplementary Table
write.csv(gata1s_module,'~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/GATA1s_leuk_module_757genes_from908genes.csv')


table(gata1s_module$tam_vs_mempT21_group)
gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)

## Select top genes only
# 30% reproduce the plot I used in the paper, but 20% is what I said in the method...
# 250802 --> for revision at NatComms, decided to stick to the method of using 20%
#gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 20/100,]
gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 20/100,]
table(gata1s_module$group)


moduleList = split(gata1s_module$ensID,gata1s_module$group)
module_type = 'GATA1s_topGenes'
moduleList[['all']] = gata1s_module$ensID[gata1s_module$group %in% c("TAM.MLDS.down", "TAM.MLDS.up")]
names(moduleList)
length(moduleList[[4]])
moduleList = moduleList[unlist(lapply(moduleList,function(x){length(x)>=5}))]

allScore = compute_bulk_module_scores(moduleList=moduleList,
                                      module_type=module_type,gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T,
                                      tpm_fp = NULL,
                                      oxford_meta_fp = NULL,
                                      oxford_tpm_fp = NULL,
                                      plot_filter = TRUE,
                                      min_genes=5)

figS6_gata1s_moduleScore_inBulkSamples(plotDir=plotDir,prefix='Sup.FigS6x_GATA1s_topGenes_moduleScore_bulkSamples_newBulkMdat_1124_cf20perc',title=NULL,allScore,module_type = module_type)



##------- MLDS - Leukaemia module
mlds_module = read.csv('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/TableS6_MLDS_gene_module.csv',skip = 2)
table(mlds_module$direction)
mlds_vs_goodTAM = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24/DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')
mlds_vs_goodTAM = mlds_vs_goodTAM[abs(mlds_vs_goodTAM$cellFrac.diff) >= 0.2,]
moduleList = split(mlds_vs_goodTAM$ensID,mlds_vs_goodTAM$direction)
moduleList[['all']] = mlds_vs_goodTAM$ensID
module_type = 'MLDS_topGenes'

allScore = compute_bulk_module_scores(moduleList=moduleList,
                                      module_type=module_type,gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T,
                                      tpm_fp = NULL,
                                      oxford_meta_fp = NULL,
                                      oxford_tpm_fp = NULL,
                                      plot_filter = TRUE,
                                      min_genes=2) # to keep the MLDS_down set which has 2 genes



figS6_gata1s_moduleScore_inBulkSamples(plotDir=plotDir,prefix='Sup.FigSxx_MLDS_topGenes_moduleScore_bulkSamples',title=NULL,allScore,module_type = module_type)






##------- L076 Copy Number profile ----------
chr_to_plot = 'chr9'
cnv_data_fp = list.files('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots',pattern = 'L076_.*_purpleCNprofile_rawData.tsv',full.names = T)

purple_CNV = do.call(rbind,lapply(seq_along(cnv_data_fp),function(i){
  d = read.delim(cnv_data_fp[i],sep = '\t')
  d$sample = gsub('.*L076_|_purple.*$','',basename(cnv_data_fp[i]))
  return(d)
}))


purple_CNV = purple_CNV[purple_CNV$chromosome == chr_to_plot,]
purple_CNV = purple_CNV[purple_CNV$end <= 47e6,]
purple_CNV$x_start = purple_CNV$start
purple_CNV$x_end = purple_CNV$end
purple_CNV$sample = factor(purple_CNV$sample,c('PD62331c','PD62331a','PD64665a','PD66167a'))

chromLabel = purple_CNV
chromLabel = chromLabel %>% group_by(chromosome) %>% summarise(chromLen = max(end),
                                                               chromStart_xpos = min(x_start)) %>% 
  mutate(pos = chromStart_xpos + 0.5*chromLen,
         chromosome = gsub('chr','',as.character(chromosome)))
chromLabel$chromosome[chromLabel$chromosome %in% c('20','22')] = ''
chromLabel$chromosome[chromLabel$chromosome == '9'] = '9p'

plotFun_PURPLE_CNprofile_chr9 = function(noFrame=FALSE,noPlot=FALSE){
  p=ggplot(purple_CNV)+
    geom_rect(data=purple_CNV[purple_CNV$chr_odd,],aes(xmin=x_start,xmax=x_end,ymin=0,ymax=5),fill=grey(0.9))+
    
    # Mark CDKN2A gene
    geom_vline(xintercept = 22005554,linewidth=0.6,alpha=0.7,col='blue')+
    #geom_vline(xintercept = 21974795,linewidth=0.5,alpha=0.7,col='blue')+
    #geom_vline(xintercept = 21995301,linewidth=0.5,alpha=0.7,col='blue')+
    
    #geom_segment(data=tmp,aes(x=chromStart_xpos,xend=chromStart_xpos,y=0,yend=5),color=grey(0),lty=2,size=0.2)+
    geom_segment(aes(x=x_start,xend=x_end,y=round(copyNumber),yend=round(copyNumber)),size=1)+
    geom_segment(aes(x=x_start,xend=x_end,y=round(minorAlleleCopyNumber),yend=round(minorAlleleCopyNumber)),colour='red',size=1)+
    
    
    facet_grid(sample~.)+
    theme_classic()+
    scale_y_continuous(breaks = c(0,1,2,3,4,5),labels = c(0,1,2,3,4,5),limits = c(0,6))+
    scale_x_continuous(breaks = seq(min(purple_CNV$start),max(purple_CNV$end),1e7),
                       labels = round(seq(min(purple_CNV$start),max(purple_CNV$end),1e7)/1e6,2))+
    geom_text(data = chromLabel,aes(x=pos,y=5.4,label = chromosome),size = 3)+
    theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
          strip.background = element_blank(),
          #axis.text.x = element_blank(),
          axis.text = element_text(colour = 'black'),
          #axis.ticks.x = element_blank(),
          axis.ticks.y=element_line(colour = 'black')) + 
    xlab('Genomic position (Mb)') + ylab('')
  print(p)
}
saveFig(file.path(plotDir,paste0('Sup.Figxx_L076_purpleCNprofile_chr9p_CDKN2A')),plotFun_PURPLE_CNprofile_chr9,rawData=purple_CNV,width = 3.5,height = 5,res = 500)


## Quantify the stepwise transcriptional changes ##
diploid_to_mlds = read.csv('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/TableS4_MLDSimprints_in_T21.fLiver.MEMP.csv')
nGenes_diploid_to_mlds = sum(diploid_to_mlds$is_DEG..MLDS.vs.diploid.FL.MEP.)
nGenes_T21_leuk = 37  

gata1s_module = read.csv('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
table(gata1s_module$tam_vs_mempT21,gata1s_module$direction)
nGenes_T21_to_mlds = nrow(gata1s_module)
nGenes_gata1_leuk = sum(gata1s_module$tam_vs_mempT21_group != 'notDE')

nGenes_mlds_leuk = 198

signal_step1= 100*nGenes_T21_leuk/nGenes_diploid_to_mlds
signal_step2= (100*nGenes_gata1_leuk/nGenes_T21_to_mlds)*(100-signal_step1)/100+signal_step1


df = data.frame(step=c('step1','step1','step2','step2','step3','step3'),
                type=c('background','signal','background','signal','background','signal'),
                x_start = c(0,0,
                            signal_step1,signal_step1,
                            signal_step2,signal_step2),
                x_end = c(100,signal_step1,
                          100,signal_step2,
                          100,100),
                y_start = rev(c(1,1,2,2,3,3)),
                y_end = rev(c(1,1,2,2,3,3)))

plotFun_stepwise_transcriptional_changes = function(noFrame=FALSE,noPlot=FALSE){
  p = ggplot(df)+
    geom_segment(aes(x=x_start,xend=x_end,y=y_start,yend=y_end,color=type),linewidth=4)+
    geom_vline(xintercept = c(signal_step1,signal_step2),linetype=2,col=grey(0.6))+
    scale_color_manual(values = c('background' = 'black','signal'='red'))+
    scale_y_continuous(breaks = c(1,2,3),labels = rev(c('Step 1 - T21','Step 2 - GATA1s','Step 3 - MLDS')))+
    theme_classic()+ylab('')+xlab('')+
    theme(panel.border = element_blank(),
          axis.line = element_blank(),
          strip.background=element_rect(linewidth=0),
          strip.text = element_text(size = 9,colour = 'black'),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'none',
          axis.text = element_text(colour = 'black')
    )
  
  print(p)
}

saveFig(file.path(plotDir,'abstract_stepwise_transcriptional_changes'),plotFun_stepwise_transcriptional_changes,rawData=df,width = 5,height = 2,res = 500)



## Supplementary table s5: mark flow-cytometry markers
gene_list = read.csv('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/TableS6_MLDS_gene_module.csv',skip = 2)
csm = read.csv('~/lustre_mt22/generalResources/proteomic_databases/curated_transmembrane_cellsurface_proteins.csv')
csm_classification = csm[grepl('cellphoneDB',csm$source),]

# Select the relevant columns
csm_classification_long <- csm_classification[!is.na(csm_classification$gene_name),] %>%
  dplyr::select(gene_name,
         cpdb_transmembrane, cpdb_peripheral, cpdb_secreted,
         cpdb_receptor, cpdb_integrin) %>%
  mutate(across(-gene_name, as.character)) %>%  # convert all non-gene_name columns to character
  pivot_longer(
    cols = -gene_name,
    names_to = "type",
    values_to = "value"
  )

other_membrane_proteins = csm[!is.na(csm$hpa_Subcellular.main.location) & 
                                csm$hpa_Subcellular.main.location == 'Plasma membrane' & 
                                !is.na(csm$gene_name) & 
                                !csm$gene_name %in% csm_classification_long$gene_name,] %>% 
  mutate(type = 'plasma_membrane',value=T) %>% 
  dplyr::select(gene_name,type,value)
other_membrane_proteins = other_membrane_proteins[!duplicated(other_membrane_proteins),]

csm_classification_long = rbind(csm_classification_long,other_membrane_proteins)

# Check result
head(csm_classification_long)
table(csm_classification_long$value)

csm_classification_long = csm_classification_long[csm_classification_long$value == T & 
                                                    !is.na(csm_classification_long$gene_name),]
csm_classification_long = csm_classification_long[!duplicated(csm_classification_long),]
dim(csm_classification_long)
n_distinct(csm_classification_long$gene_name)
csm_classification_long = csm_classification_long %>% group_by(gene_name) %>% 
  summarise(type = paste0(type,collapse = '::'))
table(csm_classification_long$type)

gene_list = gene_list %>% 
  mutate(marker_categories = dplyr::case_when(
  is_DEG == F ~ '-',
  is_DEG == T & gene_symbol %in% c('COL23A1','P2RX1','TSPAN4') ~ 'transmembrane',
  is_DEG == T & gene_symbol %in% csm_classification_long$gene_name[csm_classification_long$type == 'cpdb_transmembrane::cpdb_secreted::cpdb_receptor'] ~ 'transmembrane_receptor_secreted',
  is_DEG == T & gene_symbol %in% csm_classification_long$gene_name[csm_classification_long$type == 'cpdb_transmembrane::cpdb_secreted'] ~ 'transmembrane_secreted',
  is_DEG == T & gene_symbol %in% csm_classification_long$gene_name[grepl('secreted',csm_classification_long$type)] ~ 'secreted',
  is_DEG == T & gene_symbol %in% csm_classification_long$gene_name[grepl('transmembrane',csm_classification_long$type)] ~ 'transmembrane',
  is_DEG == T & gene_symbol %in% csm_classification_long$gene_name[grepl('receptor',csm_classification_long$type)] ~ 'receptor',
  is_DEG == T & gene_symbol %in% csm_classification_long$gene_name[grepl('plasma',csm_classification_long$type)] ~ 'plasma_membrane',
  .default='others'
))
table(gene_list$marker_categories,gene_list$is_DEG)
#table(gene_list$marker_categories,gene_list$is_cell_surface_markers_or_transmembrane)
table(csm_classification_long$gene_name %in% csm$gene_name)


unsuitable but a membrane protein
View(gene_list$gene_symbol[gene_list$flow_cytometry_suitability=='unsuitable_for_flow_cytometry' & 
                             !gene_list$marker_categories %in% c('-','others')])



# gene_list$is_cell_surface_markers_or_transmembrane[gene_list$is_DEG == T] = (gene_list$gene_symbol[gene_list$is_DEG == T] %in% csm$gene_name)
# gene_list$is_transmembrane = '-'
# gene_list$is_transmembrane[gene_list$is_DEG == T] = (gene_list$gene_symbol[gene_list$is_DEG == T] %in% csm$gene_name[!is.na(csm$cpdb_transmembrane) & csm$cpdb_transmembrane==T])
# a = gene_list[gene_list$is_DEG == T,c('gene_symbol','direction','is_cell_surface_markers_or_transmembrane','is_transmembrane')]
# table(a$cell_surface_markers,a$transmembrane)
# View(a[a$direction == 'upregulated_in_MLDS' & a$cell_surface_markers == T,])

good_flow_marker = c(
  "CD81", "HLA-A", "HLA-B", "HLA-C", "B2M",
  "ITGAM", "BST2", "CD40", "CD36",
  "LAIR1", "ITGA2B",
  
  # Down-regulated
  'ICAM1', 'TNFRSF4', 'CD58', 
  
  # not ideal 
  "FCGRT","KLRC2", "CLEC2B", 
  "SLITRK6","SLC5A3","PTGER3","ATP2B1",
  
  # manual curation "others category"
  'COL23A1','P2RX1',"CERCAM","TSPAN4"
)

gene_list = gene_list %>% 
  mutate(flow_cytometry_suitability = dplyr::case_when(
    direction == 'not_DE' ~ '-',
    gene_list$gene_symbol %in% c(good_flow_marker) ~ marker_categories,
    .default = 'unsuitable_for_flow_cytometry')
    )

write.csv(gene_list,'~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/TableS6_MLDS_gene_module_w.flowMarkers_2509.csv')


table(gene_list$direction,gene_list$flow_cytometry_suitability)
table(gene_list$flow_cytometry_suitability,gene_list$marker_categories,gene_list$direction)
table(gene_list$flow_cytometry_suitability)
View(gene_list[!gene_list$flow_cytometry_suitability %in% c('-','unsuitable_for_flow_cytometry'),])
View(gene_list[gene_list$flow_cytometry_suitability %in% c('unsuitable_for_flow_cytometry') & 
                  !gene_list$marker_categories %in% c('-','others'),])

gene_list$gene_symbol[!gene_list$marker_categories %in% c('-','others')]



table(gene_list$marker_categories[gene_list$gene_symbol %in% good_flow_marker])
# Bad options for flow
"GABRG2","CLU",'ANGPT1'


# of these DE genes, seperate them into 
# Category A: Surface marker (flow cytometry-ready)
# Category B: Intracellular protein (can be measured by intracellular flow after fixation/permeabilization)
# Category C: Not flow-detectable / unsuitable (e.g., transcription factor, nuclear protein, no antibody for flow)

## 141 upregulated_in_MLDS

length(a$gene_symbol[a$direction == 'upregulated_in_MLDS' & a$is_cell_surface_markers_or_transmembrane == T])



"LTK",
groupB_intracellular_flow_with_permeabilisation = c(
  "SLC5A3", 
  "CPA3",  "CPT1A",
  "MINPP1", "CPVL", "CLU", "PTGER3", "DPP7",
  "CPQ", "ITM2C", "SNORC", "TRDN", "CACNA2D1",
  "MYCT1"
)

"TRPC1", ion channel
KCNK1,
groupC_not_suitable = c(
  "GABRG2", "GREB1", "PLSCR4", 
  "IFI27", "TBC1D32", "SYNGR1", "MAGI1",
  "RAB6B", "NFIA", "NFIX",
  # vesicular / intracellular membrane protein
  "TMX4", "SPINT2", "TAP1",
)

table(a$direction,a$cell_surface_markers)
a$gene_symbol[a$direction == 'upregulated_in_MLDS' & a$cell_surface_markers == F]


## Down-regulated genes
a$gene_symbol[a$direction == 'downregulated_in_MLDS' & a$cell_surface_markers == T]

groupB_intracellular_flow_with_permeabilisation_down = c('ABCG1', 'ATP2B1', 'SLC7A5', 'VIM', 'TMBIM1', 'HACD1', 'SCD', 'ANTKMT', 'ZDHHC18')

a$gene_symbol[a$direction == 'downregulated_in_MLDS' & a$cell_surface_markers == F]
groupA = c()

