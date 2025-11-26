## Generates the Figures for paper ##

#################
# Materials map #
#################

## Main Figures ##






## Supplementary Figures ##



##----    Set working directory  -----##
setwd('~/lustre_mt22/Aneuploidy/')


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
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')
source('~/lustre_mt22/generalScripts/utils/sc_utils.R')

plotDir='~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/Plots'
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}


##----------------------------------##
##    Set general color schemes   ####
##----------------------------------##
keep_CC3 = F
geno_cols = c('Diploid' = grey(0.7),
              'T21' = '#93221E',
              'T18' = '#3d5dad',
              'T22' = '#679551',
              'T13' = '#526691',
              'MX' = '#b18db8',
              'Triploid' = '#e07d26')




##--------------------------------------##
##    Figure 1 - Dataset overview     ####
##--------------------------------------##
#akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver_2n/jan24/liver_liverREFmerged_clean_processed_annotated_noUnknowns_noPublished2n_0124.RDS'
#akLiv_mdat_fp = gsub('_0124.RDS','_mdat_0124.csv',akLiv_srat_fp)
#akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_annotated_0424.RDS'
#akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver_2n/jan24/liver_liverREFmerged_clean_processed_annotated_noUnknowns_noPublished2n_0124.RDS'
#mdat = read.csv('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_tmp_0324_mdat.csv')
akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'
akLiv_mdat_fp = gsub('_0824.RDS','_0824_mdat.csv',akLiv_srat_fp)

fig1a_fLiverCohort = function(){
  if(file.exists(file.path(plotDir,'Figxx_fLiver_dataset_rawData.tsv'))){
    mdat = read.delim(file.path(plotDir,'Figxx_fLiver_dataset_rawData.tsv'),sep = '\t',header = T)
  }else{
    ## fLiver dataset ##
    mdat = read.csv(akLiv_mdat_fp)
    
    mdat$Genotype[mdat$Genotype == 'diploid'] = 'Diploid'
    mdat$Genotype[mdat$Genotype == 'complete_trisomy'] = 'Triploid'
    mdat$Genotype[mdat$donorID %in% c('Hsb36','Hsb37')] = 'T21'
    
  }
  
  mdat$Genotype = factor(mdat$Genotype,c('Diploid','T21','T18','T22','MX','Triploid'))
  
  ## Write table of number of cells per donorID
  fLiver_dataset = as.data.frame(table(mdat$donorID,mdat$gestationalAge,mdat$Genotype,mdat$Sex))
  colnames(fLiver_dataset) = c('donorID','gestationalAge','Genotype','Sex','nCell')
  write.csv(fLiver_dataset[fLiver_dataset$nCell >0,],file.path(plotDir,'..','TableS1_fLiver_dataset.csv'),row.names = F)
  
  ## Write table of number of cells per channel ID
  projectMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'Aneuploidy_mani')
  projectMani = projectMani[!is.na(projectMani$sangerSampleID),]
  mdat$sangerSampleID = projectMani$sangerSampleID[match(mdat$orig.ident,projectMani$chanelID)]
  mdat$sangerSampleID[is.na(mdat$sangerSampleID) & grepl('^MY',mdat$orig.ident)] = mdat$orig.ident[is.na(mdat$sangerSampleID) & grepl('^MY',mdat$orig.ident)]
  fLiver_dataset = as.data.frame(table(mdat$sangerSampleID,mdat$sorting_strategy,mdat$assay,mdat$donorID,mdat$gestationalAge,mdat$Genotype,mdat$Sex))
  colnames(fLiver_dataset) = c('sangerSampleID','sortingStrategy','assay','donorID','gestationalAge','Genotype','Sex','nCell')
  
  dd = fLiver_dataset[fLiver_dataset$nCell >0,]
  dd$tissue = 'Foetal liver'
  dd = dd[order(dd$donorID),c('tissue','donorID','Genotype','gestationalAge','Sex','sangerSampleID','sortingStrategy','assay','nCell')]
  
  # change assay
  dd$assay = as.character(dd$assay)
  dd$assay[dd$assay == 'GEX5p'] = "scRNA 5' v2"
  dd$assay = factor(dd$assay)
  # change genotype
  dd$Genotype = as.character(dd$Genotype)
  dd$Genotype[dd$Genotype != 'Triploid'] = gsub('^T','Trisomy ',dd$Genotype[dd$Genotype != 'Triploid'])
  dd$Genotype[dd$Genotype == 'MX'] = 'Mono X'
  dd$Genotype = factor(dd$Genotype,c('Diploid','Trisomy 21','Trisomy 18','Trisomy 22','Mono X','Triploid'))
  
  # order donorID
  dd = dd[order(dd$gestationalAge),]
  dd = dd[order(dd$Genotype),]
  
  # Add PDID
  projectMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'Aneuploidy_mani')
  projectMani = projectMani[projectMani$assay == 'WGS',]
  # PDID generated but did not sequence
  projectMani = projectMani[!projectMani$donorID %in% c('Hsb36','Hsb37','Hsb38') & projectMani$Status == 'completed' & projectMani$Tissue != 'Placenta',]
  projectMani$donorID[projectMani$donorID == 'Hsb22/16223'] = 'Hsb22'
  projectMani$donorID[projectMani$donorID == 'Hsb21/16165'] = 'Hsb21'
  
  dd$PDID = projectMani$PDID[match(dd$donorID,projectMani$donorID)]
  dd$PDID[is.na(dd$PDID)] = '-'
  dd$sangerSampleID = gsub('\\.','_',as.character(dd$sangerSampleID))
  dd$cellRanger_version = ifelse(dd$assay == "scRNA 3' v3.1",'cellranger_302','cellranger_700')
  dd$refGenome_version = ifelse(dd$assay == "scRNA 3' v3.1",'GRCh38 1.2.0','GRCh38 2020-A')
  
  write.csv(dd,file.path(plotDir,'..','TableS1_fLiver_dataset_nCellBy10XSample.csv'),row.names = F)
  
  
  ## Make heatmap of samples
  fLiver_dataset = as.data.frame(table(mdat$donorID,mdat$gestationalAge))
  colnames(fLiver_dataset) = c('donorID','gestationalAge','nCell')
  
  fLiver_dataset$sampleAvailable = ifelse(fLiver_dataset$nCell > 0,1,0)
  fLiver_dataset = pivot_wider(fLiver_dataset,id_cols = 'gestationalAge',names_from = 'donorID',values_from = 'sampleAvailable') 
  fLiver_dataset = column_to_rownames(fLiver_dataset,'gestationalAge')
  fLiver_dataset = (as.matrix(fLiver_dataset))
  fLiver_dataset = fLiver_dataset[,c('Hsb32','Hsb31','Hsb35','Hsb40',
                                     'Hsb33','Hsb37','Hsb38','15724','Hsb36','Hsb34','15877','Hsb39',
                                     '16049','15756','15806','Hsb22',
                                     '15733',
                                     '15905','Hsb21','15680')]
  library(ComplexHeatmap)
  # 
  # geno_cols = c('Diploid' = grey(0.7),
  #               'T21' = '#b18db8',#7a4b82',#dea6af',
  #               'T18' = '#3d5dad',#'#89cff0',#pal37H[17],#'#3B87C7',
  #               'T22' = '#679551',
  #               'T13' = '#526691',
  #               'MX' = '#e07d26',
  #               'Triploid' = '#93221E')
  
  
  
  plotFun_fLiver_dataset = function(noFrame=FALSE,noPlot=FALSE){
    botAnno = HeatmapAnnotation(df=data.frame(Sex = mdat$Sex[match(colnames(fLiver_dataset),mdat$donorID)],
                                              Genotype = mdat$Genotype[match(colnames(fLiver_dataset),mdat$donorID)]),
                                annotation_name_side = 'right',
                                col = list(Genotype = geno_cols,
                                           #Sex = c('female'=colAlpha('#D13063',0.5),'male'='#3477B6')),
                                           Sex = c('female'=grey(0.9),'male'=grey(0.4))))
    
    par(mar=c(0.1,0.1,1,0.1))
    hm = Heatmap(as.matrix(fLiver_dataset),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
                 cluster_rows = F,cluster_columns = F,row_names_side = "left",show_column_dend = F,
                 column_split = mdat$Genotype[match(colnames(fLiver_dataset),mdat$donorID)],
                 border = T,rect_gp = gpar(col = "black", lwd = 1.2),
                 bottom_annotation = botAnno)
    draw(hm)
  }
  saveFig(file.path(plotDir,'Figxx_fLiver_dataset'),plotFun_fLiver_dataset,rawData=mdat,width = 7.8,height = 3.6,res = 500,useDingbats = F)
  
  
  
  
  plotFun_fLiver_dataset_vertical = function(noFrame=FALSE,noPlot=FALSE){
    rowAnno = rowAnnotation(df=data.frame(Genotype = mdat$Genotype[match(colnames(fLiver_dataset),mdat$donorID)]
                                          #Sex = mdat$Sex[match(colnames(fLiver_dataset),mdat$donorID)]
    ),
    col = list(Genotype = geno_cols,
               #Sex = c('female'=colAlpha('#D13063',0.5),'male'='#3477B6')),
               Sex = c('female'=grey(0.9),'male'=grey(0.4))))
    
    par(mar=c(0.1,0.1,1,0.1))
    hm = Heatmap(t(as.matrix(fLiver_dataset)),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
                 cluster_rows = F,cluster_columns = F,row_names_side = "right",show_column_dend = F,show_row_dend = F,
                 split = mdat$Genotype[match(colnames(fLiver_dataset),mdat$donorID)],
                 column_names_side = 'top',row_names_gp = gpar(fontsize=8.5),row_title_gp = gpar(fontsize=12),
                 column_title_gp = gpar(fontsize=15),
                 border = T,rect_gp = gpar(col = "black", lwd = 1.2),
                 left_annotation = rowAnno)
    draw(hm)
  }
  
  
  saveFig(file.path(plotDir,'Fig1B_fLiver_dataset_vertical'),plotFun_fLiver_dataset_vertical,rawData=fLiver_dataset,width = 3.7,height = 8.2,res = 500,useDingbats = F)
}



mlds_srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS'
mlds_mdat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv'

fig1b_MLDScohort = function(){
  ## MLDS dataset ##
  #mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_clean_annotated_noUnknowns_jan24_mdat.csv')
  #mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
  mdat = read.csv(mlds_mdat_fp,row.names = 1) 
  mdat$annot = mdat$annot_aug24
  mdat$finalAnn_broad = mdat$annot_aug24
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
  mlds_dataset2$Sex[mlds_dataset2$Sex == 'F'] = 'female'
  mlds_dataset2$Sex[mlds_dataset2$Sex == 'M'] = 'male'
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
  projectMani = read_excel('~/lustre_mt22/projectManifest.xlsx','MLDS_GOSH')
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
  
  
  write.csv(mlds_dataset2,file.path(plotDir,'..','TableS1_MLDS_dataset.csv'),row.names = F)
  
  
  mlds_dataset = as.data.frame(table(mdat$donorID,mdat$timePoint))
  colnames(mlds_dataset) = c('donorID','timePoint','nCell')
  mlds_dataset = mlds_dataset[order(mlds_dataset$donorID),]
  
  mlds_dataset$sampleAvailable = ifelse(mlds_dataset$nCell > 0,1,0)
  mlds_dataset = pivot_wider(mlds_dataset,id_cols = 'timePoint',names_from = 'donorID',values_from = 'sampleAvailable') 
  mlds_dataset = column_to_rownames(mlds_dataset,'timePoint')
  mlds_dataset = (as.matrix(mlds_dataset))
  
  library(ComplexHeatmap)
  outcomeCols = c('Remission'='#2D4372','Refractory'=col25[2],'Relapse'=col25[5],'unknown'=grey(0.6))
  mtx = mlds_dataset[,colnames(mlds_dataset) %in% mdat$donorID[mdat$disease == 'MLDS']]
  mtx = mtx[rownames(mtx) != 'Recurrent',]
  if(keep_CC3){
    mtx = mtx[,c('CC4','CC5','L091','L042','L040','L039','L019','L041','L178','CC3','L076','L038')]  
  }else{
    mtx = mtx[,c('CC4','CC5','L091','L042','L040','L039','L019','L041','L178','L076','L038')]
  }
  
  
  plotFun_MLDS_dataset = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    botAnno = HeatmapAnnotation(clinicalOutcome = mdat$clinicalOutcome[match(colnames(mtx),mdat$donorID)],
                                annotation_name_side = 'right',
                                col = list(clinicalOutcome = outcomeCols))
    
    
    hm = Heatmap(as.matrix(mtx),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
                 cluster_rows = F,cluster_columns = F,row_names_side = "left",show_column_dend = F,
                 column_split = mdat$disease[match(colnames(mtx),mdat$donorID)],
                 row_names_gp = gpar(fontsize = 10,family='arial'),
                 border = T,rect_gp = gpar(col = "black", lwd = 1.2),
                 bottom_annotation = botAnno)
    draw(hm)
  }
  
  saveFig(file.path(plotDir,'Figxx_MLDS_dataset'),plotFun_MLDS_dataset,rawData=mtx,width = 5.6,height = 2.8,res = 500,useDingbats = F)
  
  plotFun_MLDS_dataset_vertical = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    rowAnno = rowAnnotation(clinicalOutcome = mdat$clinicalOutcome[match(colnames(mtx),mdat$donorID)],
                            col = list(clinicalOutcome = outcomeCols))
    
    hm = Heatmap(t(as.matrix(mtx)),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
                 cluster_rows = F,cluster_columns = F,row_names_side = "right",show_column_dend = F,
                 row_names_gp = gpar(fontsize = 10),column_names_side = 'top',
                 border = T,rect_gp = gpar(col = "black", lwd = 1.2),
                 right_annotation = rowAnno)
    draw(hm)
  }
  
  saveFig(file.path(plotDir,'Figxx_MLDS_dataset_vertical'),plotFun_MLDS_dataset_vertical,rawData=mtx,width = 3.25,height = 5,res = 500,useDingbats = F)
  
  
  mtx = mlds_dataset[,colnames(mlds_dataset) %in% mdat$donorID[mdat$disease == 'TAM']]
  mtx = mtx[rownames(mtx) %in% c('D (PBMC)','Recurrent'),]
  mtx = mtx[,c('CC1','CC2','L075','L114','L182','CC6','CC7','CC8','L156')]
  plotFun_TAM_dataset = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    botAnno = HeatmapAnnotation(clinicalOutcome = mdat$clinicalOutcome[match(colnames(mtx),mdat$donorID)],
                                annotation_name_side = 'right',
                                col = list(clinicalOutcome = outcomeCols))
    
    hm = Heatmap(as.matrix(mtx),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
                 cluster_rows = F,cluster_columns = F,row_names_side = "left",show_column_dend = F,
                 column_split = mdat$disease[match(colnames(mtx),mdat$donorID)],
                 row_names_gp = gpar(fontsize = 10,family='Arial'),
                 column_names_gp = gpar(family='Arial'),
                 border = T,rect_gp = gpar(col = "black", lwd = 1.2),
                 bottom_annotation = botAnno)
    draw(hm)
  }
  saveFig(file.path(plotDir,'Figxx_TAM_dataset'),plotFun_TAM_dataset,rawData=mdat,width = 4.8,height = 1.6,res = 500,useDingbats = F)
  
  
  
  plotFun_TAM_dataset_vertical = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    rowAnno = rowAnnotation(clinicalOutcome = mdat$clinicalOutcome[match(colnames(mtx),mdat$donorID)],
                            col = list(clinicalOutcome = outcomeCols))
    
    hm = Heatmap(t(as.matrix(mtx)),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
                 cluster_rows = F,cluster_columns = F,row_names_side = "right",show_column_dend = F,
                 #split = mdat$disease[match(colnames(mtx),mdat$donorID)],
                 row_names_gp = gpar(fontsize = 10),column_names_side = 'top',
                 border = T,rect_gp = gpar(col = "black", lwd = 1.2),
                 right_annotation = rowAnno)
    draw(hm)
  }
  
  
  saveFig(file.path(plotDir,'Figxx_TAM_dataset'),plotFun_TAM_dataset_vertical,rawData=mtx,width = 2.6,height = 4.7,res = 500,useDingbats = F)
  
}



otherLeuk_srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2408.RDS'
otherLeuk_mdat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/otherLeuk/otherLeuk_clean_annotated_noUnknown_2408_mdat.csv'

fig1c_otherLeuk_cohort = function(){
  ## Import big.srat metadata
  mdat = read.csv(otherLeuk_mdat_fp,row.names = 1)
  
  ## Incorporate AH_LJ annotation
  #ahlj_annot_myeloid = read.csv('/lustre/scratch126/casm/team274sb/project_folders/GOSH_Leuk/sample_metadata/Patient_profiles/Myeloid_other/AHLJ_annot_oct2024_concatenated.csv')
  ahlj_annot = read.csv('/lustre/scratch126/casm/team274sb/project_folders/GOSH_Leuk/sample_metadata/Patient_profiles/AHLJ_annot_oct2024_overall.csv')
  ahlj_annot$orig.ident = gsub('::.*$','',ahlj_annot$CellID)
  
  mdat$cellID_ahlj = gsub('_','::',mdat$cellID)
  mdat$cellID_ahlj = gsub('^RNA','ALeuk_RNA',mdat$cellID_ahlj)
  mdat$cellID_ahlj = gsub('^Leuk','SB_Leuk',mdat$cellID_ahlj)
  mdat$cellID_ahlj = gsub('^NB','CG_SB_NB',mdat$cellID_ahlj)
  table(mdat$orig.ident,mdat$cellID_ahlj %in% ahlj_annot$CellID)
  table(ahlj_annot$orig.ident %in% gsub('::.*$','',mdat$cellID_ahlj))
  table(gsub('::.*$','',mdat$cellID_ahlj) %in% ahlj_annot$orig.ident)
  
  mdat$annot_ahlj = ahlj_annot$AHLJ_annot_oct2024[match(mdat$cellID_ahlj,ahlj_annot$CellID)]
  
  a = table(mdat$annot,mdat$annot_ahlj)
  a = data.frame(table(mdat$annot,mdat$annot_ahlj))
  a=a[a$Freq > 0,]
  View(a)
  
  ## Finalise annotation
  mdat$annot[mdat$donorID == 'L010' & mdat$annot == 'T_cells'] = 'unsure'
  mdat$broadLineage[mdat$donorID == 'L010' & mdat$annot == 'unsure'] = 'unsure'
  
  
  ##---- Remove irrelevant samples -----##
  ## Remove P9_TP1
  mdat = mdat[mdat$timePoint != 'TP1',]
  
  mdat$disease[mdat$donorID == 'L010'] = 'AMKL'
  mdat$disease[mdat$donorID == 'L069'] = 'AEL'
  table(mdat$donorID,mdat$disease)
  table(mdat$donorID[mdat$disease == 'pBALL'],mdat$timePoint[mdat$disease == 'pBALL'])
  table(mdat$donorID[mdat$disease == 'pAML'],mdat$timePoint[mdat$disease == 'pAML'])
  table(mdat$broadLineage,mdat$timePoint)
  dim(mdat)
  
  ##------ Redefine clinical data (outcome and timepoint and tissue) ----##
  mdat$Genotype[mdat$Genotype == 'diploid'] = 'Diploid'
  mdat$Genotype[mdat$Genotype == 'T21'] = 'Trisomy 21'
  table(mdat$clinicalOutcome,mdat$donorID,mdat$disease)
  mdat$disease = factor(mdat$disease,(c('MDS','AMKL','AEL','pAML','pBALL','LPD')))
  mdat$clinicalOutcome[mdat$donorID == 'L080'] = 'Relapse'
  mdat$clinicalOutcome[mdat$donorID == 'P9'] = 'Remission'
  mdat$clinicalOutcome[mdat$clinicalOutcome == '?'] = 'Unknown'
  mdat$clinicalOutcome[grepl('^Refractory|^refractory',mdat$clinicalOutcome)] = 'Refractory'
  mdat$clinicalOutcome[grepl('^Remission',mdat$clinicalOutcome)] = 'Remission'
  mdat$clinicalOutcome[grepl('^Relapse',mdat$clinicalOutcome)] = 'Relapse'
  table(mdat$clinicalOutcome,mdat$donorID)
  
  table(mdat$timePoint,mdat$donorID)
  mdat$timePoint[mdat$timePoint == 'RelapseD0'] = 'Relapse Diagnostic'
  
  mdat$donorID = factor(mdat$donorID,c('L067','L010','P9','L069',
                                       'L051','L044','L058','L100',
                                       'L001','L080','L007','L002','L014','L016','L050','L068',
                                       'L027','L062'))
  mdat$age_yrs[mdat$age_yrs == '1Y,0m'] = '1Y'
  
  mdat$tissue[mdat$orig.ident %in% c('Leuk13760341','Leuk13760342','Leuk13760340','Leuk13645527')] = 'BM'
  mdat$tissue[mdat$orig.ident %in% c('NB14406184','NB14406185')] = 'Brain'
  
  mdat$blastPerc[mdat$disease == 'LPD'] = '?'
  
  ##---- Remove low quality cells, calculate Leukaemia vs Normal cells -----##
  mdat$group = ifelse(mdat$broadLineage == 'Tumour','Tumour',
                      ifelse(mdat$broadLineage != 'unsure','Normal','unsure'))
  nLeuk = as.data.frame(table(mdat$donorID,mdat$orig.ident,mdat$group))
  colnames(nLeuk) = c('donorID','channelID','group','Freq')
  nLeuk = nLeuk[nLeuk$Freq >0,]
  
  
  
  
  ## Select for relevant columns only
  # Write supplementary table
  pLeuk_dataset = as.data.frame(table(mdat$donorID,mdat$Genotype,mdat$sex,mdat$age_yrs,mdat$disease,mdat$orig.ident,mdat$timePoint,mdat$tissue,mdat$clinicalOutcome,mdat$blastPerc))
  colnames(pLeuk_dataset) = c('donorID','Genotype','Sex','Age (yrs)','Disease','channelID','timePoint','Tissue','Clinical outcome','Blast percentage','nCell')
  pLeuk_dataset = pLeuk_dataset[pLeuk_dataset$nCell > 0,]
  pLeuk_dataset$Sex = as.character(pLeuk_dataset$Sex)
  pLeuk_dataset$Sex[pLeuk_dataset$Sex == 'F'] = 'female'
  pLeuk_dataset$Sex[pLeuk_dataset$Sex == 'M'] = 'male'
  pLeuk_dataset$timePoint = gsub(' (.*)','',pLeuk_dataset$timePoint)
  pLeuk_dataset$Tissue = as.character(pLeuk_dataset$Tissue)
  pLeuk_dataset$Tissue[pLeuk_dataset$Tissue == 'BM'] = 'Bone marrow'
  pLeuk_dataset$Tissue[pLeuk_dataset$Tissue %in% c('Blood','PBMC')] = 'Peripheral blood'
  
  pLeuk_dataset$timePoint[pLeuk_dataset$timePoint == 'Relapse'] = 'Relapse diagnostic'
  #pLeuk_dataset$timePoint[pLeuk_dataset$timePoint == 'Diagnostic'] = 'Day 0'
  #pLeuk_dataset$timePoint[pLeuk_dataset$timePoint == 'TP1'] = 'Timepoint 1 - post Chemo'
  pLeuk_dataset$timePoint = factor(pLeuk_dataset$timePoint,c('Diagnostic','Relapse diagnostic','Timepoint 1 - post Chemo'))
  pLeuk_dataset = pLeuk_dataset[order(pLeuk_dataset$timePoint),]
  pLeuk_dataset = pLeuk_dataset[order(pLeuk_dataset$donorID),]
  pLeuk_dataset = pLeuk_dataset[order(pLeuk_dataset$Disease),]
  
  
  ## Fix disease
  pLeuk_dataset$Disease = as.character(pLeuk_dataset$Disease)
  pLeuk_dataset$Disease = gsub('^p','',pLeuk_dataset$Disease)
  pLeuk_dataset$Disease = factor(pLeuk_dataset$Disease,c('MDS','AMKL','AEL','AML','BALL','LPD'))
  
  ##--- Add leuk / normal cell count ----##
  leuk = nLeuk[nLeuk$group == 'Tumour',]
  norm = nLeuk[nLeuk$group == 'Normal',]
  pLeuk_dataset$nLeuk = leuk$Freq[match(pLeuk_dataset$channelID,leuk$channelID)]
  pLeuk_dataset$nLeuk[is.na(pLeuk_dataset$nLeuk)] = 0
  pLeuk_dataset$nNorm = norm$Freq[match(pLeuk_dataset$channelID,norm$channelID)]
  pLeuk_dataset$nNorm[is.na(pLeuk_dataset$nNorm)] = 0
  which((pLeuk_dataset$nLeuk + pLeuk_dataset$nNorm) != pLeuk_dataset$nCell)
  all((pLeuk_dataset$nNorm + pLeuk_dataset$nLeuk) == pLeuk_dataset$nCell)
  table(pLeuk_dataset$donorID)
  
  sum(pLeuk_dataset$nCell)
  sum(pLeuk_dataset$nLeuk)
  sum(pLeuk_dataset$nNorm)
  
  
  
  
  write.csv(pLeuk_dataset,file.path(plotDir,'..','TableS1_otherLeuk_dataset.csv'),row.names = F)
  
  
  # Table for heatmap plot
  pLeuk_dataset = as.data.frame(table(mdat$donorID,mdat$disease,mdat$timePoint))
  colnames(pLeuk_dataset) = c('donorID','Disease','timePoint','nCell')
  pLeuk_dataset$sampleAvailable = ifelse(pLeuk_dataset$nCell > 0,1,0)
  pLeuk_dataset$id = paste0(pLeuk_dataset$donorID,':',pLeuk_dataset$Disease)
  mdat$id = paste0(mdat$donorID,':',mdat$disease)
  pLeuk_dataset = pivot_wider(pLeuk_dataset,id_cols = 'timePoint',names_from = 'id',values_from = 'sampleAvailable') 
  pLeuk_dataset = column_to_rownames(pLeuk_dataset,'timePoint')
  pLeuk_dataset = (as.matrix(pLeuk_dataset))
  pLeuk_dataset = pLeuk_dataset[,colSums(pLeuk_dataset) > 0]
  
  library(ComplexHeatmap)
  outcomeCols = c('Remission'='#2D4372','Refractory'=col25[2],'Relapse'=col25[5],'Unknown'=grey(0.6))
  
  plotFun_pLeuk_dataset = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    mdat$Genotype[mdat$Genotype == 'diploid'] = 'Diploid'
    mtx = pLeuk_dataset
    colnames(mtx) = gsub(':.*$','',colnames(mtx))
    botAnno = HeatmapAnnotation(df = data.frame(#Genotype = mdat$Genotype[match(colnames(mtx),mdat$donorID)],
      clinicalOutcome = mdat$clinicalOutcome[match(colnames(mtx),mdat$donorID)]),
      annotation_name_side = 'left',
      col = list(clinicalOutcome = outcomeCols))
    rowAnno = rowAnnotation(tissue = gsub('.*\\(|).*$','',rownames(mtx)),
                            col = list(tissue = c('BM'='#FBE8D8','PBMC'='#D3A292','Brain'='#5A7969','CSF'='#9EAFB6')))
    
    hm = Heatmap(as.matrix(mtx),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
                 cluster_rows = F,cluster_columns = F,row_names_side = "left",show_column_dend = F,
                 show_row_names = T,
                 column_split = mdat$disease[match(colnames(pLeuk_dataset),mdat$id)],
                 row_names_gp = gpar(fontsize = 10),
                 row_title_rot = 0,
                 #right_annotation = rowAnno,
                 border = T,rect_gp = gpar(col = "black", lwd = 1.2),column_gap = unit(0.4,'cm'),
                 #row_split = gsub(' \\(.*$','',rownames(mtx)),
                 bottom_annotation = botAnno)
    draw(hm)
  }
  
  saveFig(file.path(plotDir,'FigXX_pLeuk_dataset_v2'),plotFun_pLeuk_dataset,rawData=mdat,width = 8.3,height = 1.57,res = 500,useDingbats = F)
  
  # plotFun_pLeuk_dataset_vertical = function(noFrame=FALSE,noPlot=FALSE){
  #   par(mar=c(0.1,0.1,1,0.1))
  #   mtx = pLeuk_dataset
  #   
  #   rowAnno = rowAnnotation(clinicalOutcome = mdat$clinicalOutcome[match(colnames(mtx),mdat$donorID)],
  #                           col = list(clinicalOutcome = outcomeCols))
  #   
  #   hm = Heatmap(t(as.matrix(mtx)),name = 'sample',col = c("1"=grey(0.2),"0"=grey(1)),
  #                cluster_rows = F,cluster_columns = F,row_names_side = "right",show_column_dend = F,
  #                row_names_gp = gpar(fontsize = 10),column_names_side = 'top',
  #                border = T,rect_gp = gpar(col = "black", lwd = 1.2),
  #                right_annotation = rowAnno)
  #   draw(hm)
  # }
  # 
  # saveFig(file.path(plotDir,'Figxx_pLeuk_dataset_vertical'),plotFun_pLeuk_dataset_vertical,rawData=mdat,width = 3.25,height = 4.3,res = 500,useDingbats = F)
  
}




##---------------------------------------------##
##    Figure 2 - Foetal liver AK dataset     ####
##---------------------------------------------##
#akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver_2n/jan24/liver_liverREFmerged_clean_processed_annotated_noUnknowns_noPublished2n_0124.RDS'
#akLiv_mdat_fp = gsub('_0124.RDS','_mdat_0124.csv',akLiv_srat_fp)
#akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/mar24/liver_liverREFmerged_clean_processed_annotated_0424.RDS'
akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'
akLiv_mdat_fp = gsub('_0824.RDS','_0824_mdat.csv',akLiv_srat_fp)




celltype_cols = c(
  'HSC_MPP'='#EDBD74',"MEMP_MEP"='#d91623',"CMP_GMP"='#4E8920',"LMPP_ELP"='#e6306a',
  "earlyMK"='gold1','MK'='yellow4',
  'Mast.cell'=pal34H[23],
  'EE'='#7ec0ee','ME'='#3f7a91','LE'='#2e4472',
  "promyelocyte"=pal34H[31],"myelocyte"='#de8d02', #orange
  "pDC"='#8a4108',
  'proMono'=pal37H[31],'Monocyte'=pal37H[30],"DC2"='#e07f4f',"DC1"='#273e9c','Kupffer.cell'='#f0f030','Macrophage'=pal34H[19], # blue
  "pro.B.cell"='#d66b9f',"pre.B.cell"='#9e6c88',"B.cell"='#702963', # pink / purple
  "ILC.precursor"=grey(0.8),"T.cell"=grey(0.6),"NK_T"=grey(0.5), # greys
  "Hepatocyte"='#c3afcc',"Fibroblast"=pal37H[26],"Endo"='#e6a5d5',"Mesenchyme"=pal37H[28],"Neuron"='#ac8847',#green
  'Cholangiocytes'=pal37H[21],'Mesothelial_cells'=pal37H[32],'Trophoblast'='#b0ee95'
)



celltype_cols_old = c(#col25[6],pal34H[c(2)],pal34H[c(14)],pal37H[8], #yellow
  'HSC_MPP'='#FFD700',"MEMP_MEP"='#F0E68C',"CMP_GMP"='#FFC72C',"LMPP_ELP"='#FEBE10',
  "earlyMK"='#E0D4EB','MK'='#9C6DA5',
  #col25[c(10)],pal34H[23],#purple
  'Mast.cell'='#FFB7CE',# baby pink
  #pal34H[c(12,17,19)], # red
  'EE'='#E9967A','ME'='#E44D2E','LE'='#9e1b32',
  #'#FFAF91','#FB785E','#BD422E',
  #pal34H[c(26,25,27)], #blue
  "MOP"=pal37H[22],'proMono'=pal37H[37],'Monocyte'=pal37H[19],'Macrophage'=pal37H[17],'Kupffer.cell'=pal34H[25], # blue
  "DC1"=pal34H[33],"DC2"=pal34H[1],"pDC"=pal34H[30],#brown
  "promyelocyte"=pal34H[31],"myelocyte"=col25[5], #orange
  
  #pal37H[c(11,13,12)], # pink/purple
  #pal34H[c(11,24,22,20,34)],
  "pro.B.cell"='#B5338A',"pre.B.cell"='#8D4585',"B.cell"='#702963', # pink / purple
  "ILC.precursor"=grey(0.8),"T.cell"=grey(0.6),"NK_T"=grey(0.5), # greys
  "Hepatocyte"=pal37H[26],"Fibroblast"=pal37H[23],"Endo"=pal37H[30],"Mesenchyme"=pal37H[28],"NPC"=pal37H[33],#green
  'Cholangiocytes'=pal37H[21],'Mesothelial_cells'=pal37H[32],'Trophoblast'='#90EE90'
)


lineage_cols = c('#F0E68C',
                 '#9C6DA5',
                 '#FFB7CE',
                 '#E44D2E',
                 pal37H[17],
                 pal34H[1],
                 col25[5],
                 '#8D4585',
                 grey(0.6),
                 pal37H[32])

fLiver_broadLin_cols_v4 = c('HSC_MPP' = '#FFC72C',"Ery/MegK/Mast lineage" = '#D3A292',"Myeloid lineage" = '#1F4072',
                            'B lineage'='#702963','T/NK lineage'=grey(0.8),'Stromal'='#678551','others'='#90EE90')



fig2a_2nAK_Liv_UMAP = function(){
  
  if(file.exists(file.path(plotDir,'Fig1b_2nAKLiv_geno_UMAP_rawData.tsv'))){
    dd = read.delim(file.path(plotDir,'Fig1b_2nAKLiv_geno_UMAP_rawData.tsv'),header = T,sep = '\t')
    
  }else{
    mdat = read.csv(akLiv_mdat_fp)
    
    dd = mdat[,c("cellID","donorID","broadLineage",'finalAnn_broad','Genotype','UMAP_1','UMAP_2')]
    dd = dd[dd$finalAnn_broad != 'doublets',]
    
    
    set.seed(2397)
    
    # dd = rbind(dd[!dd$finalAnn_broad %in% c('EE','ME','LE','Hepatocyte'),],
    #            dd[dd$finalAnn_broad %in% c('EE','ME','LE','Hepatocyte'),][sample(1:nrow(dd),150000),])
    dd=dd[sample(1:nrow(dd),150000),]
    # ggplot(dd,aes(UMAP_1,UMAP_2,col=Genotype))+
    #   geom_point(size=0.01)+
    #   scale_color_manual(values = col25)+
    #   theme_classic()
  }
  
  
  plotFun_celltype = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0,0,0.8,0))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','Foetal Liver'),
         frame.plot=F)
    
    if(!noPlot){
      celltype_cols_1 = colAlpha(celltype_cols,alphas = 0.7)
      #names(celltype_cols_1) = levels(dd$finalAnn_broad)
      names(celltype_cols_1) = names(celltype_cols)
      
      dd1 = rbind(dd[dd$finalAnn_broad %in% c('MEMP_MEP'),][sample(1:nrow(dd[dd$finalAnn_broad %in% c('MEMP_MEP'),]),1000),],
                  dd[!dd$finalAnn_broad %in% c('MEMP_MEP','earlyMK','MK'),])
      dd1 = dd1[sample(1:nrow(dd1),nrow(dd1)),]
      dd2 = dd[dd$finalAnn_broad %in% c('MEMP_MEP') & !dd$cellID %in% dd1$cellID,]
      dd3 = dd[dd$finalAnn_broad %in% c('earlyMK','MK'),]
      
      points(dd1$UMAP_1,dd1$UMAP_2,
             col = celltype_cols_1[as.character(dd1$finalAnn_broad)],
             pch = 19,
             cex=0.07)
      points(dd2$UMAP_1,dd2$UMAP_2,
             col = celltype_cols_1[as.character(dd2$finalAnn_broad)],
             pch = 19,
             cex=0.07)
      points(dd3$UMAP_1,dd3$UMAP_2,
             col = celltype_cols_1[as.character(dd3$finalAnn_broad)],
             pch = 19,
             cex=0.07)
    }
    
    if(!noFrame){
      #Add coloured labels
      
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ finalAnn_broad,data=dd,FUN=mean)
      # mids[mids$broadLineage == 'HSC & prog.',]$UMAP_1 =  mids[mids$broadLineage == 'HSC & prog.',]$UMAP_1 + 2
      # mids[mids$broadLineage == 'HSC & prog.',]$UMAP_2 =  mids[mids$broadLineage == 'HSC & prog.',]$UMAP_2 + 2
      # mids[mids$broadLineage == 'Stromal',]$UMAP_1 =  mids[mids$broadLineage == 'Stromal',]$UMAP_1 - 2
      # 
      
      # mids$label = as.numeric(factor(mids$broadLineage,levels = c('HSC & prog.','Megakaryocytes','Mast.cell','Erythroblasts',
      #                                                             'Monocyte/Macrophage','Dendritic cells','Myelocytes',
      #                                                             'B lineage','T/NK lineage','Stromal','others')))
      mids$label = mids$finalAnn_broad
      # mids$legend = factor(paste0(mids$label,' - ',mids$broadLineage),levels = c('1 - HSC & prog.',
      #                                                                            '2 - Megakaryocytes',
      #                                                                            '3 - Mast.cell',
      #                                                                            '4 - Erythroblasts',
      #                                                                            '5 - Monocyte/Macrophage',
      #                                                                            '6 - Dendritic cells',
      #                                                                            '7 - Myelocytes',
      #                                                                            '8 - B lineage',
      #                                                                            '9 - T/NK lineage',
      #                                                                            '10 - Stromal'))#,
      #                                                                            # '11 - others'))
      
      #Position tweaks
      # mids[mids$finalAnn=='Leukocyte','UMAP_2'] = mids[mids$finalAnn=='Leukocyte','UMAP_2'] + 4.9
      # mids[mids$finalAnn=='Leukocyte','UMAP_1'] = mids[mids$finalAnn=='Leukocyte','UMAP_1'] - 0.3
      # mids[mids$finalAnn=='PTC','UMAP_2'] = mids[mids$finalAnn=='PTC','UMAP_2'] - 4.9
      # mids[mids$finalAnn=='PTC','UMAP_1'] = mids[mids$finalAnn=='PTC','UMAP_1'] -0.8
      # mids[mids$finalAnn=='Tumour','UMAP_2'] = mids[mids$finalAnn=='Tumour','UMAP_2'] + 3
      # mids[mids$finalAnn=='Tumour','UMAP_1'] = mids[mids$finalAnn=='Tumour','UMAP_1'] + 0.9
      # 
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$label,cex = 0.6,xpad = 1.3,ypad = 2.3,border = T,
                   bg=celltype_cols[mids$label],
                   #col=ccs[mids$broadLineage])
                   col='black')
      # for(i in 1:length(mids$UMAP_1)){
      #   draw.circle(mids$UMAP_1[i],mids$UMAP_2[i],
      #               #labels=mids$broadLineage,
      #               radius=0.4,col=ccs[mids$broadLineage][i])  
      # }
      # 
      
      #Add PDID barplot
      # pdid = table(dd$PDID,dd$finalAnn)
      # pdid = sweep(pdid, 2, colSums(pdid), "/")
      # p.mids = aggregate(cbind(UMAP_1,UMAP_2) ~ finalAnn,data=dd,FUN=mean)
      # for(i in unique(p.mids$finalAnn)){
      #   if(i == 'Leukocytes'){
      #     xshift = -5.6
      #     yshift = -5.3
      #   }else if(i=='PTC'){
      #     xshift = 4.3
      #     yshift = -2.8
      #   }else if(i=='Tumour'){
      #     xshift = 4.9
      #     yshift = 1.0
      #   }
      #   dat = pdid[c('PD35918','PD36793','PD37104','PD37228'),colnames(pdid) == i]*4.5
      #   xleft = p.mids[p.mids$finalAnn == i,]$UMAP_1 + xshift
      #   xright = xleft + 1
      #   ybottom = p.mids[p.mids$finalAnn == i,]$UMAP_2 + yshift + c(0,cumsum(dat)[-length(dat)])
      #   ytop= ybottom + dat
      
      #rect(xleft=xleft,
      #     xright=xright,
      #     ybottom=ybottom,
      #     ytop=ytop,
      #     col = ccs[names(dat)],lwd = 0.5,
      #     border = 'black')
      
      #legend(x=-14, y=-4.2,legend=levels(mids$legend),fill = lineage_cols[mids$broadLineage[match(levels(mids$legend),mids$legend)]],lwd = 0,cex = 0.75,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
      #legend(x=-5, y=3.2,legend=levels(mids$legend),fill = lineage_cols[mids$broadLineage[match(levels(mids$legend),mids$legend)]],lwd = 0,cex = 1,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
      
    }
    
  }
  
  saveFig(file.path(plotDir,'FigS1_2nAKLiv_celltype_UMAP'),plotFun_celltype,rawData=mdat,width = 4.1,height = 4,res = 500,useDingbats = F)
  
  
  
  
  plotFun_geno = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    dd$Genotype[dd$Genotype=='complete_trisomy'] = 'Triploid'
    dd$Genotype[dd$Genotype=='diploid'] ='Diploid'
    dd$Genotype = factor(dd$Genotype,names(geno_cols))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','Foetal Liver'),
         frame.plot=F)
    
    
    if(!noPlot){
      
      geno_cols_1 = colAlpha(geno_cols,alphas = 0.4)
      names(geno_cols_1) = levels(dd$Genotype)
      
      # points(dd$UMAP_1,dd$UMAP_2,
      #        col = geno_cols_1[dd$Genotype],
      #        pch = 19,
      #        cex=0.01)
      
      # Plot diploid first, then T21, then T18, then MX, then T22, and Triploid
      for(geno in c('Diploid','T21','T18','MX','Triploid','T22')){
        points(dd$UMAP_1[dd$Genotype == geno],dd$UMAP_2[dd$Genotype == geno],
               col = geno_cols_1[geno],
               pch = 19,
               cex=0.01)  
      }
      
      
      
      
    }
    #legend(x=-14, y=-4.2,legend=levels(dd$Genotype),fill = geno_cols[levels(dd$Genotype)],lwd = 0,cex = 0.9,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
    #legend(x=-5, y=3.2,legend=levels(mids$legend),fill = lineage_cols[mids$broadLineage[match(levels(mids$legend),mids$legend)]],lwd = 0,cex = 1,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
    
  }
  
  saveFig(file.path(plotDir,'Fig1_2nAKLiv_geno_UMAP'),plotFun_geno,rawData=mdat,width = 2.85,height = 2.65,res = 500,useDingbats = F)
}








dd$finalAnn_broad = factor(dd$finalAnn_broad,levels = c("HSC_MPP","MEMP_MEP","CMP_GMP","LMPP_ELP", # brown
                                                        "earlyMK",'MK', # purple
                                                        'Mast.cell', # slamon
                                                        'EE','ME','LE', # red
                                                        "proMono","Monocyte","Macrophage","Kupffer.cell", # blue
                                                        "DC1","DC2","pDC", # yellow
                                                        "promyelocyte","myelocyte", # orange
                                                        "pro.B.cell","pre.B.cell","B.cell", # pink / purple
                                                        "ILC.precursor","T.cell","NK_T", # greys
                                                        "Hepatocyte","Fibroblast","Endo","Mesenchyme","NPC",'Cholangiocytes','Mesothelial_cells','Trophoblast')) # green


dd$broadLineage = factor(dd$broadLineage, levels = c('HSC & prog.','Megakaryocytes','Mast.cell','Erythroblasts',
                                                     'Monocyte/Macrophage','Dendritic cells','Myelocytes',
                                                     'B lineage','T/NK lineage','Stromal'))




fig2a_2nAK_Liv_GenotypeContribution = function(){ 
  mdat = read.csv(akLiv_mdat_fp)
  
  dd = mdat[,c("cellID","donorID","broadLineage",'finalAnn_broad','Genotype','UMAP_1','UMAP_2')]
  dd$Genotype[dd$Genotype == 'complete_trisomy'] = 'Triploid'
  dd$Genotype[dd$Genotype == 'diploid'] = 'Diploid'
  
  ##--------- genotype contribution  --------##
  plotFun_genoContribution_no2n = function(noFrame=FALSE,noPlot=FALSE){
    
    dd$finalAnn_broad = factor(dd$finalAnn_broad,levels = c("HSC_MPP","MEMP_MEP","CMP_GMP","LMPP_ELP", # brown
                                                            "earlyMK",'MK', # purple
                                                            'Mast.cell', # slamon
                                                            'EE','ME','LE', # red
                                                            "MOP","proMono","Monocyte","Macrophage","Kupffer.cell", # blue
                                                            "DC1","DC2","pDC", # yellow
                                                            "promyelocyte","myelocyte", # orange
                                                            "pro.B.cell","pre.B.cell","B.cell", # pink / purple
                                                            "ILC.precursor","T.cell","NK_T", # greys
                                                            "Hepatocyte","Fibroblast","Endo","Mesenchyme","NPC",'Cholangiocytes','Mesothelial_cells')) # green
    
    
    dd$broadLineage = factor(dd$broadLineage, levels = c('HSC & prog.','Megakaryocytes','Mast.cell','Erythroblasts',
                                                         'Monocyte/Macrophage','Dendritic cells','Myelocytes',
                                                         'B lineage','T/NK lineage','Stromal'))
    
    
    dd$Genotype = factor(dd$Genotype,rev(c('Diploid','T21','T18','T22','MX','Triploid')))
    
    
    p1 = ggplot(dd,aes(x= broadLineage,fill=Genotype))+
      geom_bar(position='fill')+
      scale_fill_manual(values = geno_cols)+
      theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
      xlab('') + ylab('Genotype contribution')
    
    
    p2 = ggplot(dd,aes(y= broadLineage,fill=Genotype))+
      geom_bar(position='fill')+
      scale_fill_manual(values = geno_cols)+
      theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank()) +
      ylab('') + xlab('Genotype contribution')
    
    print(p2)
  }
  
  saveFig(file.path(plotDir,'FigS1b_2nAKLiv_genoContrByLin_barplot_vertical'),plotFun_genoContribution_no2n,rawData=dd,width = 3.4,height = 2.5,res = 500,useDingbats = F)
  
} 







# broadLin_cols_v3 = c('HSC_MPP' = '#FFD700',"Ery/MegK/Mast lineage" = '#c43c1f',"Myeloid lineage" = '#70480b',
#                      'B lineage'=pal37H[c(19)],'T/NK lineage'=grey(c(0.7)),'Stromal'=pal37H[c(32)],'others'='#90EE90')
broadLin_cols_v3 = c('HSC_MPP' = '#EDBD74',"Ery/MegK/Mast lineage" = '#3f7a91',"Myeloid lineage" = '#e07f4f',
                     'B lineage'='#702963','T/NK lineage'=grey(0.6),'Stromal'='#c3afcc','others'='#b0ee95')
fLiver_broadLin_cols_v4 = c('HSC_MPP' = '#FFC72C',"Ery/MegK/Mast lineage" = '#D3A292',"Myeloid lineage" = '#1F4072',
                            'B lineage'='#702963','T/NK lineage'=grey(0.8),'Stromal'='#678551','others'='#90EE90')

fLiver_broadLin_cols_v4 = c('HSC_MPP' = '#EDBD74',"Ery/MegK/Mast lineage" = '#74acc2',"Myeloid lineage" = '#e07f4f',
                            'B lineage'='#702963','T/NK lineage'=grey(0.6),'Stromal'='#c3afcc','others'='#b0ee95')


fig2a_2nAK_Liv_cellTypeContribution = function(){ 
  mdat = read.csv(akLiv_mdat_fp)
  #mdat = mdat[mdat$broadLineage != 'others',] # Remove lowQual + doublets
  mdat$finalAnn_broad = as.character(mdat$finalAnn)
  
  data = mdat[mdat$finalAnn != 'doublets',c("cellID","donorID","broadLineage",'finalAnn_broad','Genotype')]
  data$broadLineage2 = data$broadLineage
  data$broadLineage2[data$broadLineage2 == 'HSC & prog.'] = data$finalAnn_broad[data$broadLineage2 == 'HSC & prog.']
  data$Genotype[data$Genotype == 'complete_trisomy'] = 'Triploid'
  data$Genotype[data$Genotype == 'diploid'] = 'Diploid'
  
  
  ##---- Celltype contribution per donorID per geno
  plotFun_celltypeContribution_barPlot_perDonorID_suppFig2B = function(noFrame=FALSE,noPlot=FALSE){
    
    data$broadLineage3 = as.character(data$broadLineage2)
    data$broadLineage3[data$broadLineage3 %in% c("CMP_GMP",'Monocyte/Macrophage','Dendritic cells','Myelocytes')] = 'Myeloid lineage'
    data$broadLineage3[data$broadLineage3 %in% c("MEMP_MEP",'Megakaryocytes','Mast.cell','Erythroblasts','earlyMK')] = 'Ery/MegK/Mast lineage'
    data$broadLineage3[data$broadLineage3 %in% c("LMPP_ELP",'B lineage')] = 'B lineage'
    
    
    
    dd = data %>% group_by(donorID,broadLineage3,Genotype) %>% summarise(nCell = n_distinct(cellID)) %>% 
      group_by(donorID) %>% mutate(totalCell = sum(nCell),ctFrac = nCell/totalCell) 
    dd$donorID = factor(dd$donorID,rev(c('Hsb32','Hsb31','Hsb35','Hsb40',
                                         'Hsb33','Hsb37','Hsb38','15724','Hsb36','Hsb34','15877','Hsb39',
                                         '16049','15756','15806','Hsb22',
                                         '15733',
                                         '15905','Hsb21','15680')))
    dd$Genotype = factor(dd$Genotype,(c('Diploid','T21','T18','T22','MX','Triploid')))
    dd$broadLineage3 = factor(dd$broadLineage3,rev(c('HSC_MPP','Ery/MegK/Mast lineage','Myeloid lineage','B lineage','T/NK lineage','Stromal','others')))
    
    p1 = ggplot(dd,aes(x = donorID,y=ctFrac,fill=broadLineage3))+
      geom_col(width = 0.75)+
      facet_grid(.~Genotype,scales = 'free_x',space = 'free_x')+
      scale_fill_manual(values = fLiver_broadLin_cols_v4,name='Lineage')+
      theme_classic(base_size = 8.3)+ 
      scale_y_continuous(breaks = c(0,0.5,1))+
      theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
            axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
            legend.key.size = unit(0.45,'cm')) + 
      xlab('') + ylab('Celltype contribution')
    #print(p1)
    
    
    p2 = ggplot(dd,aes(y = donorID,x=ctFrac,fill=broadLineage3))+
      geom_col(width = 0.75)+
      facet_grid(Genotype~.,scales = 'free_y',space = 'free_y',switch = "y")+
      scale_fill_manual(values = fLiver_broadLin_cols_v4,name='Lineage')+
      theme_classic(base_size = 8.3)+ 
      scale_x_continuous(breaks = c(0,0.5,1))+
      scale_y_discrete(position = "right")+
      theme(panel.border = element_rect(fill=F),axis.line = element_blank(),
            legend.key.size = unit(0.45,'cm'),axis.text.y = element_text(size=8),
            strip.background = element_blank()) + 
      ylab('') + xlab('Celltype contribution')
    
    print(p2)  
  }
  
  saveFig(file.path(plotDir,'Supp.Fig2B_2nAKLiv_celltypeFraction_byDonorID_barPlot'),plotFun_celltypeContribution_barPlot_perDonorID_suppFig2B,rawData=dd,width = 3.8,height = 4.9,res = 500,useDingbats = F)
  
  
  
  
  
  plotFun_celltypeContribution_barPlot_v2 = function(noFrame=FALSE,noPlot=FALSE){
    
    data$broadLineage3 = as.character(data$broadLineage2)
    data$broadLineage3[data$broadLineage3 %in% c("CMP_GMP",'Monocyte/Macrophage','Dendritic cells','Myelocytes')] = 'Myeloid lineage'
    data$broadLineage3[data$broadLineage3 %in% c("MEMP_MEP",'earlyMK','Megakaryocytes','Mast.cell','Erythroblasts')] = 'Ery/MegK/Mast lineage'
    data$broadLineage3[data$broadLineage3 %in% c("LMPP_ELP",'B lineage')] = 'B lineage'
    
    
    
    ##---- individual donor level
    
    dd = data %>% filter(!broadLineage %in% c('Stromal','others')) %>% group_by(donorID,broadLineage3,Genotype) %>% summarise(nCell = n_distinct(cellID)) %>% 
      group_by(donorID) %>% mutate(totalCell = sum(nCell),ctFrac = nCell/totalCell) 
    dd$Genotype = factor(dd$Genotype,(c('Diploid','T21','T18','T22','Triploid','MX')))
    dd$broadLineage3 = factor(dd$broadLineage3,c('Ery/MegK/Mast lineage','B lineage','HSC_MPP','Myeloid lineage','T/NK lineage'))
    
    dd2 = dd %>% group_by(Genotype,broadLineage3) %>% summarise(mean_ctFrac = mean(ctFrac))
    
    library(ggbeeswarm)
    library(ggpubr)
    library(ggbreak) 
    # p1 = ggplot(dd[dd$broadLineage3 == 'B lineage',],aes(x = Genotype,y=ctFrac,fill=broadLineage3))+
    #   geom_col(data=dd2[dd2$broadLineage3 == 'B lineage',],aes(Genotype,mean_ctFrac,fill=Genotype),width = 0.7,alpha=0.8)+
    #   #scale_color_manual(values = col25)+
    #   scale_fill_manual(values = geno_cols)+
    #   geom_quasirandom(size=0.6,width = 0.2)+
    #   theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
    #   xlab('') + ylab('Celltype contribution') 
    # print(p1)
    
    if(!noPlot){
      
      p1 = ggplot(dd[dd$broadLineage3 %in% c('B lineage','Ery/MegK/Mast lineage'),],aes(x = Genotype,y=ctFrac,fill=broadLineage3))+
        geom_col(data=dd2[dd2$broadLineage3 %in% c('B lineage','Ery/MegK/Mast lineage'),],aes(Genotype,mean_ctFrac,fill=Genotype),width = 0.7,alpha=0.8)+
        #scale_color_manual(values = col25)+
        scale_fill_manual(values = geno_cols)+
        geom_quasirandom(size=0.2,width = 0.4)+
        facet_wrap(vars(broadLineage3),scales = 'free_y')+
        #scale_y_break(c(0.2,0.7), scales = 0.5)+
        theme_classic(base_size = 8.5)+ 
        theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
              strip.background = element_blank(),
              axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
              axis.text = element_text(colour = 'black'),axis.ticks = element_line(colour = 'black')) + 
        xlab('') + ylab('Fraction of haematopoietic cells') 
      print(p1)  
    }
    
    
  }
  
  saveFig(file.path(plotDir,'Fig1_2nAKLiv_mean.celltypeFractionByDonorByLineage_barPlot'),plotFun_celltypeContribution_barPlot_v2,rawData=dd,width = 4.4,height = 2.2,res = 500,useDingbats = F)
  
  
  ## Attempt to do a t-test on cell fraction
  # t.test(dd$ctFrac[dd$Genotype == 'T21' & dd$broadLineage3 == 'Ery/MegK/Mast lineage'],
  #        dd$ctFrac[dd$Genotype == 'Diploid' & dd$broadLineage3 == 'Ery/MegK/Mast lineage'],alternative = 'greater')
  wilcox.test(dd$ctFrac[dd$Genotype == 'T21' & dd$broadLineage3 == 'Ery/MegK/Mast lineage'],
              dd$ctFrac[dd$Genotype == 'Diploid' & dd$broadLineage3 == 'Ery/MegK/Mast lineage'],alternative = 'greater')
  wilcox.test(dd$ctFrac[dd$Genotype == 'T21' & dd$broadLineage3 == 'B lineage'],
              dd$ctFrac[dd$Genotype == 'Diploid' & dd$broadLineage3 == 'B lineage'],alternative = 'less')
  
  
  ## Binomial test - Triploid case
  ## expected probability of Ery/MK/Mast cells is taken as the mean proportion across all diploid samples
  ## n_trials = number of haematopoietic cells detected in Triploid case
  ## n_success = number of Ery/MK/Mast cells 
  expected_p = mean(dd$ctFrac[dd$Genotype == 'Diploid' & dd$broadLineage3 == 'Ery/MegK/Mast lineage'])
  n_trials = unique(dd$totalCell[dd$Genotype == 'Triploid'])
  k_successes = dd$nCell[dd$Genotype == 'Triploid' & dd$broadLineage3 =='Ery/MegK/Mast lineage']
  binom.test(x=k_successes, n=n_trials, p = expected_p,alternative = c("greater"))
  
  # For T21
  expected_p = mean(dd$ctFrac[dd$Genotype == 'Diploid' & dd$broadLineage3 == 'Ery/MegK/Mast lineage'])
  n_trials = sum(unique(dd$totalCell[dd$Genotype == 'T21']))
  k_successes = sum(dd$nCell[dd$Genotype == 'T21' & dd$broadLineage3 =='Ery/MegK/Mast lineage'])
  binom.test(x=k_successes, n=n_trials, p = expected_p,alternative = c("greater"))
  
  ##---- B lineage
  expected_p = mean(dd$ctFrac[dd$Genotype == 'Diploid' & dd$broadLineage3 == 'B lineage'])
  n_trials = unique(dd$totalCell[dd$Genotype == 'Triploid'])
  k_successes = dd$nCell[dd$Genotype == 'Triploid' & dd$broadLineage3 =='B lineage']
  binom.test(x=k_successes, n=n_trials, p = expected_p,alternative = c("less"))
  
  # For T21
  expected_p = mean(dd$ctFrac[dd$Genotype == 'Diploid' & dd$broadLineage3 == 'B lineage'])
  n_trials = sum(unique(dd$totalCell[dd$Genotype == 'T21']))
  k_successes = sum(dd$nCell[dd$Genotype == 'T21' & dd$broadLineage3 =='B lineage'])
  binom.test(x=k_successes, n=n_trials, p = expected_p,alternative = c("less"))
  
  
  # #--- OLD -----
  # plotFun_celltypeContribution_barPlot = function(noFrame=FALSE,noPlot=FALSE){
  #   
  #   data$broadLineage3 = as.character(data$broadLineage2)
  #   data$broadLineage3[data$broadLineage3 %in% c("CMP_GMP",'Monocyte/Macrophage','Dendritic cells','Myelocytes')] = 'Myeloid lineage'
  #   data$broadLineage3[data$broadLineage3 %in% c("MEMP_MEP",'Megakaryocytes','Mast.cell','Erythroblasts')] = 'Ery/MegK/Mast lineage'
  #   data$broadLineage3[data$broadLineage3 %in% c("LMPP_ELP",'B lineage')] = 'B lineage'
  #   
  #   
  #   
  #   
  #   dd = data %>% group_by(donorID,broadLineage3,Genotype) %>% summarise(nCell = n_distinct(cellID)) %>% 
  #     group_by(donorID) %>% mutate(totalCell = sum(nCell),ctFrac = nCell/totalCell) %>% 
  #     group_by(Genotype,broadLineage3) %>% summarise(mean_ctFrac = mean(ctFrac))
  #   
  #   
  #   dd$Genotype = factor(dd$Genotype,(c('Diploid','T21','T18','T22','Triploid','MX')))
  #   dd$broadLineage3 = factor(dd$broadLineage3,rev(c('HSC_MPP','B lineage','Ery/MegK/Mast lineage','Myeloid lineage','T/NK lineage','Stromal')))
  #   
  #   p1 = ggplot(dd,aes(x = Genotype,y=mean_ctFrac,fill=broadLineage3))+
  #     geom_col(width = 0.8)+
  #     #scale_color_manual(values = col25)+
  #     scale_fill_manual(values = broadLin_cols_v3)+
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #     xlab('') + ylab('Celltype contribution')
  #   print(p1)
  #   
  #   
  #   ## Line plot
  #   p1.2 = ggplot(dd,aes(x = Genotype,y=mean_ctFrac,col=broadLineage3))+
  #     geom_point()+
  #     geom_line(aes(x= as.numeric(Genotype)))+
  #     scale_y_log10()+
  #     scale_color_manual(values = broadLin_cols_v3)+
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #     xlab('') + ylab('Celltype contribution')
  #   #print(p1.2)
  #   
  #   
  #   
  #   ## Version 2: all Geno, without stromal cells
  #   dd_noStromal = data[!data$broadLineage %in% c('Stromal','others'),]
  #   dd_noStromal = dd_noStromal %>% group_by(donorID,broadLineage3,Genotype) %>% summarise(nCell = n_distinct(cellID)) %>% 
  #     group_by(donorID) %>% mutate(totalCell = sum(nCell),ctFrac = nCell/totalCell) %>% 
  #     group_by(Genotype,broadLineage3) %>% summarise(mean_ctFrac = mean(ctFrac))
  #   
  #   dd_noStromal$Genotype = factor(dd_noStromal$Genotype,(c('Diploid','T21','T18','T22','Triploid','MX')))
  #   dd_noStromal$broadLineage3 = factor(dd_noStromal$broadLineage3,rev(c('HSC_MPP','B lineage','Ery/MegK/Mast lineage','Myeloid lineage','T/NK lineage','Stromal')))
  #   broadLin_cols_v3 = c('HSC_MPP' = '#70480b',"Ery/MegK/Mast lineage" = '#D3A292',"Myeloid lineage" = '#38475b',
  #                        'B lineage'='#2e4d40','T/NK lineage'='#8F8681','Stromal'='#BB885E')
  #   p2 = ggplot(dd_noStromal,aes(x = Genotype,y=mean_ctFrac,fill=broadLineage3))+
  #     geom_col(width = 0.8)+
  #     #scale_color_manual(values = col25)+
  #     scale_fill_manual(values = broadLin_cols_v3)+
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #     xlab('') + ylab('Celltype contribution')
  #   print(p2)
  #   
  #   ## Line plot
  #   p2.2 = ggplot(dd_noStromal,aes(x = Genotype,y=mean_ctFrac,col=broadLineage3))+
  #     geom_point()+
  #     geom_line(aes(x= as.numeric(Genotype)))+
  #     scale_y_log10()+
  #     scale_color_manual(values = broadLin_cols_v3)+
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #     xlab('') + ylab('Celltype contribution')
  #   #print(p2.2)
  # }
  # 
  # saveFig(file.path(plotDir,'Figxx_2nAKLiv_mean.celltypeFractionByDonor_barPlot_v5'),plotFun_celltypeContribution_barPlot,rawData=dd,width = 3.4,height = 2,res = 500,useDingbats = F)
  # 
  # 
  # 
  # dd_noStromal = data[data$broadLineage != 'Stromal',]
  # dd_noStromal_t21_2n = data[data$Genotype %in% c('T21','Diploid') & grepl('Hsb',data$donorID) & data$broadLineage != 'Stromal',]
  # 
  # dd = data %>% group_by(donorID,broadLineage2,Genotype) %>% summarise(nCell = n_distinct(cellID)) %>% 
  #   group_by(donorID) %>% mutate(totalCell = sum(nCell),ctFrac = nCell/totalCell) %>% 
  #   group_by(Genotype,broadLineage2) %>% summarise(mean_ctFrac = mean(ctFrac))
  # 
  # dd_noStromal = dd_noStromal %>% group_by(donorID,broadLineage2,Genotype) %>% summarise(nCell = n_distinct(cellID)) %>% 
  #   group_by(donorID) %>% mutate(totalCell = sum(nCell),ctFrac = nCell/totalCell) %>% 
  #   group_by(Genotype,broadLineage2) %>% summarise(mean_ctFrac = mean(ctFrac))
  # 
  # dd_noStromal_t21_2n = dd_noStromal_t21_2n %>% group_by(donorID,broadLineage2,Genotype) %>% summarise(nCell = n_distinct(cellID)) %>% 
  #   group_by(donorID) %>% mutate(totalCell = sum(nCell),ctFrac = nCell/totalCell) %>% 
  #   group_by(Genotype,broadLineage2) %>% summarise(mean_ctFrac = mean(ctFrac))
  # 
  # plotFun_celltypeContribution_linePlot = function(noFrame=FALSE,noPlot=FALSE){
  #   broadLin_cols = c('HSC_MPP' = '#FFD700',"MEMP_MEP" = '#c43c1f',"CMP_GMP" = '#B5338A',"LMPP_ELP" = '#8DB5CE',
  #                     'Megakaryocytes' ='#9C6DA5','Mast.cell'='#FFB7CE','Erythroblasts'='#E9967A',
  #                     'Monocyte/Macrophage'='#AE8094','Dendritic cells'='#945a01','Myelocytes'='#DABE99',
  #                     'B lineage'=pal37H[c(19)],'T/NK lineage'=grey(c(0.6)),'Stromal'=pal37H[c(32)])
  #   
  #   dd$broadLineage2 = factor(dd$broadLineage2, levels = c('HSC_MPP',"CMP_GMP","MEMP_MEP","LMPP_ELP",'Megakaryocytes','Mast.cell','Erythroblasts',
  #                                                          'Monocyte/Macrophage','Dendritic cells','Myelocytes',
  #                                                          'B lineage','T/NK lineage','Stromal'))
  #   
  #   
  #   #dd$Genotype = factor(dd$Genotype,rev(c('MX','Diploid','T21','T18','T22','Triploid')))
  #   dd$Genotype = factor(dd$Genotype,rev(c('Diploid','T21','T18','T22','Triploid','MX')))
  #   #dd$Genotype = factor(dd$Genotype,rev(c('Diploid','T21')))
  #   
  #   
  #   dd$lin_of_interest = ifelse(dd$broadLineage2 %in% c('B lineage','LMPP_ELP','MEMP_MEP'),T,F)
  #   p1 = ggplot(dd,aes(x = Genotype,y=mean_ctFrac))+
  #     geom_point(size=1,aes(col=broadLineage2))+
  #     geom_line(aes(x=as.numeric(Genotype),col=broadLineage2))+
  #     geom_line(data=dd[dd$lin_of_interest == T,],aes(x=as.numeric(Genotype),col=broadLineage2),lwd=1.6,alpha=0.5)+
  #     scale_color_manual(values = broadLin_cols)+
  #     scale_y_log10() +
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #     xlab('') + ylab('Celltype contribution')
  #   
  #   print(p1)
  #   
  #   
  #   ## Version 2: all Geno, without stromal cells
  #   dd_noStromal$lin_of_interest = ifelse(dd_noStromal$broadLineage2 %in% c('B lineage','LMPP_ELP','MEMP_MEP'),T,F)
  #   dd_noStromal$Genotype = factor(dd_noStromal$Genotype,rev(c('MX','Diploid','T21','T18','T22','Triploid')))
  #   p2 = ggplot(dd_noStromal,aes(x = Genotype,y=mean_ctFrac))+
  #     geom_point(size=1,aes(col=broadLineage2))+
  #     geom_line(aes(x=as.numeric(Genotype),col=broadLineage2))+
  #     geom_line(data=dd_noStromal[dd_noStromal$lin_of_interest == T,],aes(x=as.numeric(Genotype),col=broadLineage2),lwd=1.6,alpha=0.5)+
  #     scale_color_manual(values = broadLin_cols)+
  #     scale_y_log10() +
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #     xlab('') + ylab('Celltype contribution')
  #   
  #   print(p2)
  #   
  #   
  #   
  #   
  #   
  #   ## Version 3: only CD45 sorted T21 + 2n, without stromal cells
  #   dd_noStromal_t21_2n$lin_of_interest = ifelse(dd_noStromal_t21_2n$broadLineage2 %in% c('B lineage','LMPP_ELP','MEMP_MEP'),T,F)
  #   dd_noStromal_t21_2n$Genotype = factor(dd_noStromal_t21_2n$Genotype,rev(c('Diploid','T21')))
  #   p3 = ggplot(dd_noStromal_t21_2n,aes(x = Genotype,y=mean_ctFrac))+
  #     geom_point(size=1,aes(col=broadLineage2))+
  #     geom_line(aes(x=as.numeric(Genotype),col=broadLineage2))+
  #     geom_line(data=dd_noStromal_t21_2n[dd_noStromal_t21_2n$lin_of_interest == T,],aes(x=as.numeric(Genotype),col=broadLineage2),lwd=1.8,alpha=0.7)+
  #     scale_color_manual(values = broadLin_cols)+
  #     scale_y_log10() +
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank()) + 
  #     xlab('') + ylab('Celltype contribution')
  #   
  #   print(p3)
  #   
  #   # ggplot(dd,aes(x = Genotype,y=log(mean_ctFrac)))+
  #   #   geom_point(size=1,aes(col=finalAnn_broad))+
  #   #   geom_line(aes(x=as.numeric(Genotype),col=finalAnn_broad))+
  #   #   scale_color_manual(values = fLiver_celltype_cols)+
  #   #   theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #   #   xlab('') + ylab('Celltype contribution')
  #   # 
  #   # p2 = ggplot(dd,aes(y= broadLineage,fill=Genotype))+
  #   #   geom_bar(position='fill')+
  #   #   scale_fill_manual(values = geno_cols)+
  #   #   theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank()) +
  #   #   ylab('') + xlab('Genotype contribution')
  #   # 
  #   # print(p2)
  #   
  # }
  # 
  # saveFig(file.path(plotDir,'Figxx_2nAKLiv_mean.celltypeFractionByDonor'),plotFun_celltypeContribution,rawData=dd,width = 4.4,height = 3,res = 500,useDingbats = F)
  # 
  # 
  # 
  # 
  # 
  # 
  # plotFun_celltypeContribution_plotForTTest = function(noFrame=FALSE,noPlot=FALSE){
  #   dd$broadLineage3 = as.character(dd$broadLineage2)
  #   dd$broadLineage3[dd$broadLineage3 %in% c("CMP_GMP",'Monocyte/Macrophage','Dendritic cells','Myelocytes')] = 'Myeloid lineage'
  #   dd$broadLineage3[dd$broadLineage3 %in% c("MEMP_MEP",'Megakaryocytes','Mast.cell','Erythroblasts')] = 'Ery/MegK/Mast lineage'
  #   dd$broadLineage3[dd$broadLineage3 %in% c("LMPP_ELP",'B lineage')] = 'B lineage'
  #   broadLin_cols_v3 = c('HSC_MPP' = '#FFD700',"Ery/MegK/Mast lineage" = '#c43c1f',"Myeloid lineage" = '#70480b',
  #                        'B lineage'=pal37H[c(19)],'T/NK lineage'=grey(c(0.7)),'Stromal'=pal37H[c(32)])
  #   dd$Genotype = factor(dd$Genotype,rev(c('Diploid','T21','T18','T22','Triploid','MX')))
  #   
  #   dd$broadLineage3 = factor(dd$broadLineage3,rev(c('HSC_MPP','B lineage','Ery/MegK/Mast lineage','Myeloid lineage','T/NK lineage','Stromal')))
  #   p1 = ggplot(dd,aes(x = Genotype,y=mean_ctFrac,fill=broadLineage3))+
  #     geom_point()+
  #     #scale_color_manual(values = col25)+
  #     scale_fill_manual(values = broadLin_cols_v3)+
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #     xlab('') + ylab('Celltype contribution')
  #   print(p1)
  #   
  #   ## Version 2: all Geno, without stromal cells
  #   dd_noStromal$Genotype = factor(dd_noStromal$Genotype,rev(c('Diploid','T21','T18','T22','Triploid','MX')))
  #   
  #   dd_noStromal$broadLineage3 = dd$broadLineage3[match(dd_noStromal$broadLineage2,dd$broadLineage2)]
  #   p2 = ggplot(dd_noStromal,aes(x = Genotype,y=mean_ctFrac,fill=broadLineage3))+
  #     geom_col(width = 0.8)+
  #     #scale_color_manual(values = col25)+
  #     scale_fill_manual(values = broadLin_cols_v3)+
  #     theme_bw(base_size = 8.5)+ theme(panel.grid = element_blank(),axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  #     xlab('') + ylab('Celltype contribution')
  #   print(p2)
  #   
  # }
  # 
  # saveFig(file.path(plotDir,'Figxx_2nAKLiv_mean.celltypeFractionByDonor_barPlot'),plotFun_celltypeContribution_barPlot,rawData=dd,width = 3.4,height = 2,res = 500,useDingbats = F)
  # 
  # 
  # #--- OLD -----
  
}  















fLiver_markers = unique(c(
  'CD34','SPINK2','HLF',
  #'CD38',#'MLLT3','PRSS57', # HSC_MPP
  #'CTSG',	'PRTN3', # CMP
  #'SERPINB1', 

  # MEMP - MEP
  'IL1B','GATA2','GATA1','KLF1','TESPA1',#'CTNNBL1',
  'HBD','CSF2RB',
  
  'FLI1','FCER1A','ITGA2B','PLEK', # MK
  'PF4','ITGA2B', #Megakaryocyte
  'PPBP','TUBB1', # Platelets
  
  'KIT','HDC','TPSAB1', # Mast.cell
  
  'EPCAM','APOC1', # early Erythroid
  'ALAS2', # mid.erythroid
  'HBB','BPGM', # late.erythroid
  
  
  # CMP-GMP
  #'CLEC11A','CALR', # not specific enough
  'AZU1','MPO','LYZ',
  # Monocytes
  'CD14',#'FCGR3A',
  #'S100A4', 'S100A6', # not specific enough
  'CD52', 'SELL', # not specific enough
  'C1QA','CD68',#'MSR1',#Macrophages
  #Kupffer
  'FABP3','CD5L','VCAM1',
  
  #'CCR2','CX3CR1',
  
  'IRF8',	 #DC.precursor
  'CLEC9A',#DC1
  'CLEC10A','CD1C', # DC2 
  'CLEC4C','IL3RA', #pDC
  
  
  #'CSF3R','FPR1','FCGR3B',
  'MNDA', # NEUTROPHILS
  'DEFA3','DEFA4','ELANE', # pro-myelocytes
  'CAMP','LCN2', #myelocytes
  #'CXCR2',
  #'CSF3R','FCGR3B',#'FUT4', # Neutrophil
  
  #'IGLL1','CD99', # Pre-pro
  'DNTT',
  #'EBF1',
  'CD19','RAG1',# pro-b-cells
  'MME','VPREB1','CD79A','CD79B',# pre-b-cells
  'TCL1A','MME',#'RAG1',
  'MS4A1',  # B-cells
  #'CD27',#plasma cell?
  
  'CD52','IL7R',# ILC precursor
  'ZBTB16',#'LTB', 
  'CD3D',#'GZMA',
  #'CD4','CD8A', #Early.lymphoid_T.lymphocyte
  #'TRDV2','TRGV9', # gamma delta T-cell
  #'SLC4A10','TRAV1-2', #MAIT t-cell
  #'PRF1', # effector T cells
  #'FOXP3',	'CDH1', # regulatory T cells
  'NKG7','KLRD1', #NK
  
  'PTPRC',
  #'ESAM','PECAM1', 
  
  # Endothelial
  'KDR', #'PTPRB', 
  #'PLVAP', # Endothelium
  #'DCN','SERPINF1',
  'COL3A1',#'COL1A1',
  #'BGN','VIM','ECM1', # endosteal fibroblast / fibroblast
  #'APOA1',
  #'SCD',
  'ALB','TTR','GSTA1', # Hepatocyte
  'SOX9','CHST9','CDH6', 'KRT19', # Cholangiocyte
  'UPK3B','MSLN','KLK11', # Mesothelial cells
  'SOX11','NNAT','DCX',# Neuron
  'PAGE4','XAGE3','PHLDA2' # Trophoblast
))

fig2b_2nAK_Liv_dotPlot = function(){
  library(viridis)
  
  # Import the seurat object
  srat = readRDS(akLiv_srat_fp)
  mdat = read.csv(akLiv_mdat_fp)
  srat$broadLineage = mdat$broadLineage[match(srat$cellID,mdat$cellID)]
  
  # Import annot 
  srat$broadLineage = factor(srat$broadLineage, levels = c('HSC & prog.','Megakaryocytes','Mast.cell','Erythroblasts',
                                                           'Monocyte/Macrophage','Dendritic cells','Myelocytes',
                                                           'B lineage','T/NK lineage','Stromal','others'))
  
  srat$broadLineage2 = factor(as.numeric(srat$broadLineage),c(1:length(levels(srat$broadLineage))))
  
  direction = 'vertical'
  if(direction == 'vertical'){
    genes = rev(fLiver_markers)
  }else{
    genes = fLiver_markers
  }
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    Idents(srat) = srat$broadLineage2
    p = DotPlot(srat,idents = unique(srat$broadLineage2),
                #cols = c("#EBFFE5", "#244b05"),
                #cols = c(colAlpha('#F1F5FA',1),'#425580'),
                cols = c(colAlpha(grey(0.95),0.8),'black'),
                #cols = c(grey(0.99), grey(0.2)),
                #group.by = 'seurat_clusters',
                #idents = unique(srat$finalAnn[srat$finalAnn != 'others']),
                features = genes)+RotatedAxis() 
    
    if(direction == 'vertical'){
      p = p + coord_flip() + 
        scale_y_discrete(position = "right")+
        theme(axis.text.x = element_text(size=11,angle = 0,vjust = 0,hjust = 0.5),
              axis.text.y = element_text(size=10),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top') + xlab('') + ylab('')
    }else if(direction == 'horizontal'){
      p = p +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top') + xlab('') + ylab('')
    }
    print(p)
  }
  
  if(direction == 'vertical'){
    saveFig(file.path(plotDir,'Fig1c_2nAK_Liv_Lin_DotPlot_vertical'),plotFun,width = 5,height = 11.7,res = 500)  
  }else{
    saveFig(file.path(plotDir,'Fig1c_2nAK_Liv_Lin_DotPlot_horizontal'),plotFun,width = 9.5,height = 4.6,res = 500)  
  }
  
}






figS1_2nAK_Liv_dotPlot_byCT = function(){
  library(viridis)
  
  # Import the seurat object
  srat = readRDS(akLiv_srat_fp)
  mdat = read.csv(akLiv_mdat_fp)
  ## Update labels
  srat$finalAnn_broad[srat$finalAnn_broad == 'Neuron'] = 'Neuronal'
  srat$finalAnn_broad[srat$finalAnn_broad == 'earlyMK'] = 'early MK'
  srat$finalAnn_broad[srat$finalAnn_broad == 'EE'] = 'early Ery'
  srat$finalAnn_broad[srat$finalAnn_broad == 'ME'] = 'mid Ery'
  srat$finalAnn_broad[srat$finalAnn_broad == 'LE'] = 'late Ery'
  srat$finalAnn_broad[srat$finalAnn_broad == 'proMono'] = 'Pro-monocyte'
  srat$finalAnn_broad[srat$finalAnn_broad == 'myelocyte'] = 'Myelocyte'
  srat$finalAnn_broad[srat$finalAnn_broad == 'Mesothelial_cells'] = 'Mesothelial cell'
  srat$finalAnn_broad[srat$finalAnn_broad == 'Cholangiocytes'] = 'Cholangiocyte'
  
  srat$finalAnn_broad = gsub('_',' / ',srat$finalAnn_broad)
  srat$finalAnn_broad = gsub('\\.',' ',srat$finalAnn_broad)
  srat$finalAnn_broad[srat$finalAnn_broad != 'pDC'] = gsub('^p','P',srat$finalAnn_broad[srat$finalAnn_broad != 'pDC'])
  
  srat$finalAnn_broad = factor(srat$finalAnn_broad,levels = c("HSC / MPP","MEMP / MEP",
                                                              "early MK",'MK', # purple
                                                              'Mast cell', # slamon
                                                              'early Ery','mid Ery','late Ery', # red
                                                              "CMP / GMP",
                                                              "Pro-monocyte","Monocyte","Macrophage","Kupffer cell", # blue
                                                              "DC1","DC2","pDC", # yellow
                                                              "Promyelocyte","Myelocyte", # orange
                                                              "LMPP / ELP", "Pro B cell","Pre B cell","B cell", # pink / purple
                                                              "ILC precursor","T cell","NK / T", # greys
                                                              "Endo","Fibroblast","Hepatocyte","Mesenchyme",'Cholangiocyte','Mesothelial cell','Neuronal',"NPC",'Trophoblast')) # green
  
  
  
  direction = 'vertical'
  if(direction == 'vertical'){
    genes = rev(fLiver_markers)
  }else{
    genes = fLiver_markers
  }
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    Idents(srat) = srat$finalAnn_broad
    p = DotPlot(srat,idents = unique(srat$finalAnn_broad),
                #cols = c("#EBFFE5", "#244b05"),
                #cols = c(colAlpha('#F1F5FA',1),'#425580'),
                cols = c(colAlpha(grey(0.95),0.8),'black'),
                #cols = c(grey(0.99), grey(0.2)),
                #group.by = 'seurat_clusters',
                #idents = unique(srat$finalAnn[srat$finalAnn != 'others']),
                features = genes)+RotatedAxis() 
    
    if(direction == 'vertical'){
      p = p + coord_flip() + 
        scale_y_discrete(position = "right")+
        theme(axis.text.x = element_text(size=11,angle = 90,vjust = 1,hjust = 0,colour = 'black'),
              axis.text.y = element_text(size=10,face='italic',hjust = 1),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top',
              axis.ticks = element_line(colour = 'black'),
              panel.background = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank()) + xlab('') + ylab('')
    }else if(direction == 'horizontal'){
      p = p +
        theme(axis.text.y = element_text(size=12),
              axis.text.x = element_text(size=10,angle = 90,vjust = 0.5,hjust = 1,face='italic',colour = 'black'),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top',
              axis.ticks = element_line(colour = 'black'),
              panel.background = element_rect(fill=F,colour = 'black',linewidth = 1),axis.line = element_blank()) + xlab('') + ylab('')
    }
    print(p)
  }
  
  if(direction == 'vertical'){
    saveFig(file.path(plotDir,'FigS1_2nAK_Liv_CT_DotPlot_vertical'),plotFun,width = 8,height = 15.5,res = 500)  
  }else{
    saveFig(file.path(plotDir,'FigS1_2nAK_Liv_CT_DotPlot_horizontal'),plotFun,width = 14,height = 7.6,res = 500)  
  }
  
}












lineageCol = c('B.cell'=colAlpha('#BBDDFF',0.01),
               #'EE' = '#d10808',
               'LE' = '#E9967A',
               'MK' = '#8d439c',
               'Mast.cell' = '#D360A3',
               'Monocyte' = pal37H[c(26)]
)


fig2D_CR_fateProbs_UMAP = function(){
  
  ## Import cellrank output
  mdat = read.csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/fLiver_2nAK_nonDS_palantired_CR_mdat_2404.csv')
  umap = read.csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/fLiver_2nAK_nonDS_palantired_CR_umap_2404.csv')
  umap$cellID = mdat$cellID
  mdat = merge(mdat[,!colnames(mdat) %in% c('UMAP_1','UMAP_2')],umap[,colnames(umap) != 'X'],by='cellID',all=T)
  
  
  plotFun_CR_fateProbs_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.5,0.5,1,0.5))
    plot(mdat$UMAP_1,mdat$UMAP_2,type='n',
         cex.main = 0.35,xaxt='n',yaxt='n',
         xlab='UMAP 1',ylab=('UMAP 2'),
         frame.plot=T)
    
    for(termState in unique(mdat$clusters_gradients)){
      d = mdat[mdat$clusters_gradients == termState,]
      #col_fun = colorRamp2(c(0.3, 0.8), c(grey(0.9),colAlpha(lineageCol[termState],0.2)))
      col_fun = colorRamp2(c(0.3,0.8), c(grey(0.9),colAlpha(lineageCol[termState],1)))
      #col_fun = colorRamp2(c(0.3, 0.6), c(colAlpha(lineageCol[termState],0.1),colAlpha(lineageCol[termState],0.001)))
      cols = col_fun(d$term_states_fwd_probs)
      
      points(d$UMAP_1,d$UMAP_2,pch=19,cex=0.05,col=cols)
    }
  }
  
  saveFig(file.path(plotDir,'Figxx_fLiver_2nAK_CRfateProb_UMAP'),plotFun_CR_fateProbs_UMAP,rawData=mdat,width = 3,height = 3,res = 500,useDingbats = F)
  #saveFig(file.path(plotDir,'Figxx_fLiver_2nT21_CRfateProb_UMAP'),plotFun_CR_fateProbs_UMAP,rawData=mdat,width = 3,height = 3,res = 500,useDingbats = F)
  
  mdat$finalAnn_broad = mdat$annot_mar24
  data = mdat[,c('cellID','Genotype','finalAnn_broad','UMAP_1','UMAP_2','umap_density_Genotype')]
  plotFun_CR_GenotypeDensity_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    plot_list = list()
    for(geno in c('diploid','T21','T18','T22','complete_trisomy','MX')){
      d = data[data$Genotype == geno,]
      if(geno == 'diploid'){
        g = 'Diploid'
      }else if(geno == 'complete_trisomy'){
        g = 'Triploid'
      }else{
        g = geno
      }
      
      p1 = ggplot(data,aes(UMAP_1,UMAP_2))+
        geom_point(size=0.001,col=grey(0.8))+
        geom_point(data=d,aes(col=`umap_density_Genotype`),size=0.1)+
        scale_color_gradient(low = brewer.pal('OrRd',n = 9)[1],high = brewer.pal('OrRd',n = 9)[9],name='UMAP density')+
        theme_classic(base_size = 5) + 
        ggtitle(g) + 
        theme(panel.border = element_rect(fill=F,linewidth = 0),
              legend.key.height = unit(0.4,'cm'),
              legend.key.width = unit(0.12,'cm'),
              plot.title = element_text(hjust = 0.5),
              axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + 
        xlab('') + ylab('') 
      
      plot_list[[geno]] = p1
    }
    print(patchwork::wrap_plots(plot_list,nrow=1))
  }
  
  saveFig(file.path(plotDir,'Figxx_fLiver_2nAK_GenotypeDensity_UMAP'),plotFun_CR_GenotypeDensity_UMAP,rawData=data,width = 10,height = 1.2,res = 500,useDingbats = F)
  
  
  plotFun_CR_GenotypeDensity_UMAP_greyBackground = function(noFrame=FALSE,noPlot=FALSE){
    
    p1 = ggplot(data,aes(UMAP_1,UMAP_2))+
      geom_point(size=0.001,col=grey(0.8))+
      #geom_point(data=d,aes(col=`umap_density_Genotype`),size=0.1)+
      #scale_color_gradient(low = brewer.pal('OrRd',n = 9)[1],high = brewer.pal('OrRd',n = 9)[9],name='UMAP density')+
      theme_classic(base_size = 5) + 
      theme(panel.border = element_rect(fill=F,linewidth = 0),
            legend.key.height = unit(0.4,'cm'),
            legend.key.width = unit(0.12,'cm'),
            plot.title = element_text(hjust = 0.5),
            axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + 
      xlab('') + ylab('') 
    print(p1)
  }
  
  saveFig(file.path(plotDir,'Figxx_fLiver_2nAK_GenotypeDensity_UMAP_greyBackground'),plotFun_CR_GenotypeDensity_UMAP_greyBackground,rawData=data,width = 2.4,height = 2.2,res = 500,useDingbats = F)
  
  
  
  
  ##--------------------------------------------------------------##
  ##    Plot fate probability for HSC in different genotypes    ####
  ##--------------------------------------------------------------##
  fateProb_cr = read.csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/9c_fLiver_trajectory/1_palantir_output/fLiver_2nAK_palantired_CR_fateProbs_2404.csv',row.names = 1)
  fateProb_cr = cbind(mdat,fateProb_cr)
  # subset to just HSC_MPP
  df = fateProb_cr[fateProb_cr$annot_mar24 == 'HSC_MPP',]
  table(df$clusters_gradients)
  df$clusters_gradients[df$clusters_gradients %in% c('LE','MK','EE')] = 'Megk/Ery'
  df$Genotype[df$Genotype == 'complete_trisomy'] = 'Triploid'
  df$Genotype[df$Genotype == 'diploid'] = 'Diploid'
  df$Genotype = factor(df$Genotype,c('Diploid','T21','T18','MX','Triploid','T22'))
  
  df = df %>% group_by(Genotype,clusters_gradients) %>% summarise(nCell = n_distinct(cellID)) %>% 
    group_by(Genotype) %>% mutate(totalCell = sum(nCell),
                                  frac = nCell/totalCell)
  
  df$clusters_gradients = factor(df$clusters_gradients,c('Megk/Ery','B.cell','Monocyte'))
  
  
  plot_CR_hscFateProb_allGeno_lineGraph = function(noFrame=FALSE,noPlot=FALSE){
    p1 = ggplot(df,aes(Genotype,frac,col=clusters_gradients,fill=clusters_gradients))+
      geom_point()+
      geom_line(aes(x=as.numeric(Genotype))) + 
      theme_classic(base_size = 12)+
      scale_fill_manual(values = col25[c(2,1,3)])+
      scale_color_manual(values = col25[c(2,1,3)]) + 
      theme(axis.text.x =  element_text(angle=90,vjust = 0.5,hjust = 1))+
      xlab('') + ylab('Fraction of HSC/MPP cells')
    print(p1)
  }
  
  
  df = df[df$Genotype %in% c('T21','Diploid'),]
  
  lineageCol = c('B.cell'='#006400',
                 'Megk/Ery' = colAlpha('#ffb703',0.7),
                 'Monocyte' = '#D360A3'
  )
  
  plot_CR_hscFateProb_2nT21_stackedBarPlot = function(noFrame=FALSE,noPlot=FALSE){
    p1 = ggplot(df,aes(Genotype,frac,fill=clusters_gradients))+
      geom_col(width = 0.6)+
      theme_classic(base_size = 14)+
      scale_y_continuous(labels = c(0,'',0.5,'',1))+
      scale_fill_manual(values = lineageCol)+
      theme(axis.text.x =  element_text(angle=90,vjust = 0.5,hjust = 1),
            panel.border = element_rect(fill=F),axis.line = element_blank())+
      xlab('') + ylab('Fraction of HSC/MPP cells')
    print(p1)
  }
  
  saveFig(file.path(plotDir,'CR_HSC.MPP_lineage_2nT21_barplot'),plotFun = plot_CR_hscFateProb_2nT21_stackedBarPlot,rawData=df,width = 3.5,height = 2.5,res = 500,useDingbats = F)
  
}








supFig2c_2nAK_Liv_LRv2_heatmap = function(){
  
  if(file.exists(file.path(plotDir,'Fig1b_2nAKLiv_geno_UMAP_rawData.tsv'))){
    dd = read.delim(file.path(plotDir,'Fig1b_2nAKLiv_geno_UMAP_rawData.tsv'),header = T,sep = '\t')
    
  }else{
    mdat = read.csv(akLiv_mdat_fp)
    
    ## Update labels
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'Neuron'] = 'Neuronal'
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'earlyMK'] = 'early MK'
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'EE'] = 'early Ery'
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'ME'] = 'mid Ery'
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'LE'] = 'late Ery'
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'proMono'] = 'Pro-monocyte'
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'myelocyte'] = 'Myelocyte'
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'Mesothelial_cells'] = 'Mesothelial cell'
    mdat$finalAnn_broad[mdat$finalAnn_broad == 'Cholangiocytes'] = 'Cholangiocyte'
    
    mdat$finalAnn_broad = gsub('_',' / ',mdat$finalAnn_broad)
    mdat$finalAnn_broad = gsub('\\.',' ',mdat$finalAnn_broad)
    mdat$finalAnn_broad[mdat$finalAnn_broad != 'pDC'] = gsub('^p','P',mdat$finalAnn_broad[mdat$finalAnn_broad != 'pDC'])
    
    mdat$group = paste0(mdat$finalAnn_broad,':',mdat$Genotype)
    
    mdat$finalAnn_broad = factor(mdat$finalAnn_broad,levels = c("HSC / MPP","MEMP / MEP",
                                                                "early MK",'MK', # purple
                                                                'Mast cell', # slamon
                                                                'early Ery','mid Ery','late Ery', # red
                                                                "CMP / GMP",
                                                                "Pro-monocyte","Monocyte","Macrophage","Kupffer cell", # blue
                                                                "DC1","DC2","pDC", # yellow
                                                                "Promyelocyte","Myelocyte", # orange
                                                                "LMPP / ELP", "Pro B cell","Pre B cell","B cell", # pink / purple
                                                                "ILC precursor","T cell","NK / T", # greys
                                                                "Endo","Fibroblast","Hepatocyte","Mesenchyme",'Cholangiocyte','Mesothelial cell','Neuronal',"NPC",'Trophoblast')) # green
    
    ## Read in LRv2 similarity matrix
    lr_output = readRDS('~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_fetalREFliver_LRwCT_output.RDS')
    lr_output = lr_output[[1]]
    # Subset to include only relevant cells    
    table(mdat$cellID %in% c(rownames(lr_output)))
    lr_output = lr_output[rownames(lr_output) %in% mdat$cellID,]
    
    ## summarise the average similarity by cell type
    # split the matrix into celltype
    meanLR_mtx = sapply(split(1:nrow(lr_output),mdat$group[match(rownames(lr_output),mdat$cellID)]),function(e){colMeans(lr_output[e,])})
    # Remove NPC
    meanLR_mtx = meanLR_mtx[,!grepl('NPC|Neuron|Cholangiocyte|Mesothelial|Trophoblast',colnames(meanLR_mtx))]
    
    ## Fix up reference classess labels
    rownames(meanLR_mtx)[rownames(meanLR_mtx) == 'Megakaryocyte'] = 'MK'
    rownames(meanLR_mtx)[rownames(meanLR_mtx) == 'Early.Erythroid'] = 'early Ery'
    rownames(meanLR_mtx)[rownames(meanLR_mtx) == 'Mid.Erythroid'] = 'mid Ery'
    rownames(meanLR_mtx)[rownames(meanLR_mtx) == 'Late.Erythroid'] = 'late Ery'
    rownames(meanLR_mtx) = gsub('_',' / ',rownames(meanLR_mtx))
    rownames(meanLR_mtx) = gsub('\\.',' ',rownames(meanLR_mtx))
    
    ref_order = c("HSC / MPP",
                  "MEMP","MK","Mast cell",
                  "early Ery", "mid Ery","late Ery",
                  "Neutrophil myeloid progenitor",
                  "Monocyte precursor","Monocyte","Mono Mac","Kupffer Cell",'VCAM1  EI macrophage',
                  "DC precursor","DC1","DC2","pDC precursor",
                  "Pre pro B cell","pro B cell","pre B cell","B cell",
                  "ILC precursor", "Early lymphoid / T lymphocyte", "NK",
                  "Hepatocyte","Endothelial cell","Fibroblast" 
    )
    
    table(ref_order[! ref_order %in% rownames(meanLR_mtx)])
    table(rownames(meanLR_mtx)[!rownames(meanLR_mtx) %in% ref_order])
    
    
    tgt_order = levels(mdat$finalAnn_broad)
    
    tgt_order = paste0(tgt_order,':',t(replicate(expr = unique(mdat$Genotype),n = length(tgt_order))))
    tgt_order = tgt_order[tgt_order %in% colnames(meanLR_mtx)]
    table(tgt_order[!tgt_order %in% colnames(meanLR_mtx)])
    table(colnames(meanLR_mtx)[!colnames(meanLR_mtx) %in% tgt_order])
    
    
    meanLR_mtx = (meanLR_mtx[ref_order,tgt_order])
  }
  
  
  direction = 'horizontal'
  plotFun_2nAK_fLiver_LRv2_clusterLevel_hm = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    logitCols = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
    cols = circlize::colorRamp2(seq(-13,10,length.out=length(logitCols)),logitCols)
    
    geno = gsub('.*:','',colnames(meanLR_mtx))
    geno[geno == 'complete_trisomy'] = 'Triploid'
    geno[geno == 'diploid'] = 'Diploid'
    geno = factor(geno,c('Diploid','T21','T18','T22','MX','Triploid'))
    
    
    if(direction == 'vertical'){
      
      botAnno = HeatmapAnnotation(Genotype = geno,
                                  col = list(Genotype = geno_cols),
                                  annotation_legend_param = list(Genotype = list(direction = "horizontal")))
      
      hm = Heatmap(meanLR_mtx,show_row_dend = F,show_column_dend = F,
                   col = cols,
                   name=paste0('Predicted Similarity','(Logit)'),
                   show_row_names = T,show_column_names = F,cluster_rows = F,
                   cluster_columns = F,cluster_row_slices = F,
                   row_names_gp = gpar(fontsize=8),
                   column_split = mdat$finalAnn_broad[match(colnames(meanLR_mtx),mdat$group)],
                   column_title_rot = 90,row_title_gp = gpar(fontsize=8),column_title_gp = gpar(fontsize=9),
                   bottom_annotation = botAnno,
                   heatmap_legend_param = list(direction = "horizontal"),
                   row_title = 'Reference classes')
    }else{
      rowAnno = rowAnnotation(Genotype = geno,
                              col = list(Genotype = geno_cols))
      
      hm = Heatmap(t(meanLR_mtx),show_row_dend = F,show_column_dend = F,
                   col = cols,
                   name=paste0('Predicted Similarity','(Logit)'),
                   show_row_names = F,show_column_names = T,cluster_rows = F,
                   cluster_columns = F,cluster_column_slices = F,
                   row_names_gp = gpar(fontsize=8),row_title_rot = 0,row_title_gp = gpar(fontsize=9),
                   split = mdat$finalAnn_broad[match(colnames(meanLR_mtx),mdat$group)],
                   column_names_gp = gpar(fontsize=9),column_title_rot = 0,column_title_gp = gpar(fontsize=9),
                   right_annotation = rowAnno,
                   heatmap_legend_param = list(direction = "vertical"),
                   column_title = 'Reference classes')
    }
    
    #geno = gsub('.*:','',rownames(meanLR_mtx))
    # geno[geno == 'complete_trisomy'] = 'Triploid'
    # geno = factor(geno,c('diploid','T21','T18','T22','MX','Triploid'))
    # rowAnno = rowAnnotation(Genotype = geno,
    #                         col = list(Genotype = geno_cols))
    # hm = Heatmap(meanLR_mtx,show_row_dend = F,show_column_dend = F,
    #              col = cols,
    #              name=paste0('Predicted\nSimilarity\n','(Logit)'),
    #              show_row_names = F,show_column_names = T,cluster_rows = F,cluster_columns = F,
    #              row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=8),
    #              split = mdat$finalAnn_broad[match(rownames(meanLR_mtx),mdat$group)],
    #              row_title_rot = 0,row_title_gp = gpar(fontsize=8),cluster_row_slices = F,
    #              right_annotation = botAnno,
    #              column_title = 'Reference classes')
    
    
    if(!noPlot){
      if(direction == 'vertical'){
        draw(hm,heatmap_legend_side = "bottom")  
      }else{
        draw(hm)
      }
      
    }
  }
  
  if(direction == 'vertical'){
    saveFig(file.path(plotDir,'Figxx_2nAKLiv_LRv2_heatmap_ver'),plotFun_2nAK_fLiver_LRv2_clusterLevel_hm,rawData=meanLR_mtx,width = 11,height = 5.5,res = 500,useDingbats = F)  
  }else{
    saveFig(file.path(plotDir,'Figxx_2nAKLiv_LRv2_heatmap_hor'),plotFun_2nAK_fLiver_LRv2_clusterLevel_hm,rawData=meanLR_mtx,width = 6,height = 11,res = 500,useDingbats = F)
  }
  
}





fig1f_2nAK_Liv_fracDEG = function(){
  
  
  if(file.exists(file.path(plotDir,'Fig1f_2nAKLiver_fractionExpressedGenesDE_v3_rawData.tsv'))){
    data = read.delim(file.path(plotDir,'Fig1f_2nAKLiver_fractionExpressedGenesDE_v3_rawData.tsv'),sep = '\t')
  }else{
    # Work generated in 3e_pseudoBulk_DESeq2_withCyclCells.R (or 3.2_DEanalysis_plots.R)
    data = read.csv('~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/liver_nDEGs_perCTperGeno_summary.csv')
  }
  
  
  plotFun_degFrac_perCT_perGeno = function(noFrame=FALSE,noPlot=FALSE){
    
    ## Define lineage for facet
    data$group = ifelse(data$ct %in% c('HSC_MPP'),'HSC_MPP',
                        ifelse(data$ct %in% c('MEMP_MEP','MK','Mast.cell','EE','ME','LE'),'Meg/Ery/Mast',
                               ifelse(data$ct %in% c('CMP_GMP','Monocyte','Macrophage','Kupffer.cell','myelocyte','DC2','pDC'),'Myeloid',
                                      ifelse(data$ct %in% c('B.cell.prog','B.cell','NK.T'),'Lymphoid',
                                             ifelse(data$ct %in% c('Hepatocyte','Fibroblast','Endo'),'Stromal','others')))))
    data$group = factor(data$group,c('HSC_MPP','Meg/Ery/Mast','Myeloid','Lymphoid','Stromal'))
    ## Define cell type order
    data$ct = factor(as.character(data$ct),c('HSC_MPP',
                                             'MEMP_MEP','MK','Mast.cell','EE','ME','LE',
                                             'CMP_GMP','Monocyte','Macrophage','Kupffer.cell','myelocyte','DC2','pDC',
                                             'Endo','Fibroblast','Hepatocyte',
                                             'B.cell.prog','B.cell',"NK.T"
    ))
    
    # For the stick
    data$ct_numericalID = NA
    for(g in unique(data$group)){
      d = data[data$group == g,]
      d$ct = factor(d$ct,levels = levels(d$ct)[levels(d$ct) %in% d$ct])
      d$ct_numericalID = as.numeric(d$ct)
      data$ct_numericalID[data$group == g] = d$ct_numericalID[match(data$ct[data$group == g],d$ct)]
    }
    
    data2 = data %>% group_by(group,ct,ct_numericalID) %>% summarise(max_frac = max(frac_total_nDEG))
    
    p3 = ggplot(data,aes(ct,log10(frac_total_nDEG)))+
      facet_grid(.~group,scales = 'free_x',space = 'free_x')+
      geom_point(aes(col=geno,size=log10(nCell)),alpha=0.5)+
      geom_segment(data=data2,aes(x = ct_numericalID,xend=ct_numericalID,y=min(log10(data$frac_total_nDEG)),yend=log10(max_frac)),col=grey(0.7))+
      scale_y_continuous(labels = c(10^(-2),10^(-3),10^(-4)),breaks = c(-2,-3,-4))+
      scale_color_manual(values = geno_cols)+
      theme_classic(base_size = 15)+
      theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
      theme(panel.border = element_rect(fill=F,linewidth = 1.5),axis.line = element_blank())+
      xlab('Cell type') + ylab('Fraction of DEGs') 
    print(p3)
  }
  
  saveFig(file.path(plotDir,'Fig1f_2nAKLiver_fractionExpressedGenesDE_v3_log10yScale'),plotFun_degFrac_perCT_perGeno,rawData=data,width = 11,height = 5.5,res = 500,useDingbats = F)
}












##----------------------------------------##
##    Figure 3 - TAM - MLDS dataset     ####
##----------------------------------------##
mlds_srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS'
mlds_mdat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv'

mlds_linCol = c('B lineage' = '#7bb9ba',
                'Dendritic cells'=grey(0.7),
                'Erythroblasts'='#FDD7B7',
                'HSC & prog.' = '#946e04',
                'Megakaryocytes' = '#7f0996',
                'Monocyte/Macrophage' = '#5492BA',#'#8ca8ba',
                'Neutrophil'='#DFCDE4',
                'T/NK lineage' = '#B6C484',
                'Tumour'=brewer.pal(6,'Greys')[c(5,4,3)],
                # TAM  - darkBlue
                'TAM'='#294478')
# Aug2024
mlds_linCol = c('B lineage' = '#b28ba1',
                'Dendritic cells'='#78a9c9',
                'Erythroblasts'='#ecd0cf',
                'HSC & prog.' = '#EDBD74',
                'Megakaryocytes' = 'yellow4',
                'Monocyte/Macrophage' = '#cad5a6',#'#8ca8ba',
                'Neutrophil'='#DFCDE4',
                'T/NK lineage' = grey(0.8),
                'Tumour'=brewer.pal(6,'Greys')[c(5,4,3)],
                'TAM' = '#b5d5ef','TAM_Relapse' = '#0A8546',
                'Tumour'='#463A2F','Tumour_Refractory'='#B20000','Tumour_postChemo'='#8f4aa8',
                'Tumour_Relapse'='#fcef11','Tumour_Relapse2'='#CC4D00')

#
ccs = c('TAM_D' = '#b5d5ef',#'TAM_D' = 'dodgerblue2',
        'TAM_TP1' = '#0A8546',#'#005900', #'TAM_TP1' = '#7A4B82',
        'Tumour'='#463A2F','Tumour_Refractory'='#B20000',
        'Tumour_postChemo'='#8f4aa8','Tumour_Relapse'='#fcef11',#'Tumour_Relapse'='#E1BA00',
        'Tumour_Relapse2'='#CC4D00','Normal' = grey(0.8))

col_jan24 = c(
  #B lineage
  '#DABE99',#'#9a7c54','#635547',
  #Dendritic cells
  '#ED9300',#'#ECBE00','#E3DE00', 
  #Erythroblasts
  brewer.pal(6,'OrRd')[c(2)],
  #HSC & prog.
  brewer.pal(7,'OrRd')[c(5)],
  
  #Megakaryocytes
  '#F397C0',#'#EF5A9D','#F6BFCB',
  
  #Monocyte/Macrophage
  #'#005579',
  #'#8DB5CE',#'#C9EBFB',
  '#8ca8ba',
  
  #Neutrophil
  '#DFCDE4',#'#532C8A','#8870ad',
  
  #T/NK lineage
  #'#647a4f','#8EC792',
  '#CDE088',
  
  # TAM  - darkBlue
  '#294478',
  #Tumour
  brewer.pal(6,'Greys')[c(5,4,3)],
  #Tumour?
  brewer.pal(7,'Greens')[c(7)],
  '#E48600',brewer.pal(6,'PuRd')[c(2)],
  '#DD5F21', '#E2691A', '#EEA800'
)



fig3a_MLDS_UMAP = function(){
  
  if(file.exists(file.path(plotDir,'Fig3a_MLDS_TumNorm_UMAP_rawData.tsv')) & skipIfExists){
    dd = read.delim(file.path(plotDir,'Fig3a_MLDS_TumNorm_UMAP_rawData.tsv'),sep = '\t',header = T)
    
  }else{
    # Import the seurat object
    # mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_clean_annotated_noUnknowns_jan24.RDS')
    # broadAnno = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24_mdat.csv')
    # mlds$broadLineage = broadAnno$broadLineage[match(mlds$cellID,broadAnno$cellID)]
    
    mdat = read.csv(mlds_mdat_fp)
    mdat$annot = mdat$annot_aug24_new
    mdat$finalAnn_broad = mdat$annot_aug24_new
    mdat$broadLineage[mdat$finalAnn_broad == 'Tumour' & mdat$broadLineage != 'Tumour'] = 'Tumour'
    mdat$broadLineage[mdat$finalAnn_broad %in% c('Tum_MK?','Tumour_WT')] = 'Tumour_unsure'
    
    if(!keep_CC3){
      # Remove CC3 and L041 Diagnostic samples
      mdat = mdat[!(mdat$donorID %in% c('CC3','L041') & mdat$timePoint == 'Diagnostic'),]
    }
    
    dd = mdat
    #dd = cbind(mdat[,c("cellID","donorID","broadLineage",'timePoint','disease')],mlds@reductions$umap@cell.embeddings)
    
    dd$broadLineage_TP = paste0(dd$broadLineage,':',dd$timePoint)
    dd = dd[!dd$broadLineage %in% c('doublets','unsure_others','others','Tumour?','lowQual'),]
    dd$broadLineage[dd$timePoint == 'TP1' & dd$broadLineage == 'Tumour' & dd$donorID == 'L038'] = 'Tumour_Refractory'
    dd$broadLineage[dd$timePoint %in% c('TP1','TP2','TP4') & dd$broadLineage == 'Tumour' & dd$disease == 'MLDS' & dd$donorID != 'L038'] = 'Tumour_postChemo'
    dd$broadLineage[dd$timePoint == 'D.Relapse' & dd$broadLineage == 'Tumour' & dd$disease == 'MLDS'] = 'Tumour_Relapse'
    dd$broadLineage[dd$timePoint == 'D.Relapse' & dd$broadLineage == 'Tumour' & dd$disease == 'TAM'] = 'TAM_Relapse'
    dd$broadLineage[dd$timePoint == 'D.Relapse2' & dd$broadLineage == 'Tumour'] = 'Tumour_Relapse2'
    dd$broadLineage[dd$timePoint == 'Diagnostic' & dd$broadLineage == 'Tumour' & dd$disease == 'TAM'] = 'TAM'
    
    
  }
  # Prepare source data
  data = mdat[!mdat$broadLineage %in% c('doublets','unsure_others','others','Tumour?','lowQual'),c('cellID','UMAP_1','UMAP_2','annot','orig.ident','donorID','timePoint','broadLineage','clinicalOutcome')]
  data$broadLineage = dd$broadLineage[match(data$cellID,dd$cellID)]
  
  ## Downsample the object for plotting
  dd = dd %>% group_by(broadLineage) %>% mutate(nCell = n())
  dd1 = dd[dd$nCell <= 10000,]
  dd2 = dd[dd$nCell > 10000,]
  set.seed(1234)
  dd2 = dd2 %>% group_by(broadLineage) %>% mutate(id = 1:n(),
                                                  selected = ifelse(id %in% sample(1:n(),10000),T,F))
  dd = rbind(dd1,dd2[dd2$selected==T,!colnames(dd2) %in% c('id','selected')])
  
  plotFun_byTumNorm_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    dd$tumNormCat = ifelse(dd$broadLineage == 'TAM' & dd$timePoint == 'Diagnostic','TAM_D',
                           ifelse(dd$broadLineage == 'TAM_Relapse' & dd$timePoint != 'Diagnostic','TAM_TP1',
                                  ifelse(dd$broadLineage %in% c('Tumour','Tumour_Refractory','Tumour_postChemo','Tumour_Relapse','Tumour_Relapse2'),dd$broadLineage,'Normal')))
    ccs = c('TAM_D' = '#b5d5ef',#'TAM_D' = 'dodgerblue2',
            'TAM_TP1' = '#0A8546',#'#005900', #'TAM_TP1' = '#7A4B82',
            'Tumour'='#463A2F','Tumour_Refractory'='#B20000',
            'Tumour_postChemo'='#8f4aa8','Tumour_Relapse'='#fcef11',#'Tumour_Relapse'='#E1BA00',
            'Tumour_Relapse2'='#CC4D00','Normal' = grey(0.8))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','TAM - MLDS'),
         frame.plot=F)
    
    if(!noPlot){
      # points(dd$UMAP_1,dd$UMAP_2,
      #        col = colAlpha(ccs[dd$tumNormCat],alphas = 0.75),
      #        pch = 19,
      #        cex=0.03)
      points(dd$UMAP_1[!dd$tumNormCat %in% c('Tumour_postChemo','Tumour_Relapse','Tumour_Refractory','TAM_TP1')],dd$UMAP_2[!dd$tumNormCat %in% c('Tumour_postChemo','Tumour_Relapse','Tumour_Refractory','TAM_TP1')],
             col = colAlpha(ccs[dd$tumNormCat[!dd$tumNormCat %in% c('Tumour_postChemo','Tumour_Relapse','Tumour_Refractory','TAM_TP1')]],alphas = 0.75),
             pch = 19,
             cex=0.03)
      
      points(dd$UMAP_1[dd$tumNormCat %in% c('Tumour_postChemo','Tumour_Relapse','Tumour_Refractory','TAM_TP1')],dd$UMAP_2[dd$tumNormCat %in% c('Tumour_postChemo','Tumour_Relapse','Tumour_Refractory','TAM_TP1')],
             col = colAlpha(ccs[dd$tumNormCat[dd$tumNormCat %in% c('Tumour_postChemo','Tumour_Relapse','Tumour_Refractory','TAM_TP1')]],alphas = 0.7),
             pch = 19,
             cex=0.03)
      
    }
    
    if(noPlot & !noFrame){
      #Add coloured labels
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ broadLineage,data=dd,FUN=mean)
      mids$label = as.character(factor(mids$broadLineage,levels = c('HSC & prog.','Megakaryocytes','Erythroblasts',
                                                                    'Monocyte/Macrophage','Dendritic cells','Neutrophil',
                                                                    'B lineage','T/NK lineage','TAM','TAM_Relapse','Tumour','Tumour_postChemo','Tumour_Refractory','Tumour_Relapse','Tumour_Relapse2','Tumour_unsure')))
      
      plotrix::boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$label,cex = 0.9,xpad = 1.9,ypad = 1.8,border = F,
                   col='black')
    }
  }
  
  saveFig(file.path(plotDir,'Fig2_MLDS_TumNorm_UMAP'),plotFun_byTumNorm_UMAP,rawData=data,width = 3.3,height = 3,res = 500)
  
  
  
  plotFun_byDiseae_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    dd$disease_category = ifelse(dd$broadLineage %in% c('TAM','TAM_Relapse'),'TAM',
                                 ifelse(dd$broadLineage %in% c('Tumour','Tumour_postChemo','Tumour_Refractory','Tumour_Relapse','Tumour_Relapse2','Tumour_unsure'),'ML-DS','Normal'))
    
    ccs = c('TAM' = '#b5d5ef',#'TAM_D' = 'dodgerblue2',
            'ML-DS'='#463A2F',
            'Normal' = grey(0.8))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','TAM - MLDS'),
         frame.plot=F)
    
    if(!noPlot){
      # points(dd$UMAP_1,dd$UMAP_2,
      #        col = colAlpha(ccs[dd$tumNormCat],alphas = 0.75),
      #        pch = 19,
      #        cex=0.03)
      points(dd$UMAP_1,dd$UMAP_2,
             col = colAlpha(ccs[dd$disease_category],alphas = 0.75),
             pch = 19,
             cex=0.03)
    }

  }
  
  saveFig(file.path(plotDir,'SuppFig2_MLDS_byDisease_UMAP'),plotFun_byDiseae_UMAP,rawData=data,width = 3.3,height = 3,res = 500)
  
  
  
  
  
  plotFun_byLin_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    #dd$broadLineage[dd$tumNormCat != 'Normal'] = dd$tumNormCat[dd$tumNormCat != 'Normal']
    ccs2 = c(mlds_linCol)
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         #xlim=c(-13,17),
         #ylim=c(-13,17), 
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','MLDS Bone Marrow'),
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = colAlpha(ccs2[dd$broadLineage],alphas = 0.75),
             pch = 19,
             cex=0.03)
      
      # #Add coloured labels
      # mids = aggregate(cbind(UMAP_1,UMAP_2) ~ broadLineage,data=dd,FUN=mean)
      # mids$label = as.numeric(factor(mids$broadLineage,levels = c('HSC & prog.','Megakaryocytes','Erythroblasts',
      #                                                             'Monocyte/Macrophage','Dendritic cells','Neutrophil',
      #                                                             'B lineage','T/NK lineage','Tumour','Tumour?')))
      # mids$legend = factor(paste0(mids$label,' - ',mids$broadLineage),levels = c('1 - HSC & prog.',
      #                                                                            '2 - Erythroblasts',
      #                                                                            '3 - Megakaryocytes',
      #                                                                            '4 - Monocyte/Macrophage',
      #                                                                            '5 - Dendritic cells',
      #                                                                            '6 - Neutrophil',
      #                                                                            '7 - B lineage',
      #                                                                            '8 - T/NK lineage',
      #                                                                            '9 - Tumour',
      #                                                                            '10 - Tumour?'))
      # 
      # #Position tweaks
      # # mids[mids$finalAnn=='Leukocyte','UMAP_2'] = mids[mids$finalAnn=='Leukocyte','UMAP_2'] + 4.9
      # # mids[mids$finalAnn=='Leukocyte','UMAP_1'] = mids[mids$finalAnn=='Leukocyte','UMAP_1'] - 0.3
      # # mids[mids$finalAnn=='PTC','UMAP_2'] = mids[mids$finalAnn=='PTC','UMAP_2'] - 4.9
      # # mids[mids$finalAnn=='PTC','UMAP_1'] = mids[mids$finalAnn=='PTC','UMAP_1'] -0.8
      # # mids[mids$finalAnn=='Tumour','UMAP_2'] = mids[mids$finalAnn=='Tumour','UMAP_2'] + 3
      # # mids[mids$finalAnn=='Tumour','UMAP_1'] = mids[mids$finalAnn=='Tumour','UMAP_1'] + 0.9
      # # 
      # ccs2 = col_lab[1:n_distinct(dd$broadLineage_TP)]
      # 
      # 
      # names(ccs2) = mids$broadLineage
      # # boxed.labels(mids$UMAP_1,mids$UMAP_2,
      # #              labels=mids$label,cex = 0.9,xpad = 1.9,ypad = 1.8,border = T,
      # #              bg=ccs2[mids$broadLineage],
      # #              #col=ccs[mids$broadLineage])
      # #              col='black')
      # # for(i in 1:length(mids$UMAP_1)){
      # #   draw.circle(mids$UMAP_1[i],mids$UMAP_2[i],
      # #               #labels=mids$broadLineage,
      # #               radius=0.4,col=ccs[mids$broadLineage][i])  
      # # }
      # # 
      
    }
  }
  
  saveFig(file.path(plotDir,'Supp.Fig3_MLDS_Lin_UMAP'),plotFun_byLin_UMAP,rawData=data,width = 3.3,height = 3,res = 500)
  
  
  plotFun_byPDID_UMAP = function(noFrame=FALSE,noPlot=FALSE){
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
    # ccs_donorID_2 = sample(col25[!col25 %in% ccs_donorID],c(n_distinct(dd$donorID) - length(ccs_donorID)))
    # names(ccs_donorID_2)= unique(dd$donorID[!dd$donorID %in% names(ccs_donorID)])
    # ccs_donorID = c(ccs_donorID,ccs_donorID_2)
    table(duplicated(ccs_donorID))
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','TAM - MLDS'),
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = colAlpha(ccs_donorID[dd$donorID],alphas = 1),
             pch = 19,
             cex=0.03)
    }
    
    if(noPlot & !noFrame){
      #Add coloured labels
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ donorID,data=dd[dd$broadLineage == 'Tumour',],FUN=mean)
      mids$label = mids$donorID
      
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$label,cex = 0.9,xpad = 1.9,ypad = 1.8,border = T,
                   col='black',bg = ccs_donorID[mids$label])
    }
  }
  
  saveFig(file.path(plotDir,'SuppFigXX_MLDS_donorID_UMAP'),plotFun_byPDID_UMAP,rawData=data,width = 3.3,height = 3,res = 500)
  
}


DimPlot(mlds,group.by = 'donorID',cols = col25,label = T,repel = T,label.box = T,label.size = 2) +
  NoAxes() + theme(panel.border = element_rect(fill=F,color = 'black',linewidth = 1))





## figure 2C - UMAP coloured by GATA1s 

fig2a_MLDS_GATA1s_status = function(){
  if(file.exists(file.path(plotDir,'Fig2c_MLDS_GATA1s_UMAP_rawData.tsv')) & skipIfExists){
    df = read.delim(file.path(plotDir,'Fig2c_MLDS_GATA1s_UMAP_rawData.tsv'),sep = '\t',header = T)
  }else{
    
    # Import the seurat object
    # mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_clean_annotated_noUnknowns_jan24.RDS')
    # broadAnno = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_dsBALL_MDS_combined_feb24_mdat.csv')
    # mlds$broadLineage = broadAnno$broadLineage[match(mlds$cellID,broadAnno$cellID)]
    mdat = read.csv(mlds_mdat_fp)
    mdat$broadLineage[mdat$annot_aug24 == 'Tumour'] = 'Tumour'
    mdat$broadLineage[mdat$annot_aug24 %in% c('Tumour_WT','Tum_MK?')] = 'Tumour_unsure'
    
    if(!keep_CC3){
      # Remove CC3 and L041 Diagnostic samples
      mdat = mdat[!(mdat$donorID %in% c('CC3','L041') & mdat$timePoint == 'Diagnostic'),]
    }
    mdat$annot = as.character(mdat$annot_aug24)
    dd = mdat
    
    dd = dd[,c("cellID","donorID",'disease','annot',"broadLineage",'timePoint','tissue','GATA1_status','UMAP_1','UMAP_2')]
    #colnames(dd)[colnames(dd) == 'GATA1s_status2'] = 'GATA1s_status'
    dd$GATA1s_status = dd$GATA1_status
    dd$GATA1s_status[dd$donorID %in% c('CC1','L114') & dd$GATA1s_status == 'uninformative'] = 'uninformative'
    #dd$GATA1s_status[dd$donorID %in% c('CC1','L114') & dd$GATA1s_status == 'uninformative'] = 'TBD'
    
    dd = dd[!dd$broadLineage %in% c('doublets','unsure_others','others','lowQual'),]
    dd$GATA1s_status[dd$GATA1s_status %in% c('WT','GATA1s_WT')] = 'GATA1 wild type'
    dd$GATA1s_status[dd$GATA1s_status %in% c('Mut','GATA1s_mutant')] = 'GATA1s mutation'
    dd$GATA1s_status[dd$GATA1s_status %in% c('unsure','GATA1s_unsure')] = 'Not informative'
    dd$GATA1s_status[dd$GATA1s_status %in% c('noCov')] = 'No coverage at mutation site'
    dd$GATA1s_status[dd$GATA1s_status %in% c('noGATA1expr','no_GATA1_expr')] = 'No GATA1 expression'
    dd$GATA1s_status2 = dd$GATA1s_status
    dd$GATA1s_status2[dd$GATA1s_status2 %in% c('Not informative','No coverage at mutation site','uninformative')] = 'Uninformative'
    
    dd_summary = as.data.frame(table(dd$GATA1s_status2,dd$broadLineage == 'Tumour',dd$donorID))
    dd_summary = dd_summary[dd_summary$Freq > 0,]
    colnames(dd_summary) = c('GATA1_status','isTumour','donorID','nCell')
    dd_summary = dd_summary[dd_summary$isTumour == T & 
                              !dd_summary$donorID %in% c('CC1','L041'),]
    dd_summary = dd_summary %>% group_by(donorID) %>% mutate(totalTumourCell = sum(nCell),
                                                             perc = 100*nCell / totalTumourCell)
    100*73545/sum(dd_summary$nCell)
    sum(dd_summary$nCell[dd_summary$GATA1_status == 'GATA1s mutation'])
    sum(dd_summary$nCell[dd_summary$GATA1_status == 'Uninformative'])
    range(dd_summary$perc[dd_summary$GATA1_status == 'GATA1s mutation'])
    
    
  }
  
  #005579 - WT
  ccs = c("GATA1 wild type"='#0666b5',"GATA1s mutation"='#A92821',"Uninformative"=colAlpha(grey(0.8),0.3),"No GATA1 expression"=colAlpha(grey(0),0.15),"TBD"=colAlpha('#a17f45',0.4))
  table(is.na(ccs))
  
  version = 'v2'
  plotFun_GATA1s_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         #xlim=c(-13,17),
         #ylim=c(-13,17), 
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','MLDS Bone Marrow'),
         frame.plot=F)
    
    if(!noPlot){
      if(version == 'v1'){
        points(dd$UMAP_1[dd$GATA1s_status2 != 'GATA1s mutation'],dd$UMAP_2[dd$GATA1s_status2 != 'GATA1s mutation'],
               col = ccs[dd$GATA1s_status2[dd$GATA1s_status2 != 'GATA1s mutation']],
               pch = 19,
               cex=0.01)
        points(dd$UMAP_1[dd$GATA1s_status2 == 'GATA1s mutation'],dd$UMAP_2[dd$GATA1s_status2 == 'GATA1s mutation'],
               col = ccs[dd$GATA1s_status2[dd$GATA1s_status2 == 'GATA1s mutation']],
               pch = 19,
               cex=0.01)  
      }else{
        points(dd$UMAP_1,dd$UMAP_2,
               col = ccs[dd$GATA1s_status2],
               pch = 19,
               cex=0.01)  
      }
      
      
    }
    #legend(x=-7, y=-7,legend=unique(dd$GATA1s_status2),fill = ccs[match(unique(dd$GATA1s_status2),names(ccs))],lwd = 0,cex = 0.4,lty = NA,xjust = 0,seg.len=0.01,box.lwd = 0.0,bty = 'n')
    
  }
  
  if(version == 'v1'){
    saveFig(file.path(plotDir,'Fig2_MLDS_GATA1s_UMAP_v1'),plotFun_GATA1s_UMAP,rawData=dd[,!colnames(dd) %in% c('GATA1s_status','GATA1_status')],width = 3.3,height = 3.0,res = 500)  
  }else{
    saveFig(file.path(plotDir,'Fig2_MLDS_GATA1s_UMAP_v2'),plotFun_GATA1s_UMAP,rawData=dd[,!colnames(dd) %in% c('GATA1s_status','GATA1_status')],width = 3.3,height = 3.0,res = 500)
  }
  
  
  
  
  
  
  
  
  
  
  dd$broadAnno = ifelse(dd$broadLineage %in% c('Tumour','Tumour?','Megakaryocytes','Erythroblasts'), dd$broadLineage,
                        ifelse(dd$annot %in% c('MEP'),'MEP','others'))
  dd$donorID[dd$donorID == 'L076' & dd$tissue == 'Blood'] = 'L076_PB'
  df = dd %>% group_by(donorID,disease,broadAnno,tissue,timePoint) %>% mutate(totalCells = n()) %>% 
    group_by(donorID,disease,broadAnno,timePoint,tissue,totalCells,GATA1s_status2) %>% summarise(n = n()) %>% mutate(fracCells = n/totalCells)
  
  table(df$donorID,df$tissue)
  df$broadAnno = factor(df$broadAnno,levels = c('Tumour','Tumour?','Megakaryocytes','MEP','Erythroblasts','others'))
  dd = df[df$broadAnno != 'Tumour?',]
  dd$broadAnno2 = as.character(dd$broadAnno)
  dd$broadAnno2[dd$broadAnno2 == 'Megakaryocytes']  = 'MK'
  dd$broadAnno2[dd$broadAnno2 == 'Erythroblasts']  = 'Ery'
  dd$broadAnno2[dd$broadAnno2 == 'Tumour']  = 'Blast'
  dd$broadAnno2 = factor(dd$broadAnno2,levels = c('Blast','MEP','MK','Ery','others'))
  
  dd$GATA1s_status2 = factor(dd$GATA1s_status2,rev(c('GATA1s mutation','GATA1 wild type','Uninformative','No GATA1 expression','TBD')))
  dd$timePoint[dd$timePoint == 'D.Relapse'] = 'Relapse diagnostic'
  dd$timePoint[dd$timePoint == 'D.Relapse2'] = 'Relapse2 diagnostic'
  dd$timePoint = factor(dd$timePoint,c('Diagnostic','TP1','TP2','TP4','Relapse diagnostic','Relapse2 diagnostic'))
  dd$disease = factor(dd$disease,c('TAM','MLDS'))
  
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    
    p1 = ggplot(dd,aes(broadAnno2,fracCells,fill=GATA1s_status2))+
      geom_col()+
      scale_fill_manual(values = rev(c('#005579','#A92821',#colAlpha('#9E2718',0.4),
                                       grey(0.65),grey(0.8))),name = 'GATA1 status')+
      facet_grid(vars(donorID),vars(timePoint))+
      scale_y_continuous(breaks=seq(0,1,by=0.5))+#,labels = c(0,0.5,1)
      theme_classic(base_size = 17)+
      theme(axis.text.x = element_text(angle=90,size=14,vjust = 0.5,hjust = 1),
            axis.text.y = element_text(size=11),
            legend.position = 'left')+
      ylab('Fraction of cells')+xlab('')
    
    p2 = ggplot(dd,aes(broadAnno2,fracCells,fill=GATA1s_status2))+
      geom_col(width=0.7)+
      scale_fill_manual(values = ccs,name = 'GATA1 status')+
      facet_grid(timePoint ~ disease+donorID)+
      scale_y_continuous(breaks=seq(0,1,by=0.5),labels = c(0,0.5,1))+
      theme_classic(base_size = 17)+
      theme(axis.text.x = element_text(angle=90,size=12,vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text.y = element_text(size=11,colour = 'black'),
            panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(colour='black'),
            axis.ticks = element_line(colour = 'black'),
            legend.position = 'bottom')+
      ylab('Fraction of cells')+xlab('')
    print(p2)
  }
  saveFig(file.path(plotDir,'SuppFigxx_MLDS_GATA1s_status_barplot_hor'),plotFun,rawData=df,width = 22,height = 10,res = 500)
  
  #saveFig(file.path(plotDir,'Figxx_MLDS_GATA1s_status_barplot'),plotFun,rawData=df,width = 7.7,height = 7.65,res = 500)
  
}





normalCell_markers = unique(c(
  'SPINK2', # HSC_MPP
  'CTSG',	'PRTN3', # CMP
  'AZU1','MPO',
  'CD14',#'FCGR3A',# Monocytes
  #'CD68',#'MSR1',#Macrophages
  #Kupffer
  #'FABP3','CD5L','VCAM1',
  
  
  'IRF8',	 #DC.precursor
  'CLEC9A','CLEC10A','CD1C',#DC1
  'CLEC4C','IL3RA', #pDC


  'CSF3R',
  'MNDA', # NEUTROPHILS
  'DEFA3',#'DEFA4', # pro-myelocytes
  'CAMP','LCN2', #myelocytes
  #'CXCR2',
  'FCGR3B',#'FUT4', # Neutrophil
  #'ANPEP','FCGR3A',
  #'HLA-DRA','HLA-DRB1',
  #'IGLL1','CD99', # Pre-pro
  'DNTT',
  #'EBF1',
  'CD19',#'RAG1',# pro-b-cells
  'MME','VPREB1','CD79A','CD79B',# pre-b-cells
  'TCL1A','MME',
  'MS4A1',  # B-cells
  #'CD27',#plasma cell?
  
  'CD52','IL7R',# ILC precursor
  #'CD3D',#'GZMA',
  'NKG7','KLRD1' #NK
))

# Jack's Flow markers
# positive markers from Jack's flow
positive_markers = unique(c('CD34','CD38',
                            # MEMP - MEP
                            'IL1B','GATA2','GATA1','KLF1','TESPA1',#'CTNNBL1',
                            'HBD','CSF2RB','CTNNBL1',#MEMP
                            'FLI1','FCER1A','ITGA2B', # earlyMK
                            'KIT','HDC','CPA3',
                            'CD4','CD33','CD7','ZBTB16',#'LTB'
                            
                            
                            'PLEK','PF4', #Megakaryocyte
                            'GP1BA','ITGB3', #Megakaryocyte
                            'PPBP','TUBB1', # Platelets
                            'TPSAB1', # Mast.cell
                            'EPCAM','APOC1', # early Erythroid
                            'ALAS2', # mid.erythroid
                            'HBB','BPGM' # late.erythroid
                            
                            
)) 

negative_markers = c('CD19','CD2','MME','MS4A1', # Negative
                     'CD5','CD3D','CD3E','CD8A'# Negative
)

# Laura's Markers - both positive and negative markers
lauras_markers = c('CD34','NCAM1','KIT','ANPEP','CD33','CD7','CD4','ITGA2B','GP1BA','MPL','IL3RA','CD36','ITGB3','TFRC',
                   'MPO','ICAM3','CD40','GYPA')

# Curated markers ##
FeaturePlot(tgt.srat.Tumour,features = c('GATA1','GATA2','IKZF1','MYB','MYC','CD34','KIT','ANPEP','CD33','CD7','NCAM1','MPO'))



figSup3c_MLDS_dotPlot = function(){
  ## Import MLDS object
  mlds = readRDS(mlds_srat_fp)
  mdat = read.csv(mlds_mdat_fp)
  mdat$broadLineage[mdat$annot_aug24 == 'Tumour'] = 'Tumour'
  mlds$broadLineage = mdat$broadLineage[match(mlds$cellID,mdat$cellID)]
  mlds$broadLineage[mlds$broadLineage != 'Tumour' & mlds$annot_aug24 == 'Tumour' ] = 'Tumour'
  mlds$broadLineage[mlds$annot_aug24 != 'Tumour' & mlds$broadLineage == 'Tumour'] = 'Tumour_unsure' # this category includes tumour_WT, Tum_MK?, tumour_unsure_postTreatment
  
  #mlds$broadLineage[mlds$broadLineage =='Tumour' & mlds$donorID == 'L041'] = 'others'
  mlds$timePoint = ifelse(grepl('Diagnostic',mlds$timePoint),gsub('Diagnostic','D',mlds$timePoint),
                          ifelse(grepl('Relapse',mlds$timePoint),gsub('Relapse','R',mlds$timePoint),mlds$timePoint))
  mlds$timePoint[mlds$timePoint == 'D.R'] = 'RD'
  mlds$timePoint[mlds$timePoint == 'D.R2'] = 'R2D'
  
  mlds$broadLineage[mlds$broadLineage == 'Tumour'] = paste0(mlds$disease[mlds$broadLineage == 'Tumour'],':',mlds$donorID[mlds$broadLineage == 'Tumour'],':',mlds$timePoint[mlds$broadLineage == 'Tumour'],':',mlds$tissue[mlds$broadLineage == 'Tumour'])
  mlds$broadLineage[grepl('MLDS|TAM',mlds$broadLineage) & !mlds$timePoint %in% c('Diagnostic','D.Relapse','D.Relapse2','D','RD','R2D') & mlds$donorID != 'L038'] = 'Blasts post-treatment'
  mlds$broadLineage[mlds$broadLineage == 'Tumour_unsure'] = 'Blast_unsure'
  nCells_perGroup = table(mlds$broadLineage)
  mlds$nCell = nCells_perGroup[match(mlds$broadLineage,names(nCells_perGroup))]
  mlds$broadLineage = paste0(mlds$broadLineage,' (n=',prettyNum(mlds$nCell,big.mark = ','),')')
  mlds$broadLineage2 = as.character(mlds$broadLineage)
  mlds$broadLineage2 = gsub('TAM:|MLDS:|:D|:BM|:Blood','',mlds$broadLineage2)
  
  # Import annot 
  mlds$broadLineage2 = factor(mlds$broadLineage2, levels = c(unique(mlds$broadLineage2[grepl('TAM',mlds$broadLineage)])[order(unique(mlds$broadLineage2[grepl('TAM',mlds$broadLineage)]))],
                                                           unique(mlds$broadLineage2[grepl('^MLDS',mlds$broadLineage)])[order(unique(mlds$broadLineage2[grepl('MLDS',mlds$broadLineage)]))],
                                                           unique(mlds$broadLineage2[grepl('Blast',mlds$broadLineage)])[order(unique(mlds$broadLineage2[grepl('Blast',mlds$broadLineage)]))],
                                                           #unique(mlds$broadLineage2[grepl('Tumour_unsure',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('HSC & prog',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('Megakaryocytes',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('Mast',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('Erythroblasts',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('Monocyte',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('Dendritic',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('Myelocytes',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('Neutrophil',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('B ',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('T/NK',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('Stromal',mlds$broadLineage)]),
                                                           unique(mlds$broadLineage2[grepl('others',mlds$broadLineage)]),
                                                           
                                                           'HSC & prog.','Megakaryocytes','Mast.cell','Erythroblasts',
                                                           'Monocyte/Macrophage','Dendritic cells','Myelocytes','Neutrophil',
                                                           'B lineage','T/NK lineage','Stromal','others'))
  
  table(is.na(mlds$broadLineage2))
  
  #mlds$broadLineage2 = factor(as.numeric(mlds$broadLineage),c(1:length(levels(mlds$broadLineage))))
  
  direction = 'horizontal'
  if(direction == 'vertical'){
    genes = rev(unique(c(positive_markers,normalCell_markers,negative_markers)))
  }else{
    genes = unique(c(positive_markers,normalCell_markers,negative_markers))
  }
  
  if(keep_CC3){
    idents_toPlot = unique(mlds$broadLineage2[mlds$broadLineage != 'others'])  
  }else{
    idents_toPlot = unique(mlds$broadLineage2[!grepl('others|Blast_unsure|CC3',mlds$broadLineage)])
  }
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    Idents(mlds) = mlds$broadLineage2
    #Idents(mlds) = mlds$annot_aug24
    p = DotPlot(mlds,idents = idents_toPlot,
                #idents = c('Tumour','MEP','HSC_MPP','MK','EE','EE_1','ME','LE'),
                #cols = c("#EBFFE5", "#244b05"),
                #cols = c(colAlpha('#F1F5FA',1),'#425580'),
                cols = c(colAlpha(grey(0.95),0.8),'black'),
                #cols = c(grey(0.99), grey(0.2)),
                #group.by = 'seurat_clusters',
                #idents = unique(srat$finalAnn[srat$finalAnn != 'others']),
                features = genes)+RotatedAxis() 
    
    if(direction == 'vertical'){
      p = p + coord_flip() + 
        scale_y_discrete(position = "right")+
        theme(axis.text.x = element_text(size=11,angle = 0,vjust = 0,hjust = 0.5),
              axis.text.y = element_text(size=10),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top') + xlab('') + ylab('')
    }else if(direction == 'horizontal'){
      p = p +
        theme(axis.text.y = element_text(size=11),
              axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top') + xlab('') + ylab('')
    }
    print(p)
  }
  
  if(direction == 'vertical'){
    saveFig(file.path(plotDir,'FigSuppXX_MLDS_Lin_DotPlot_vertical'),plotFun,width = 5,height = 11.7,res = 500)  
  }else{
    saveFig(file.path(plotDir,'FigSuppXX_MLDS_Lin_DotPlot_horizontal'),plotFun,width = 14,height = 8,res = 500)  
  }
  
}




supFig3D_MLDS_LRv1_heatmap = function(){
  ref_dataset = '2n_FL'
  
  if(file.exists(file.path(plotDir,'Fig1b_2nAKLiv_geno_UMAP_rawData.tsv'))){
    #dd = read.delim(file.path(plotDir,'Fig1b_2nAKLiv_geno_UMAP_rawData.tsv'),header = T,sep = '\t')
    print('whoops')
  }else{
    mdat = read.csv(mlds_mdat_fp)
    mdat$broadLineage[mdat$broadLineage != 'Tumour' & mdat$annot_aug24 == 'Tumour'] = 'Tumour'
    mdat$broadLineage[mdat$broadLineage %in% c('Tumour_unsure','lowQual')] = 'others'
    mdat$broadLineage[mdat$broadLineage =='Tumour' & mdat$donorID == 'L041'] = 'others'
    
    mdat$group = ifelse(mdat$broadLineage == 'Tumour',paste0(mdat$disease,':',mdat$donorID,':',mdat$timePoint,':',mdat$tissue),
                        ifelse(mdat$broadLineage == 'others','others',mdat$annot_aug24))
    mdat$group[mdat$group %in% c('T_gd','T_reg','T_CD4','T_CD8','T_cells')] = 'T cell'
    mdat$group[mdat$group %in% c('activated_neutrophil','Neutrophil')] = 'Neutrophil'
    mdat$group[mdat$group %in% c('Mono_CD14','Mono_CD16')] = 'Monocyte'
    mdat$group[mdat$group %in% c('naive.B')] = 'B cell'
    mdat$group[grepl('\\.',mdat$group) & !grepl('Leuk',mdat$group)] = gsub('\\.',' ',mdat$group[grepl('\\.',mdat$group) & !grepl('Leuk',mdat$group)])
    mdat$group = gsub('D Relapse','RD',mdat$group)
    
    
    
    ## Read in LRv1 similarity matrix
    #lr_output = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/LRv1_published_FLref/MLDS_2nFLref_maxCells_70perc_raw_LR_outputs.RDS')
    lr_output = readRDS('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/2_LRv1_FLref/MLDS_2nFLref_maxCells_70perc_raw_LR_outputs.RDS')
    lr_output = lr_output[[1]]
    
    
    if(ref_dataset == 'publishedFL'){
      ref_order = c("HSC_MPP",
                    "MEMP","Megakaryocyte","Mast cell",
                    "Early Erythroid", "Mid Erythroid","Late Erythroid",
                    "Neutrophil-myeloid progenitor",
                    "Monocyte precursor","Monocyte","Mono-Mac","Kupffer Cell",'VCAM1..EI.macrophage',
                    "DC precursor","DC1","DC2","pDC precursor",
                    "Pre pro B cell","pro-B cell","pre-B cell","B cell",
                    "ILC precursor", "Early lymphoid_T lymphocyte", "NK",
                    "Hepatocyte","Endothelial cell","Fibroblast",'SCPs' 
      )
    }else if(ref_dataset == '2n_FL'){
      ref_order = c("HSC_MPP","MEMP_MEP",'CMP_GMP','LMPP_ELP',
                    'earlyMK',"MK","Mast.cell",
                    "EE", "ME","LE",
                    'promyelocyte','myelocyte','proMono','Monocyte','Macrophage','Kupffer.cell',
                    "DC precursor","DC1","DC2","pDC",
                    "pro.B.cell","pre.B.cell","B.cell",
                    "ILC.precursor", "T.cell", "NK_T",
                    "Hepatocyte","Endo",'Mesenchyme','Fibroblast',
                    'Cholangiocytes','Mesothelial_cells','Neuron','SCPs' 
      )
    }
    
    
    colnames(lr_output)[!colnames(lr_output) %in% ref_order]
    ref_order = ref_order[ref_order %in% colnames(lr_output)]
    table(ref_order %in% colnames(lr_output))
    
    rownames(lr_output) = gsub('D Relapse','RD',rownames(lr_output))
    tgt_order = c(rownames(lr_output)[grepl('ref_',rownames(lr_output))],
                  rownames(lr_output)[grepl('TAM',rownames(lr_output))][order(rownames(lr_output)[grepl('TAM',rownames(lr_output))])],
                  rownames(lr_output)[grepl('MLDS',rownames(lr_output))][order(rownames(lr_output)[grepl('MLDS',rownames(lr_output))])],
                  "HSC_MPP","MEP",'CMP_GMP',
                  "MK",
                  "EE", "ME","LE",
                  'Myelocytes','Neutrophil','Monocyte','Macrophage',
                  "DC1","DC2","pDC",
                  "pro B cell","pre B cell","B cell",'Plasma cell',
                  'T cell',"NK_T",'NK',
                  'others'
    )
    
    
    table(rownames(lr_output) %in% tgt_order)
    rownames(lr_output)[!rownames(lr_output) %in% tgt_order]
    tgt_order = tgt_order[tgt_order %in% rownames(lr_output)]
    
    mdat$broadLineage[mdat$broadLineage == 'Tumour'] = mdat$disease[mdat$broadLineage == 'Tumour']
    ## Assign "unsure - post treatment" blasts as "unsure"
    mdat$broadLineage[grepl('TP',mdat$group) & mdat$donorID != 'L038'] = 'unsure_MLDS'
    
    mdat$broadLineage = factor(mdat$broadLineage,levels = unique(mdat$broadLineage[match(tgt_order[!grepl('ref',tgt_order)],mdat$group)]))
    
    
    table(tgt_order[!tgt_order %in% rownames(lr_output)])
    table(rownames(lr_output)[!rownames(lr_output) %in% tgt_order])
    
    
    lr_output = lr_output[tgt_order,rev(ref_order)]
    lr_output = lr_output[!grepl('ref_|others',rownames(lr_output)),] %>% t()
  }
  
  
  
  plotFun_MLDS_LRv1_clusterLevel_hm = function(noFrame=FALSE,noPlot=FALSE){
    library(ComplexHeatmap)
    par(mar=c(0.1,0.1,1,0.1))
    logitCols = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
    cols = circlize::colorRamp2(seq(-3,3,length.out=length(logitCols)),logitCols)
    
    mtx = lr_output
    query_cell_names = data.frame(og = colnames(mtx),
                                  new = colnames(mtx))
    query_cell_names$new = gsub('TAM:|MLDS:','',query_cell_names$new)
    query_cell_names$new = gsub('Diagnostic','D',query_cell_names$new)
    query_cell_names$new = gsub('Relapse','R',query_cell_names$new)
    query_cell_names$new = gsub('D R','RD',query_cell_names$new)
    query_cell_names$new = gsub('Blood','PB',query_cell_names$new)
    query_cell_names$new = gsub('_',' / ',query_cell_names$new)
    query_cell_names$new = gsub('^pr','Pr',query_cell_names$new)
    query_cell_names$new = gsub('EE','early Ery',query_cell_names$new)
    query_cell_names$new = gsub('ME$','mid Ery',query_cell_names$new)
    query_cell_names$new = gsub('LE','late Ery',query_cell_names$new)
    
    mdat$group = query_cell_names$new[match(mdat$group,query_cell_names$og)]
    colnames(mtx) = query_cell_names$new[match(colnames(mtx),query_cell_names$og)]
    #o = order(colnames(mtx))
    #mtx = mtx[,o]
    #lr_output = lr_output[,o]
    
    rownames(mtx)[rownames(mtx) == 'Mesothelial_cells'] = 'Mesothelial cell'
    rownames(mtx)[rownames(mtx) == 'SCPs'] = 'SCP'
    rownames(mtx)[rownames(mtx) == 'Cholangiocytes'] = 'Cholangiocyte'
    rownames(mtx)[rownames(mtx) == 'Endo'] = 'Endothelium'
    rownames(mtx) = gsub('_',' / ',rownames(mtx))
    rownames(mtx) = gsub('\\.',' ',rownames(mtx))
    rownames(mtx) = gsub('^pr','Pr',rownames(mtx))
    rownames(mtx) = gsub('myelocyte','Myelocyte',rownames(mtx))
    rownames(mtx) = gsub('roM','ro-m',rownames(mtx))
    rownames(mtx) = gsub('earlyMK','early MK',rownames(mtx))
    rownames(mtx) = gsub('EE','early Ery',rownames(mtx))
    rownames(mtx) = gsub('ME$','mid Ery',rownames(mtx))
    rownames(mtx) = gsub('LE','late Ery',rownames(mtx))
    
    #colnames(mtx)[grepl(':',colnames(mtx))] = paste0(sapply(split(colnames(mtx)[grepl(':',colnames(mtx))],':')))
    # mdat$group2 = ifelse(mdat$annot_aug24 == 'Tumour',paste0(mdat$donorID,':',mdat$timePoint,':',mdat$tissue),as.character(mdat$broadLineage))
    # mdat$group2 = gsub('Diagnostic','D',mdat$group2)
    # mdat$group2 = gsub('Relapse','R',mdat$group2)
    # mdat$group2 = gsub('D R','RD',mdat$group2)
    # mdat$group2 = gsub('Blood','PB',mdat$group2)
    
    mdat$clinicalOutcome[mdat$donorID %in% c('L156','L038')] = 'Refractory'
    mdat$clinicalOutcome[mdat$donorID %in% c('L076','CC3')] = 'Relapse'
    mdat$clinicalOutcome[mdat$donorID %in% c('L019','L039','L040','L042','L091','L041','CC4','CC5','CC1','CC2','L178',
                                             'L075','L114','L182')] = 'Remission'
    mdat$clinicalOutcome[mdat$clinicalOutcome == '?'] = 'unknown'
    
    mdat$clinicalOutcome_group = ifelse(grepl(':',mdat$group), mdat$clinicalOutcome, '')
    
    
    outcomeCols = c('Remission'='#2D4372','Refractory'=col25[2],'Relapse'=col25[5],'unknown'=grey(0.6))
    botAnno = HeatmapAnnotation(df = data.frame(clinicalOutcome = mdat$clinicalOutcome_group[match(colnames(mtx),mdat$group)]),
                                annotation_name_side = 'left',
                                col = list(clinicalOutcome = outcomeCols))
    hm = Heatmap(mtx,show_row_dend = F,show_column_dend = F,
                 col = cols,
                 name=paste0('Predicted Similarity','(Logit)'),
                 show_row_names = T,show_column_names = T,cluster_rows = F,
                 cluster_columns = F,cluster_row_slices = F,
                 row_names_gp = gpar(fontsize=10),row_title_gp = gpar(fontsize=14),
                 column_split = mdat$broadLineage[match(colnames(mtx),mdat$group)],
                 column_title_rot = 90,
                 column_names_gp = gpar(fontsize=11),column_title_gp = gpar(fontsize=12),
                 bottom_annotation = botAnno,
                 #heatmap_legend_param = list(direction = "horizontal"),
                 row_title = 'Reference classes')
    #geno = gsub('.*:','',rownames(meanLR_mtx))
    # geno[geno == 'complete_trisomy'] = 'Triploid'
    # geno = factor(geno,c('diploid','T21','T18','T22','MX','Triploid'))
    # rowAnno = rowAnnotation(Genotype = geno,
    #                         col = list(Genotype = geno_cols))
    # hm = Heatmap(meanLR_mtx,show_row_dend = F,show_column_dend = F,
    #              col = cols,
    #              name=paste0('Predicted\nSimilarity\n','(Logit)'),
    #              show_row_names = F,show_column_names = T,cluster_rows = F,cluster_columns = F,
    #              row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=8),
    #              split = mdat$finalAnn_broad[match(rownames(meanLR_mtx),mdat$group)],
    #              row_title_rot = 0,row_title_gp = gpar(fontsize=8),cluster_row_slices = F,
    #              right_annotation = botAnno,
    #              column_title = 'Reference classes')
    
    
    if(!noPlot){
      draw(hm,heatmap_legend_side = "right")
    }
  }
  
  saveFig(file.path(plotDir,'FigSup3_MLDS_LRv1_hm'),plotFun_MLDS_LRv1_clusterLevel_hm,rawData=mtx,width = 13.5,height = 6.3,res = 500,useDingbats = F)
}


## figure supplementary - L067 MDS - UMAP coloured by broad lineage + GATA1s
figSuppXX_MDS_UMAP = function(){
  
  if(file.exists(file.path(plotDir,'Fig3a_MLDS_TumNorm_UMAP_rawData.tsv')) & skipIfExists){
    dd = read.delim(file.path(plotDir,'Fig3a_MLDS_TumNorm_UMAP_rawData.tsv'),sep = '\t',header = T)
    
  }else{
    # Import the seurat object
    #mds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MDS/MDS_clean_annotated_tmp.RDS')
    dd = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MDS/MDS_clean_annotated_2404_mdat.csv')
    dd = dd[!dd$broadLineage %in% c('doublets','unsure_others','others','lowQual'),]
    
  }
  
  
  plotFun_byLin_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    ccs2 = c(mlds_linCol,'Tumour'='#463A2F')
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         #xlim=c(-13,17),
         #ylim=c(-13,17), 
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','MDS Bone Marrow'),
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = colAlpha(ccs2[dd$broadLineage],alphas = 1),
             pch = 19,
             cex=0.3)
    }
    
    if(noPlot & !noFrame){
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ broadLineage,data=dd,FUN=mean)
      
      mids$label = mids$broadLineage
      
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$label,cex = 0.6,xpad = 1.3,ypad = 2.3,border = T,
                   bg=ccs2[mids$label],
                   col='black')
    }
  }
  
  saveFig(file.path(plotDir,'FigSuppXX_MDS_Lin_UMAP'),plotFun_byLin_UMAP,rawData=dd,width = 2,height = 2,res = 500)
  
  
  
  
  
  
  
  
  ##------    GATA1s status  --------##
  dd$GATA1_status[dd$GATA1_status %in% c('WT','GATA1s_WT')] = 'GATA1 wild type'
  dd$GATA1_status[dd$GATA1_status %in% c('Mut','GATA1s_mutant')] = 'GATA1s mutation'
  dd$GATA1_status[dd$GATA1_status %in% c('noGATA1expr','no_GATA1_expr')] = 'No GATA1 expression'
  dd$GATA1_status[dd$GATA1_status %in% c('Not informative','No coverage at mutation site','uninformative')] = 'Uninformative'
  ccs = c("GATA1 wild type"='#005579',"GATA1s mutation"='#A92821',"Uninformative"=grey(0.65),"No GATA1 expression"=grey(0.8),"TBD"=colAlpha('#a17f45',0.4))
  table(is.na(ccs))
  
  plotFun_GATA1s_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         #xlim=c(-13,17),
         #ylim=c(-13,17), 
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','MLDS Bone Marrow'),
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1[dd$GATA1_status != 'GATA1s mutation'],dd$UMAP_2[dd$GATA1_status != 'GATA1s mutation'],
             col = ccs[dd$GATA1_status[dd$GATA1_status != 'GATA1s mutation']],
             pch = 19,
             cex=0.1)
      points(dd$UMAP_1[dd$GATA1_status == 'GATA1s mutation'],dd$UMAP_2[dd$GATA1_status == 'GATA1s mutation'],
             col = ccs[dd$GATA1_status[dd$GATA1_status == 'GATA1s mutation']],
             pch = 19,
             cex=0.1)
      
      
    }
  }
  
  saveFig(file.path(plotDir,'FigSuppXX_MDS_GATA1s_UMAP'),plotFun_GATA1s_UMAP,rawData=dd,width = 2,height = 2,res = 500)
  
}




supFig3D_MLDS_LRv2_heatmap = function(){
  
  if(file.exists(file.path(plotDir,'Fig1b_2nAKLiv_geno_UMAP_rawData.tsv'))){
    dd = read.delim(file.path(plotDir,'Fig1b_2nAKLiv_geno_UMAP_rawData.tsv'),header = T,sep = '\t')
    
  }else{
    mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
    mdat$broadLineage[mdat$broadLineage %in% c('Tumour_unsure','lowQual')] = 'others'
    mdat$broadLineage[mdat$broadLineage =='Tumour' & mdat$donorID == 'L041'] = 'others'
    
    mdat$broadLineage[mdat$broadLineage == 'Tumour'] = paste0('Leuk:',mdat$disease[mdat$broadLineage == 'Tumour'],':',mdat$donorID[mdat$broadLineage == 'Tumour'])
    
    mdat$group = ifelse(mdat$broadLineage == 'Tumour',paste0('Leuk:',mdat$disease,':',mdat$donorID),
                        ifelse(mdat$broadLineage == 'others','others',mdat$annot_mar24))
    mdat$group[mdat$group %in% c('T_gd','T_reg','T_CD4','T_CD8')] = 'T.cell'
    mdat$group[mdat$group %in% c('activated_neutrophil','Neutrophil')] = 'Neutrophil'
    mdat$group[mdat$group %in% c('Mono_CD14','Mono_CD16')] = 'Monocyte'
    mdat$group[mdat$group %in% c('naive.B')] = 'B.cell'
    
    
    ## Read in LRv2 similarity matrix
    lr_output = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24_publishedAnno/MLDS_fetalREFliver_LRwCT_output.RDS')
    lr_output = lr_output[[1]]
    
    # Subset to include only relevant cells    
    
    table(mdat$cellID %in% c(rownames(lr_output)))
    lr_output = lr_output[rownames(lr_output) %in% mdat$cellID,]
    
    ## summarise the average similarity by cell type
    # split the matrix into celltype
    #meanLR_mtx = sapply(split(1:nrow(lr_output),mdat$group[match(rownames(lr_output),mdat$cellID)]),function(e){colMeans(lr_output[e,])})
    meanLR_mtx = sapply(split(1:nrow(lr_output),mdat$group[match(rownames(lr_output),mdat$cellID)]),function(e){apply(lr_output[e,],2,median)})
    
    ref_order = c("HSC_MPP",
                  "MEMP","Megakaryocyte","Mast.cell",
                  "Early.Erythroid", "Mid.Erythroid","Late.Erythroid",
                  "Neutrophil.myeloid.progenitor",
                  "Monocyte.precursor","Monocyte","Mono.Mac","Kupffer.Cell",'VCAM1..EI.macrophage',
                  "DC.precursor","DC1","DC2","pDC.precursor",
                  "Pre.pro.B.cell","pro.B.cell","pre.B.cell","B.cell",
                  "ILC.precursor", "Early.lymphoid_T.lymphocyte", "NK",
                  "Hepatocyte","Endothelial.cell","Fibroblast" 
    )
    table(ref_order[! ref_order %in% rownames(meanLR_mtx)])
    table(rownames(meanLR_mtx)[!rownames(meanLR_mtx) %in% ref_order])
    
    
    tgt_order = c(unique(mdat$group[grepl('MLDS',mdat$group)]),
                  unique(mdat$group[grepl('TAM',mdat$group)]),
                  "HSC_MPP",
                  "MEMP_MEP",'MEP',
                  "earlyMK","MK",
                  "Mast.cell",
                  "EE", "ME","LE",
                  "CMP_GMP","promyelocyte","Myelocytes",
                  "proMono","Monocyte","Macrophage","Kupffer.cell",'Neutrophil',
                  "DC1","DC2","pDC",
                  "LMPP_ELP","pro.B.cell","pre.B.cell","B.cell",'Plasma.cell',
                  "ILC.precursor", "T.cell","NK_T", 'NK',
                  "Hepatocyte","Endo","Fibroblast",'Mesenchyme','NPC','others'
    )
    
    
    table(colnames(meanLR_mtx) %in% tgt_order)
    colnames(meanLR_mtx)[!colnames(meanLR_mtx) %in% tgt_order]
    tgt_order = tgt_order[tgt_order %in% colnames(meanLR_mtx)]
    
    mdat$group = factor(mdat$group,levels = tgt_order)
    
    
    table(tgt_order[!tgt_order %in% colnames(meanLR_mtx)])
    table(colnames(meanLR_mtx)[!colnames(meanLR_mtx) %in% tgt_order])
    
    
    meanLR_mtx = (meanLR_mtx[ref_order,tgt_order])
  }
  
  
  
  plotFun_2nAK_fLiver_LRv2_clusterLevel_hm = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    logitCols = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')
    cols = circlize::colorRamp2(seq(-13,10,length.out=length(logitCols)),logitCols)
    
    
    rowAnno = rowAnnotation(Genotype = geno,
                            col = list(Genotype = geno_cols))
    
    hm = Heatmap(meanLR_mtx,show_row_dend = F,show_column_dend = F,
                 col = cols,
                 name=paste0('Predicted Similarity','(Logit)'),
                 show_row_names = T,show_column_names = T,cluster_rows = F,
                 cluster_columns = F,cluster_row_slices = F,
                 row_names_gp = gpar(fontsize=8),
                 column_split = mdat$broadLineage[match(colnames(meanLR_mtx),mdat$group)],
                 column_title_rot = 90,row_title_gp = gpar(fontsize=8),column_title_gp = gpar(fontsize=9),
                 heatmap_legend_param = list(direction = "horizontal"),
                 row_title = 'Reference classes')
    #geno = gsub('.*:','',rownames(meanLR_mtx))
    # geno[geno == 'complete_trisomy'] = 'Triploid'
    # geno = factor(geno,c('diploid','T21','T18','T22','MX','Triploid'))
    # rowAnno = rowAnnotation(Genotype = geno,
    #                         col = list(Genotype = geno_cols))
    # hm = Heatmap(meanLR_mtx,show_row_dend = F,show_column_dend = F,
    #              col = cols,
    #              name=paste0('Predicted\nSimilarity\n','(Logit)'),
    #              show_row_names = F,show_column_names = T,cluster_rows = F,cluster_columns = F,
    #              row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=8),
    #              split = mdat$finalAnn_broad[match(rownames(meanLR_mtx),mdat$group)],
    #              row_title_rot = 0,row_title_gp = gpar(fontsize=8),cluster_row_slices = F,
    #              right_annotation = botAnno,
    #              column_title = 'Reference classes')
    
    
    if(!noPlot){
      draw(hm,heatmap_legend_side = "bottom")
    }
  }
  
  saveFig(file.path(plotDir,'Figxx_2nAKLiv_celltype_UMAP'),plotFun_2nAK_fLiver_LRv2_clusterLevel_hm,rawData=meanLR_mtx,width = 11,height = 5.5,res = 500,useDingbats = F)
}













# # Compute the statistics for Figure 3b -----------------------------------------
# # Old version, the median value does not match the plot in figure 3b
# fig3b_data = read.delim('~/lustre_mt22/Aneuploidy/manuscriptDraft_0424/Plots/Fig3x_medianScorePerGroup_GATA1s_topGenes_rawData.tsv',sep = '\t')
# ggplot(fig3b_data,aes(groupID,medianScore))+
#   geom_point()+
#   facet_grid(.~ group + dataset_ID ,space = 'free_x',scales = 'free_x')+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))


# Figure 4C --------------------------------------------------------------------
mlds_srat_fp <- "~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns.RDS"
mlds_mdat_fp <- "~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns_mdat_2508.csv"
mlds = readRDS(mlds_srat_fp)

l076_fp <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L076/L076_sratObj.RDS"
l038_fp <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L038/L038_sratObj.RDS"
l038_d_mdat_fp <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L038/L038_Diagnostic_AIresult_mdat.csv"
l038_tp1_mdat_fp <- "~/ML-DS/Results/06_MLDS_refractory_relapse/L038/L038_TP1_AIresult_mdat.csv"

##------------------------------------------------##
##    Quantify expression levels of key genes   ####
##------------------------------------------------##
gene =  c('H1F0','ADAMTS1','FHL2','CXCL8','EPS8',#'IFI27','GJA1',
          #'SMAD3',
          'EPS15','MEIS2',#'ARID1B', # up in R2D_vs_D(BM)
          'CD82',
          'LDLR','CD81',
          'COL18A1')
#'CD34','EPCAM','NR4A1')#,
#'CD34','ZFP36','CD7','SAAL1')
#
tmp = mlds@meta.data
cnt = as.data.frame(t(mlds@assays$RNA@data[gene,]))
cnt$cellID = rownames(cnt)
cnt = pivot_longer(cnt,cols = gene,names_to = 'gene',values_to = 'norm_count')
tmp = cbind(tmp[match(cnt$cellID,tmp$cellID),],cnt[,colnames(cnt)!='cellID'])


tmp$group = ifelse(tmp$group_tmp != 'others' & tmp$disease == 'MLDS' & tmp$timePoint == 'Diagnostic' & !tmp$donorID %in% c('L038','L076'),'responsive_MLDS',
                   ifelse(tmp$group_tmp != 'others' & tmp$disease == 'MLDS' & tmp$donorID %in% c('L038','L076'),as.character(tmp$group_tmp),
                          ifelse(tmp$group_tmp != 'others' & tmp$disease == 'TAM' & tmp$timePoint == 'Diagnostic','TAM',as.character(tmp$group_tmp))))


tmp = tmp %>% dplyr::mutate(group = dplyr::case_when(
  !annot %in% c("Tum_MK?",'Tumour','Tumour_WT','unsure_Tumour') ~ 'normal',
  annot == 'Tumour' & disease == 'TAM' & timePoint == 'Diagnostic' ~ 'TAM',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & !donorID %in% c('L038','L076') ~ 'responsive_MLDS',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & donorID %in% c('L038') & cellID %in% l038_d_mdat$cellID[l038_d_mdat$group == 'clone1_D'] ~ 'L038_D_clone1',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & donorID %in% c('L038') & cellID %in% l038_d_mdat$cellID[l038_d_mdat$group == 'clone2_D'] ~ 'L038_D_clone2',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'TP1' & donorID %in% c('L038') & cellID %in% l038$cellID[l038$group == 'clone2']  ~ 'L038_TP1_clone2',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'TP1' & donorID %in% c('L038') & cellID %in% l038$cellID[l038$group == 'clone2_ery'] ~ 'clone2_ery',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & donorID %in% c('L076') & tissue == 'Blood' ~ 'L076_Blood',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'Diagnostic' & donorID %in% c('L076') & tissue == 'BM' ~ 'L076_BM',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'D.Relapse' & donorID %in% c('L076') & tissue == 'BM' ~ 'L076_D.Relapse',
  annot == 'Tumour' & disease == 'MLDS' & timePoint == 'D.Relapse2' & donorID %in% c('L076') & tissue == 'BM' ~ 'L076_D.Relapse2',
  .default = 'others'
),
group = factor(group,c('normal','TAM','responsive_MLDS',
                       'L038_D_clone1','L038_D_clone2','L038_TP1_clone2','clone2_ery',
                       'L076_Blood','L076_BM','L076_D.Relapse','L076_D.Relapse2','others'))) %>% 
  dplyr::filter(group != 'others')

table(!tmp$group %in% c('normal','TAM','responsive_MLDS',
                        'L038_D_clone1','L038_D_clone2','L038_TP1_clone2','clone2_ery',
                        'L038_Diagnostic','L038_TP1','L076_Blood','L076_BM','L076_D.Relapse','L076_D.Relapse2'))                  

# tmp$group[tmp$cellID %in% l038$cellID[l038$group == 'clone1']] = 'L038_D_clone1'
# tmp$group[tmp$cellID %in% l038$cellID[l038$group == 'clone2']] = 'L038_TP1_clone2'
# tmp$group[tmp$cellID %in% l038$cellID[l038$group == 'clone2_D']] = 'L038_D_clone2'
# tmp = tmp[tmp$group !='L038_TP1',]

tmp$group2 = ifelse(!tmp$donorID %in% c('L038','L076','L156'),'Remission',
                    ifelse(tmp$donorID == 'L076','Relapse',
                           ifelse(tmp$donorID == 'L156','Recurrent',
                                  ifelse(tmp$donorID == 'L038','Refractory','others'))))
tmp$group2[tmp$group == 'normal'] = 'Normal'
table(tmp$group2)
table(tmp$group)
tmp$group2 = factor(tmp$group2,c('Normal','Remission','Recurrent','Relapse','Refractory'))


tmp = tmp[tmp$gene %in% gene,]
tmp$gene = factor(tmp$gene,gene)
library(ggbeeswarm)
ccs = c('Normal'=grey(0.8),'Remission'='#2D4372','Refractory'=col25[2],'Relapse'=col25[5],'Recurrent'=col25[5])

#tmp = read.delim(file.path(plotDir,'Fig4E_badMLDS_markers_rawData.tsv'),sep = '\t')

plot_normalisedExpression = function(noFrame=FALSE,noPlot=FALSE){
  p = ggplot(tmp[tmp$gene %in% gene,],aes(group,norm_count,fill=group2))+
    geom_boxplot(outliers=F,outlier.shape = NA,colour='black',width=0.7,linewidth=0.3)+
    scale_fill_manual(values = ccs)+
    facet_grid(gene~group2,scales = 'free',space = 'free_x')+
    #geom_quasirandom(size=0.1,alpha=0.1)+
    theme_classic(base_size = 9)+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          strip.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(colour='black'),
          panel.spacing.y = unit(0.1,'cm'),
          axis.ticks = element_line(colour='black'),axis.title = element_text(colour='black'),
          panel.border = element_rect(fill=F,colour='black'))+xlab('')+ylab('Normalized expression')
  print(p)
}

#saveFig(file.path(plotDir,paste0('Fig4E_badMLDS_markers')),plot_normalisedExpression,rawData=tmp,width = 4.5,height = 13,res = 500,useDingbats = F)
saveFig(file.path(plotDir,paste0('Fig4E_badMLDS_markers')),plot_normalisedExpression,rawData=NULL,width = 4,height = 6,res = 500,useDingbats = F)
saveFig(file.path(plotDir,paste0('Fig4E_badMLDS_markers_hor')),plot_normalisedExpression,rawData=tmp,width = 14,height = 4.5,res = 500,useDingbats = F)
saveFig(file.path(plotDir,paste0('Fig4E_badMLDS_markers_hor')),plot_normalisedExpression,rawData=tmp,width = 10,height = 4,res = 500,useDingbats = F)

##------------------------------------##
##        Supplementary Tables      ####
##------------------------------------##
#Define genomic coordinates
library(GenomicFeatures)
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

geneMap = read.delim('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data/Donor15680/liver/fLiver_MY_200531_10043298/filtered_feature_bc_matrix/features.tsv.gz',header = F)
colnames(geneMap) = c('ensID','geneSym','GEX')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))




## Table S1: List of samples used in the study (sc.Leuk + sc.fetalLiver + bulk.Samples) ####
shared_columns = c('sc_or_bulk','dataset','donorID','Genotype','Sex','age_category','Disease','channelID','timePoint','Tissue','Clinical.outcome','blastPerc',
                   'sortingStrategy',"assay","reference_genome","cellranger_version",'nCell','nLeuk',"nNorm","PDID")


##--- sc.TAM.MLDS
scTAM.MLDS = read.csv('~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/TableS1_MLDS_dataset.csv')
scTAM.MLDS$sc_or_bulk = 'single-cell RNA-seq'
scTAM.MLDS$dataset = 'TAM / ML-DS'
scTAM.MLDS$age_category = ifelse(is.na(as.numeric(scTAM.MLDS$Age..year.)),'unknown',
                               ifelse(as.numeric(scTAM.MLDS$Age..year.) < 0.5,'< 6 months','6 months - 5 years'))
scTAM.MLDS$channelID = gsub('\\.','_',scTAM.MLDS$channelID)
scTAM.MLDS$channelID = ifelse(grepl('^RNA\\d+',scTAM.MLDS$channelID),paste0('ALeuk_',scTAM.MLDS$channelID),
                              ifelse(grepl('^Leuk\\d+',scTAM.MLDS$channelID),paste0('SB_',scTAM.MLDS$channelID),scTAM.MLDS$channelID))
table(scTAM.MLDS$channelID)
scTAM.MLDS$timePoint[scTAM.MLDS$timePoint == "Relapse diagnostic"] = "Relapse 1 diagnosis"
scTAM.MLDS$assay = "scRNA 5' v2"
scTAM.MLDS$reference_genome = "GRCh38 2020-A"
scTAM.MLDS$cellranger_version = "cellranger_700"
scTAM.MLDS$PDID[is.na(scTAM.MLDS$PDID)] = '-'
scTAM.MLDS$sortingStrategy = 'none'
table(shared_columns %in% colnames(scTAM.MLDS))
shared_columns[!shared_columns %in% colnames(scTAM.MLDS)]

##--- sc.Leukaemia
scLeuk = read.csv('~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/TableS1_otherLeuk_dataset.csv')
scLeuk$sc_or_bulk = 'single-cell RNA-seq'
scLeuk$dataset = 'Other leukaemia'
scLeuk$age_category = NA
scLeuk$age_category[scLeuk$Age..yrs. %in% c('0m','infant','6m','4m')] = '< 1 year'
scLeuk$age_category[scLeuk$Age..yrs. %in% c('1Y','2y,1m','2y,2m','2y,9m','3y,5m','4Y,2m')] = '1-5 years'
scLeuk$age_category[scLeuk$Age..yrs. %in% c('5Y,2m','5y,4m','6y,2m','7y','7y,9m','8y,11m','11y','11y, 6m')] = '5-12 years'
table(is.na(scLeuk$age_category))
table(scLeuk$age_category,scLeuk$Age..yrs.)
scLeuk$timePoint[scLeuk$timePoint == "Relapse diagnostic"] = 'Relapse 1 diagnosis'
scLeuk$channelID = gsub('\\.','_',scLeuk$channelID)
scLeuk$channelID = ifelse(grepl('^RNA\\d+',scLeuk$channelID),paste0('ALeuk_',scLeuk$channelID),
                              ifelse(grepl('^Leuk\\d+',scLeuk$channelID),paste0('SB_',scLeuk$channelID),
                                     ifelse(grepl('^NB\\d+',scLeuk$channelID),paste0('CG_SB_',scLeuk$channelID),scLeuk$channelID)))
table(scLeuk$channelID)

scLeuk$assay = "scRNA 5' v2"
scLeuk$reference_genome = "GRCh38 2020-A"
scLeuk$cellranger_version = "cellranger_700"
scLeuk$blastPerc = scLeuk$Blast.percentage
scLeuk$PDID = '-'
scLeuk$sortingStrategy = 'none'

table(shared_columns %in% colnames(scLeuk))
shared_columns[!shared_columns %in% colnames(scLeuk)]

##--- sc.fetal.Livers
sc.fLivers = read.csv('~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/TableS1_fLiver_dataset_nCellBy10XSample.csv')
sc.fLivers$sc_or_bulk = 'single-cell RNA-seq'
sc.fLivers$dataset = 'Fetal liver'
sc.fLivers$tissue = gsub('Foetal','Fetal',sc.fLivers$tissue)
sc.fLivers$age_category = sc.fLivers$gestationalAge
unique(sc.fLivers$sangerSampleID)
sc.fLivers$channelID = sc.fLivers$sangerSampleID
sc.fLivers$reference_genome = sc.fLivers$refGenome_version
sc.fLivers$cellranger_version = sc.fLivers$cellRanger_version
sc.fLivers$Tissue = sc.fLivers$tissue

table(shared_columns %in% colnames(sc.fLivers))
shared_columns[!shared_columns %in% colnames(sc.fLivers)]
for(col in unique(shared_columns[!shared_columns %in% colnames(sc.fLivers)])){
  sc.fLivers[[col]]  = '-'
}


##--- bulk samples
bulk = readxl::read_excel('~/lustre_mt22/Hennings_bulkRNAseq_metadata.xlsx')
bulk = bulk[!bulk$`included for analysis?` %in% c("Removed (Konstantin email)","Removed (PDX)","Removed (unclassified)"),]
bulk$sc_or_bulk = 'bulk RNA-seq'
bulk$dataset = bulk$Group
bulk$dataset[bulk$dataset == 'TAM_MLDS'] = 'TAM / ML-DS (bulk)'
bulk$dataset[bulk$dataset == 'Other Leukaemia'] = 'Other leukaemia (bulk)'
bulk$dataset[bulk$dataset == 'Nomal_HSPCs'] = 'Nomal HSPCs (bulk)'
bulk$Tissue = bulk$Source
bulk$Tissue[bulk$Tissue == 'cord blood'] = 'Cord blood'
bulk$Tissue[bulk$Tissue == 'fetal liver cells'] = 'Fetal liver cells'
bulk$Tissue[bulk$Tissue == 'patient sample'] = 'Patient sample'
bulk$Tissue[bulk$Tissue == 'peripheral blood'] = 'Peripheral blood'
bulk$donorID = bulk$Sample
bulk$donorID = gsub('t_\\d*_\\d*','t',bulk$donorID)
bulk$donorID = sapply(strsplit(bulk$donorID,'_'),'[',3)
bulk$donorID[grepl('hFL',bulk$Sample)] = gsub('_5_2','_5',gsub('_4_2','_4',gsub('_2_2','_2',gsub('healthy_.*_hFL-','',bulk$Sample[grepl('hFL',bulk$Sample)]))))
bulk$donorID[grepl('healthy_.*_CB_',bulk$Sample)] = gsub('healthy_.*_CB_','',bulk$Sample[grepl('healthy_.*_CB_',bulk$Sample)])

bulk$Genotype = ifelse(bulk$dataset == 'TAM / ML-DS (bulk)','Trisomy 21','Diploid')
bulk$Disease = ifelse(bulk$dataset != 'Nomal HSPCs (bulk)', 
                      ifelse(grepl('t\\(.*',bulk$Subgroup),'MLL',bulk$Subgroup),'-')
bulk$Disease[bulk$Disease == 'TMD'] = 'TAM'
bulk$Disease[bulk$Disease == 'MLDS'] = 'ML-DS'
bulk$bulk_sampleID = bulk$Sample
bulk$timePoint = ifelse(bulk$dataset != 'Nomal HSPCs (bulk)', 'Diagnostic','-')
bulk$Clinical.outcome = ifelse(bulk$Disease == 'TAM', bulk$Category,'-')
bulk$Clinical.outcome[bulk$Clinical.outcome == "TAM_MLDS"] = "Progressive TAM"
bulk$Clinical.outcome[bulk$Clinical.outcome == "TAM_conventional"] = "Conventional TAM"
bulk$Clinical.outcome[bulk$Clinical.outcome == "TAM_earlyDeath"] = "Early death TAM"
bulk$Clinical.outcome[bulk$Clinical.outcome == "TAM_progressive"] = "Progressive TAM"
bulk$age_category = 'unknown'
bulk$channelID = bulk$Sample
bulk$assay = 'bulk RNA-seq'
bulk$reference_genome = 'GRCh38'

## Determine sex
library(edgeR)
bulkSamples = import_HenningsBulkSamples(gns=gns,rm_henning_healthy=F,oxfordBulk=F,filter_lowExpr_gene=T,
                                         tpm_fp = c('~/lustre_mt22/Down_Leukemia/Down_Leukemia_Klusmann_StarSalmon_tpm.csv'))
tpm_count = bulkSamples[['tpm_count']]
# Determine sex
tpm_sex = tpm_count[geneMap$ensID[geneMap$geneSym %in% c('XIST','RPS4Y1')],]
rownames(tpm_sex) = geneMap$geneSym[match(rownames(tpm_sex),geneMap$ensID)]
sex = apply(tpm_sex,2,function(s){
  names(s)=rownames(tpm_sex)
  if(s['XIST'] > s['RPS4Y1']){
    return('F')
  }else{
    return('M')
  }
})
bulk$Sex = sex[match(gsub('-','.',bulk$Sample),names(sex))]
bulk$Sex[bulk$Sex == 'M'] = 'male'
bulk$Sex[bulk$Sex == 'F'] = 'female'
bulk$timePoint[bulk$Sample %in% c('Patient_TMD_2724_2','Patient_t_10_11_3412_2','Patient_t_10_11_3697_2')] = 'unknown'
bulk$sortingStrategy = 'FACS-sorted'

table(shared_columns %in% colnames(bulk))
shared_columns[!shared_columns %in% colnames(bulk)]
for(col in unique(shared_columns[!shared_columns %in% colnames(bulk)])){
  bulk[[col]]  = '-'
}


##--- combine all 4 tables
df = do.call(rbind,list(scTAM.MLDS[,shared_columns],
                        scLeuk[,shared_columns],
                        sc.fLivers[,shared_columns],
                        bulk[bulk$Disease == '-',shared_columns],
                        bulk[bulk$Disease == 'TAM',shared_columns],
                        bulk[bulk$Disease == 'ML-DS',shared_columns],
                        bulk[bulk$Disease == 'AMKL',shared_columns],
                        bulk[bulk$Disease == 'MLL',shared_columns]))
df$Sex[df$Sex == 'male'] = 'Male'
df$Sex[df$Sex == 'female'] = 'Female'
df$timePoint = gsub('iagnostic','iagnosis',df$timePoint)
df$blastPerc[df$blastPerc %in% c('?','-')] = 'unknown'
df$Clinical.outcome[df$Disease !='-' & df$Clinical.outcome %in% c('-','TAM_unknown')] = 'unknown'
df$assay[df$assay == "scRNA 3' v3.1"] = "scRNA 3' v3"

colnames(df) = c("single-cell or bulk RNA-seq",	"Dataset",	"Donor ID",	"Genotype",	"Sex",	"Age category", "Disease",	
"Sample ID",	"Timepoint",	"Tissue",	"Clinical outcome",	"%blasts by flow cytometry", "Sorting strategy"	, "10X assay",	"Reference genome version",	"CellRanger version",	
"Number of cells post QC",	"Number of cancer cells post QC",	"Number of normal cells post QC",	"WGS ID")
for(i in 1: ncol(df)){
  print(colnames(df)[i])
  #print(table(df[[i]]))
  #print(sum(is.na(df[[i]])))
  print(sum(df[[i]] == '-'))

}
table(df$`WGS ID`,df$Dataset)
table(df$`Age category`,df$Dataset)
table(df$`CellRanger version`,df$Dataset)
table(df$`10X assay`,df$`CellRanger version`)
table(df$`Donor ID`[is.na(df$`WGS ID`)])
table(df$Dataset[df$`Donor ID` == '-'])
table(is.na(df))

for(d in unique(df$Dataset)){
  print(d)
  print(n_distinct(df$`Donor ID`[df$Dataset == d]))
}

for(d in unique(df$Disease[df$Dataset == 'Other leukaemia (bulk)'])){
  print(d)
  print(n_distinct(df$`Donor ID`[df$Disease == d & df$Dataset == 'Other leukaemia (bulk)']))
}
View(df[df$Dataset == 'Other leukaemia (bulk)' & df$Disease=='MLL',])


write.csv(df,file.path(plotDir,'../TableS1_leukaemia_and_fLiver_datasets_overview.csv'),row.names = F)





## Table S3: Differential Gene Expression analysis - Foetal liver ####
resultDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jul24/without_published2n/genoAssay_withX'
tissue = 'liver'

# Import log2FC output
allGenes = import_pbDEGresults(outDir = resultDir, allGenes = T,tissue='liver')
allGenes$geneSym = geneMap$geneSym[match(allGenes$ensID,geneMap$ensID)]
allGenes$chr = geneMap$chr[match(allGenes$ensID,geneMap$ensID)]
allGenes$geno[allGenes$geno == 'complete'] = '3n'
allGenes$comp = paste0(allGenes$geno,'_',allGenes$ct,'_',allGenes$geneSym)

## Mark which genes are considered DEG
allDEGs = import_pbDEGresults(outDir = resultDir,tissue = 'liver')
allDEGs$comp = paste0(allDEGs$geno,'_',allDEGs$ct,'_',allDEGs$geneSym)
table(allGenes$comp %in% allDEGs$comp,allGenes$geno)

allGenes$is_DEG = F
allGenes$is_DEG[allGenes$comp %in% allDEGs$comp] = T
table(allGenes$geno[allGenes$is_DEG & !grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',allGenes$geneSym)])


## Prepare to write the final table
allGenes = allGenes[,c('ct','geno','ensID','geneSym','chr','baseMean','log2FoldChange','lfcSE','pvalue','padj','is_DEG')]
allGenes$geno[allGenes$geno == '3n'] = 'Triploid'
## Re-write cell type names
allGenes$ct[allGenes$ct == 'EE'] = 'early Ery'
allGenes$ct[allGenes$ct == 'ME'] = 'mid Ery'
allGenes$ct[allGenes$ct == 'LE'] = 'late Ery'
allGenes$ct[allGenes$ct == 'Endo'] = 'Endothelium'
allGenes$ct[allGenes$ct == 'Mesothelial_cells'] = 'Mesothelial cell'
allGenes$ct[allGenes$ct == 'Cholangiocytes'] = 'Cholangiocyte'
allGenes$ct[allGenes$ct == 'NK.T'] = 'NK / T'
allGenes$ct = gsub('\\.',' ',allGenes$ct)
allGenes$ct = gsub('_',' / ',allGenes$ct)

colnames(allGenes) = c('cell_type','genotype','ensembl_ID','gene_symbol','chromosome','baseMean','log2FoldChange','lfcSE','pvalue','padj','is_DEG')

write.csv(allGenes,file.path(plotDir,'../TableS3_fLiver_DEanalysis_allGenes.csv'),row.names = F)






## Table S4: Differential Gene Expression analysis - T21 leukaemic gene module ####

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


## Import MLDS-vs-FL_diploid result
outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/3_MLDSimprint_in_fLiverT21/jul24/'
setwd(outDir)
## list of DEGs
mlds_vs_2n_deg = read.csv('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset_topDEGs.csv')
## All genes 
mlds_vs_2n_res = readRDS('DESeq2_dMLDS.vs.2nfLiverMEMP_dataset_res.RDS')
mlds_vs_2n_res = as.data.frame(mlds_vs_2n_res[['res']])
mlds_vs_2n_res$geneSym = geneMap$geneSym[match(rownames(mlds_vs_2n_res),geneMap$ensID)]
mlds_vs_2n_res$chr = geneMap$chr[match(rownames(mlds_vs_2n_res),geneMap$ensID)]
mlds_vs_2n_res$ensID = rownames(mlds_vs_2n_res)
mlds_vs_2n_res$is_DEG = (mlds_vs_2n_res$geneSym %in% mlds_vs_2n_deg$geneSym)
mlds_vs_2n_res$direction = ifelse(mlds_vs_2n_res$is_DEG & mlds_vs_2n_res$log2FoldChange > 0,'upregulated_in_MLDS',
                                  ifelse(mlds_vs_2n_res$is_DEG & mlds_vs_2n_res$log2FoldChange < 0,'downregulated_in_MLDS','not_DE'))

## Import FL_T21-vs-FL_diploid result
t21_vs_2n_result = read.csv('DESeq2_T21.vs.diploid.fLiver.MEMP_geno.assay_results_allGenes.csv',row.names = 1)
mlds_vs_2n_deg = mlds_vs_2n_deg %>% filter(!is.na(padj)) %>%  filter(padj < 0.05)
mlds_vs_2n_deg$DE = ifelse((mlds_vs_2n_deg$cellFrac_g1 >= max_cellfrac | mlds_vs_2n_deg$cellFrac_g2 >= max_cellfrac),T,F)
mlds_degs = mlds_vs_2n_deg[mlds_vs_2n_deg$DE == T,]

## subset for MLDS DEGs
max_cellfrac = 0.1
mlds_degs_in_T21 = t21_vs_2n_result[t21_vs_2n_result$geneSym %in% mlds_degs$geneSym,]
mlds_degs_in_T21_filtered = rbind(t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_up'] & t21_vs_2n_result$log2FoldChange > 0 & !is.na(t21_vs_2n_result$pvalue),],
                             t21_vs_2n_result[abs(t21_vs_2n_result$cellFrac_max) >= max_cellfrac & t21_vs_2n_result$geneSym %in% mlds_degs$geneSym[mlds_degs$direction == 'MLDS_down'] & t21_vs_2n_result$log2FoldChange < 0 & !is.na(t21_vs_2n_result$pvalue),])
mlds_imprints_in_T21 = read.csv('MLDS_imprint_in_T21.fLiver.MEMP_2411.csv',row.names = 1)

mlds_degs_in_T21$is_DEG = ifelse(mlds_degs_in_T21$cellFrac_max < max_cellfrac,'lowly_expressed',
                                            ifelse(!mlds_degs_in_T21$geneSym %in% mlds_degs_in_T21_filtered$geneSym,'opposite_direction',
                                                   ifelse(mlds_degs_in_T21$geneSym %in% mlds_imprints_in_T21$geneSym,'DE','not_significant')))
mlds_degs_in_T21$padj[mlds_degs_in_T21$is_DEG != 'DE'] = NA
mlds_degs_in_T21$padj[mlds_degs_in_T21$is_DEG == 'DE'] = mlds_imprints_in_T21$padj[match(mlds_degs_in_T21$geneSym[mlds_degs_in_T21$is_DEG == 'DE'],mlds_imprints_in_T21$geneSym)]

## Combine the two tables
mlds_vs_2n_res = mlds_vs_2n_res[order(mlds_vs_2n_res$padj,decreasing = F),]
mlds_vs_2n_res$direction = factor(mlds_vs_2n_res$direction,c('upregulated_in_MLDS','downregulated_in_MLDS'))
mlds_vs_2n_res = mlds_vs_2n_res[order(mlds_vs_2n_res$direction),]
mlds_vs_2n_res = rbind(mlds_vs_2n_res[mlds_vs_2n_res$is_DEG,],mlds_vs_2n_res[!mlds_vs_2n_res$is_DEG,])

# Prepare the table 
mlds_degs_in_T21 = mlds_degs_in_T21[,colnames(mlds_vs_2n_res)[colnames(mlds_vs_2n_res) !='direction']]
colnames(mlds_degs_in_T21) = paste0(colnames(mlds_degs_in_T21),'_T21')
mlds_imprints_in_T21 = cbind(mlds_vs_2n_res,mlds_degs_in_T21[match(mlds_vs_2n_res$geneSym,mlds_degs_in_T21$geneSym_T21),])
mlds_imprints_in_T21 = mlds_imprints_in_T21[,c("ensID","geneSym","chr","baseMean","log2FoldChange","lfcSE","pvalue","padj","is_DEG","direction",
                                               "ensID_T21","geneSym_T21","chr_T21","baseMean_T21","log2FoldChange_T21","lfcSE_T21","pvalue_T21","padj_T21","is_DEG_T21")]
colnames(mlds_imprints_in_T21) = c('ensembl_ID (MLDS-vs-diploid.FL.MEP)','gene_symbol (MLDS-vs-diploid.FL.MEP)','chromosome (MLDS-vs-diploid.FL.MEP)','baseMean (MLDS-vs-diploid.FL.MEP)','log2FoldChange (MLDS-vs-diploid.FL.MEP)','lfcSE (MLDS-vs-diploid.FL.MEP)','pvalue (MLDS-vs-diploid.FL.MEP)','padj (MLDS-vs-diploid.FL.MEP)','is_DEG (MLDS-vs-diploid.FL.MEP)','direction (MLDS-vs-diploid.FL.MEP)',
                                   'ensembl_ID (T21.FL.MEP-vs-diploid.FL.MEP)','gene_symbol (T21.FL.MEP-vs-diploid.FL.MEP)','chromosome (T21.FL.MEP-vs-diploid.FL.MEP)','baseMean (T21.FL.MEP-vs-diploid.FL.MEP)','log2FoldChange (T21.FL.MEP-vs-diploid.FL.MEP)','lfcSE (T21.FL.MEP-vs-diploid.FL.MEP)','pvalue (T21.FL.MEP-vs-diploid.FL.MEP)','padj (T21.FL.MEP-vs-diploid.FL.MEP)','is_DEG (T21.FL.MEP-vs-diploid.FL.MEP)')


write.csv(mlds_imprints_in_T21,file.path(plotDir,'../TableS4_MLDSimprints_in_T21.fLiver.MEMP.csv'),row.names = F)




## Table S5: Differential Gene Expression analysis - GATA1s gene module ####
outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/oct24'
setwd(outDir)

## Import MLDS-vs-bonemarrow_MEP result

#-- DEG
mlds.vs.gMEP_deg = read.csv('DESeq2_dMLDS.vs.goshMEP_topDEGs.csv')
#-- all genes
mlds.vs.gMEP_res = readRDS('DESeq2_dMLDS.vs.goshMEP_res.RDS')
mlds.vs.gMEP_res = as.data.frame(mlds.vs.gMEP_res[['res']])
mlds.vs.gMEP_res$geneSym = geneMap$geneSym[match(rownames(mlds.vs.gMEP_res),geneMap$ensID)]
mlds.vs.gMEP_res$chr = geneMap$chr[match(rownames(mlds.vs.gMEP_res),geneMap$ensID)]
mlds.vs.gMEP_res$ensID = rownames(mlds.vs.gMEP_res)
mlds.vs.gMEP_res$is_DEG = (mlds.vs.gMEP_res$geneSym %in% mlds.vs.gMEP_deg$geneSym)
mlds.vs.gMEP_res$direction = ifelse(mlds.vs.gMEP_res$is_DEG & mlds.vs.gMEP_res$log2FoldChange > 0,'upregulated_in_MLDS',
                                  ifelse(mlds.vs.gMEP_res$is_DEG & mlds.vs.gMEP_res$log2FoldChange < 0,'downregulated_in_MLDS','not_DE'))


## Import TAM-vs-bonemarrow_MEP result
#-- all genes
tam.vs.gMEP_results = read.csv('DESeq2_goodTAM.vs.goshMEP_result_allGene.csv')
#-- subset to just DEG in MLDS
mlds_degs_in_TAM = tam.vs.gMEP_results[tam.vs.gMEP_results$geneSym %in% mlds.vs.gMEP_deg$geneSym,]
#-- subset to just DEG in MLDS + meet criteria
mlds_degs_in_TAM_filtered = rbind(tam.vs.gMEP_results[abs(tam.vs.gMEP_results$cellFrac_max) >= max_cellfrac & tam.vs.gMEP_results$geneSym %in% mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_up'] & tam.vs.gMEP_results$log2FoldChange > 0 & !is.na(tam.vs.gMEP_results$pvalue),],
                                  tam.vs.gMEP_results[abs(tam.vs.gMEP_results$cellFrac_max) >= max_cellfrac & tam.vs.gMEP_results$geneSym %in% mlds.vs.gMEP_deg$geneSym[mlds.vs.gMEP_deg$direction == 'MLDS_down'] & tam.vs.gMEP_results$log2FoldChange < 0 & !is.na(tam.vs.gMEP_results$pvalue),])
#-- final gene module
gata1s_module = read.csv('DESeq2_goodTAM.vs.goshMEP_MLDSgenes.csv',row.names = 1)

mlds_degs_in_TAM$is_DEG = ifelse(mlds_degs_in_TAM$cellFrac_max < max_cellfrac,'lowly_expressed',
                                 ifelse(!mlds_degs_in_TAM$geneSym %in% mlds_degs_in_TAM_filtered$geneSym,'opposite_direction',
                                        ifelse(mlds_degs_in_TAM$geneSym %in% gata1s_module$geneSym,'DE','not_significant')))
table(mlds_degs_in_TAM$is_DEG)
mlds_degs_in_TAM$padj[mlds_degs_in_TAM$is_DEG != 'DE'] = NA
mlds_degs_in_TAM$padj[mlds_degs_in_TAM$is_DEG == 'DE'] = gata1s_module$padj[match(mlds_degs_in_TAM$geneSym[mlds_degs_in_TAM$is_DEG == 'DE'],gata1s_module$geneSym)]

## Combine the two tables
mlds.vs.gMEP_res = mlds.vs.gMEP_res[order(mlds.vs.gMEP_res$padj,decreasing = F),]
mlds.vs.gMEP_res$direction = factor(mlds.vs.gMEP_res$direction,c('upregulated_in_MLDS','downregulated_in_MLDS','not_DE'))
mlds.vs.gMEP_res = mlds.vs.gMEP_res[order(mlds.vs.gMEP_res$direction),]
mlds.vs.gMEP_res = rbind(mlds.vs.gMEP_res[mlds.vs.gMEP_res$is_DEG,],mlds.vs.gMEP_res[!mlds.vs.gMEP_res$is_DEG,])
table(mlds.vs.gMEP_res$direction,mlds.vs.gMEP_res$is_DEG)
dim(mlds.vs.gMEP_res)

# Prepare the table 
mlds_degs_in_TAM = mlds_degs_in_TAM[,colnames(mlds.vs.gMEP_res)[colnames(mlds.vs.gMEP_res) !='direction']]
colnames(mlds_degs_in_TAM) = paste0(colnames(mlds_degs_in_TAM),'_TAM')
gata1s_module = cbind(mlds.vs.gMEP_res,mlds_degs_in_TAM[match(mlds.vs.gMEP_res$geneSym,mlds_degs_in_TAM$geneSym_TAM),])
gata1s_module = gata1s_module[,c("ensID","geneSym","chr","baseMean","log2FoldChange","lfcSE","pvalue","padj","is_DEG","direction",
                                 "ensID_TAM","geneSym_TAM","chr_TAM","baseMean_TAM","log2FoldChange_TAM","lfcSE_TAM","pvalue_TAM","padj_TAM","is_DEG_TAM")]
colnames(gata1s_module) = c('ensembl_ID (MLDS-vs-T21.MEP)','gene_symbol (MLDS-vs-T21.MEP)','chromosome (MLDS-vs-T21.MEP)','baseMean (MLDS-vs-T21.MEP)','log2FoldChange (MLDS-vs-T21.MEP)','lfcSE (MLDS-vs-T21.MEP)','pvalue (MLDS-vs-T21.MEP)','padj (MLDS-vs-T21.MEP)','is_DEG (MLDS-vs-T21.MEP)','direction (MLDS-vs-T21.MEP)',
                            'ensembl_ID (TAM-vs-T21.MEP)','gene_symbol (TAM-vs-T21.MEP)','chromosome (TAM-vs-T21.MEP)','baseMean (TAM-vs-T21.MEP)','log2FoldChange (TAM-vs-T21.MEP)','lfcSE (TAM-vs-T21.MEP)','pvalue (TAM-vs-T21.MEP)','padj (TAM-vs-T21.MEP)','is_DEG (TAM-vs-T21.MEP)')


write.csv(gata1s_module,file.path(plotDir,'../TableS5_GATA1s_gene_module.csv'),row.names = F)





## Table S6: Differential Gene Expression analysis - MLDS gene module ####
outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/5_goodTAM_vs_MLDS/jul24'
setwd(outDir)

## Import MLDS-vs-TAM result

#-- DEGs between TAM and MLDS only
tam.vs.all.mlds_deg = read.csv('DESeq2_goodTAM.vs.all.dMLDS_topDEGs_2410.csv')

#-- all genes
tam.vs.all.mlds_res = readRDS('DESeq2_goodTAM.vs.all.dMLDS_res_2410.RDS')
tam.vs.all.mlds_res = as.data.frame(tam.vs.all.mlds_res[['res']])
tam.vs.all.mlds_res$ensID = rownames(tam.vs.all.mlds_res)
tam.vs.all.mlds_res$geneSym = geneMap$geneSym[match(tam.vs.all.mlds_res$ensID,geneMap$ensID)]
tam.vs.all.mlds_res$chr = geneMap$chr[match(tam.vs.all.mlds_res$ensID,geneMap$ensID)]

tam.vs.all.mlds_res$is_DEG = (tam.vs.all.mlds_res$geneSym %in% tam.vs.all.mlds_deg$geneSym)
tam.vs.all.mlds_res$direction = ifelse(tam.vs.all.mlds_res$is_DEG == T & tam.vs.all.mlds_res$log2FoldChange > 0, 'upregulated_in_MLDS',
                                       ifelse(tam.vs.all.mlds_res$is_DEG == T & tam.vs.all.mlds_res$log2FoldChange < 0, 'downregulated_in_MLDS','not_DE'))

table(tam.vs.all.mlds_res$direction,tam.vs.all.mlds_res$is_DEG)
dim(tam.vs.all.mlds_res)

# Prepare the table 
tam.vs.all.mlds_res$direction = factor(tam.vs.all.mlds_res$direction, c('upregulated_in_MLDS','downregulated_in_MLDS','not_DE'))
tam.vs.all.mlds_res = tam.vs.all.mlds_res[order(tam.vs.all.mlds_res$padj,decreasing = F),]
tam.vs.all.mlds_res = tam.vs.all.mlds_res[order(tam.vs.all.mlds_res$direction),]
tam.vs.all.mlds_res = rbind(tam.vs.all.mlds_res[tam.vs.all.mlds_res$is_DEG,],tam.vs.all.mlds_res[!tam.vs.all.mlds_res$is_DEG,])

tam.vs.all.mlds_res = tam.vs.all.mlds_res[,c("ensID","geneSym","chr","baseMean","log2FoldChange","lfcSE","pvalue","padj","is_DEG","direction")]
colnames(tam.vs.all.mlds_res) = c('ensembl_ID','gene_symbol','chromosome','baseMean','log2FoldChange','lfcSE','pvalue','padj','is_DEG','direction')


write.csv(tam.vs.all.mlds_res,file.path(plotDir,'../TableS6_MLDS_gene_module.csv'),row.names = F)
#tam.vs.all.mlds_res = read.csv(file.path(plotDir,'../TableS6_MLDS_gene_module.csv'))




cosmicGenes = read.delim('~/lustre_mt22/generalResources/COSMIC_v100_202408/Cosmic_CancerGeneCensus_v100_GRCh38.tsv.gz',sep = '\t')
cosmicMutation = read.delim('~/lustre_mt22/generalResources/COSMIC_v100_202408/Cosmic_MutantCensus_v100_GRCh38.tsv.gz',sep = '\t')
cosmicMutation_tally = cosmicMutation %>% group_by(GENE_SYMBOL) %>% summarise(n_mutation = n(),n_mutationID = n_distinct(MUTATION_ID))

## Table S7: List of somatic mutation in L076 ####
mutations = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L076/L076_4samples_finalSomaticVariants_jbrowseFiltered.csv')
## Only keep the intersection between Blood and BM - diagnosis
mutations = mutations[mutations$somaticVar_type %in% c('shared','Blood_D:D','RD:R2D','RD','R2D'),]
mutations = mutations[,c('Chrom','Pos','Ref','Alt','Gene','Type','Impact','AAchange','PD62331c','PD62331a', 'PD64665a', 'PD66167a','somaticVar_type')]															
mutations$COSMIC_gene_tier = ifelse(mutations$Gene %in% cosmicGenes$GENE_SYMBOL,cosmicGenes$TIER[match(mutations$Gene,cosmicGenes$GENE_SYMBOL)],'-')
mutations$COSMIC_mutation_count = ifelse(mutations$Gene %in% cosmicMutation_tally$GENE_SYMBOL,cosmicMutation_tally$n_mutation[match(mutations$Gene,cosmicMutation_tally$GENE_SYMBOL)],'-')

mutations$mutation_group = as.character(mutations$somaticVar_type)
mutations$mutation_group[mutations$mutation_group == 'shared'] = 'shared across all samples'
mutations$mutation_group[mutations$mutation_group == 'Blood_D:D'] = 'unique to diagnostic clone'
mutations$mutation_group[mutations$mutation_group == 'RD:R2D'] = 'shared between relapse 1 and relapse 2 clones'
mutations$mutation_group[mutations$mutation_group == 'RD'] = 'unique to relapse 1 clone'
mutations$mutation_group[mutations$mutation_group == 'R2D'] = 'unique to relapse 2 clone'
mutations$somaticVar_type = factor(mutations$somaticVar_type,c('shared','Blood_D:D','RD:R2D','RD','R2D'))
mutations$Chrom = factor(mutations$Chrom,paste0('chr',c(1:22,'X','Y')))
mutations = mutations[order(mutations$Pos,decreasing = F),]
mutations = mutations[order(mutations$Chrom),]
mutations = mutations[order(mutations$somaticVar_type),]

colnames(mutations) = c('Chrom','Pos','Ref','Alt','Gene','Type','Impact','AAchange','PD62331c_D.Blood_VAF','PD62331a_D.BM_VAF', 'PD64665a_Relapse1_VAF', 'PD66167a_Relapse2_VAF','somaticVar_type','COSMIC_gene_tier','COSMIC_mutation_count','mutation_group')
mutations = mutations[,c('Chrom','Pos','Ref','Alt','Gene','Type','Impact','AAchange','COSMIC_gene_tier','COSMIC_mutation_count','mutation_group','PD62331c_D.Blood_VAF','PD62331a_D.BM_VAF', 'PD64665a_Relapse1_VAF', 'PD66167a_Relapse2_VAF')]

write.csv(mutations,file.path(plotDir,'../TableS7_L076_RelapseMLDS_somaticMutationCatalogue.csv'),row.names = F)



## Table S8: List of somatic mutation in L038 ####
mutations = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/L038_D.vs.TP1_finalSomaticVariants.csv')
## Only keep the intersection between Blood and BM - diagnosis
mutations = mutations[,c('Chrom','Pos','Ref','Alt','Gene','Type','Impact','AAchange','PD60301a','PD61846a','somaticVar_type')]															
mutations$COSMIC_gene_tier = ifelse(mutations$Gene %in% cosmicGenes$GENE_SYMBOL,cosmicGenes$TIER[match(mutations$Gene,cosmicGenes$GENE_SYMBOL)],'-')
mutations$COSMIC_mutation_count = ifelse(mutations$Gene %in% cosmicMutation_tally$GENE_SYMBOL,cosmicMutation_tally$n_mutation[match(mutations$Gene,cosmicMutation_tally$GENE_SYMBOL)],'-')

mutations$mutation_group = as.character(mutations$somaticVar_type)
mutations$mutation_group[mutations$mutation_group == 'shared'] = 'shared across all samples'
mutations$mutation_group[mutations$mutation_group %in% c('unique_D:CN_region','unique_D:shearwater_Failed')] = 'unique to diagnostic major clone'
mutations$mutation_group[mutations$mutation_group == 'unique_TP1:shearwater_Failed'] = 'unique to refractory clone'
mutations$somaticVar_type = factor(mutations$somaticVar_type,c('shared','unique_D:CN_region','unique_D:shearwater_Failed','unique_TP1:shearwater_Failed'))
mutations$Chrom = factor(mutations$Chrom,paste0('chr',c(1:22,'X','Y')))
mutations = mutations[order(mutations$Pos,decreasing = F),]
mutations = mutations[order(mutations$Chrom),]
mutations = mutations[order(mutations$somaticVar_type),]

colnames(mutations) = c('Chrom','Pos','Ref','Alt','Gene','Type','Impact','AAchange','PD60301a_Diagnosis_VAF','PD61846a_TP1_VAF','somaticVar_type','COSMIC_gene_tier','COSMIC_mutation_count','mutation_group')
mutations = mutations[,c('Chrom','Pos','Ref','Alt','Gene','Type','Impact','AAchange','COSMIC_gene_tier','COSMIC_mutation_count','mutation_group','PD60301a_Diagnosis_VAF','PD61846a_TP1_VAF')]

write.csv(mutations,file.path(plotDir,'../TableS8_L038_RefractoryMLDS_somaticMutationCatalogue.csv'),row.names = F)

## Table S9: Differential Gene Expression analysis - relapse / refractory MLDS gene module ####

## Combination of:
# 1. DEGs between L076_R2 vs L076_D(bm)
# 2. DEGs between L038_TP1 vs L038_D (no TP53 loss)
# 3. DEGs between (L038_TP1 and L076_R2) vs responsive diagnostic ML-DS
outDir = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/7_bad.vs.good_MLDS'
setwd(outDir)
## get L076 comparison
l076_markers = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/7_bad.vs.good_MLDS/L076_relapse_FindMarkers.csv',row.names = 1)
l076_markers$direction = ifelse(l076_markers$avg_log2FC > 0 & l076_markers$comp == 'R2D_vs_D','D_down',
                           ifelse(l076_markers$avg_log2FC > 0 & l076_markers$comp == 'RD_vs_D','D_down',
                                  ifelse(l076_markers$avg_log2FC < 0 & l076_markers$comp == 'R2D_vs_D','D_up',
                                         ifelse(l076_markers$avg_log2FC < 0 & l076_markers$comp == 'RD_vs_D','D_up',
                                                'others'))))

min_l2FC = 0.5
max_cellfrac = 0.1
#l076_markers = l076_markers[abs(l076_markers$avg_log2FC) >= min_l2FC & pmax(l076_markers$pct.1,l076_markers$pct.2) >= max_cellfrac,]

l076_markers = l076_markers %>% group_by(geneSym,direction) %>% mutate(nComp = n_distinct(comp))
View(l076_markers)

l076_markers$ensID = geneMap$ensID[match(l076_markers$geneSym,geneMap$geneSym)]
l076_markers$chr = geneMap$chr[match(l076_markers$geneSym,geneMap$geneSym)]
l076_markers = l076_markers[,c("ensID","geneSym","chr","avg_log2FC","pct.1","pct.2","p_val","p_val_adj","direction",'comp','nComp')]
l076_markers$is_DEG = ifelse(abs(l076_markers$avg_log2FC) < min_l2FC,'low_l2FC',
                             ifelse(pmax(l076_markers$pct.1,l076_markers$pct.2) < max_cellfrac,'low_cellFrac_expressed','DEG'))

colnames(l076_markers) = c('ensembl_ID','gene_symbol','chromosome','log2FoldChange',
                           'fraction of cells expressed - group 1',
                           'fraction of cells expressed - group 2',
                           'pvalue','padj','direction','comparison','number_of_comparison','is_DEG')

l076_markers$comparison = paste0('L076_',l076_markers$comparison)

## get L038 comparison
l038_tp1.vs.clone1 = read.csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/7_bad.vs.good_MLDS/L038_tp1.vs.clone1.D_markers.csv',row.names = 1)
l038_tp1.vs.clone1$direction[l038_tp1.vs.clone1$direction == 'clone2_up'] = 'D_down'
l038_tp1.vs.clone1$direction[l038_tp1.vs.clone1$direction == 'clone2_down'] = 'D_up'
l038_tp1.vs.clone1$comp = 'L038_clone2.TP_vs_clone1.D'
l038_tp1.vs.clone1$nComp = 1
l038_tp1.vs.clone1 = l038_tp1.vs.clone1[,c("ensID","geneSym","chr","avg_log2FC","pct.1","pct.2","p_val","p_val_adj","direction",'comp','nComp')]
l038_tp1.vs.clone1$is_DEG = ifelse(abs(l038_tp1.vs.clone1$avg_log2FC) < min_l2FC,'low_l2FC',
                             ifelse(pmax(l038_tp1.vs.clone1$pct.1,l038_tp1.vs.clone1$pct.2) < max_cellfrac,'low_cellFrac_expressed','DEG'))

colnames(l038_tp1.vs.clone1) = c('ensembl_ID','gene_symbol','chromosome','log2FoldChange',
                           'fraction of cells expressed - group 1',
                           'fraction of cells expressed - group 2',
                           'pvalue','padj','direction','comparison','number_of_comparison','is_DEG')

bad_MLDS_markers = rbind(l076_markers,l038_tp1.vs.clone1)
bad_MLDS_markers = bad_MLDS_markers %>% group_by(gene_symbol,direction) %>% mutate('number_of_comparison' = n_distinct(comparison))
table(bad_MLDS_markers$is_DEG,bad_MLDS_markers$number_of_comparison,bad_MLDS_markers$direction)

write.csv(bad_MLDS_markers,file.path(plotDir,'../TableS9_L076.L038_relapse.refractory_vs_diagnostic_blasts_FindMarkers.csv'),row.names = F)
bad_MLDS_markers = read.csv(file.path(plotDir,'../TableS9_L076.L038_relapse.refractory_vs_diagnostic_blasts_FindMarkers.csv'))





##------------------------##
##      EGA requests    ####
##------------------------##
sc_dataset = read.csv(file.path(plotDir,'../TableS1_leukaemia_and_fLiver_datasets_overview.csv'))
sc_dataset=sc_dataset[sc_dataset$single.cell.or.bulk.RNA.seq != 'bulk RNA-seq',]
sc_dataset$supplierName = NA
sc_dataset$sequencescape_study = NA

## Loop through the sequencescape manifests and add supplier Sample name
for(study in c(6404,6918,7507)){
  d = file.path('~/lustre_mt22/Aneuploidy/sequenceScape_manifest',study)
  files = list.files(d,full.names = T)
  for(f in files){
    dd = read_excel(f,skip = 8)
    if(any(sc_dataset$Sample.ID %in% dd$`SANGER SAMPLE ID`)){
      sc_dataset$supplierName[sc_dataset$Sample.ID %in% dd$`SANGER SAMPLE ID`]  = dd$`SUPPLIER SAMPLE NAME`[match(sc_dataset$Sample.ID[sc_dataset$Sample.ID %in% dd$`SANGER SAMPLE ID`],dd$`SANGER SAMPLE ID`)]
      sc_dataset$sequencescape_study[sc_dataset$Sample.ID %in% dd$`SANGER SAMPLE ID`] = study
    }
    
  }
}

## Add the unique case
sc_dataset$supplierName[sc_dataset$Sample.ID == '4602STDY7920965'] = 'CongAML_diagnosis'
sc_dataset$sequencescape_study[sc_dataset$Sample.ID == '4602STDY7920965'] = '4602'

sc_dataset$supplierName[sc_dataset$Sample.ID == 'CG_SB_NB14406184'] = 'Y41-TUM-0-FT-1a'
sc_dataset$supplierName[sc_dataset$Sample.ID == 'CG_SB_NB14406185'] = 'Y41-TUM-0-FT-1b'
sc_dataset$sequencescape_study[sc_dataset$Sample.ID %in% c('CG_SB_NB14406184','CG_SB_NB14406185')] = '5676'

table(is.na(sc_dataset$supplierName))
table(is.na(sc_dataset$sequencescape_study))
table(sc_dataset$Dataset,is.na(sc_dataset$supplierName))
View(sc_dataset[is.na(sc_dataset$supplierName),])

sc_dataset = sc_dataset[,c("single.cell.or.bulk.RNA.seq","Sample.ID","supplierName","sequencescape_study")]
sc_dataset$Project = 'MLDS_sc'

##--- DNA data
dna_dataset = read.csv(file.path(plotDir,'../TableS1_leukaemia_and_fLiver_datasets_overview.csv'))
dna_dataset=dna_dataset[dna_dataset$single.cell.or.bulk.RNA.seq != 'bulk RNA-seq' & dna_dataset$WGS.ID != '-',]
dna_dataset$supplierName = dna_dataset$WGS.ID
dna_dataset$sequencescape_study = NA
## Loop through the sequencescape manifests and add supplier Sample name
for(study in c(6404,6918,7507,7153,6466)){
  d = file.path('~/lustre_mt22/Aneuploidy/sequenceScape_manifest',study)
  files = list.files(d,full.names = T)
  for(f in files){
    dd = read_excel(f,skip = 8)
    if(any(dna_dataset$WGS.ID %in% dd$`SUPPLIER SAMPLE NAME`)){
      dna_dataset$sequencescape_study[dna_dataset$WGS.ID %in% dd$`SUPPLIER SAMPLE NAME`] = study
    }
    
  }
}

dna_dataset$sequencescape_study[dna_dataset$supplierName == 'PD53663c'] = '6142'
table(is.na(dna_dataset$sequencescape_study))
View(dna_dataset[is.na(dna_dataset$sequencescape_study),])

dna_dataset = dna_dataset[,c("single.cell.or.bulk.RNA.seq","supplierName","sequencescape_study")]
dna_dataset$Project = 'MLDS_WGS'
dna_dataset = dna_dataset[!duplicated(dna_dataset),]
dna_dataset$Sample.ID = dna_dataset$supplierName
dna_dataset$single.cell.or.bulk.RNA.seq = 'WGS'

mlds_data_for_EGA = rbind(sc_dataset,dna_dataset[,colnames(sc_dataset)])
colnames(mlds_data_for_EGA) = c('Data_type','Sample_ID','Supplier_Sample_Name','sequencescape_study','Project')
mlds_data_for_EGA = mlds_data_for_EGA[,c('Project','Data_type','Sample_ID','Supplier_Sample_Name','sequencescape_study')]
dim(mlds_data_for_EGA)
write.csv(mlds_data_for_EGA,'~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/data_info_for_EGA.csv')



##------------------------------------------------------------------------------------------------##
##      SUPPLEMENTARY FIGURES   ####
##------------------------------------------------------------------------------------------------##

##----- Sup.Fig X
## Plot PURPLE CN profile for L076 and L038
patientID = 'L076'
if(patientID == 'L076'){
  samples_manifest = data.frame(PDID=c('PD62331c','PD62331a','PD64665a','PD66167a'),
                                projectid=rep(3484,4),
                                tissue=c('Blood_D','BM_D','BM_RD','BM_R2D'))
  
}else if(patientID == 'L038'){
  samples_manifest = data.frame(PDID=c('PD60301a','PD61846a'),
                                projectid=c(3030,3030),
                                timepoint=c('Diagnostic','TP1'))
  
}





figSuppXX_PURPLE_CNprofile = function(){
  
  for(s in unique(samples_manifest$PDID)){
    purple_fp=paste0('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_PURPLE/',patientID,'/',s,'/purple_output/without_GRIDDS/',s,'.purple.cnv.somatic.tsv')
    purple_CNV = read.delim(purple_fp,sep = '\t')
    ## calculate cumulative x-position
    purple_CNV$chromosome = factor(purple_CNV$chromosome,paste0('chr',c(1:22,'X','Y')))
    purple_CNV = purple_CNV[order(purple_CNV$start,decreasing = F),]
    purple_CNV = purple_CNV[order(purple_CNV$chromosome),]
    purple_CNV$x_start = NA
    purple_CNV$x_end = NA
    for(i in 1:nrow(purple_CNV)){
      if(purple_CNV$chromosome[i] == 'chr1'){
        purple_CNV$x_start[i] = purple_CNV$start[i]
        purple_CNV$x_end[i] = purple_CNV$end[i]
      }else{
        purple_CNV$x_start[i] = purple_CNV$x_end[i-1] +1
        purple_CNV$x_end[i] = purple_CNV$x_start[i] + (purple_CNV$end[i]-purple_CNV$start[i])
      }
    }
    
    purple_CNV$chr_odd = (purple_CNV$chromosome %in% paste0('chr',c(seq(1,22,2),'X')))
    purple_CNV$copyNumber[purple_CNV$chromosome == 'chr21'] = 3
    chromLabel = purple_CNV
    chromLabel = chromLabel %>% group_by(chromosome) %>% summarise(chromLen = max(end),
                                                                   chromStart_xpos = min(x_start)) %>% 
      mutate(pos = chromStart_xpos + 0.5*chromLen,
             chromosome = gsub('chr','',as.character(chromosome)))
    chromLabel$chromosome[chromLabel$chromosome %in% c('20','22')] = ''
    
    
    tmp = rbind(chromLabel,chromLabel[chromLabel$chromosome=='Y',])
    tmp$chromStart_xpos[tmp$chromosome == 'Y'][2] = max(purple_CNV$x_end)
    
    plotFun_PURPLE_CNprofile = function(noFrame=FALSE,noPlot=FALSE){
      p = ggplot(purple_CNV)+
        geom_rect(data=purple_CNV[purple_CNV$chr_odd,],aes(xmin=x_start,xmax=x_end,ymin=0,ymax=5),fill=grey(0.9))+
        geom_segment(data=tmp,aes(x=chromStart_xpos,xend=chromStart_xpos,y=0,yend=5),color=grey(0),lty=2,size=0.2)+
        geom_segment(aes(x=x_start,xend=x_end,y=round(copyNumber),yend=round(copyNumber)),size=1)+
        geom_segment(aes(x=x_start,xend=x_end,y=round(minorAlleleCopyNumber),yend=round(minorAlleleCopyNumber)),colour='red',size=1)+
        theme_classic()+
        scale_y_continuous(breaks = c(0,1,2,3,4,5),labels = c(0,1,2,3,4,5),limits = c(0,6))+
        geom_text(data = chromLabel,aes(x=pos,y=5.4,label = chromosome),size = 3)+
        theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
              strip.background = element_blank(),
              axis.text.x = element_blank(),
              axis.text = element_text(colour = 'black'),axis.ticks.x = element_blank(),
              axis.ticks.y=element_line(colour = 'black')) + 
        xlab('Genomic position') + ylab('')
      
      print(p)
    }
    saveFig(file.path(plotDir,paste0('Sup.Figxx_',patientID,'_',s,'_purpleCNprofile')),plotFun_PURPLE_CNprofile,rawData=purple_CNV,width = 8,height = 2.5,res = 500)
  }  
}





##-------- Plot scRNA-seq CN profile - L076
l076 = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/L076_sratObj.RDS')
mDat = l076@meta.data
mDat$group = ifelse(!mDat[['annot_aug24']] %in% c('Tumour','Tum_MK?','Tumour_WT','unsure_Tumour','MK'),'Normal',
                    ifelse(mDat[['annot_aug24']] == 'Tumour',paste0(mDat$timePoint,'_',mDat$tissue),mDat$annot_aug24))
df = plot_BAF_byCellClusters(mDat=mDat, cellID_column='cellID',group = 'group', normalGroups = 'Normal',
                             outDir='~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L076/L076/tumourDNA_onlyAnalysis/',
                             minRead = 3,
                             patientID = patientID,PDID = PDID,
                             phCnts_fp=NULL,dropImprinted=dropImprinted,segs=segs,
                             tgtChrs=c(1:22))

chromInfo = read.delim('/lustre/scratch125/casm/team274sb/mt22/generalResources/chrom_abspos_kb.txt',sep = '\t')
chromInfo$chr = paste0('chr',chromInfo$chrom)
chromInfo = chromInfo[chromInfo$arm == 'q',]

df$chr_end_abspos = 1e3*chromInfo$abspos_kb[match(df$chr,paste0('chr',chromInfo$chrom))]
df$chrom = as.numeric(gsub('chr','',df$chr))
df$chr_start_abspos = ifelse(df$chr == 'chr1',1,
                             1e3*chromInfo$abspos_kb[match(df$chrom-1,as.numeric(chromInfo$chrom))])
  
df$abs_pos = df$pos + df$chr_start_abspos
ccs = c('major_allele' = 'red',
        'minor_allele' = 'black',
        'uninformative' = grey(0.7))  
dd = df[df$chr %in% c('chr5','chr8','chr13','chr21') &
          df$totCount > 5 & df$clusterID %in% c('Diagnostic_Blood (n=1040)','Diagnostic_BM (n=2008)','D.Relapse_BM (n=693)','D.Relapse2_BM (n=12252)','Normal (n=9463)'),]
dd$clusterID  = factor(dd$clusterID, c('Diagnostic_Blood (n=1040)','Diagnostic_BM (n=2008)','D.Relapse_BM (n=693)','D.Relapse2_BM (n=12252)','Normal (n=9463)'))

plotFun_scRNAseq_CNprofile = function(noFrame=FALSE,noPlot=FALSE){
  p = ggplot(dd,aes(pos,altFreq))+
    geom_hline(yintercept = 0.5,lty=1,lwd=0.3,col='black')+
    #geom_segment(data=chromInfo[chromInfo$chr %in% dd$chr,],aes(x=1e3*abspos_kb,xend=1e3*abspos_kb,y=0,yend=1),color='black',lty=2,size=0.2)+
    geom_point(size=0.0001,aes(col=phasingAssign),alpha=0.5)+
    #scale_color_manual(values = ccs)+
    scale_color_manual(values = c('#FF4D00','black',grey(0.7)))+
    facet_grid(clusterID~chr,scales = 'free_x')+
    theme_classic()+
    scale_y_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1),limits = c(0,1))+
    #geom_text(data = chromLabel,aes(x=pos,y=5.4,label = chromosome),size = 3)+
    
    theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
          strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_text(colour = 'black'),axis.ticks.x = element_blank(),
          axis.ticks.y=element_line(colour = 'black'),
          legend.position = 'none') + 
    xlab('Genomic position') + ylab('Aggregated alternate-allele frequency at hetSNPs')
  
  print(p)
}
saveFig(file.path(plotDir,paste0('Sup.Figxx_',patientID,'_scRNAseq.CNprofile')),plotFun_scRNAseq_CNprofile,rawData=dd,width = 5,height = 4,res = 500,useDingbats = F)



