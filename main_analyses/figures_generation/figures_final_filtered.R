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

plotDir='~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots'
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}


##----------------------------------##
##    Set general color schemes   ####
##----------------------------------##

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
akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'
akLiv_mdat_fp = gsub('_0824.RDS','_0824_mdat.csv',akLiv_srat_fp)

fig1a_fLiverCohort = function(){
  library(ComplexHeatmap)
  
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
  write.csv(fLiver_dataset[fLiver_dataset$nCell >0,],'~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/TableS1_fLiver_dataset.csv',row.names = F)
  
  ## Write table of number of cells per channel ID
  projectMani = read_excel('~/lustre_mt22/projectManifest.xlsx',sheet = 'Aneuploidy_mani')
  projectMani = projectMani[!is.na(projectMani$sangerSampleID),]
  mdat$sangerSampleID = projectMani$sangerSampleID[match(mdat$orig.ident,projectMani$chanelID)]
  mdat$sangerSampleID[is.na(mdat$sangerSampleID) & grepl('^MY',mdat$orig.ident)] = mdat$orig.ident[is.na(mdat$sangerSampleID) & grepl('^MY',mdat$orig.ident)]
  fLiver_dataset = as.data.frame(table(mdat$sangerSampleID,mdat$sorting_strategy,mdat$assay,mdat$donorID,mdat$gestationalAge,mdat$Genotype,mdat$Sex))
  colnames(fLiver_dataset) = c('sangerSampleID','sortingStrategy','assay','donorID','gestationalAge','Genotype','Sex','nCell')
  write.csv(fLiver_dataset[fLiver_dataset$nCell >0,],'~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/TableS1_fLiver_dataset_nCellBy10XSample.csv',row.names = F)
  
  
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
  mdat = read.csv(mlds_mdat_fp,row.names = 1) 
  
  mdat$disease = factor(mdat$disease,rev(c('TAM','MLDS')))
  mdat$clinicalOutcome[mdat$donorID %in% c('L156','L038')] = 'Refractory'
  mdat$clinicalOutcome[mdat$donorID %in% c('L076','CC3')] = 'Relapse'
  mdat$clinicalOutcome[mdat$donorID %in% c('L019','L039','L040','L042','L091','L041','CC4','CC5','CC1','CC2','L178',
                                           'L075','L114','L182')] = 'Remission'
  mdat$clinicalOutcome[mdat$clinicalOutcome == '?'] = 'unknown'
  table(mdat$clinicalOutcome,mdat$donorID)
  
  mdat$donorID = factor(mdat$donorID,c('L041','L040','L042','L091','L039','L019','L076','L038','L075','CC2','CC1','L156','CC3','CC4','CC5','CC6','CC7','CC8',
                                       'L114','L178','L182'))
  
  
  mdat$timePoint = as.character(mdat$timePoint)
  mdat$timePoint[mdat$donorID == 'L156' & mdat$orig.ident %in% c('MY.200531.14784986','MY.200531.14784987')] = 'Recurrent'
  mdat$timePoint[mdat$timePoint == 'Diagnostic' & mdat$donorID == 'L076' & mdat$tissue == 'Blood'] = 'D (PBMC)'
  mdat$timePoint[mdat$timePoint == 'Diagnostic' & mdat$disease == 'TAM'] = 'D (PBMC)'
  mdat$timePoint[mdat$timePoint == 'Diagnostic'] = 'D (BM)'
  mdat$timePoint = factor(mdat$timePoint,(c('D (BM)','D (PBMC)','TP1','TP2','TP4','D.Relapse','D.Relapse2','Recurrent')))
  
  mdat$group = ifelse(mdat$broadLineage == 'Tumour','Tumour',ifelse(mdat$broadLineage %in% c('lowQual','Tumour_unsure'),'lowQual','normal'))
  nLeuk = as.data.frame(table(mdat$donorID,mdat$orig.ident,mdat$group))
  colnames(nLeuk) = c('donorID','channelID','group','Freq')
  nLeuk = nLeuk[nLeuk$Freq >0,]
  
  mdat = mdat[mdat$group != 'lowQual',]
  mlds_dataset = as.data.frame(table(mdat$donorID,mdat$orig.ident,mdat$timePoint,mdat$disease,mdat$tissue,mdat$clinicalOutcome,mdat$blastPerc,mdat$Genotype,mdat$sex,mdat$age_yrs))
  colnames(mlds_dataset) = c('donorID','channelID','timePoint','Disease','Tissue','Clinical outcome','blastPerc','Genotype','Sex','Age (year)','nCell')
  mlds_dataset = mlds_dataset[order(mlds_dataset$donorID),]
  
  
  
  
  mlds_dataset2 = mlds_dataset[mlds_dataset$nCell > 0,]
  mlds_dataset2$Sex = as.character(mlds_dataset2$Sex)
  mlds_dataset2$Sex[mlds_dataset2$Sex == 'F'] = 'female'
  mlds_dataset2$Sex[mlds_dataset2$Sex == 'M'] = 'male'
  mlds_dataset2$timePoint = gsub(' (.*)','',mlds_dataset2$timePoint)
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
  table(mlds_dataset2$donorID)
  
  
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$timePoint,decreasing = F),]
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$donorID),]
  mlds_dataset2 = mlds_dataset2[order(mlds_dataset2$Disease),]
  mlds_dataset2$blastPerc = as.character(mlds_dataset2$blastPerc)
  mlds_dataset2$blastPerc[mlds_dataset2$blastPerc != '?'] = as.numeric(as.character(mlds_dataset2$blastPerc[mlds_dataset2$blastPerc != '?']))
  
  write.csv(mlds_dataset2,file.path(plotDir,'../TableS1_MLDS_dataset.csv'),row.names = F)
  
  
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
  mtx = mtx[,c('CC4','CC5','L091','L042','L040','L039','L019','L041','L178','CC3','L076','L038')]
  
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
  mdat$annot[mdat$donorID == 'L010' & mdat$annot == 'T_cells'] = 'unsure'
  mdat$broadLineage[mdat$donorID == 'L010' & mdat$annot == 'unsure'] = 'unsure'
  
  
  ##---- Remove irrelevant samples -----##
  ## Remove P9_TP1
  mdat = mdat[mdat$timePoint != 'TP1',]
  
  mdat$disease[mdat$donorID == 'L010'] = 'pAML'
  table(mdat$donorID,mdat$disease)
  table(mdat$donorID[mdat$disease == 'pBALL'],mdat$timePoint[mdat$disease == 'pBALL'])
  table(mdat$donorID[mdat$disease == 'pAML'],mdat$timePoint[mdat$disease == 'pAML'])
  table(mdat$finalAnn_broad,mdat$timePoint)
  
  ##------ Redefine clinical data (outcome and timepoint and tissue) ----##
  table(mdat$clinicalOutcome,mdat$donorID,mdat$disease)
  mdat$disease = factor(mdat$disease,(c('MDS','AMKL','pAML','pBALL','LPD')))
  mdat$clinicalOutcome[mdat$clinicalOutcome == '?'] = 'Unknown'
  mdat$clinicalOutcome[grepl('Refractory|refractory',mdat$clinicalOutcome)] = 'Refractory'
  mdat$clinicalOutcome[grepl('emission',mdat$clinicalOutcome)] = 'Remission'
  
  table(mdat$timePoint,mdat$donorID)
  mdat$timePoint[mdat$timePoint == 'RelapseD0'] = 'Relapse Diagnostic'
  
  mdat$donorID = factor(mdat$donorID,c('L067','P9',
                                       'L010','L051','L058','L044','L100','L069',
                                       'L001','L080','L007','L014','L016','L050','L002','L068',
                                       'L027','L062'))
  mdat$age_yrs[mdat$age_yrs == '1Y,0m'] = '1Y'
  
  mdat$tissue[mdat$orig.ident %in% c('Leuk13760341','Leuk13760342','Leuk13760340','Leuk13645527')] = 'BM'
  mdat$tissue[mdat$orig.ident %in% c('NB14406184','NB14406185')] = 'Brain'
  
  ##---- Remove low quality cells, calculate Leukaemia vs Normal cells -----##
  mdat$group = ifelse(mdat$broadLineage == 'Tumour','Tumour',ifelse(mdat$broadLineage != 'unsure','Normal','unsure'))
  nLeuk = as.data.frame(table(mdat$donorID,mdat$orig.ident,mdat$group))
  colnames(nLeuk) = c('donorID','channelID','group','Freq')
  nLeuk = nLeuk[nLeuk$Freq >0,]
  
  
  
  
  ## Select for relevant columns only
  # Write supplementary table
  pLeuk_dataset = as.data.frame(table(mdat$donorID,mdat$orig.ident,mdat$timePoint,mdat$disease,mdat$tissue,mdat$clinicalOutcome,mdat$blastPerc,mdat$Genotype,mdat$sex,mdat$age_yrs))
  colnames(pLeuk_dataset) = c('donorID','channelID','timePoint','Disease','Tissue','Clinical outcome','Blast percentage','Genotype','Sex','Age (yrs)','nCell')
  pLeuk_dataset = pLeuk_dataset[pLeuk_dataset$nCell > 0,]
  pLeuk_dataset$Sex = as.character(pLeuk_dataset$Sex)
  pLeuk_dataset$Sex[pLeuk_dataset$Sex == 'F'] = 'female'
  pLeuk_dataset$Sex[pLeuk_dataset$Sex == 'M'] = 'male'
  pLeuk_dataset$timePoint = gsub(' (.*)','',pLeuk_dataset$timePoint)
  pLeuk_dataset$Tissue = as.character(pLeuk_dataset$Tissue)
  pLeuk_dataset$Tissue[pLeuk_dataset$Tissue == 'BM'] = 'Bone marrow'
  pLeuk_dataset$Tissue[pLeuk_dataset$Tissue == 'Blood'] = 'Peripheral blood'
  
  pLeuk_dataset$timePoint[pLeuk_dataset$timePoint == 'Relapse'] = 'Relapse day 0'
  pLeuk_dataset$timePoint[pLeuk_dataset$timePoint == 'Diagnostic'] = 'Day 0'
  #pLeuk_dataset$timePoint[pLeuk_dataset$timePoint == 'TP1'] = 'Timepoint 1 - post Chemo'
  pLeuk_dataset$timePoint = factor(pLeuk_dataset$timePoint,c('Day 0','Relapse day 0','Timepoint 1 - post Chemo'))
  pLeuk_dataset = pLeuk_dataset[order(pLeuk_dataset$timePoint),]
  pLeuk_dataset = pLeuk_dataset[order(pLeuk_dataset$donorID),]
  pLeuk_dataset = pLeuk_dataset[order(pLeuk_dataset$Disease),]
  
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
  
  
  write.csv(pLeuk_dataset,'~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/TableS1_otherLeuk_dataset.csv',row.names = F)
  
  
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
##    Figure 1 - Foetal liver AK dataset     ####
##---------------------------------------------##
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





fLiver_broadLin_cols_v4 = c('HSC_MPP' = '#EDBD74',"Ery/MegK/Mast lineage" = '#74acc2',"Myeloid lineage" = '#e07f4f',
                            'B lineage'='#702963','T/NK lineage'=grey(0.6),'Stromal'='#c3afcc','others'='#b0ee95')


fig2a_2nAK_Liv_cellTypeContribution = function(){ 
  mdat = read.csv(akLiv_mdat_fp)
  
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
    
    dd = data %>% filter(!broadLineage %in% c('Stromal','others','Trophoblast')) %>% group_by(donorID,broadLineage3,Genotype) %>% summarise(nCell = n_distinct(cellID)) %>% 
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
  
  
  ## Two-proportion z-test
  ##---- Ery/MK/Mast
  # For T21
  prop.test(x=c(sum(dd$nCell[dd$Genotype == 'T21' & dd$broadLineage3 =='Ery/MegK/Mast lineage']),
                sum(dd$nCell[dd$Genotype == 'Diploid' & dd$broadLineage3 =='Ery/MegK/Mast lineage'])),
            n=c(sum(unique(dd$totalCell[dd$Genotype == 'T21'])),
                sum(unique(dd$totalCell[dd$Genotype == 'Diploid']))), 
            p = NULL, alternative = "greater",
            correct = FALSE)
  
  # For Triploid
  prop.test(x=c(sum(dd$nCell[dd$Genotype == 'Triploid' & dd$broadLineage3 =='Ery/MegK/Mast lineage']),
                sum(dd$nCell[dd$Genotype == 'Diploid' & dd$broadLineage3 =='Ery/MegK/Mast lineage'])),
            n=c(sum(unique(dd$totalCell[dd$Genotype == 'Triploid'])),
                sum(unique(dd$totalCell[dd$Genotype == 'Diploid']))), 
            p = NULL, alternative = "greater",
            correct = FALSE)
  
  
  ##---- B lineage
  # For T21
  prop.test(x=c(sum(dd$nCell[dd$Genotype == 'T21' & dd$broadLineage3 =='B lineage']),
                sum(dd$nCell[dd$Genotype == 'Diploid' & dd$broadLineage3 =='B lineage'])),
            n=c(sum(unique(dd$totalCell[dd$Genotype == 'T21'])),
                sum(unique(dd$totalCell[dd$Genotype == 'Diploid']))), 
            p = NULL, alternative = "less",
            correct = FALSE)
  
  # For Triploid
  prop.test(x=c(sum(dd$nCell[dd$Genotype == 'Triploid' & dd$broadLineage3 =='B lineage']),
                sum(dd$nCell[dd$Genotype == 'Diploid' & dd$broadLineage3 =='B lineage'])),
            n=c(sum(unique(dd$totalCell[dd$Genotype == 'Triploid'])),
                sum(unique(dd$totalCell[dd$Genotype == 'Diploid']))), 
            p = NULL, alternative = "less",
            correct = FALSE)
  
}  


##---------------------------------------------##
##    Figure 3 - other leukaemia dataset     ####
##---------------------------------------------##

fig3_otherLeuk_UMAP = function(){
  
  if(file.exists(file.path(plotDir,'Fig3_otherLeuk_TumNorm_UMAP_rawData.tsv')) & skipIfExists){
    dd = read.delim(file.path(plotDir,'Fig3a_otherLeuk_TumNorm_UMAP_rawData.tsv'),sep = '\t',header = T)
    
  }else{
    # Import the seurat object
    mdat = read.csv(otherLeuk_mdat_fp)
    mdat$annot[mdat$donorID == 'L010' & mdat$annot == 'T_cells'] = 'unsure'
    mdat$broadLineage[mdat$donorID == 'L010' & mdat$annot == 'unsure'] = 'unsure'
    mdat$disease[mdat$donorID == 'L010'] = 'pAML'
    
    ##---- Remove irrelevant samples -----##
    ## Remove P9_TP1
    mdat = mdat[mdat$timePoint != 'TP1',]
    mdat$broadLineage[mdat$donorID == 'L010' & mdat$annot == 'unsure'] = 'unsure'
    
    
    dd = mdat
    dd = dd[!dd$broadLineage %in% c('doublets','unsure_others','Tumour?','lowQual'),]
    
    
  }
  # Prepare source data
  data = mdat[!mdat$broadLineage %in% c('doublets','unsure_others','Tumour?','Tumour_unsure','lowQual'),c('cellID','UMAP_1','UMAP_2','annot','orig.ident','donorID','timePoint','broadLineage','clinicalOutcome')]
  data$broadLineage = dd$broadLineage[match(data$cellID,dd$cellID)]
  
  ## Downsample the object for plotting
  dd = dd %>% group_by(broadLineage,donorID) %>% mutate(nCell = n())
  dd1 = dd[dd$nCell <= 5000,]
  dd2 = dd[dd$nCell > 5000,]
  set.seed(1234)
  dd2 = dd2 %>% group_by(broadLineage,donorID) %>% mutate(id = 1:n(),
                                                  selected = ifelse(id %in% sample(1:n(),5000),T,F))
  dd = rbind(dd1,dd2[dd2$selected==T,!colnames(dd2) %in% c('id','selected')])
  
  plotFun_byTumNorm_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    dd$tumNormCat = ifelse(dd$broadLineage == 'Tumour' & dd$timePoint == 'Diagnostic' & dd$disease == 'pBALL','Tumour_D0_pBALL',
                           ifelse(dd$broadLineage == 'Tumour' & dd$timePoint == 'Diagnostic' & dd$disease == 'pAML','Tumour_D0_pAML',
                                  ifelse(dd$broadLineage == 'Tumour' & dd$timePoint == 'Diagnostic' & dd$disease == 'AMKL','Tumour_D0_AMKL',
                                         ifelse(dd$broadLineage == 'Tumour' & dd$timePoint == 'Diagnostic' & dd$disease == 'MDS','Tumour_D0_MDS',
                                                ifelse(dd$broadLineage == 'Tumour' & dd$timePoint == 'Diagnostic' & dd$disease == 'LPD','Tumour_D0_LPD',
                                                       ifelse(dd$broadLineage == 'Tumour'  & dd$timePoint != 'Diagnostic','Tumour_RD0','Normal'))))))
                           
    ccs = c('Tumour_D0_pBALL' = '#756d4d','Tumour_D0_pAML'='#0078b3',
            'Tumour_D0_AMKL'='#E1BA00','Tumour_D0_MDS'='#B20000','Tumour_D0_LPD'='#d9a1a0',
            'Tumour_RD0'='#373024','Normal' = grey(0.75),
            'TAM_TP1' = '#7A4B82')
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','Other Leukaemia'),
         frame.plot=F)
    
    if(!noPlot){
      # points(dd$UMAP_1,dd$UMAP_2,
      #        col = colAlpha(ccs[dd$tumNormCat],alphas = 0.75),
      #        pch = 19,
      #        cex=0.03)
      points(dd$UMAP_1[!grepl('Tumour',dd$tumNormCat)],dd$UMAP_2[!grepl('Tumour',dd$tumNormCat)],
             col = colAlpha(ccs[dd$tumNormCat[!grepl('Tumour',dd$tumNormCat)]],alphas = 0.75),
             pch = 19,
             cex=0.03)
      
      points(dd$UMAP_1[grepl('Tumour',dd$tumNormCat)],dd$UMAP_2[grepl('Tumour',dd$tumNormCat)],
             col = colAlpha(ccs[dd$tumNormCat[grepl('Tumour',dd$tumNormCat)]],alphas = 0.8),
             pch = 19,
             cex=0.03)
      points(dd$UMAP_1[grepl('MDS',dd$tumNormCat)],dd$UMAP_2[grepl('MDS',dd$tumNormCat)],
             col = colAlpha(ccs[dd$tumNormCat[grepl('MDS',dd$tumNormCat)]],alphas = 0.8),
             pch = 19,
             cex=0.1)
      
    }
    
    if(noPlot & !noFrame){
      #Add coloured labels
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ broadLineage,data=dd,FUN=mean)
      mids$label = as.character(factor(mids$broadLineage,levels = c('HSC & prog.','Megakaryocytes','Erythroblasts',
                                                                    'Monocyte/Macrophage','Dendritic cells','Neutrophil',
                                                                    'B lineage','T/NK lineage','TAM','TAM_Relapse','Tumour','Tumour_postChemo','Tumour_Refractory','Tumour_Relapse','Tumour_Relapse2','Tumour_unsure')))
      
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$label,cex = 0.9,xpad = 1.9,ypad = 1.8,border = F,
                   col='black',bg = ccs[mids$broadLineage])
    }
  }
  
  saveFig(file.path(plotDir,'Fig3A_otherLeuk_TumNorm_UMAP'),plotFun_byTumNorm_UMAP,rawData=data,width = 3.3,height = 3,res = 500)
  
  
  plotFun_byPDID_UMAP = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    ccs_donorID = sample(col25,n_distinct(dd$donorID))
    names(ccs_donorID)= unique(dd$donorID)
    
    plot(dd$UMAP_1,dd$UMAP_2,
         las=1,
         type='n',
         #xlim=c(-13,17),
         #ylim=c(-13,17), 
         cex.main = 0.85,xaxt='n',yaxt='n',
         xlab='',ylab='',
         main=ifelse(noFrame,'','Other Leukaemia'),
         frame.plot=F)
    
    if(!noPlot){
      points(dd$UMAP_1,dd$UMAP_2,
             col = colAlpha(ccs_donorID[dd$donorID],alphas = 1),
             pch = 19,
             cex=0.03)
    }
    if(noPlot & !noFrame){
      #Add coloured labels
      mids = aggregate(cbind(UMAP_1,UMAP_2) ~ donorID,data=dd,FUN=mean)
      mids$label = mids$donorID
      
      boxed.labels(mids$UMAP_1,mids$UMAP_2,
                   labels=mids$label,cex = 0.9,xpad = 1.9,ypad = 1.8,border = T,
                   col='black',bg = ccs_donorID[mids$label])
    }
  }
  
  saveFig(file.path(plotDir,'SuppFigXX_otherLeuk_donorID_UMAP'),plotFun_byPDID_UMAP,rawData=data,width = 3.3,height = 3,res = 500)
  
}











figS3_GATA1s_moduleExpr_inFLiver = function(){
  
  gata1s_module = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/4_GATA1s_module/goodTAM/jul24/DESeq2_goshMEP.goodTAM.MLDS_topGeneModule.csv',row.names = 1)
  table(gata1s_module$tam_vs_mempT21_group)
  gata1s_module$group = ifelse(gata1s_module$tam_vs_mempT21_group == 'notDE',gata1s_module$direction,gata1s_module$tam_vs_mempT21_group)
  ## Select top genes only
  #gata1s_module = gata1s_module[abs(gata1s_module$cellFrac.diff) >= 20/100 & !grepl('^LINC\\d+|AC\\d+',gata1s_module$geneSym),]
  
  ##---- In fLiver
  akLiv_srat_fp = '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0824.RDS'
  akLiv_mdat_fp = gsub('_0824.RDS','_0824_mdat.csv',akLiv_srat_fp)
  
  fLiver = readRDS(akLiv_srat_fp)
  
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
  mtx = avgExpr[rownames(avgExpr) %in% topGenes$geneSym[topGenes$group %in% c('TAM.MLDS.down','TAM.MLDS.up')],grepl('diploid',colnames(avgExpr)) & !grepl('NPC|Endo|Fibroblast|Hepa|Mesen|doublets|Cholangiocytes|Mesothelial|Neuron|Trophoblast',colnames(avgExpr))]
  #mtx = avgExpr[rownames(avgExpr) %in% topGenes$geneSym[topGenes$group %in% c('TAM.MLDS.down','TAM.MLDS.up')],grepl('diploid',colnames(avgExpr)) & !grepl('NPC|doublets',colnames(avgExpr))]
  mtx=mtx[rowSums(mtx) > 0,]
  
  # Remove lowly expressed genes?
  mtx = mtx[apply(mtx,1,max) > 0.5,]
  
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
               row_names_gp = gpar(fontsize=6),column_names_gp = gpar(fontsize=6),column_title_gp = gpar(fontsize=10),column_title_rot = 90,km=10,
               #split = topGenes$group[match(rownames(mtx),topGenes$geneSym)],row_title_rot = 0,
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
  
  gata1s_module$fLiver_geneGroup = ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[[c('10')]],row_order(ht)[[c('2')]])],'Ery',
                                          ifelse(gata1s_module$geneSym %in% rownames(mtx)[row_order(ht)[['3']]],'MK',
                                                 ifelse(gata1s_module$geneSym %in% rownames(mtx)[row_order(ht)[['1']]],'Prog',
                                                        ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[[c('7')]],row_order(ht)[[c('8')]])],'NK.T',
                                                               ifelse(gata1s_module$geneSym %in% rownames(mtx)[row_order(ht)[['4']]],'Mast',
                                                                      ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['5']],row_order(ht)[['6']])],'Myeloid',
                                                                             ifelse(gata1s_module$geneSym %in% rownames(mtx)[c(row_order(ht)[['9']])],'B',
                                                                                    ifelse(!gata1s_module$geneSym %in% rownames(mtx),'-','others'))))))))
  
  table(gata1s_module$fLiver_geneGroup,gata1s_module$group)
  gata1s_module$fLiver_geneGroup = factor(gata1s_module$fLiver_geneGroup,c('Prog','MK','Mast','Ery','Myeloid','NK.T','B','others','-'))
  
  fLiver$broadLineage = factor(fLiver$broadLineage,levels = c( "HSC & prog.","Megakaryocytes","Mast.cell","Erythroblasts","Monocyte/Macrophage", "Dendritic cells","Myelocytes","T/NK lineage","B lineage","Stromal"))
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
  mtx = avgExpr[rownames(avgExpr) %in% topGenes$geneSym[topGenes$group %in% c('TAM.MLDS.down','TAM.MLDS.up')],grepl('diploid',colnames(avgExpr)) & !grepl('NPC|Endo|Fibroblast|Hepa|Mesen|doublets|Cholangiocytes|Mesothelial|Neuron|Trophoblast',colnames(avgExpr))]
  mtx=mtx[rowSums(mtx) > 0,]
  # Remove lowly expressed genes?
  mtx = mtx[apply(mtx,1,max) > 0.5,]
  
  mtx2 = mtx
  colnames(mtx2) = gsub('\\.diploid$','',colnames(mtx2))
  mtx2 = mtx2[,c("HSC_MPP","MEMP_MEP","CMP_GMP","LMPP_ELP","earlyMK","MK","Mast.cell",'EE','ME','LE','Mono.Mac','Dendritic.cells','myelocyte',"ILC.precursor",'NK_T',"early.B.cell",'B.cell')]
  mtx = mtx[,paste0(colnames(mtx2),'.diploid')]
  
  genes_toMark = gata1s_module$geneSym[gata1s_module$isTF==T | gata1s_module$isCSM==T | gata1s_module$isCosmic==T | gata1s_module$chr == 'chr21']
  #rownames(mtx2) = ifelse(rownames(mtx2) %in% genes_toMark,rownames(mtx2),'')
  
  
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
                                          avg_l2FC = gene_order$log2FoldChange[match(rownames(mtx2),gene_order$geneSym)]),
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
  
  saveFig(file.path(plotDir,paste0('FigSupp4_GATA1s_moduleExpr_in2nFLiver_heatmap_vertical')),plotFun_GATA1s_moduleExpr_inFLiver_vertical,rawData=data_forPlot,width = 6.5,height = 12,res = 500)  
  
  
  
  
  
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

