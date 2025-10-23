#----------------------------##
#          Libraries       ####
#----------------------------##
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(ggbeeswarm)
library(Seurat)
library(RColorBrewer)
source('/lustre/scratch126/casm/team274sb/mt22/CN_methods/scripts/finalScripts/R/misc.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')
source('~/lustre_mt22/generalScripts/utils/sc_utils.R')

plotDir='~/lustre_mt22/Aneuploidy/manuscriptDraft_0724/Plots/'


##----  1. Import the projection results -----##
#results = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/6_MK_lineage/jan24/fLiver2n_allcells/3_projection_output/fLiver2n_MKlin_traj1200_240315_combinedBranches.REF_wholeTrans_PCA_40PCs_50binByPseudotime.noDS_projection_output.RDS')
#results = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/3_projection_output/fLiver2n_MKlin_traj1200_240406_combinedBranches.REF_wholeTrans_PCA_40PCs_50binByPseudotime.noDS_projection_output.RDS')
results = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/3_projection_output/fLiver2n_MKlin_traj1200_240412_combinedBranches.REF_wholeTrans_PCA_40PCs_50binByPseudotime.noDS_projection_output.RDS')
results = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/HVG/3_projection_output/fLiver2n_MKlin_traj1200_240412_combinedBranches.REF_HVG_PCA_40PCs_50binByPseudotime.noDS_projection_output.RDS')
results = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/jul24/wholeTranscriptome/3_projection_output/fLiver2n_MKlin_traj1200_240823_combinedBranches.REF_wholeTrans_PCA_20PCs_50binByPseudotime.noDS_projection_output.RDS')
processed_output = results[['processed_output']]
processed_output$branchAssigned = gsub('^.*/|:all','',processed_output$binName)

##----  2. Quantify the distribution of cells across pseudotime -----##
col_fun = colorRamp2(c(0,0.1, 0.25, 0.7), c('white',grey(0.7),grey(0.5),'black'))

data = processed_output
#data$group = data$celltype

## Define groups to include in heatmap:
# 1. Positive controls: fLiver:2n_test
# 2. Negative controls: fLiver:2n:B.cell and fAdr:diploid:SCPs
# 3. Cancer cells: MLDS and TAM (Diagnostic samples only, except for L038)
d = data[grepl('fLiver:2n_train|fLiver:2n_test|fLiver:2n:B.cell|fAdr|TAM:.*:Tumour$|MLDS:L.*:Diagnostic:Tumour$|MLDS:L076_BM:D.Relapse:Tumour$|MLDS:CC.*:Diagnostic:Tumour$|Tumour_CNA|L038:TP1:Tumour$',data$group),]
d = d[!grepl('unsure',d$group),]
d$disease = gsub(':.*$','',d$group)

# Group pseudotime bins into bigger bin - merging 3 adjacent bins into 1
d$pseudobin = as.numeric(gsub('/.*$','',d$binName))
d$pseudobin_group = floor(d$pseudobin / 3.1)
table(d$pseudobin,d$pseudobin_group)

## Quantify the ratio of cells in mid vs late time point in TAM vs MLDS ##
d$pseudobin_group2 = ifelse(d$pseudobin_group > 6,'late','early')
d$group2 = d$group
d$group2[d$group2 == 'MLDS:L038:TP1:Tumour_CNA']='MLDS:L038:TP1:Tumour'
d$group2[d$group2 == 'MLDS:L038:Diagnostic:Tumour_CNA']='MLDS:L038:Diagnostic:Tumour'

## add donorID info to fLiver positive control group
d$group2[grepl('fLiver',d$group2)] = paste0(d$group2[grepl('fLiver',d$group2)],':',d$donorID[grepl('fLiver',d$group2)])
dd = d[!grepl('fLiver_train',d$group2),] %>% group_by(group2,pseudobin_group2) %>% summarise(nCell = n()) %>% mutate(totalCell = sum(nCell),frac = nCell/totalCell)
dd$group3 = gsub(':.*$','',dd$group2)

##------- bar plot: Proportion of cells mapped to early vs late pseudobins
ggplot(dd,aes(x=reorder(group2,frac),frac,fill=pseudobin_group2))+
  geom_col()+
  facet_grid(.~group3,scales = 'free_x',space = 'free_x')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1,size=7))


## scatter plot: ratio of late:early proportion of cells
dd = pivot_wider(dd,id_cols = c('group2','group3','totalCell'),names_from = 'pseudobin_group2',values_from = 'nCell')
dd$donorID = d$donorID[match(dd$group2,d$group2)]
dd$early2late = dd$early / dd$late
dd$late2early = dd$late / dd$early

dd$celltype = ifelse(dd$group3 %in% c('fAdr','fLiver'),sapply(strsplit(dd$group2,split=':'),'[',3),dd$group3)
dd$celltype = factor(dd$celltype,c('HSC_MPP','MEMP_MEP','earlyMK','MK','EE','ME','LE','Mast.cell','B.cell','SCPs','TAM','MLDS'))
dd$group_toCol = ifelse(dd$group3 %in% c('TAM','MLDS'),dd$donorID,dd$group3)
dd$group_toCol[dd$group2 == 'MLDS:L038:TP1:Tumour'] = 'L038:TP1'
dd$group_toCol[dd$group2 == 'MLDS:L038:Diagnostic:Tumour'] = 'L038:D'
dd$group_toCol[dd$group2 == 'MLDS:L076_BM:Diagnostic:Tumour'] = 'L076_BM'
dd$group_toCol[dd$group2 == 'MLDS:L076:Diagnostic:Tumour'] = 'L076_Blood'
dd$group_toCol[dd$group2 == 'MLDS:L076_BM:D.Relapse:Tumour'] = 'L076_BM.R'

dd$group_toCol = factor(dd$group_toCol,c('fAdr','fLiver',unique(dd$group_toCol[!dd$group_toCol %in% c('fAdr','fLiver')])))
dd$facet_col = ifelse(dd$group3 %in% c('MLDS','TAM'),'TAM / ML-DS',gsub(':.*$','',dd$group2))
dd$clinicalOutcome = ifelse(dd$donorID %in% c('L038','L076','L156'),'bad','good')
donorID_cols = c('L019'=col25[1],
                 'L039'=col25[4],
                 'L040'=col25[3],
                 'L042'='#8DB5CE',
                 'L091'='#ab7149',
                 'L076_Blood'=col25[15],
                 'L076_BM'=col25[16],
                 'L038:D'=col25[2],
                 'L038:TP1'=col25[24],
                 'L075'='#ff891c',
                 'L156'='#532C8A',
                 'CC1'='#139992',
                 'CC2'='#C594BF'
                 )
plotFun_cellProjection_late2early = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(dd[!grepl('B.cell|fAdr',dd$group2),],aes(celltype,late2early,col=group_toCol,shape=clinicalOutcome,size=clinicalOutcome))+
    geom_hline(yintercept = 1,lty=2,linewidth=0.5,alpha=0.6)+
    geom_quasirandom(width = 0.4)+
    scale_size_manual(values = c(2.1,1.5))+
    scale_shape_manual(values = c(17,19))+
    facet_grid(.~facet_col,space = 'free_x',scales = 'free_x')+
    #scale_color_manual(values = c(grey(0.6),donorID_cols,col25))+
    scale_color_manual(values = c(grey(0.6),col25))+
    scale_y_log10()+
    theme_classic(base_size = 12)+
    theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1),
          panel.border = element_rect(fill=F,linewidth = 1),axis.line = element_blank(),
          strip.background = element_rect(linewidth = 1))+
    xlab('')
  print(p1)  
}
#saveFig(file.path(plotDir,'Fig3x_cellProjection_late2early_scatterPlot'),plotFun_cellProjection_late2early,rawData=dd,width = 5.2,height = 3.4,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'Fig3x_cellProjection_late2early_scatterPlot_v2'),plotFun_cellProjection_late2early,rawData=dd,width = 6.5,height = 3.4,res = 500)


plotFun_cellProjection_late2early_hor = function(noFrame=FALSE,noPlot=FALSE){
  dd$celltype = as.character(dd$celltype)
  dd$celltype[dd$celltype %in% c('MK','EE','Mast.cell','ME','LE')] = 'MK/Ery/Mast'
  dd$celltype = factor(dd$celltype,rev(c('HSC_MPP','MEMP_MEP','earlyMK','MK/Ery/Mast','B.cell','SCPs','TAM','MLDS')))
  dd$group_toCol = dd$group3
  p1 = ggplot(dd[!grepl('B.cell|fAdr',dd$group2),],aes(y=celltype,x=late2early,col=group_toCol
                                                       #shape=clinicalOutcome,size=clinicalOutcome
                                                       ))+
    geom_vline(xintercept = 1,lty=2,linewidth=0.5,alpha=0.6)+
    geom_quasirandom(width = 0.4,groupOnX = F)+
    #geom_jitter(height = 0.2)+
    scale_size_manual(values = c(2.1,1.5))+
    scale_shape_manual(values = c(17,19))+
    facet_grid(facet_col~.,space = 'free_y',scales = 'free_y')+
    scale_color_manual(values = c(grey(0.7),'brown','#3B87C7'))+
    scale_x_log10()+
    theme_classic(base_size = 12)+
    theme(axis.text.y = element_text(),
          panel.border = element_rect(fill=F,linewidth = 1),axis.line = element_blank(),
          strip.background = element_rect(linewidth = 1))+
    ylab('') + xlab('late / early cell fraction')
  print(p1)  
}
#saveFig(file.path(plotDir,'Fig3x_cellProjection_late2early_scatterPlot'),plotFun_cellProjection_late2early,rawData=dd,width = 5.2,height = 3.4,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'Fig4_cellProjection_late2early_scatterPlot_hor'),plotFun_cellProjection_late2early_hor,rawData=dd,width = 6,height = 3,res = 500)


##------- Correlation between min_ED and correlation per Celltype_group  ----##
ggplot(processed_output,aes(max_ED,correlation_score))+
  geom_point(size=0.3)+
  geom_hline(yintercept = 0.5)+
  facet_wrap(vars(group))+
  theme_classic()+
  theme(panel.border = element_rect(fill=F),axis.line = element_blank())

##------- Heatmap: of distribution of cells across pseudotime  ----##
maxScoreVariable = 'cor'
use_corr.weighted_distribution = F

d2 = d
d2$group[grepl('L038:Diagnostic',d2$group)] = 'MLDS:L038:Diagnostic:Tumour'
d2$group[grepl('L038:TP1',d2$group)] = 'MLDS:L038:TP1:Tumour'

if(maxScoreVariable == 'cor'){
  # Calculate unexplained fraction from the correlation
  d2$unexplained = 1-d2$max_correlation
  
  
  # for(b in c('EE','MK','Mast.cell')){
  #   d2 = d %>% filter(branch %in% c('prog',b)) %>% group_by(group,pseudobin_group) %>% summarise(nCell = n()) %>% 
  #     group_by(group) %>% mutate(totalCell = sum(nCell),frac=nCell/totalCell)
  #   d2 = d2[order(d2$pseudobin_group,decreasing = F),]
  # }
  
  
  if(use_corr.weighted_distribution){
    ## Weighted cell distribution
    d3 = d2 %>% group_by(group,pseudobin_group) %>% summarise(nCell_real = n(),
                                                              nCell_weighted = sum(max_correlation),
                                                              nCell_unexplained = sum(unexplained)) %>%
      group_by(group) %>% mutate(totalCell = sum(nCell_real),frac=nCell_weighted/totalCell,
                                 frac_unexplained = nCell_unexplained/totalCell,
                                 sum_unexplained = sum(frac_unexplained))
    unexplainedSignal = d3 %>% group_by(group) %>% summarise(unexplained = unique(sum_unexplained))
    
  }else{
    ## Un-weighted cell distribution
    # Calculate mean correlation across all cells in the group
    matchScore_perGroup = d2 %>% group_by(group) %>% summarise(meanScore = mean(max_correlation))
    
  }
}else if(maxScoreVariable == 'ed'){
  # Calculate mean ED across all cells in the group
  matchScore_perGroup = d2 %>% group_by(group) %>% summarise(meanScore = mean(max_ED))
}

d2 = d2 %>% group_by(group,pseudobin_group) %>% summarise(nCell_real = n()) %>% 
  group_by(group) %>% mutate(totalCell = sum(nCell_real),frac=nCell_real/totalCell)
d2 = d2[order(d2$pseudobin_group,decreasing = F),]


d2 = pivot_wider(d2,names_from = 'pseudobin_group',values_from = 'frac',id_cols = c('group')) %>% as.data.frame()
d2[is.na(d2)] = 0
rownames(d2) = d2$group
#d2 = cbind(d2,unexplainedSignal$unexplained[match(rownames(d2),unexplainedSignal$group)])
#colnames(d2)[12] = 'unExplained'


##--- Plotting the heatmap ----##

pR2Cols  = c('#fff5eb','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#a63603','#7f2704')
bCols_unExplained = circlize::colorRamp2(seq(0.2,0.5,length.out=length(pR2Cols)),pR2Cols)
bCols_meanCor = circlize::colorRamp2(seq(0.5,0.8,length.out=length(pR2Cols)),pR2Cols)
bCols_meanED = circlize::colorRamp2(seq(30,90,length.out=length(pR2Cols)),pR2Cols)
if(maxScoreVariable == 'cor'){
  bCols_matchScore = bCols_meanCor
}else{
  bCols_matchScore = bCols_meanED
}

col_fun = colorRamp2(c(0,0.05,0.1, 0.2, 0.3), c('white',grey(0.6),grey(0.45),grey(0.3),'black'))


## Plotting reference / control heatmap
mtx = as.matrix(d2[grepl('fLiver',d2$group),-1])
rowAnno = rowAnnotation(df = data.frame(dataset = gsub(':.*$','',rownames(mtx)),
                                        unexplained = unexplainedSignal$unexplained[match(rownames(mtx),unexplainedSignal$group)]),
                        col = list(dataset = c('TAM' = col25[1],'MLDS'=col25[2],'fLiver'=col25[3]),
                                   unexplained = bCols))
Heatmap((mtx),show_column_names = T,show_row_names = T,
        col = col_fun,
        show_row_dend = T,show_column_dend = F,cluster_rows = T,cluster_columns = F,
        right_annotation=rowAnno,km=4,column_names_rot = 0)


## Plotting TAM/MLDS
mtx = as.matrix(d2[!grepl('fLiver',d2$group),-1])


rowAnno = rowAnnotation(df = data.frame(dataset = gsub(':.*$','',rownames(mtx)),
                                        unexplained = unexplainedSignal$unexplained[match(rownames(mtx),unexplainedSignal$group)]),
                        col = list(dataset = c('TAM' = col25[1],'MLDS'=col25[2],'fLiver'=col25[3]),
                                   unexplained = bCols))
Heatmap((mtx),show_column_names = T,show_row_names = T,
        col = col_fun,
        show_row_dend = T,show_column_dend = F,cluster_rows = T,cluster_columns = F,
        right_annotation=rowAnno,km=4,column_names_rot = 0)


## Plotting reference + Cancer heatmap
mtx = as.matrix(d2[!grepl('fLiver:2n_train',d2$group),-1])
mtx = mtx[match(c(rownames(mtx)[grepl('HSC',rownames(mtx))],
                  rownames(mtx)[grepl('MEMP',rownames(mtx))],
                  rownames(mtx)[grepl('earlyMK',rownames(mtx))],
                  rownames(mtx)[grepl(':MK',rownames(mtx))],
                  rownames(mtx)[grepl('EE',rownames(mtx))],
                  #rownames(mtx)[grepl('ME$',rownames(mtx))],
                  #rownames(mtx)[grepl('LE',rownames(mtx))],
                  rownames(mtx)[grepl('Mast.cell',rownames(mtx))],
                  rownames(mtx)[grepl('MLDS',rownames(mtx))],
                  rownames(mtx)[grepl('TAM',rownames(mtx))],
                  rownames(mtx)[grepl('B\\.cell',rownames(mtx))],
                  #rownames(mtx)[grepl('NK_T',rownames(mtx))],
                  rownames(mtx)[grepl('fAdr',rownames(mtx))]
),rownames(mtx)),]

rowSplit_group = ifelse(grepl('B.cell|NK_T|fAdr',rownames(mtx)),'neg_ref',
                        ifelse(grepl('fLiver',rownames(mtx)),'pos_ref','cancer'))
rowSplit_group = factor(rowSplit_group,c('pos_ref','cancer','neg_ref'))
rowAnno = rowAnnotation(df = data.frame(dataset = gsub(':.*$','',rownames(mtx)),
                                        #unexplained = unexplainedSignal$unexplained[match(rownames(mtx),unexplainedSignal$group)],
                                        matchScore = matchScore_perGroup$meanScore[match(rownames(mtx),matchScore_perGroup$group)]),
                        annotation_legend_param = list(matchScore = list(direction = "horizontal")),
                        col = list(dataset = c('TAM' = col25[1],'MLDS'=col25[2],'fLiver'=col25[3],'fAdr'=col25[4]),
                                   #unexplained = bCols_unExplained,
                                   matchScore = bCols_matchScore
                                   ))
mtx2 = mtx
rownames(mtx2) = gsub('fLiver:2n_test:|MLDS:|TAM:|:D:Tumour|:Tumour|fAdr:diploid:|fLiver:2n:','',rownames(mtx2))
rownames(mtx2) = gsub('Diagnostic','D',rownames(mtx2))

plotFun_cellProjection_heatmap = function(noFrame=FALSE,noPlot=FALSE){
  hm = Heatmap((mtx2),show_column_names = F,show_row_names = T,
               col = col_fun,name = 'Frac. of cells',
               show_row_dend = T,show_column_dend = F,
               cluster_rows = F,cluster_columns = F,
               row_names_side = 'left',row_names_gp = gpar(fontsize=10),
               split=rowSplit_group,
               left_annotation=rowAnno,#km=4,
               heatmap_legend_param = list(direction = "horizontal"),
               column_names_rot = 0)
  draw(hm,heatmap_legend_side = "bottom")
}

saveFig(file.path(plotDir,'Sup.FigX_cellProjection_heatmap'),plotFun_cellProjection_heatmap,rawData=mtx2,width = 5.4,height = 7.8,res = 500,useDingbats = T)



##--- Plotting correlation distribution ----##
dd = data[data$group %in% rownames(mtx),]
d = dd %>% group_by(group) %>% mutate(nCell_perGroupID = n())
d$genoAlpha_group = ifelse(d$nCell_perGroupID > 500,'low',
                           ifelse(d$nCell_perGroupID < 100,'high','mid'))
#d$genoAlpha_group[grepl('MLDS_TAM:|MLDS_MLDS:',d$groupID)] = 'mid'
d$genoAlpha_group = factor(d$genoAlpha_group,c('low','mid','high'))
d$group = factor(d$group,rev(rownames(mtx)))
d$ref_cancer = ifelse(grepl('B.cell|NK_T|SCP|NPC',d$group),'neg_ref',
                      ifelse(grepl('fLiver',d$group),'pos_ref','cancer'))

# Downsample data to plot
d1 = d[d$nCell_perGroupID <= 1000,]
d2 = d[d$nCell_perGroupID > 1000,]
set.seed(1234)
d2 = d2 %>% group_by(group) %>% mutate(id = 1:n(),
                                       selected = ifelse(id %in% sample(1:n(),1000),T,F))
d = rbind(d1,d2[d2$selected == T,!colnames(d2) %in% c('id','selected')])
data_toPlot = d
groupID_order = levels(data_toPlot$group)

## Calculate boxplot quantiles
quant_df = data.frame()
for(i in 1:n_distinct(dd$group)){
  d = dd[dd$group == unique(dd$group)[i],]
  quant = quantile(d$max_correlation) %>% t() %>% as.data.frame()
  quant$group = unique(dd$group)[i]
  quant_df = rbind(quant_df,quant)
}
quant_df$ref_cancer = data_toPlot$ref_cancer[match(quant_df$group,data_toPlot$group)]

quant_df$groupID_level = NA
for(g in unique(quant_df$ref_cancer)){
  quant_df$groupID_level[quant_df$ref_cancer == g] = as.numeric(factor(quant_df$group[quant_df$ref_cancer == g],
                                                                       levels = groupID_order[groupID_order %in% quant_df$group[quant_df$ref_cancer == g]]))    
}



library(ggbeeswarm)
data_toPlot$ref_cancer = factor(data_toPlot$ref_cancer,(c('pos_ref','cancer','neg_ref')))
quant_df$ref_cancer = factor(quant_df$ref_cancer,(c('pos_ref','cancer','neg_ref')))
alphaValue = c(0.035,0.06,0.1)

plotFun_cellProjection_correlation_distribution = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(data_toPlot,aes(x=max_correlation,y=group))+
    geom_vline(xintercept = c(0.53,0.83),col='red')+
    geom_jitter(aes(alpha = genoAlpha_group),width = 0,size=0.01)+
    scale_alpha_manual(values = alphaValue)+
    geom_segment(data=quant_df,aes(y=groupID_level-0.3,x=`50%`,yend=groupID_level+0.3,xend=`50%`),lwd=0.5)+
    geom_segment(data=quant_df,aes(y=groupID_level,x=`25%`,yend=groupID_level,xend=`75%`),lwd=0.5)+
    facet_grid(ref_cancer~.,scales = 'free_y',space = 'free')+
    xlim(0,1)+
    scale_x_continuous(breaks = c(0.0,0.5,1))+
    #scale_x_log10()+
    xlab('')+ylab('MLDS module score')+
    theme_classic(base_size = 7)+
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size = 7),
          strip.text = element_text(size = 7),strip.background = element_blank(),
          panel.border = element_rect(fill=F),axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = 'bottom',legend.text = element_text(size=5)) 
  
  print(p1)
}

saveFig(file.path(plotDir,'Fig3x_cellProjection_corDist'),plotFun_cellProjection_correlation_distribution,rawData=mtx2,width = 3.3,height = 7,res = 500,useDingbats = T)





##--- Plotting raw projection along pseudotime ----##

## Example: fLiver_MK
d = data[data$group %in% c('fLiver:2n_test:MK','fLiver:2n_test:HSC_MPP','MLDS:L091:Diagnostic:Tumour'),]
# d = data[data$group %in% c('fLiver:2n_test:MK','fLiver:2n_test:HSC_MPP','MLDS:CC5:Diagnostic:Tumour'),]
# d = data[data$group %in% c('fLiver:2n_test:MK','fLiver:2n_test:HSC_MPP','TAM:CC8:Tumour'),]
d$main_refCT = factor(d$main_refCT,c('HSC_MPP','MEMP_MEP','earlyMK','MK','EE','Mast.cell'))
celltype_col = c('HSC_MPP'=col25[1],
                 'MEMP_MEP'=col25[4],
                 'earlyMK'="#2d8e91",
                 'MK'= col25[3],
                 'EE' = '#FBBE92',
                 'Mast.cell' = grey(0.6))
plotFun_cellProjection_example = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(d[d$main_refCT %in% c('MEMP_MEP'),],aes(avg_pseudotime,max_correlation,col = main_refCT))+
    geom_quasirandom(size=0.2,width = 0.02,alpha=0.3)+
    geom_quasirandom(data = d[!d$main_refCT %in% c('MEMP_MEP'),],size=0.3,width = 0.01,alpha=0.6)+
    facet_wrap(vars(group),nrow = 1)+
    scale_color_manual(values = celltype_col)+
    xlim(0,1)+ylim(0.3,1)+
    scale_y_continuous(breaks = c(0.4,1),position = "right")+
    scale_x_continuous(breaks = c(0,0.5,1))+
    theme_classic()+
    theme(axis.text.x = element_text(size=8),
          axis.text.y = element_text(size = 8),
          strip.text = element_text(size = 7),strip.background = element_blank(),
          panel.border = element_rect(fill=F),axis.line = element_blank(),
          legend.position = 'bottom',legend.text = element_text(size=8)) 
  
  print(p1)
}
saveFig(file.path(plotDir,'Fig3x_cellProjection_example'),plotFun_cellProjection_example,rawData=d,width = 7.2,height = 2.3,res = 500,useDingbats = T)


plotFun_cellProjection_example_Frame = function(noFrame=FALSE,noPlot=FALSE){
  p1 = ggplot(d,aes(avg_pseudotime,max_correlation,col = main_refCT))+
    #geom_quasirandom(size=0.2,width = 0.02,alpha=0.3)+
    #geom_quasirandom(data = d[!d$main_refCT %in% c('MEMP_MEP'),],size=0.3,width = 0.01,alpha=0.6)+
    facet_grid(.~group)+
    scale_color_manual(values = celltype_col)+
    xlim(0,1)+ylim(0.3,1)+
    scale_y_continuous(breaks = c(0.4,1),labels = c(0.4,1))+
    scale_x_continuous(breaks = c(0,0.5,1))+
    theme_classic()+
    theme(axis.text.x = element_text(size=9),
          axis.text.y = element_text(size = 9),
          strip.text = element_text(size = 7),strip.background = element_blank(),
          panel.border = element_rect(fill=F),axis.line = element_blank(),
          legend.position = 'bottom',legend.text = element_text(size=8)) 
  
  print(p1)
}

saveFig(file.path(plotDir,'Fig3x_cellProjection_example_Frame'),plotFun_cellProjection_example_Frame,rawData=d,width = 7.6,height = 1.7,res = 500,useDingbats = T)




##----    4. FeaturePlot of key markers in Cancer cells ---------####
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
processed_output$branchAssigned = gsub('.*/|:all','',processed_output$binName)
table(processed_output$donorID[grepl('Tumour',processed_output$group)],processed_output$branchAssigned[grepl('Tumour',processed_output$group)])

mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
mlds_mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405_mdat.csv')
table(mlds$annot == mlds_mdat$annot[match(mlds$cellID,mlds_mdat$cellID)])
mlds$broadLineage = mlds_mdat$broadLineage[match(mlds$cellID,mlds_mdat$cellID)]
mlds$timePoint = mlds_mdat$timePoint[match(mlds$cellID,mlds_mdat$cellID)]

##------------- L091 ---------------####
#   Highly similar to normal cells    ##
l091_allCell = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'L091' & !mlds$broadLineage %in% c('lowQual','Tumour_unsure')])
l091_allCell = standard_clustering(l091_allCell)
l091_allCell = RunUMAP(l091_allCell, dims=seq(75))
l091_allCell = FindClusters(l091_allCell, resolution = 2)
DimPlot(l091_allCell,group.by = 'seurat_clusters',label = T,label.box = T,repel = T) + NoLegend()
DimPlot(l091_allCell,group.by = 'GATA1_status',label = T,label.box = T,repel = T,cols=col25) + NoLegend()
DimPlot(l091,cells.highlight = l091_allCell$cellID[l091_allCell$seurat_clusters == 20])
DimPlot(mlds,cells.highlight = l091_allCell$cellID[l091_allCell$annot == 'Tumour'])
DimPlot(l091_allCell,cells.highlight = l091_allCell$cellID[l091_allCell$seurat_clusters %in% c(15,20) & l091_allCell$annot == 'Tumour'])
l091_allCell$annot[l091_allCell$seurat_clusters == 25] = 'unsure_EE'
l091_allCell$annot[l091_allCell$seurat_clusters %in% c(15,20) & l091_allCell$annot == 'Tumour'] = 'unsure_Tumour'
l091_allCell = subset(l091_allCell,subset = cellID %in% l091_allCell$cellID[!l091_allCell$annot %in% c('unsure_Tumour','unsure_EE')])

## Add back to the MLDS metadata table
mlds_mdat$annot_apr24 = as.character(mlds_mdat$annot)
mlds_mdat$annot_apr24[mlds_mdat$cellID %in% l091_allCell$cellID] = l091_allCell$annot[match(mlds_mdat$cellID[mlds_mdat$cellID %in% l091_allCell$cellID],l091_allCell$cellID)]



dd_allCells = cbind(l091_allCell@meta.data,l091_allCell@reductions$umap@cell.embeddings)
dd_allCells = dd_allCells[!dd_allCells$annot %in% c('unsure_Tumour','unsure_EE'),]

ccs = c('TAM_D' = '#d9a1a0','TAM_TP1' = '#9d6771','Tumour'='#463A2F','Tumour_TP1'='#CC4D00','Tumour_Relapse'='#E1BA00','Normal' = grey(0.8))
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
plotFun_donor_allCell_UMAP = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  dd = dd_allCells
  donorID = unique(dd$donorID)
  
  dd$tumNormCat = ifelse(dd$broadLineage == 'TAM' & dd$timePoint == 'Diagnostic','TAM_D',
                         ifelse(dd$broadLineage == 'TAM' & dd$timePoint != 'Diagnostic','TAM_TP1',
                                ifelse(dd$broadLineage %in% c('Tumour','Tumour_TP1','Tumour_Relapse'),dd$broadLineage,'Normal')))
  
  dd$broadLineage[dd$tumNormCat != 'Normal'] = dd$tumNormCat[dd$tumNormCat != 'Normal']
  
  ccs2 = c(mlds_linCol,ccs)
  
  p1 = ggplot(dd_allCells,aes(UMAP_1,UMAP_2,col=broadLineage))+
    geom_point(size=0.005)+
    xlab('') + ylab('')+
    scale_color_manual(values = ccs2)+
    ggtitle(donorID)+
    theme_classic(base_size = 4)+theme(panel.border = element_rect(fill=F,color='black'),
                                       axis.line = element_blank(),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank())
  print(p1)
}
saveFig(file.path(plotDir,'Fig3f_L091_allCells_UMAP'),plotFun_donor_allCell_UMAP,rawData = dd_allCells,width = 3,height = 2.25,res = 500,useDingbats = T)













dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('WT','GATA1s_WT')] = 'GATA1 wild type'
dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('Mut','GATA1s_mutant')] = 'GATA1s mutation'
dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('noGATA1expr','no_GATA1_expr')] = 'No GATA1 expression'
dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('Not informative','No coverage at mutation site','uninformative')] = 'Uninformative'

ccs = c("GATA1 wild type"='#005579',"GATA1s mutation"='#A92821',"Uninformative"=grey(0.65),"No GATA1 expression"=grey(0.8),"TBD"=colAlpha('#a17f45',0.4))

plotFun_GATA1s_UMAP = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  dd = dd_allCells
  #[dd$GATA1_status != 'GATA1s mutation',]
  p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=GATA1_status))+
    geom_point(size=0.005)+
    #geom_point(data = dd[dd$GATA1_status == 'GATA1s mutation',])+
    xlab('') + ylab('')+
    scale_color_manual(values = ccs)+
    ggtitle(donorID)+
    theme_classic(base_size = 4)+theme(panel.border = element_rect(fill=F,color='black'),
                                       axis.line = element_blank(),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank())
  print(p1)
  
}

saveFig(file.path(plotDir,'Fig3f_L091_allCells_GATA1status_UMAP'),plotFun_GATA1s_UMAP,rawData = dd_allCells,width = 3,height = 2.2,res = 500,useDingbats = T)





l091 = subset(l091_allCell,subset = cellID %in% l091_allCell$cellID[l091_allCell$annot == 'Tumour'])
l091 = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'L091' & mlds$annot == 'Tumour' & 
                                                      !mlds$cellID %in% mlds_mdat$cellID[mlds_mdat$annot_apr24 %in% c('unsure_Tumour','unsure_EE')]])

l091 = subset(l091,subset = cellID %in% l091$cellID[l091$seurat_clusters != 14])
l091 = standard_clustering(l091)
l091 = RunUMAP(l091, dims=1:20, min.dist = 0.8,n.neighbors = 80,seed.use = 23)
DimPlot(l091,group.by = 'seurat_clusters',label = T,label.box = T,repel = T) + NoLegend()
DimPlot(l091,group.by = 'timePoint',label = T,label.box = T,repel = T)
DimPlot(l091_allCell,cells.highlight = l091$cellID[l091$seurat_clusters==13])

FeaturePlot(l091,'nFeature_RNA')
l091$branch = gsub('^.*\\/|:all','',processed_output$binName[match(l091$cellID,processed_output$cellID)])
l091$projection_pseudotime = processed_output$avg_pseudotime[match(l091$cellID,processed_output$cellID)]
l091$projection_bin = ifelse(l091$projection_pseudotime < 0.4,'early',
                             ifelse(l091$projection_pseudotime >= 0.4 & l091$projection_pseudotime < 0.6,'mid1',
                                    ifelse(l091$projection_pseudotime >= 0.6 & l091$projection_pseudotime < 0.8,'mid2','late')))
l091$projection_bin = factor(l091$projection_bin,levels = c('early','mid1','mid2','late'))


DimPlot(l091,group.by = 'projection_bin',cols = c(grey(0.9),grey(0.7),grey(0.45),grey(0))) + ggtitle('L091 leukaemic cells')+theme(legend.position = 'bottom')
DimPlot(l091,group.by = 'branch',cols = col25) + ggtitle('L091 leukaemic cells')+theme(legend.position = 'bottom')
DotPlot(l091,group.by = 'projection_bin',
        features = c('CD34','CD38',
                     'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1', 'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
                     'PF4','PPBP',
                     'ALAS2','HBA1','BPGM',
                     'CSF2RB','HDC','TPSAB1','KIT','CD63','ENPP3','CPA3')) + RotatedAxis() + ggtitle('L091 Tumour')


DotPlot(REF.srat,group.by = 'finalAnn_broad',
        features = c('CD34','CD38','SPINK2',
                     'SERPINB1', 'GATA1',	'GATA2', 'TESPA1', 'CTNNBL1', 'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
                     'PF4','PPBP',
                     'ALAS2','HBA1','BPGM',
                     'CSF2RB','HDC','TPSAB1','KIT','CD63','ENPP3','CPA3')) + RotatedAxis() + ggtitle('L091 Tumour')


## Determine ratio of GATA1 and GATA2
donorID = 'L091'
mtx = t(as.matrix(mlds@assays$RNA@data[c('GATA1','GATA2','CPA3','PLEK'),mlds$cellID[mlds$annot == 'Tumour' & mlds$donorID == donorID]]))
mtx = as.data.frame(mtx)
mtx$GATA1_over_GATA2 = mtx$GATA1 / mtx$GATA2

dd = read.delim('~/lustre_mt22/Aneuploidy/manuscriptDraft_0424/Plots/Fig3f_L091_branch_UMAP_rawData.tsv',sep = '\t')
dd = read.delim('~/lustre_mt22/Aneuploidy/manuscriptDraft_0424/Plots/FigSuppXX_L075_branch_UMAP_rawData.tsv',sep = '\t')
dd = cbind(l091@meta.data,l091@reductions$umap@cell.embeddings)
dd = merge(dd,mtx,by.x='cellID',by.y=0,all=T)
dd$ratio = mtx$GATA1_over_GATA2[match(dd$cellID,rownames(mtx))]
dd$ratio[dd$GATA1 == 0 & dd$GATA2 == 0] = NA
dd$ratio[!(dd$GATA1 == 0 & dd$GATA2 == 0)] = dd$GATA1_over_GATA2[!(dd$GATA1 == 0 & dd$GATA2 == 0)]
dd$ratio[dd$GATA1 >0 & dd$GATA2 == 0] = 7
dd$ratio[dd$GATA1 == 0 & dd$GATA2 > 0] = 1/7
#dd$ratio[dd$ratio > 2] = 2
dd$diff = dd$GATA1 - dd$GATA2


ggplot(dd,aes(projection_pseudotime,PLEK))+
  geom_point()+
  geom_smooth()+
  facet_wrap(vars(branch),ncol = 1)

plotFun_donor_GATA1.to.GATA2_UMAP = function(noFrame=FALSE,noPlot=FALSE){
  donorID = unique(dd$donorID)
  p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=log2(ratio)))+
    geom_point(size=0.01)+
    xlab('') + ylab('')+
    scale_color_gradientn(colours = c('#006762','#f3f7c4','#8a1631'),values = c(0,0.5,1),na.value = 'grey')+
    ggtitle(donorID)+
    theme_classic(base_size = 4)+theme(panel.border = element_rect(fill=F,color='black'),
                                       axis.line = element_blank(),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank(),legend.key.width = unit(2,'mm'))
  print(p1)
}
saveFig(file.path(plotDir,'Fig4_L091_GATA1.to.GATA2_UMAP'),plotFun_donor_GATA1.to.GATA2_UMAP,rawData = dd,width = 2,height = 1.75,res = 500,useDingbats = F)


dd = as.data.frame(l091@reductions$umap@cell.embeddings)
dd$cellID = l091$cellID
dd$projection_bin = l091$projection_bin
dd$branch = l091$branch

pseudoBin_col = c('early'=grey(0.9),
                  'mid1'=grey(0.7),
                  'mid2'=grey(0.4),
                  'late'=grey(0.1))

plotFun_donor_cancerUMAP_pseudotime = function(noFrame=FALSE,noPlot=FALSE){
  donorID = unique(dd$donorID)
  p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=projection_bin))+
    geom_point(size=0.01)+
    xlab('') + ylab('')+
    scale_color_manual(values = pseudoBin_col)+
    ggtitle(donorID)+
    theme_classic(base_size = 4)+theme(panel.border = element_rect(fill=F,color='black'),
                          axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())
  print(p1)
}
saveFig(file.path(plotDir,'Fig4_L091_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = F)




celltype_col = c('HSC_MPP'=col25[1],
                 'MEMP_MEP'=col25[4],
                 'earlyMK'="#2d8e91",
                 'MK'= col25[3],
                 'EE' = '#FBBE92',
                 'Mast.cell' = grey(0.6))
plotFun_donor_cancerUMAP_branch = function(noFrame=FALSE,noPlot=FALSE){
  donorID = unique(dd$donorID)
  p1 = ggplot(dd,aes(UMAP_1,UMAP_2,col=branch))+
    geom_point(size=0.01)+
    xlab('') + ylab('')+
    scale_color_manual(values = celltype_col)+
    ggtitle(donorID)+
    theme_classic(base_size = 4)+theme(panel.border = element_rect(fill=F,color='black'),
                                       axis.line = element_blank(),
                                       axis.text = element_blank(),
                                       axis.ticks = element_blank())
  print(p1)
}
saveFig(file.path(plotDir,'Fig3f_L091_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)


s = l091
plotFun_markerExpr_donor = function(noFrame=FALSE,noPlot=FALSE){
  p = FeaturePlot(s,features = c('CD34','CD38','GATA2', # General early progenitors
                                    'SERPINB1', 'GATA1','CSF2RB','ITGA2B','FCER1A',# MEP
                                    'KLF1', 'ALAS2',# EE
                                 #'GYPA','EPCAM','APOC1','PF4',
                                 #'TPSB2','MS4A2','IL1RL1' # Mast cells
                                    'FLI1', 'PLEK',# MK
                                    'HDC','CD63','CPA3','TPSAB1' # Mast
  ),pt.size = 0.2,ncol = 5,combine = F) 
    
  
  if(!noPlot){
    for(i in 1:length(p)) {
      p[[i]] <- p[[i]]+ NoLegend() + theme(panel.border = element_rect(fill=F,color = 'black'),
                                           axis.text = element_blank(),
                                           axis.ticks = element_blank(),
                                           axis.line = element_blank(),
                                           title = element_text(size=6),
                                           plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
        scale_color_gradientn(colours = c(grey(0.85),'#759BD6','#2D4372'),limits=c(0,3))+
        xlab('') + ylab('')
    }
    
    print(cowplot::plot_grid(plotlist = p,nrow = 4))  
  }
  
  if(noPlot & !noFrame){
    print(p[[1]] + theme(legend.position = 'bottom')+
            scale_color_gradientn(colours = c(grey(0.85),'#759BD6','#2D4372'),limits=c(0,3)))
  }
  
}
  
saveFig(file.path(plotDir,'Fig3f_L091_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)

plotFun_markerExpr_donor_legend = function(noFrame=FALSE,noPlot=FALSE){
  p = FeaturePlot(s,features = c('CD34','CD38','GATA2', # General early progenitors
                                    'SERPINB1', 'GATA1','CSF2RB','ITGA2B','FCER1A',# MEP
                                    'KLF1', 'ALAS2',# EE
                                    'FLI1', 'PLEK',# MK
                                    'HDC','CD63','CPA3','TPSAB1' # Mast
  ),pt.size = 0.2,ncol = 5,combine = F)
    
  
  
  print(p[[1]] + theme(legend.position = 'bottom')+
          scale_color_gradientn(colours = c(grey(0.85),'#759BD6','#2D4372'),limits=c(0,3)))
}

saveFig(file.path(plotDir,'Fig3f_L091_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)










##--------- L038 ---------------#####
#   Highly similar to normal cells    ##
l038_allCells = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'L038' & !mlds$broadLineage %in% c('lowQual','Tumour_unsure')])
l038_allCells = standard_clustering(l038_allCells)
FeaturePlot(l038_allCells,c('nFeature_RNA','nCount_RNA','percent.mt'))
DimPlot(l038_allCells,group.by = 'seurat_clusters',label = T,label.box = T,repel = T) + NoLegend()
DimPlot(l038_allCells,group.by = 'GATA1_status',label = T,label.box = T,repel = T,cols=c(col25,pal34H)) + NoLegend()

DimPlot(mlds,cells.highlight = l038_allCells$cellID[l038_allCells$annot != 'Tumour' & l038_allCells$seurat_clusters == 9])
DimPlot(l091,cells.highlight = l091_allCell$cellID[l091_allCell$seurat_clusters == 20])

dd_allCells = cbind(l038_allCells@meta.data,l038_allCells@reductions$umap@cell.embeddings)
dd_allCells = dd_allCells[!dd_allCells$annot %in% c('unsure_Tumour','unsure_EE'),]

ccs = c('TAM_D' = '#d9a1a0','TAM_TP1' = '#9d6771','Tumour'='#463A2F','Tumour_TP1'='#CC4D00','Tumour_Relapse'='#E1BA00','Normal' = grey(0.8))
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

saveFig(file.path(plotDir,'Fig3f_L038_allCells_UMAP'),plotFun_donor_allCell_UMAP,rawData = dd_allCells,width = 3,height = 2.25,res = 500,useDingbats = F)




dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('WT','GATA1s_WT')] = 'GATA1 wild type'
dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('Mut','GATA1s_mutant')] = 'GATA1s mutation'
dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('noGATA1expr','no_GATA1_expr')] = 'No GATA1 expression'
dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('Not informative','No coverage at mutation site','uninformative')] = 'Uninformative'

ccs = c("GATA1 wild type"='#005579',"GATA1s mutation"='#A92821',"Uninformative"=grey(0.65),"No GATA1 expression"=grey(0.8),"TBD"=colAlpha('#a17f45',0.4))

saveFig(file.path(plotDir,'Fig3f_L038_allCells_GATA1status_UMAP'),plotFun_GATA1s_UMAP,rawData = dd_allCells,width = 3,height = 2.2,res = 500,useDingbats = F)


##----- Sub-cluster L038 Leukaemic blasts  ------##
l038 = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'L038' & mlds$annot == 'Tumour'])
l038 = standard_clustering(l038)

DimPlot(l038,group.by = 'timePoint',label = T,label.box = T,repel = T)
FeaturePlot(l091,'nFeature_RNA')
l038$branch = processed_output$branchAssigned[match(l038$cellID,processed_output$cellID)]
l038$projection_pseudotime = processed_output$avg_pseudotime[match(l038$cellID,processed_output$cellID)]
l038$projection_bin = ifelse(l038$projection_pseudotime < 0.375,'early',
                             ifelse(l038$projection_pseudotime >= 0.375 & l038$projection_pseudotime < 0.625,'mid',
                                    ifelse(l038$projection_pseudotime >= 0.625 & l038$projection_pseudotime < 0.875,'late1','late2')))
l038$projection_bin = factor(l038$projection_bin,levels = c('early','mid','late1','late2'))


DimPlot(l038,group.by = 'projection_bin',cols = c(grey(0.9),grey(0.7),grey(0.45),grey(0))) + ggtitle('l038 leukaemic cells')+theme(legend.position = 'bottom')
DimPlot(l038,group.by = 'branch',cols = col25) + ggtitle('l038 leukaemic cells')+theme(legend.position = 'bottom')
DotPlot(l038,group.by = 'projection_bin',
        features = c('CD34','CD38',
                     'SERPINB1', 'GATA1',	'GATA2', 'TESPA1',	'CTNNBL1', 'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
                     'PF4','PPBP',
                     'ALAS2','HBA1','BPGM',
                     'CSF2RB','HDC','TPSAB1','KIT','CD63','ENPP3','CPA3')) + RotatedAxis() + ggtitle('l038 Tumour')


DotPlot(REF.srat,group.by = 'finalAnn_broad',
        features = c('CD34','CD38','SPINK2',
                     'SERPINB1', 'GATA1',	'GATA2', 'TESPA1', 'CTNNBL1', 'FCER1A', 'ITGA2B', 'HBD','KLF1','PLEK', # MEP
                     'PF4','PPBP',
                     'ALAS2','HBA1','BPGM',
                     'CSF2RB','HDC','TPSAB1','KIT','CD63','ENPP3','CPA3')) + RotatedAxis() + ggtitle('L038 Tumour')

dd = as.data.frame(l038@reductions$umap@cell.embeddings)
dd$cellID = l038$cellID
dd$projection_bin = l038$projection_bin
dd$branch = l038$branch
dd$donorID = 'L038'


plotFun_donor_cancerUMAP_pseudotime()
plotFun_donor_cancerUMAP_branch()
s=l038
plotFun_markerExpr_donor()
plotFun_markerExpr_donor_legend

saveFig(file.path(plotDir,'FigSuppXX_L038_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_L038_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_L038_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_L038_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)




ccs = c("GATA1 wild type"='#005579',"GATA1s mutation"='#A92821',"Uninformative"=grey(0.65),"No GATA1 expression"=grey(0.8),"TBD"=colAlpha('#a17f45',0.4))


mlds_mdat$annot_apr24 = as.character(mlds_mdat$annot)
'CC5'
donor = 'CC8'
for(donor in unique(mlds$donorID)){
  s_allCells = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == donor & !mlds$broadLineage %in% c('lowQual','Tumour_unsure')])
  s_allCells = standard_clustering(s_allCells)
  
  FeaturePlot(s_allCells,c('nFeature_RNA','nCount_RNA','percent.mt','PF4'))
  DimPlot(s_allCells,group.by = 'seurat_clusters',label = T,label.box = T,repel = T) + NoLegend()
  DimPlot(s_allCells,group.by = 'GATA1_status',label = T,label.box = T,repel = T,cols=c(col25,pal34H)) + NoLegend()
  
  DimPlot(mlds,cells.highlight = s_allCells$cellID[s_allCells$seurat_clusters == 3])
  
  ## Add back to the MLDS metadata table
  if(donor %in% c('L039')){
    mlds_mdat$annot_apr24[mlds_mdat$cellID %in% s_allCells$cellID[s_allCells$annot == 'Tum_MK?' & s_allCells$seurat_clusters == 26]] = 'Tumour'
    s_allCells$annot[s_allCells$annot == 'Tum_MK?' & s_allCells$seurat_clusters == 26] = 'Tumour'
  }else if(donor %in% c('L019')){
    mlds_mdat$annot_apr24[mlds_mdat$cellID %in% s_allCells$cellID[s_allCells$annot == 'Tum_MK?']] = 'Tumour'
    s_allCells$annot[s_allCells$annot == 'Tum_MK?'] = 'Tumour'
  }else if(donor %in% c('L040')){
    s_allCells$annot[s_allCells$annot == 'MEP' & s_allCells$seurat_clusters == 23] = 'Tumour'
    mlds_mdat$annot_apr24[mlds_mdat$annot_apr24 == 'MEP' & mlds_mdat$cellID %in% s_allCells$cellID[s_allCells$seurat_clusters==23]] = 'Tumour'
  }
  
  dd_allCells = cbind(s_allCells@meta.data,s_allCells@reductions$umap@cell.embeddings)
  dd_allCells = dd_allCells[!dd_allCells$annot %in% c('unsure_Tumour','unsure_EE'),]
  dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('WT','GATA1s_WT')] = 'GATA1 wild type'
  dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('Mut','GATA1s_mutant')] = 'GATA1s mutation'
  dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('noGATA1expr','no_GATA1_expr')] = 'No GATA1 expression'
  dd_allCells$GATA1_status[dd_allCells$GATA1_status %in% c('Not informative','No coverage at mutation site','uninformative')] = 'Uninformative'
  
  # UMAP colored by broadLinage (normal + cancer)
  saveFig(file.path(plotDir,paste0('FigSupXX_',donor,'_allCells_UMAP')),plotFun_donor_allCell_UMAP,rawData = dd_allCells,width = 3,height = 2.25,res = 500,useDingbats = F)
  # UMAP colored by GATA1s status (normal + cancer)
  saveFig(file.path(plotDir,paste0('FigSupXX_',donor,'_allCells_GATA1status_UMAP')),plotFun_GATA1s_UMAP,rawData = dd_allCells,width = 3,height = 2.2,res = 500,useDingbats = F)
  
}



##--------- L075 ---------------#####
l075 = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'L075' & mlds$annot == 'Tumour'])
l075 = standard_clustering(l075)

DimPlot(l075,group.by = 'branch',label = T,label.box = T,repel = T)
FeaturePlot(l075,'nFeature_RNA')

l075$branch = processed_output$branchAssigned[match(l075$cellID,processed_output$cellID)]
l075$projection_pseudotime = processed_output$avg_pseudotime[match(l075$cellID,processed_output$cellID)]
l075$projection_bin = ifelse(l075$projection_pseudotime < 0.375,'early',
                             ifelse(l075$projection_pseudotime >= 0.375 & l075$projection_pseudotime < 0.625,'mid',
                                    ifelse(l075$projection_pseudotime >= 0.625 & l075$projection_pseudotime < 0.875,'late1','late2')))
l075$projection_bin = factor(l075$projection_bin,levels = c('early','mid','late1','late2'))

dd = as.data.frame(l075@reductions$umap@cell.embeddings)
dd$cellID = l075$cellID
dd$projection_bin = l075$projection_bin
dd$branch = l075$branch
dd$donorID = 'L075'

plotFun_donor_cancerUMAP_pseudotime()
plotFun_donor_cancerUMAP_branch()
s = l075
plotFun_markerExpr_donor()
plotFun_markerExpr_donor_legend()

saveFig(file.path(plotDir,'FigSuppXX_L075_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_L075_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_L075_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_L075_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)








##--------- CC1 ---------------#####
cc1 = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'CC1' & mlds$annot == 'Tumour'])
cc1 = standard_clustering(cc1)

DimPlot(cc1,group.by = 'GATA1_status',label = T,label.box = T,repel = T)
FeaturePlot(cc1,'nFeature_RNA')

cc1$branch = processed_output$branchAssigned[match(cc1$cellID,processed_output$cellID)]
cc1$projection_pseudotime = processed_output$avg_pseudotime[match(cc1$cellID,processed_output$cellID)]
cc1$projection_bin = ifelse(cc1$projection_pseudotime < 0.375,'early',
                             ifelse(cc1$projection_pseudotime >= 0.375 & cc1$projection_pseudotime < 0.625,'mid',
                                    ifelse(cc1$projection_pseudotime >= 0.625 & cc1$projection_pseudotime < 0.875,'late1','late2')))
cc1$projection_bin = factor(cc1$projection_bin,levels = c('early','mid','late1','late2'))

dd = as.data.frame(cc1@reductions$umap@cell.embeddings)
dd$cellID = cc1$cellID
dd$projection_bin = cc1$projection_bin
dd$branch = cc1$branch
dd$donorID = 'CC1'

plotFun_donor_cancerUMAP_pseudotime()
plotFun_donor_cancerUMAP_branch()
s = cc1
plotFun_markerExpr_donor()
plotFun_markerExpr_donor_legend()

saveFig(file.path(plotDir,'FigSuppXX_CC1_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC1_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC1_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC1_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)










##--------- CC2 ---------------#####
cc2 = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'CC2' & mlds$annot == 'Tumour'])
cc2 = standard_clustering(cc2)

DimPlot(cc2,group.by = 'seurat_clusters',label = T,label.box = T,repel = T)
FeaturePlot(cc1,'nFeature_RNA')

cc2$branch = processed_output$branchAssigned[match(cc2$cellID,processed_output$cellID)]
cc2$projection_pseudotime = processed_output$avg_pseudotime[match(cc2$cellID,processed_output$cellID)]
cc2$projection_bin = ifelse(cc2$projection_pseudotime < 0.375,'early',
                            ifelse(cc2$projection_pseudotime >= 0.375 & cc2$projection_pseudotime < 0.625,'mid',
                                   ifelse(cc2$projection_pseudotime >= 0.625 & cc2$projection_pseudotime < 0.875,'late1','late2')))
cc2$projection_bin = factor(cc2$projection_bin,levels = c('early','mid','late1','late2'))

dd = as.data.frame(cc2@reductions$umap@cell.embeddings)
dd$cellID = cc2$cellID
dd$projection_bin = cc2$projection_bin
dd$branch = cc2$branch
dd$donorID = 'CC2'

plotFun_donor_cancerUMAP_pseudotime()
plotFun_donor_cancerUMAP_branch()
s = cc2
plotFun_markerExpr_donor()
plotFun_markerExpr_donor_legend()

saveFig(file.path(plotDir,'FigSuppXX_CC2_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC2_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC2_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC2_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)



## quick cluster markers
library(SoupX)
qm = quickMarkers(cc2@assays$RNA@counts,cc2$seurat_clusters)






##--------- CC6 ---------------#####
cc6 = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'CC6' & mlds$annot == 'Tumour'])
cc6 = standard_clustering(cc6)

DimPlot(cc6,group.by = 'seurat_clusters',label = T,label.box = T,repel = T)
FeaturePlot(cc1,'nFeature_RNA')

cc6$branch = processed_output$branchAssigned[match(cc6$cellID,processed_output$cellID)]
cc6$projection_pseudotime = processed_output$avg_pseudotime[match(cc6$cellID,processed_output$cellID)]
cc6$projection_bin = ifelse(cc6$projection_pseudotime < 0.375,'early',
                            ifelse(cc6$projection_pseudotime >= 0.375 & cc6$projection_pseudotime < 0.625,'mid',
                                   ifelse(cc6$projection_pseudotime >= 0.625 & cc6$projection_pseudotime < 0.875,'late1','late2')))
cc6$projection_bin = factor(cc6$projection_bin,levels = c('early','mid','late1','late2'))

dd = as.data.frame(cc6@reductions$umap@cell.embeddings)
dd$cellID = cc6$cellID
dd$projection_bin = cc6$projection_bin
dd$branch = cc6$branch
dd$donorID = 'CC6'

plotFun_donor_cancerUMAP_pseudotime()
plotFun_donor_cancerUMAP_branch()
s = cc6
plotFun_markerExpr_donor()
plotFun_markerExpr_donor_legend()

saveFig(file.path(plotDir,'FigSuppXX_CC6_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC6_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC6_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC6_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)



## quick cluster markers
library(SoupX)
DimPlot(cc6,cells.highlight = cc6$cellID[cc6$seurat_clusters %in% c(10,1,7,2,6)])
cc6$group = ifelse(cc6$seurat_clusters %in% c(13,0,3,11,4,14),'g1',
                   ifelse(cc6$seurat_clusters %in% c(10,1,7,2,6),'g2',as.character(cc6$seurat_clusters)))
qm = quickMarkers(cc6@assays$RNA@counts,cc6$group)






##--------- CC7 ---------------#####
cc7 = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID == 'CC7' & mlds$annot == 'Tumour'])
cc7 = standard_clustering(cc7)

DimPlot(cc7,group.by = 'seurat_clusters',label = T,label.box = T,repel = T)
FeaturePlot(cc7,c('nFeature_RNA','percent.mt'))

cc7$branch = processed_output$branchAssigned[match(cc7$cellID,processed_output$cellID)]
cc7$projection_pseudotime = processed_output$avg_pseudotime[match(cc7$cellID,processed_output$cellID)]
cc7$projection_bin = ifelse(cc7$projection_pseudotime < 0.375,'early',
                            ifelse(cc7$projection_pseudotime >= 0.375 & cc7$projection_pseudotime < 0.625,'mid',
                                   ifelse(cc7$projection_pseudotime >= 0.625 & cc7$projection_pseudotime < 0.875,'late1','late2')))
cc7$projection_bin = factor(cc7$projection_bin,levels = c('early','mid','late1','late2'))

dd = as.data.frame(cc7@reductions$umap@cell.embeddings)
dd$cellID = cc7$cellID
dd$projection_bin = cc7$projection_bin
dd$branch = cc7$branch
dd$donorID = 'CC7'

plotFun_donor_cancerUMAP_pseudotime()
plotFun_donor_cancerUMAP_branch()
s = cc7
plotFun_markerExpr_donor()
plotFun_markerExpr_donor_legend()

saveFig(file.path(plotDir,'FigSuppXX_CC7_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC7_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC7_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_CC7_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)



##--------- TAM ---------------#####
tam = subset(mlds,subset = cellID %in% mlds$cellID[mlds$disease == 'TAM' & mlds$annot == 'Tumour'])
tam = standard_clustering(tam)

DimPlot(tam,group.by = 'donorID',label = T,label.box = T,repel = T)
FeaturePlot(tam,c('nFeature_RNA','percent.mt'))

tam$branch = processed_output$branchAssigned[match(tam$cellID,processed_output$cellID)]
tam$projection_pseudotime = processed_output$avg_pseudotime[match(tam$cellID,processed_output$cellID)]
tam$projection_bin = ifelse(tam$projection_pseudotime < 0.375,'early',
                            ifelse(tam$projection_pseudotime >= 0.375 & tam$projection_pseudotime < 0.625,'mid',
                                   ifelse(tam$projection_pseudotime >= 0.625 & tam$projection_pseudotime < 0.875,'late1','late2')))
tam$projection_bin = factor(tam$projection_bin,levels = c('early','mid','late1','late2'))

dd = as.data.frame(tam@reductions$umap@cell.embeddings)
dd$cellID = tam$cellID
dd$projection_bin = tam$projection_bin
dd$branch = tam$branch
dd$donorID = tam$donorID[match(dd$cellID,tam$cellID)]
dd$donorID = 'TAM'


plotFun_donor_cancerUMAP_pseudotime()
plotFun_donor_cancerUMAP_branch()
s = tam
plotFun_markerExpr_donor()
plotFun_markerExpr_donor_legend()


saveFig(file.path(plotDir,'FigSuppXX_TAM_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_TAM_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_TAM_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_TAM_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)







##--------- ML-DS ---------------#####
mlds.cancer = subset(mlds,subset = cellID %in% mlds$cellID[mlds$disease == 'MLDS' & mlds$annot == 'Tumour'])
mlds.cancer = standard_clustering(mlds.cancer)

DimPlot(mlds.cancer,group.by = 'donorID',label = T,label.box = T,repel = T)
FeaturePlot(mlds.cancer,c('nFeature_RNA','percent.mt'))

mlds.cancer$branch = processed_output$branchAssigned[match(mlds.cancer$cellID,processed_output$cellID)]
mlds.cancer$projection_pseudotime = processed_output$avg_pseudotime[match(mlds.cancer$cellID,processed_output$cellID)]
mlds.cancer$projection_bin = ifelse(mlds.cancer$projection_pseudotime < 0.375,'early',
                            ifelse(mlds.cancer$projection_pseudotime >= 0.375 & mlds.cancer$projection_pseudotime < 0.625,'mid',
                                   ifelse(mlds.cancer$projection_pseudotime >= 0.625 & mlds.cancer$projection_pseudotime < 0.875,'late1','late2')))
mlds.cancer$projection_bin = factor(mlds.cancer$projection_bin,levels = c('early','mid','late1','late2'))

dd = as.data.frame(mlds.cancer@reductions$umap@cell.embeddings)
dd$cellID = mlds.cancer$cellID
dd$projection_bin = mlds.cancer$projection_bin
dd$branch = mlds.cancer$branch
dd$donorID = mlds.cancer$donorID[match(dd$cellID,mlds.cancer$cellID)]
dd$donorID = 'MLDS cancer'


plotFun_donor_cancerUMAP_pseudotime()
plotFun_donor_cancerUMAP_branch()
s = mlds.cancer
plotFun_markerExpr_donor()
plotFun_markerExpr_donor_legend()


saveFig(file.path(plotDir,'FigSuppXX_MLDScancer_pseudotime_UMAP'),plotFun_donor_cancerUMAP_pseudotime,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_MLDScancer_branch_UMAP'),plotFun_donor_cancerUMAP_branch,rawData = dd,width = 2,height = 1.55,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_MLDScancer_featurePlot'),plotFun_markerExpr_donor,width = 4.4,height = 5,res = 500,useDingbats = T)
saveFig(file.path(plotDir,'FigSuppXX_MLDScancer_featurePlot_legend'),plotFun_markerExpr_donor_legend,width = 2,height = 2,res = 500,useDingbats = T)

