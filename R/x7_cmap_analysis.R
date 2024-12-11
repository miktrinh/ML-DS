## Process Connectivity Map results

signature = c('MLDS','GATA1')
signature = 'GATA1'

if(signature == 'GATA1'){
  outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/CLUE/GATA1_2410_full_TAMupOnly/'  
  gct_fp = file.path(outDir,'my_analysis.sig_queryl1k_tool.6704fa740abb5900137a886a/arfs/TAG/query_result.gct')
}else{
  outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/CLUE/MLDS_2410_v2_MLDSupOnly/'
  gct_fp = file.path(outDir,'my_analysis.sig_queryl1k_tool.6706561826e57d00131f063d/arfs/TAG/query_result.gct')
}


if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

setwd(outDir)

top_n = 20
cellline = read.csv('../cmap_cellline.csv')
cellline$group = ifelse(cellline$Culture.Properties == 'Suspension','Leukaemia',cellline$Type)
cellline = rbind(cellline,
                 data.frame(Cell.Line = c('MCF10A','HEK293','JURKAT','NEU','PHH','HUVEC','SKB','HPTEC'),
                            Culture.Properties = '?',
                            Description = '?',
                            Type = 'non-cancer',
                            group = 'non-cancer'))

cellline = rbind(cellline,
                 data.frame(Cell.Line = c('HUH7','NPC','MDAMB231','HELA','HCT116','U2OS','SKBR3','SW480','H1299','HS578T','NCIH508','SHSY5Y','NCIH2073','NCIH596','BT20'),
                            Culture.Properties = '?',
                            Description = '?',
                            Type = 'Cancer',
                            group = 'Cancer'))

cellline = rbind(cellline,
                 data.frame(Cell.Line = c('U937','HAP1','ASC','SKL','HL60','NOMO1'),
                            Culture.Properties = '?',
                            Description = '?',
                            Type = 'Leukaemia',
                            group = 'Leukaemia'))


cellline$group[cellline$group == 'Cancer'] = 'Solid cancer'










## Import cmap results
dat.gct <- read.delim(gct_fp,skip = 2)

## Extract compound only
cp = dat.gct[dat.gct$pert_type == 'trt_cp',]
cp$norm_cs = as.numeric(cp$norm_cs)

## extract top 20 similar and dissimilar entries
cp = cp[cp$moa != '-666',]
cp_toKeep = cp[cp$moa != '-666',]
cp_toKeep$norm_cs_ranking = order(cp_toKeep$norm_cs)
cp_toKeep_dissimilar = cp_toKeep[cp_toKeep$norm_cs_ranking %in% c(1:top_n),]
cp_toKeep_dissimilar$group = 'dissimilar'
cp_toKeep_similar = cp_toKeep[cp_toKeep$norm_cs_ranking %in% c(max(cp_toKeep$norm_cs_ranking):(max(cp_toKeep$norm_cs_ranking)-top_n)),]
cp_toKeep_similar$group = 'similar'
cp_toKeep = rbind(cp_toKeep_similar,cp_toKeep_dissimilar)

cp.sub = cp[cp$pert_iname %in% cp_toKeep$pert_iname,]
mtx = pivot_wider(cp.sub[,c('pert_iname','cell_iname','norm_cs')],id_cols = 'pert_iname',values_from = 'norm_cs', names_from = 'cell_iname')
mtx = column_to_rownames(mtx,'pert_iname')

for(i in 1:nrow(mtx)){
  for(j in 1:ncol(mtx)){
    mtx[i,j] = paste(unlist(mtx[i,j]),collapse = ':')
    if(!grepl(':',mtx[i,j])){
      mtx[i,j] = as.numeric(mtx[i,j])
    }else{
      # both -ve or +ve
      values = as.numeric(strsplit(mtx[i,j][[1]],split = ':')[[1]])
      if(all(values > 0)){
        mtx[i,j] = min(values)
      }else if(all(values < 0)){
        mtx[i,j] = max(values)
      }else{
        mtx[i,j] = NA
      }
    }
  }
}

#mtx2 = do.call(cbind,apply(mtx,2,unlist))
mtx2 = apply(mtx,2,unlist)
mtx2 = apply(mtx2, 2, as.numeric)
rownames(mtx2) = rownames(mtx)

library(ComplexHeatmap)
mtx2 = mtx2[match((cp_toKeep$pert_iname[order(cp_toKeep$norm_cs_ranking,decreasing = F)]),rownames(mtx2)),]
rownames(mtx2) = paste0(rownames(mtx2),' (',cp_toKeep$moa[match(rownames(mtx2),cp_toKeep$pert_iname)],')')
cp_toKeep$pert_name = paste0(cp_toKeep$pert_iname,' (',cp_toKeep$moa,')')
table(rownames(mtx2) %in% cp_toKeep$pert_name)
rownames(mtx2)[!rownames(mtx2) %in% cp_toKeep$pert_name]


## Make heatmap

colSplit = cellline$group[match(colnames(mtx2),cellline$Cell.Line)]
colSplit[is.na(colSplit)] = 'unknown'
colSplit = factor(colSplit,c('Leukaemia','Solid cancer','non-cancer','unknown'))
logitCols = c('#00008b','#adcae6','white','#ffd1d1','#d80000')

if(signature=='GATA1'){
  pdf('GATA1_heatmap.pdf',width = 12,height = 8.5)
  hm = Heatmap(mtx2,name='Similarity to GATA1 up-reg module',na_col = grey(1),show_column_dend = F,show_row_dend = F,
               col = circlize::colorRamp2(seq(-2,2,length.out=length(logitCols)),logitCols),
               row_names_gp = gpar(fontsize=9),cluster_rows = F,cluster_columns = F,
               column_names_gp = gpar(fontsize=9),border = T,
               split = cp_toKeep$group[match(rownames(mtx2),cp_toKeep$pert_name)],
               column_split = colSplit,  heatmap_legend_param = list(direction = "horizontal"))
  
  draw(hm,heatmap_legend_side = "bottom")
  dev.off()  
}else if (signature == 'MLDS'){
  pdf('MLDS_heatmap.pdf',width = 12,height = 7.5)
  hm = Heatmap(mtx2,name='Similarity to MLDS up-reg module',na_col = grey(1),show_column_dend = F,show_row_dend = F,
               col = circlize::colorRamp2(seq(-2,2,length.out=length(logitCols)),logitCols),
               row_names_gp = gpar(fontsize=9),cluster_rows = F,cluster_columns = F,
               column_names_gp = gpar(fontsize=9),border = T,
               split = cp_toKeep$group[match(rownames(mtx2),cp_toKeep$pert_name)],
               column_split = colSplit,  heatmap_legend_param = list(direction = "horizontal"))
  
  draw(hm,heatmap_legend_side = "bottom")
  dev.off() 
}








