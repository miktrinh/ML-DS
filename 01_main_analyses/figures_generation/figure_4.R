# Supplementary plots 
library(dplyr)


## Import L038 seurat objects
l038_d$group = ifelse(l038_d$seurat_clusters %in% c(19,23,12,0,33),'clone2_D',
                      ifelse(l038_d$annot_aug24 == 'Tumour','clone1_D','normal_D'))
l038_d$group = ifelse(l038_d$seurat_clusters %in% c(23,12,19,0,33) | l038_d$annot_aug24 == 'Tumour',paste0(as.character(l038_d$seurat_clusters),'_D'),'normal_D')

DimPlot(l038_d,group.by = 'group',label = T,repel = T,label.box = T)


l038_tp1 = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x_alleleIntegrator/L038/L038_TP1.RDS')
DimPlot(l038_tp1,group.by = 'seurat_clusters',label = T,repel = T,label.box = T)
l038_tp1$group = ifelse(l038_tp1$seurat_clusters %in% c(1,19,9,13,15,10,24,14,25) | l038_tp1$annot_aug24 == 'Tumour','clone2_TP1',
                        ifelse(l038_tp1$seurat_clusters %in% c(11),as.character(l038_tp1$seurat_clusters),'normal_TP1'))

l038_tp1$group = ifelse(l038_tp1$seurat_clusters %in% c(19,13,15,1,10,24,14,25) | l038_tp1$annot_aug24 == 'Tumour',paste0(as.character(l038_tp1$seurat_clusters),'_TP1'),'normal_TP1')

DimPlot(l038_tp1,group.by = 'group',label = T,repel = T,label.box = T)
DimPlot(l038_tp1,cells.highlight = l038_tp1$cellID[l038_tp1$annot_aug24 %in% c('Mono_CD14','T_cells')])
DimPlot(l038_d,cells.highlight = l038_d$cellID[l038_d$annot_aug24 %in% c('Mono_CD14','T_cells')])

mdat$group = l038_d$group[match(mdat$cellID,l038_d$cellID)]
mdat$group[mdat$cellID %in% l038_tp1$cellID] = l038_tp1$group[match(mdat$cellID[mdat$cellID %in% l038_tp1$cellID],l038_tp1$cellID)]



mdat = l038_srat@meta.data
passCellIDs = mdat$cellID
clusterIDs = setNames(mdat$cellID,ifelse(!is.na(mdat$group),as.character(mdat$group),
                                         ifelse(mdat$annot == 'Tumour',paste0(mdat$annot,':',mdat$timePoint),mdat$annot)))
normIDs = setNames(mdat$cellID[mdat$annot %in% c('Mono_CD14','T_cells')],mdat$annot[mdat$annot %in% c('Mono_CD14','T_cells')])


phCnts = readRDS(file.path(outDir,paste0(PDID,'.PD61846a_phCnts.RDS')))

phCnts$matCount_og = phCnts$matCount
phCnts$patCount_og = phCnts$patCount

table(phCnts$informative)
phCnts$matCount[is.na(phCnts$matCount)] = ifelse(phCnts$altCountTum[is.na(phCnts$matCount)]>phCnts$refCountTum[is.na(phCnts$matCount)],phCnts$altCount[is.na(phCnts$matCount)],phCnts$refCount[is.na(phCnts$matCount)])
phCnts$patCount[is.na(phCnts$patCount)] = ifelse(phCnts$altCountTum[is.na(phCnts$patCount)]<phCnts$refCountTum[is.na(phCnts$patCount)],phCnts$refCount[is.na(phCnts$patCount)],phCnts$refCount[is.na(phCnts$patCount)])



# phCnts$matCount = phCnts$matCount_og
# phCnts$patCount = phCnts$patCount_og
# phCnts$matCount[is.na(phCnts$matCount)] = phCnts$altCount[is.na(phCnts$matCount)]
# phCnts$patCount[is.na(phCnts$patCount)] = phCnts$refCount[is.na(phCnts$patCount)]


gCnts = filterCells(phCnts,dropUninformative = F,clusterIDs=clusterIDs,normIDs=normIDs,regionsToKeep = c('Exonic','Genic','Intronic'))
mdat = l038_srat@meta.data
gCnts.sub = gCnts[gCnts$cellID %in% mdat$cellID]

gCnts.segs = gCnts.sub
gCnts.segs$regionID = names(gCnts.segs)


#clusterIDs = setNames(mdat$cellID,mdat$group)
#clusterIDs = mdat$cellID[mdat$cellID %in% gCnts.segs$cellID]
l038_tumour$group_1 = mdat$group[match(l038_tumour$cellID,mdat$cellID)]
l038_tumour$group_2 = l038_tumour$group_1
l038_tumour$group_2[l038_tumour$group_1 == 'L038_TP1_others'] = 'TP1_wCNA'
l038_tumour$group_2[l038_tumour$seurat_clusters %in% c(10,11,8,17)] = 'TP1_wCNA'
l038_tumour$group_2[l038_tumour$group_1 == 'L038_D_others' & l038_tumour$seurat_clusters %in% c(0,15,4,9)] = 'D_wCNA'
l038_tumour$group_2[l038_tumour$group_1 == 'L038_D_others' & l038_tumour$seurat_clusters %in% c(1,12)] = 'D_noCNA'
l038_tumour$group_2[l038_tumour$seurat_clusters == 0 & l038_tumour$timePoint == 'Diagnostic'] = 'D_wCNA'



# l038_tumour$group_2[l038_tumour$group_2 == 'L038_D_others'] = paste0(l038_tumour$group_2[l038_tumour$group_2 == 'L038_D_others'],':',as.character(l038_tumour$seurat_clusters[l038_tumour$group_2 == 'L038_D_others']))
# l038_tumour$group_2[l038_tumour$group_1 == 'D_noCNA' & l038_tumour$seurat_clusters %in% c(0,4,9,15)] = paste0('D_noCNA:',as.character(l038_tumour$seurat_clusters[l038_tumour$group_1 == 'D_noCNA' & l038_tumour$seurat_clusters %in% c(0,4,9,15)]))
# l038_tumour$group_2[l038_tumour$group_2 == 'D_noCNA' & l038_tumour$seurat_clusters %in% c(14)] = paste0('D_noCNA:',as.character(l038_tumour$seurat_clusters[l038_tumour$group_2 == 'D_noCNA' & l038_tumour$seurat_clusters %in% c(14)]))
# l038_tumour$group_2[l038_tumour$group_1 == 'D_wCNA' & l038_tumour$seurat_clusters %in% c(7,14)] = paste0('D_wCNA:',as.character(l038_tumour$group_2[l038_tumour$group_1 == 'D_wCNA' & l038_tumour$seurat_clusters %in% c(7,14)]))
# 
# l038_tumour$group_2[l038_tumour$group_2 == 'TP1_noCNA' & l038_tumour$seurat_clusters %in% c(0,6)] = paste0('TP1_noCNA:',as.character(l038_tumour$seurat_clusters[l038_tumour$group_2 == 'TP1_noCNA' & l038_tumour$seurat_clusters %in% c(0,6)]))
# l038_tumour$group_2[l038_tumour$group_2 == 'TP1_wCNA' & l038_tumour$seurat_clusters %in% c(5,14,18,21)] = paste0('TP1_wCNA:',as.character(l038_tumour$seurat_clusters[l038_tumour$group_2 == 'TP1_wCNA' & l038_tumour$seurat_clusters %in% c(5,14,18,21)]))


DimPlot(l038_tumour,group.by = 'group_2',cols = col25)
DimPlot(l038_tumour,cells.highlight = l038_tumour$cellID[l038_tumour$group_2 == 'D_wCNA:D_wCNA'])

names(clusterIDs) = ifelse(clusterIDs %in% l038_tumour$cellID,as.character(l038_tumour$group_2[match(clusterIDs,l038_tumour$cellID)]),mdat$group[match(clusterIDs,mdat$cellID)])
gCnts.segs$clusterID = names(clusterIDs)[match(gCnts.segs$cellID,clusterIDs)]


clusterCnts = aggregateByLists(gCnts.segs, assays = c("altCount", "refCount"), gCnts.segs$clusterID, gCnts.segs$regionID)
clusterCnts.segs = clusterCnts
clusterCnts.segs$pos = as.numeric(gsub('.*:|_.*$','',clusterCnts.segs$regionID))
clusterCnts.segs$chr = gsub(':.*$','',clusterCnts.segs$regionID)
clusterCnts.segs$altFreq = clusterCnts.segs$altCount / (clusterCnts.segs$altCount + clusterCnts.segs$refCount)
clusterCnts.segs$totCount = clusterCnts.segs$altCount + clusterCnts.segs$refCount
clusterCnts.segs$cov = ifelse(clusterCnts.segs$totCount < 5,'<5',
                              ifelse(clusterCnts.segs$totCount < 10,'<10',
                                     ifelse(clusterCnts.segs$totCount < 20,'<20','>=20')))
clusterCnts.segs$cov = factor(clusterCnts.segs$cov,c('<5','<10','<20','>=20'))
clusterCnts.segs$timePoint = gsub('^.*:','',clusterCnts.segs$cellID)
clusterCnts.segs$celltype = gsub(':.*$','',clusterCnts.segs$cellID)

chr_toPlot = paste0('chr',c(5,17))

clusterCnts.segs$chr = factor(clusterCnts.segs$chr,paste0('chr',c(1:22,'X')))
#clusterCnts.segs$cellID = factor(clusterCnts.segs$cellID,c('D_noCNA','D_wCNA','TP1_noCNA','TP1_wCNA','L038_D_others','L038_TP1_others','Norm_others'))
#clusterCnts.segs$cellID = factor(clusterCnts.segs$cellID,c(unique(clusterCnts.segs$cellID[grepl('D$',clusterCnts.segs$cellID)]),unique(clusterCnts.segs$cellID[grepl('TP1$',clusterCnts.segs$cellID)])))
clusterCnts.segs$cellID = factor(clusterCnts.segs$cellID,c('11','clone1_D','clone2_D','clone2_TP1','normal_D','normal_TP1'))
clusterCnts.segs$altIsMum = gCnts$altIsMum[match(clusterCnts.segs$regionID,gCnts$regionID)]
clusterCnts.segs$altIsMum[is.na(clusterCnts.segs$altIsMum)] = 'uninformative'
plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/Plots'


for(chrom in chr_toPlot){
  # if(file.exists(file.path(plotDir,paste0('Fig5D_L038Tum_rawBAF_',chrom,'_rawData.tsv')))){
  #   df = read.delim(file.path(plotDir,paste0('Fig5D_L038Tum_rawBAF_',chrom,'_rawData.tsv')),sep = '\t')
  # }
  
  df = clusterCnts.segs[clusterCnts.segs$totCount > 10 & clusterCnts.segs$chr %in% chrom & clusterCnts.segs$cellID != '11',]
  #df = clusterCnts.segs[clusterCnts.segs$totCount > 10 & clusterCnts.segs$chr %in% chrom,]
  plotFun_rawBAF_L038Tum = function(noFrame=FALSE,noPlot=FALSE){
    par(mar=c(0.1,0.1,1,0.1))
    
    if(chrom == 'chr5'){
      x_max = max(df$pos)/1e6
    }else if (chrom == 'chr17'){
      x_max = 23
    }
    
    if(!noPlot & !noFrame){
      p = ggplot(df,aes(pos/1e6,altFreq))+
        #geom_point(aes(col=cov),size=0.3,alpha=1)+
        geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+
        geom_point(aes(col=altIsMum),size=0.05,alpha=0.5)+
        #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
        facet_grid(cellID  ~ chr,scales = 'free_x')+
        geom_hline(yintercept = 0.5,lty=1,lwd=0.3,col='black')+
        #scale_color_manual(values = c('black','red',grey(0.7)))+
        scale_color_manual(values = c('black','#FF4D00',grey(0.7)))+
        #scale_color_manual(values = c(col25[1:2],grey(0.85)))+
        #geom_vline(xintercept = 48800e3/1e6,col='red')+
        #geom_vline(xintercept = 35e6/1e6,col='red')+
        scale_y_continuous(breaks = c(0.0,0.5,1.0),labels = c(0,0.5,1))+
        theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(fill=F),
          axis.ticks = element_line(colour='black',linewidth = 0.3),
          axis.line = element_blank(),axis.text = element_text(size=8,colour = 'black'))+
        xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
      
      
      if(chrom == 'chr17'){
        p = p + 
          geom_rect(aes(xmin=7661779/1e6, xmax=7687550/1e6, ymin=-0.02, ymax=1.02), col='darkblue')
      }
      # png('tmp.png',width = 1000,height = 5000)
      # png('tmp.png',width = 500,height = 1500)
      print(p)
      #dev.off()    
    }
    
    if(!noPlot & noFrame){
      p = ggplot(df,aes(pos/1e6,altFreq))+
        #geom_point(aes(col=cov),size=0.3,alpha=1)+
        geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+
        geom_point(aes(col=altIsMum),size=0.05,alpha=0.5)+
        #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
        facet_grid(cellID  ~ chr,scales = 'free_x')+
        geom_hline(yintercept = 0.5,lty=1,lwd=0.3,col='black')+
        #scale_color_manual(values = c('black','red',grey(0.7)))+
        scale_color_manual(values = c('black','#FF4D00',grey(0.7)))+
        #scale_color_manual(values = c(col25[1:2],grey(0.85)))+
        #geom_vline(xintercept = 48800e3/1e6,col='red')+
        #geom_vline(xintercept = 35e6/1e6,col='red')+
        scale_y_continuous(breaks = c(0.0,0.5,1.0),labels = c(0,0.5,1))+
        theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
          strip.background = element_blank(),
          #panel.border = element_rect(fill=F),
          axis.ticks = element_line(colour='black',linewidth = 0.3),
          axis.line = element_blank(),axis.text = element_text(size=8,colour = 'black'))+
        xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
      
      if(chrom == 'chr17'){
        p = p + 
          geom_rect(aes(xmin=7661779/1e6, xmax=7687550/1e6, ymin=-0.02, ymax=1.02), col='darkblue')
      }
      
      print(p)
    }
    
    
    if(noPlot & !noFrame){
      df2 = do.call(rbind,lapply(levels(df$cellID),function(i){
        w = which(df$cellID == i)
        return(df[w[1:2],])
      }))
      if(chrom == 'chr17'){
        df2 = rbind(df2,
                    df[df$chr == chrom & df$pos == max(df$pos[df$chr == chrom]),][1,])
      }
      p = ggplot(df2[!is.na(df2$cellID),])+
        #geom_point(aes(col=cov),size=0.3,alpha=1)+
        geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+
        geom_point(aes(pos/1e6,altFreq,col=altIsMum),size=0.05,alpha=0.5)+
        #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
        facet_grid(cellID ~ chr,scales = 'free_x')+
        geom_hline(yintercept = 0.5,lty=1,lwd=0.3,col='black')+
        scale_color_manual(values = c('black','#FF4D00',grey(0.7)))+
        #scale_color_manual(values = c('black','red',grey(0.7)))+
        #scale_color_manual(values = c(col25[1:2],grey(0.85)))+
        #geom_vline(xintercept = 48800e3/1e6,col='red')+
        #geom_vline(xintercept = 35e6/1e6,col='red')+
        scale_y_continuous(breaks = c(0.0,0.5,1.0),labels = c(0,0.5,1))+
        theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(fill=F),
          axis.ticks = element_line(colour='black',linewidth = 0.3),
          axis.line = element_blank(),axis.text = element_text(size=8,colour = 'black'))+
        xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
      
      if(chrom == 'chr17'){
        p = p + 
          geom_rect(aes(xmin=7661779/1e6, xmax=7687550/1e6, ymin=-0.02, ymax=1.02), col='darkblue')
      }
      
      print(p)
    }
    
  }
# Figure 4E
# Import l038 srat object
l038_srat = readRDS(file.path("~/ML-DS/Results/06_MLDS_refractory_relapse/L038", "L038_sratObj.RDS"))


mlds_srat_fp <- "~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns.RDS"
mlds_mdat_fp <- "~/ML-DS/Results/03_MLDS_annotation/MLDS_clean_annotated_2408_noUnknowns_mdat_2508.csv"


##------ Generate Copy Number plots  ---------------
## Import MLDS metadata
mlds_mdat = read.csv(mlds_mdat_fp,row.names = 1)
mlds_mdat$annot <- mlds_mdat$annot_aug24_new

## Add alleleIntegrator results
l038_aiRes = read.csv(file.path('~/ML-DS/Results/06_MLDS_refractory_relapse/L038','L038_AI_res.csv'))
library(dplyr)

l038_aiRes <- l038_aiRes %>%
  mutate(
    AIres = AI_output,
    AIres = case_when(
      AIres %in% c("?", "Uncalled") ~ "Uninformative",
      AIres %in% c("normFrac") ~ "without CNA",
      AIres %in% c("abbFrac") ~ "with CNA",
      .default = AIres
    ),
    group = case_when(
      AIres == "without CNA" & timePoint == "Diagnostic" & donorID == "L038" & annot == "Tumour" ~ "D_noCNA",
      AIres == "with CNA"    & timePoint == "Diagnostic" & donorID == "L038" & annot == "Tumour" ~ "D_wCNA",
      AIres == "with CNA"    & timePoint == "TP1"        & donorID == "L038" & annot == "Tumour" ~ "TP1_wCNA",
      AIres == "without CNA" & timePoint == "TP1"        & donorID == "L038" & annot == "Tumour" ~ "TP1_noCNA",
      annot == "Tumour"      & donorID == "L038" ~ "L038_others",
      annot == "Tumour" ~ "Tum_others",
      .default = "Norm_others"
    ),
    group = case_when(
      group == "L038_others" & timePoint == "TP1" ~ "L038_TP1_others",
      group == "L038_others" & timePoint == "Diagnostic" ~ "L038_D_others",
      .default = group
    )
  )


## Plot raw aggregate BAF per cell for cells/clusters of interest ####
# when not droping uninformative SNP, the code struggle with "droping low coverage" snp
# this is because the code use phCnts$matCount + phCnts$patCount as total coverage. but Mat/Pat count is only available for informative snps....
# so, need to manually run this bit of the code to use phCnts$totCount = phCnts$altCount + phCnts$refCount
passCellIDs = mdat$cellID
clusterIDs = setNames(mdat$cellID,ifelse(!is.na(mdat$group),as.character(mdat$group),
                                         ifelse(mdat$annot == 'Tumour',paste0(mdat$annot,':',mdat$timePoint),mdat$annot)))
normIDs = setNames(mdat$cellID[mdat$annot %in% c('Mono_CD14','T_cells')],mdat$annot[mdat$annot %in% c('Mono_CD14','T_cells')])


phCnts$matCount_og = phCnts$matCount
phCnts$patCount_og = phCnts$patCount

table(phCnts$informative)
phCnts$matCount[is.na(phCnts$matCount)] = ifelse(phCnts$altCountTum[is.na(phCnts$matCount)]>phCnts$refCountTum[is.na(phCnts$matCount)],phCnts$altCount[is.na(phCnts$matCount)],phCnts$refCount[is.na(phCnts$matCount)])
phCnts$patCount[is.na(phCnts$patCount)] = ifelse(phCnts$altCountTum[is.na(phCnts$patCount)]<phCnts$refCountTum[is.na(phCnts$patCount)],phCnts$refCount[is.na(phCnts$patCount)],phCnts$refCount[is.na(phCnts$patCount)])



# phCnts$matCount = phCnts$matCount_og
# phCnts$patCount = phCnts$patCount_og
# phCnts$matCount[is.na(phCnts$matCount)] = phCnts$altCount[is.na(phCnts$matCount)]
# phCnts$patCount[is.na(phCnts$patCount)] = phCnts$refCount[is.na(phCnts$patCount)]


gCnts = filterCells(phCnts,dropUninformative = F,clusterIDs=clusterIDs,normIDs=normIDs,regionsToKeep = c('Exonic','Genic','Intronic'))

gCnts.sub = gCnts[gCnts$cellID %in% mdat$cellID]

gCnts.segs = gCnts.sub
gCnts.segs$regionID = names(gCnts.segs)


#clusterIDs = setNames(mdat$cellID,mdat$group)
#clusterIDs = mdat$cellID[mdat$cellID %in% gCnts.segs$cellID]
l038_tumour$group_1 = mdat$group[match(l038_tumour$cellID,mdat$cellID)]
l038_tumour$group_2 = l038_tumour$group_1
l038_tumour$group_2[l038_tumour$group_1 == 'L038_TP1_others'] = 'TP1_wCNA'
l038_tumour$group_2[l038_tumour$seurat_clusters %in% c(10,11,8,17)] = 'TP1_wCNA'
l038_tumour$group_2[l038_tumour$group_1 == 'L038_D_others' & l038_tumour$seurat_clusters %in% c(0,15,4,9)] = 'D_wCNA'
l038_tumour$group_2[l038_tumour$group_1 == 'L038_D_others' & l038_tumour$seurat_clusters %in% c(1,12)] = 'D_noCNA'
l038_tumour$group_2[l038_tumour$seurat_clusters == 0 & l038_tumour$timePoint == 'Diagnostic'] = 'D_wCNA'



# l038_tumour$group_2[l038_tumour$group_2 == 'L038_D_others'] = paste0(l038_tumour$group_2[l038_tumour$group_2 == 'L038_D_others'],':',as.character(l038_tumour$seurat_clusters[l038_tumour$group_2 == 'L038_D_others']))
# l038_tumour$group_2[l038_tumour$group_1 == 'D_noCNA' & l038_tumour$seurat_clusters %in% c(0,4,9,15)] = paste0('D_noCNA:',as.character(l038_tumour$seurat_clusters[l038_tumour$group_1 == 'D_noCNA' & l038_tumour$seurat_clusters %in% c(0,4,9,15)]))
# l038_tumour$group_2[l038_tumour$group_2 == 'D_noCNA' & l038_tumour$seurat_clusters %in% c(14)] = paste0('D_noCNA:',as.character(l038_tumour$seurat_clusters[l038_tumour$group_2 == 'D_noCNA' & l038_tumour$seurat_clusters %in% c(14)]))
# l038_tumour$group_2[l038_tumour$group_1 == 'D_wCNA' & l038_tumour$seurat_clusters %in% c(7,14)] = paste0('D_wCNA:',as.character(l038_tumour$group_2[l038_tumour$group_1 == 'D_wCNA' & l038_tumour$seurat_clusters %in% c(7,14)]))
#
# l038_tumour$group_2[l038_tumour$group_2 == 'TP1_noCNA' & l038_tumour$seurat_clusters %in% c(0,6)] = paste0('TP1_noCNA:',as.character(l038_tumour$seurat_clusters[l038_tumour$group_2 == 'TP1_noCNA' & l038_tumour$seurat_clusters %in% c(0,6)]))
# l038_tumour$group_2[l038_tumour$group_2 == 'TP1_wCNA' & l038_tumour$seurat_clusters %in% c(5,14,18,21)] = paste0('TP1_wCNA:',as.character(l038_tumour$seurat_clusters[l038_tumour$group_2 == 'TP1_wCNA' & l038_tumour$seurat_clusters %in% c(5,14,18,21)]))


DimPlot(l038_tumour,group.by = 'group_2',cols = col25)
DimPlot(l038_tumour,cells.highlight = l038_tumour$cellID[l038_tumour$group_2 == 'D_wCNA:D_wCNA'])

names(clusterIDs) = ifelse(clusterIDs %in% l038_tumour$cellID,as.character(l038_tumour$group_2[match(clusterIDs,l038_tumour$cellID)]),mdat$group[match(clusterIDs,mdat$cellID)])
gCnts.segs$clusterID = names(clusterIDs)[match(gCnts.segs$cellID,clusterIDs)]


clusterCnts = aggregateByLists(gCnts.segs, assays = c("altCount", "refCount"), gCnts.segs$clusterID, gCnts.segs$regionID)
clusterCnts.segs = clusterCnts
clusterCnts.segs$pos = as.numeric(gsub('.*:|_.*$','',clusterCnts.segs$regionID))
clusterCnts.segs$chr = gsub(':.*$','',clusterCnts.segs$regionID)
clusterCnts.segs$altFreq = clusterCnts.segs$altCount / (clusterCnts.segs$altCount + clusterCnts.segs$refCount)
clusterCnts.segs$totCount = clusterCnts.segs$altCount + clusterCnts.segs$refCount
clusterCnts.segs$cov = ifelse(clusterCnts.segs$totCount < 5,'<5',
                              ifelse(clusterCnts.segs$totCount < 10,'<10',
                                     ifelse(clusterCnts.segs$totCount < 20,'<20','>=20')))
clusterCnts.segs$cov = factor(clusterCnts.segs$cov,c('<5','<10','<20','>=20'))
clusterCnts.segs$timePoint = gsub('^.*:','',clusterCnts.segs$cellID)
clusterCnts.segs$celltype = gsub(':.*$','',clusterCnts.segs$cellID)

chr_toPlot = paste0('chr',c(5,17))

clusterCnts.segs$chr = factor(clusterCnts.segs$chr,paste0('chr',c(1:22,'X')))
#clusterCnts.segs$cellID = factor(clusterCnts.segs$cellID,c('D_noCNA','D_wCNA','TP1_noCNA','TP1_wCNA','L038_D_others','L038_TP1_others','Norm_others'))
#clusterCnts.segs$cellID = factor(clusterCnts.segs$cellID,c(unique(clusterCnts.segs$cellID[grepl('D$',clusterCnts.segs$cellID)]),unique(clusterCnts.segs$cellID[grepl('TP1$',clusterCnts.segs$cellID)])))
clusterCnts.segs$cellID = factor(clusterCnts.segs$cellID,c('11','clone1_D','clone2_D','clone2_TP1','normal_D','normal_TP1'))
clusterCnts.segs$altIsMum = gCnts$altIsMum[match(clusterCnts.segs$regionID,gCnts$regionID)]
clusterCnts.segs$altIsMum[is.na(clusterCnts.segs$altIsMum)] = 'uninformative'
plotDir = '~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/Plots'



#dev.off()


head(clusterCnts.segs)
View(clusterCnts.segs[clusterCnts.segs$totCount > 3 &
                        clusterCnts.segs$chr %in% chr_toPlot &
                        grepl('T:',clusterCnts.segs$cellID),] %>% group_by(chr,cellID) %>% summarise(n=n()))





plotFun_rawBAF_L038Tum = function(noFrame=FALSE,noPlot=FALSE){
  par(mar=c(0.1,0.1,1,0.1))
  
  if(!noPlot & !noFrame){
    print(base_plot(dd, size = 1, alpha = 1, facet_var = vars(clusterID),chr_cfg))
  }
  
  if(!noPlot & noFrame){
    p = ggplot(dd,aes(pos/1e6,altFreq))+
      geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+
      geom_point(aes(col=phasingAssign),size=0.05,alpha=0.5)+
      facet_grid(chr ~ clusterID,scales = 'free_x')+
      geom_hline(yintercept = 0.5,lty=1,lwd=0.3,col='black')+
      scale_color_manual(values = c('black','#FF4D00',grey(0.7)))+
      scale_y_continuous(breaks = c(0.0,0.5,1.0),labels = c(0,0.5,1))+
      theme_classic(base_size = 13)+
      theme(axis.title = element_text(colour='black'),
            strip.background = element_blank(),
            panel.border = element_rect(fill=F,colour = 'black'),
            axis.ticks = element_line(colour='black',linewidth = 0.3),
            axis.line = element_blank(),
            axis.text = element_text(size=8,colour = 'black'))+
      xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
    
    if(chrom == 'chr17'){
      p = p +
        geom_rect(aes(xmin=7661779/1e6, xmax=7687550/1e6, ymin=-0.02, ymax=1.02), col='darkblue')
    }
    
    print(p)
  }
  
  
  if(noPlot & !noFrame){
    df2 = do.call(rbind,lapply(levels(df$cellID),function(i){
      w = which(df$cellID == i)
      return(df[w[1:2],])
    }))
    if(chrom == 'chr17'){
      df2 = rbind(df2,
                  df[df$chr == chrom & df$pos == max(df$pos[df$chr == chrom]),][1,])
    }
    p = ggplot(df2[!is.na(df2$cellID),])+
      #geom_point(aes(col=cov),size=0.3,alpha=1)+
      geom_rect(aes(xmin=0, xmax=x_max, ymin=-0.02, ymax=1.02), fill=grey(0.9),alpha = 0.2)+
      geom_point(aes(pos/1e6,altFreq,col=altIsMum),size=0.05,alpha=0.5)+
      #scale_color_manual(values = c(grey(0.7),col25[c(3,4,5)]),name = 'Aggregated Coverage')+
      facet_grid(cellID ~ chr,scales = 'free_x')+
      geom_hline(yintercept = 0.5,lty=1,lwd=0.3,col='black')+
      scale_color_manual(values = c('black','#FF4D00',grey(0.7)))+
      #scale_color_manual(values = c('black','red',grey(0.7)))+
      #scale_color_manual(values = c(col25[1:2],grey(0.85)))+
      #geom_vline(xintercept = 48800e3/1e6,col='red')+
      #geom_vline(xintercept = 35e6/1e6,col='red')+
      scale_y_continuous(breaks = c(0.0,0.5,1.0),labels = c(0,0.5,1))+
      theme_classic(base_size = 13)+theme(#axis.text.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill=F),
        axis.ticks = element_line(colour='black',linewidth = 0.3),
        axis.line = element_blank(),axis.text = element_text(size=8,colour = 'black'))+
      xlab('Genomic position') + ylab('Aggregated Alt allele frequency')
    
    if(chrom == 'chr17'){
      p = p +
        geom_rect(aes(xmin=7661779/1e6, xmax=7687550/1e6, ymin=-0.02, ymax=1.02), col='darkblue')
    }
    
    print(p)
  }
  
}


saveFig(file.path(plotDir,paste0('Fig4D_L038Tum_rawBAF_',chrom)),plotFun_rawBAF_L038Tum,rawData=df,width = 3.7,height = 3.7,res = 500,useDingbats = F)
#saveFig(file.path(plotDir,paste0('test_',chrom)),plotFun_rawBAF_L038Tum,rawData=df,width = 4,height = 17,res = 500,useDingbats = F)
}
}
