# Extract GATA1s mutation from scRNAseq data of MLDS samples
# Aim: Genotype individual cells to help with calling leukaemic cells

# Step 1: Identify reads covering the region of GATA1s mutation - performed in bash (~/lustre_mt22/Aneuploidy/scripts/GATA1_mut_reads_v2.sh)
# Step 2: Tally up these reads to gain evidence for / against the presence of GATA1s mutation at single-cell level
# Step 3: Generate an MLDS object for 2023 (excluding 1 TAM from Gosh + 2 TAM from Hennings + new samples to be analysed after Xmas)

# Main output saved at line 249

setwd("~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/")

#############
# Libraries #
#############
# Load libraries
library(tidyverse)
library(readxl)
library(Seurat)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")



##---------------------------------##
##    1. Import MLDS object      ####
##---------------------------------##

#mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/may23/MLDS_clean_LRwCTannotated_may23.RDS')
# mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jul23/MLDS_clean_annotated_noUnknowns_jul23.RDS')
# mlds$timePoint[mlds$orig.ident == 'Leuk13645528'] = 'postChemo TP1'
# mlds$timePoint[mlds$orig.ident %in% c('Leuk13645529','Leuk13645524')] = 'postChemo TP4'

#mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/sept23/MLDS_clean_annotated_noUnknowns_sept23_tmp.RDS')
#mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/nov23/MLDS/MLDS_clean_annotated_231127.RDS')
#mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jan24/MLDS_clean_annotated_noUnknowns_jan24.RDS')
#mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/mar24/MLDS/MLDS_clean_noMTcells_annotated_240405.RDS')
mlds = readRDS('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_noMTCells.RDS')
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_mdat.csv')
mlds$annot_aug24 = mdat$annot_aug24[match(mlds$cellID,mdat$cellID)]


##-------------------------------------------##
##    2. Process GATA1 reads analysis      ####
##-------------------------------------------##

out = mdat[,c('cellID','donorID','annot_aug24')]
out$GATA1s_cov = '?'   # Number of reads mapped to GATA1 gene - which overlaps with the mutation site? (somewhere on exon 2)
out$n_GATA1s_mutant = '?'
out$n_GATA1s_wt = '?'


# tumCells_annot = data.frame()
# mlds$GATA1s_cov = '?'
# mlds$n_GATA1reads = '?'
# 
# unsure.GATAstatus.all = data.frame()
unique(mdat$donorID)
for(donorID in unique(mdat$donorID)){
  print(donorID)
  if(donorID %in% c('L041','CC1','L114')){next()}
  channels = unique(as.character(mdat$orig.ident[mdat$donorID == donorID]))
  
  for(channel in channels){
    print(channel)
    
    ##-- Which cells have GATA1 reads covering that region?
    
    # All reads covering the region of interest (mutation site) in GATA1 gene
    reads = read.delim(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/',donorID,
                                   paste0(donorID,'_',channel,'_GATA1_reads.txt')),header = F)  
    # Remove rows which are not real reads (just another column in the reads)
    reads = reads[!grepl('^UB:Z:',reads$V1),]
    # Extract cellID from read
    reads$cellID = sapply(1:nrow(reads),function(i){reads[i,which(grepl('^CB\\:Z',reads[i,]))]})
    # if(grepl('cellranger',channel)){
    #   reads$cellID = paste0(gsub('_','.',gsub('.*_MY_','MY_',channel)),'_',gsub('^CB:Z:','',reads$cellID))  
    # }else{
    #   reads$cellID = paste0(channel,'_',gsub('^CB:Z:','',reads$cellID))  
    # }
    reads$cellID = paste0(channel,'_',gsub('^CB:Z:','',reads$cellID))  
    
    message(sprintf('[%s - %s] %d / %d cells found with mutation-site GATA1 reads',donorID,channel,sum(unique(reads$cellID) %in% mdat$cellID),length(mdat$cellID[mdat$donorID == donorID & mdat$orig.ident == channel])))
    #mdat$GATA1s_cov[mdat$cellID %in% reads$cellID] = TRUE
    reads_summary = reads %>% group_by(cellID) %>% summarise(nRead= n())
    
    out$GATA1s_cov[out$cellID %in% reads_summary$cellID] = reads_summary$nRead[match(out$cellID[out$cellID %in% reads_summary$cellID],reads_summary$cellID)]
    
    
    
    ##-- GATA1 reads with actual mutation
    fp = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads',donorID,
                   paste0(donorID,'_',channel,'_GATA1mut_reads.txt'))  
    
    if(file.size(fp) == 0){
      mut = data.frame(readID = c(),
                       V25 = c(),
                       cellID = c(),
                       containsMUT=c())
    }else{
      mut = read.delim(fp,sep = '\t',header = F,strip.white = T,fill = T,col.names = paste0('V',seq(1:29)))
      mut = mut[!grepl('^UB:Z',mut$V1),]
      
      mut$cellID = sapply(1:nrow(mut),function(i){mut[i,which(grepl('^CB\\:Z',mut[i,]))]})
      mut$cellID = paste0(channel,'_',gsub('^CB:Z:','',mut$cellID))
      colnames(mut)[1] = c('readID')
      mut = mut[,c('readID','cellID')]
      mut$containsMUT = T
      # check how many cellIDs are included in the seurat object
      message(sprintf('[%s - %s] %d / %d cells found with mutant GATA1 reads',donorID,channel,sum(unique(mut$cellID) %in% mdat$cellID),length(mdat$cellID[mdat$donorID == donorID & mdat$orig.ident == channel])))
      
      # add to output
      mut_summary = mut %>% group_by(cellID) %>% summarise(nRead = n())
      out$n_GATA1s_mutant[out$cellID %in% mut_summary$cellID] = mut_summary$nRead[match(out$cellID[out$cellID %in% mut_summary$cellID],mut_summary$cellID)]
    }
    
    
    ##-- GATA1 reads withOUT mutation
    if(file.size(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads',donorID,
                           paste0(donorID,'_',channel,'_GATA1_WT_reads.txt'))) == 0){
      out$n_GATA1s_wt[out$cellID %in% reads_summary$cellID] = NA
    }else{
      wt = read.delim(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads',donorID,
                                paste0(donorID,'_',channel,'_GATA1_WT_reads.txt')),sep = '\t',header = F)  
      wt = wt[!grepl('^UB:Z',wt$V1),]
      wt$cellID = sapply(1:nrow(wt),function(i){wt[i,which(grepl('^CB\\:Z',wt[i,]))]})
      wt$cellID = paste0(channel,'_',gsub('^CB:Z:','',wt$cellID))
      
      colnames(wt)[1] = c('readID')
      wt = wt[,c('readID','cellID')]
      wt$containsMUT = F
      
      # check how many cellIDs are included in the seurat object
      message(sprintf('[%s - %s] %d / %d cells found with wild-type GATA1 reads',donorID,channel,sum(unique(wt$cellID) %in% mdat$cellID),length(mdat$cellID[mdat$donorID == donorID & mdat$orig.ident == channel])))
      
      ## check if there are duplicated readID being called between MUT vs WT
      if(sum(mut$readID %in% wt$readID) > 0){
        stop(sprintf('Duplicated readID for channel %s from patient %s. Please check!',channel,donorID))
        View(mut[mut$readID %in% wt$readID,])
        View(wt[wt$readID %in% mut$readID,])
      }
      
      # add to output
      wt_summary = wt %>% group_by(cellID) %>% summarise(nRead = n())
      out$n_GATA1s_wt[out$cellID %in% wt_summary$cellID] = wt_summary$nRead[match(out$cellID[out$cellID %in% wt_summary$cellID],wt_summary$cellID)]
      
      unsureCells = sum(as.numeric(out$n_GATA1s_mutant[out$n_GATA1s_mutant !='?' & out$n_GATA1s_wt !='?']) > 0 & as.numeric(out$n_GATA1s_wt[out$n_GATA1s_mutant !='?' & out$n_GATA1s_wt !='?']) > 0)
      message(sprintf('[%s - %s] %d cells found with both mutation and WT GATA1 reads',donorID,channel,unsureCells))
      
    }
    
    # ## summarise WT / MUT by cellID
    # df = rbind(wt,mut)
    # df = df %>% group_by(cellID,containsMUT) %>% summarise(nReads = n_distinct(readID)) %>%
    #   group_by(cellID) %>% mutate(cat = n_distinct(containsMUT),
    #                               annot = ifelse(cat == 1 & containsMUT,'Tum',
    #                                              ifelse((cat == 1 & !containsMUT),'Norm','?')))
    # 
    # # Cells with '?' are most likely those with both normal and mutated reads
    # df.unsure = df[df$cat == 2,]
    # 
    # 
    # if(nrow(df.unsure) >= 1){
    #   df.unsure$containsMUT = ifelse(df.unsure$containsMUT,'Mut','WT')
    #   df.unsure = pivot_wider(df.unsure,id_cols = c('cellID','cat','annot'),names_from = 'containsMUT',values_from = 'nReads')
    #   df.unsure$annot = ifelse(df.unsure$WT <= 5 & df.unsure$Mut > 2*df.unsure$WT,'Tum',
    #                            ifelse(df.unsure$Mut <= 5 & df.unsure$WT > 2*df.unsure$Mut,'WT','unsure'))
    #   
    #   df.unsure$channelID = channel
    #   df.unsure$donorID = donorID
    #   unsure.GATAstatus.all = rbind(unsure.GATAstatus.all,df.unsure)
    #   df$annot[df$cellID %in% df.unsure$cellID] = df.unsure$annot[match(df$cellID[df$cellID %in% df.unsure$cellID],df.unsure$cellID)]
    # }
    # 
    # 
    # df$channelID = channel
    # df$donorID = donorID
    # 
    # tumCells_annot = rbind(tumCells_annot,df)
  }
}

# Cells with GATA1s-site coverage but no WT reads --> WT = 0
out$n_GATA1s_wt[out$GATA1s_cov != '?' & out$n_GATA1s_wt == '?'] = 0
out$n_GATA1s_mutant[out$GATA1s_cov != '?' & out$n_GATA1s_mutant == '?'] = 0

##--------------  GATA expression in each cell ---------------####
# Calculate number of reads mapped to GATA1 gene - does the cell express GATA1 or not?
gata1_expr = mlds@assays$RNA@counts['GATA1',]
out$GATA1_UMI_soupXed_count = gata1_expr[match(out$cellID,names(gata1_expr))]

View(out[is.na(out$n_GATA1s_wt) & out$GATA1s_cov > 0,])
out$n_GATA1s_wt[is.na(out$n_GATA1s_wt)] = 0

out$GATA1_status = '?'
out$GATA1_status[out$GATA1s_cov != '?' & out$n_GATA1s_mutant > 0 & out$n_GATA1s_wt == 0] = 'GATA1s_mutant'
out$GATA1_status[out$GATA1s_cov != '?' & out$n_GATA1s_mutant == 0 & out$n_GATA1s_wt > 0] = 'GATA1s_WT'
out$GATA1_status[out$GATA1s_cov != '?' & out$n_GATA1s_mutant > 0 & out$n_GATA1s_wt > 0] = 'GATA1s_unsure'
out$GATA1_status[out$GATA1_status == '?' & out$GATA1_UMI_soupXed_count > 0 ] = 'uninformative'
out$GATA1_status[out$GATA1_status == '?' & out$GATA1_UMI_soupXed_count == 0 ] = 'no_GATA1_expr'

out$current_GATA1s_status = mdat$GATA1_status[match(out$cellID,mdat$cellID)]
table(out$GATA1_status,out$current_GATA1s_status)
table(out$GATA1_status,out$donorID)


##-------    Resolve the unsure cells   ------####
out$mut_vs_WT_diff = as.numeric(out$n_GATA1s_mutant) - as.numeric(out$n_GATA1s_wt)
unsure = out[out$GATA1_status == 'GATA1s_unsure',]
DimPlot(mlds,cells.highlight = unsure$cellID[(unsure$mut_vs_WT_diff) <= -10])

out$GATA1_status[out$GATA1_status == 'GATA1s_unsure' & out$mut_vs_WT_diff >= 5] = 'GATA1s_mutant' 
out$GATA1_status[out$GATA1_status == 'GATA1s_unsure' & out$mut_vs_WT_diff <= -5] = 'GATA1s_WT' 

# Some cells have GATA1s_mut reads but no coverage - this is because I didn't explicitly filter for GATA1 reads when looking for mutant reads, but I did when extracting all reads covering the region
DimPlot(mlds,cells.highlight = out$cellID[(out$GATA1_status == 'uninformative') & out$n_GATA1s_mutant > 0 & out$n_GATA1s_wt > 0])
out$GATA1_status[out$GATA1_status == 'uninformative' & out$n_GATA1s_mutant > 0] = 'GATA1s_mutant' 
out$GATA1_status[out$GATA1_status == 'uninformative' & out$n_GATA1s_wt > 0] = 'GATA1s_WT' 

# Tumour cells (from cell type annotation) with GATA1_WT
View(out[out$annot_aug24 == 'Tumour' & out$GATA1_status == 'GATA1s_WT',])
table(out[out$annot_aug24 == 'Tumour' & out$GATA1_status == 'GATA1s_WT',]$donorID)
DimPlot(mlds,cells.highlight = out$cellID[out$annot_aug24 == 'Tumour' & out$GATA1_status == 'GATA1s_WT' & out$n_GATA1s_mutant > 0])
DimPlot(mlds,cells.highlight = out$cellID[out$annot_aug24 == 'NK' & out$GATA1_status == 'GATA1s_WT'])

# Mostly L038 and CC6, both of which are females --> Perhaps random x-inactivation?
##-------    Add results back to MLDS   ------####
mlds@meta.data = mlds@meta.data[,!grepl('GATA1',colnames(mlds@meta.data))]
table(is.na(match(mdat$cellID,out$cellID)))
colnames(out)[colnames(out) == 'mut_vs_WT_diff'] = 'GATA1_mut.vs.WT_diff'
mlds@meta.data = cbind(mlds@meta.data,out[match(mdat$cellID,out$cellID),!colnames(out) %in% c('cellID','donorID','annot_aug24','current_GATA1s_status')])

#mdat = mdat[,!grepl('GATA1',colnames(mdat))]
#mdat = cbind(mdat,out[match(mdat$cellID,out$cellID),!colnames(out) %in% c('cellID','donorID','annot_aug24','current_GATA1s_status')])


DimPlot(mlds,cells.highlight = mdat$cellID[mdat$GATA1_status == 'GATA1s_mutant' & mdat$annot_aug24 != 'Tumour'])
DimPlot(mlds,group.by = 'GATA1_status',cols = col25)


# Save the output
write.csv(out, '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/GATA1mut_reads_summary_perCell_aug24.csv')
out = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/GATA1mut_reads_summary_perCell_aug24.csv')



View(table(mdat$GATA1s_status,mdat$GATA1s_cov))
View(table(mdat$GATA1s_status,mdat$sex,mdat$donorID))
View(table(mdat$GATA1s_status,mdat$finalAnn_refined_v1))


## Assign Tumour status by combining GATA1s genotyping + Cell type annotation
mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_mutant' & mlds$annot_aug24 %in% c('LE','ME','MEP','EE')] = 'Tumour'
mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_mutant' & !mlds$annot_aug24 %in% c('Tumour')] = 'doublets'

table(mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_WT'])
mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_WT' & mlds$annot_aug24 %in% c('Tumour')] = 'Tumour_WT'
mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_WT' & mlds$annot_aug24 == '?'] = 'MEP'
mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_WT' & mlds$annot_aug24 == 'MK_WT'] = 'MK'
mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_WT' & mlds$annot_aug24 == 'unsure_EE'] = 'EE'
mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_WT' & !mlds$annot_aug24 %in% c('Tumour_WT','MEP','EE','ME','LE','HSC_MPP')] = 'doublets'


## Assess Tumour with GATA1_WT 
mlds$annot_aug24[mlds$GATA1_status == 'GATA1s_WT' & mlds$annot_aug24 %in% c('Tumour','Tumour_maybe')] = 'Tumour_WT'

## Tumour_WT = doublets
                                                               
mdat = cbind(mlds@meta.data,mlds@reductions$umap@cell.embeddings)
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_mdat.csv')

#write.csv(mlds@meta.data[,!colnames(mlds@meta.data) %in% c('n_GATA1reads','annot_aug24.1')],'MLDS_GATA1mut_scGenotyping_mdat_jul24.csv')

DimPlot(mlds,group.by = 'GATA1_status',cols = c(col25[c(2,3,1)],grey(0.8),grey(0.4)),label = F,label.box = T,repel = T) + NoAxes()
DimPlot(mlds,cells.highlight = mlds$cellID[mlds$GATA1_status == 'GATA1s_mutant']) + NoAxes()
DimPlot(mlds,cells.highlight = mlds$cellID[mlds$GATA1_status == 'GATA1s_mutant' & mlds$donorID == 'L182' & mlds$timePoint == 'Diagnostic'])


FeaturePlot(mlds,'GATA1')

































##-------------------------------##
##    Further investigation    ####
##-------------------------------##


# Investigate "Tumour" with WT GATA1 reads
tumWT = mlds$finalAnn_broad[mlds$GATA1s_status == 'WT' &
                              !mlds$finalAnn_refined_v1 %in% c('LE','ME','EE','MK','MEP','HSC_MPP',
                                                               unique(mlds$finalAnn_refined_v1[grepl('Tum',mlds$finalAnn_refined_v1)]))] 

DimPlot(mlds,cells.highlight = mlds$cellID[mlds$GATA1s_status == 'WT' & mlds$finalAnn_broad == 'MK'])

## Identify some special groups of cells to have a look at
d38 = subset(mlds,subset = cellID %in% mlds$cellID[mlds$donorID=='L038'])
d38 = standard_clustering(d38)
DimPlot(d38,group.by = 'seurat_clusters',label = T,label.box = T,repel = T)+NoLegend()
DimPlot(d38,group.by = 'tumAnnot',label = T,label.box = T,repel = T)+NoLegend()
VlnPlot(d38,features = c('nCount_RNA','nFeature_RNA','percent.mt'),group.by = 'seurat_clusters',ncol = 1)
FeaturePlot(d38,features = c('nCount_RNA','nFeature_RNA','percent.mt','TPSAB1'))
FeaturePlot(mlds,features = c('PLEK','PF4','PPBP','ALAS2','TPSAB1','HDC'))

DimPlot(mlds,cells.highlight = l038_gata1s_c7)
#1. L038 GATA1s pos normClusters
l038_gata1s_c7 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 == 'Tumour' & d38$seurat_clusters == 7]
#2. L038 GATA1s neg normClusters
l038_gata1sNeg_c7 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 != 'Tumour' & d38$seurat_clusters == 7]
#3. L038 norm normClusters
l038_gata1sNeg_c19 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 != 'Tumour' & d38$seurat_clusters == 19]
l038_gata1sPos_c19 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 == 'Tumour' & d38$seurat_clusters == 19]
l038_gata1sNeg_c8 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 != 'Tumour' & d38$seurat_clusters == 8]
l038_gata1sPos_c8 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 == 'Tumour' & d38$seurat_clusters == 8]

l038_gata1sNeg_c30 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 != 'Tumour' & d38$seurat_clusters == 30]
l038_gata1sPos_c30 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 == 'Tumour' & d38$seurat_clusters == 30]

l038_gata1sNeg_c13 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 != 'Tumour' & d38$seurat_clusters == 13]
l038_gata1sPos_c13 = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 == 'Tumour' & d38$seurat_clusters == 13]

l038_gata1sNeg_tum = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 != 'Tumour' & d38$seurat_clusters %in% c(2,12,4,5,9,27,18)]
l038_gata1sPos_tum = d38$cellID[d38$donorID=='L038' & d38$finalAnn_refined_v1 == 'Tumour' & d38$seurat_clusters %in% c(2,12,4,5,9,27,18)]

mlds$specialCells = ifelse(mlds$cellID %in% c(l038_gata1s_c7,l038_gata1sNeg_c7),'c7',
                           ifelse(mlds$cellID %in% c(l038_gata1sNeg_c19,l038_gata1sPos_c19),'c19',
                                  ifelse(mlds$cellID %in% c(l038_gata1sNeg_c8,l038_gata1sPos_c8),'c8',
                                         ifelse(mlds$cellID %in% c(l038_gata1sNeg_c30,l038_gata1sPos_c30),'c30',
                                                ifelse(mlds$cellID %in% c(l038_gata1sNeg_c13,l038_gata1sPos_c13),'c13',
                                                       ifelse(mlds$cellID %in% c(l038_gata1sNeg_tum,l038_gata1sPos_tum),'tum','none'))))))




## Quantify the number of cells which can confidently be identified as tumour vs normal vs unsure


## Quantify scGenotyping method sensitivity
unsure.GATAstatus.all$sex = mlds$sex[match(unsure.GATAstatus.all$donorID,mlds$donorID)]
a = as.data.frame(table(unsure.GATAstatus.all$annot,unsure.GATAstatus.all$donorID,unsure.GATAstatus.all$sex))
View(a[a$Freq >0,])
ggplot(unsure.GATAstatus.all,aes(log10(WT),log10(Mut),col=annot))+
  geom_point(size=0.8)+
  geom_abline()+
  facet_wrap(vars(donorID))+
  theme_bw(base_size = 15)+xlab("log10 # WT reads ") + ylab("log10 # Mut reads")




DimPlot(mlds,group.by = 'finalAnn_refined_v1',label = T,label.box = T,repel = T)+NoLegend()
View(table(mlds$finalAnn_refined_v1,mlds$GATA1s_status))

mdat = mlds@meta.data
mdat$broadAnn = ifelse(mdat$finalAnn_refined_v1 %in% c('Tumour'),'Tumour',
                       ifelse(mdat$finalAnn_refined_v1 %in% c('Tum_MK?','Tum_GMP?','Tumour_noMut'),'maybe_Tumour',
                              ifelse(mdat$finalAnn_refined_v1 %in% c('MK','norm_MK'),'MK',
                                     ifelse(mdat$finalAnn_refined_v1 %in% c('MEP','EE'),as.character(mdat$finalAnn_refined_v1),
                                            ifelse(mdat$finalAnn_refined_v1 %in% c('ME','LE'),'ME_LE','others')))))
View(table(mdat$finalAnn_refined_v1,mdat$broadAnn))


mdat$nReads_cat = ifelse(mdat$n_GATA1_supporting_reads == '?','?',
                         ifelse(mdat$n_GATA1_supporting_reads == '1','1',
                                ifelse(mdat$n_GATA1_supporting_reads %in% as.character(c(2:5)),'2-5',
                                       ifelse(mdat$n_GATA1_supporting_reads %in% as.character(c(6:10)),'6-10','>10'))))
df = df %>% group_by(Cell_type,donorID,channelID,timepoint,nReads_cat) %>% summarise(n=n())


df = mdat %>% group_by(donorID,timePoint,broadAnn,GATA1s_status,nReads_cat) %>% summarise(nCell = n()) %>% 
  group_by(donorID,timePoint,broadAnn) %>% mutate(totalCells = sum(nCell),
                                                  perc = nCell / totalCells)
df$GATA1s_status = ifelse(df$nReads_cat != '?',paste0(df$GATA1s_status,':',df$nReads_cat),as.character(df$GATA1s_status))
df$GATA1s_status = factor(df$GATA1s_status,levels = c('noCov','cov_noInfo','unsure','Mut:>10','Mut:6-10','Mut:2-5','Mut:1',
                                                      'WT:>10','WT:6-10','WT:2-5','WT:1'))
df$broadAnn = factor(df$broadAnn,levels = c('Tumour','maybe_Tumour','MK','MEP','EE','ME_LE','others'))
library(RColorBrewer)
display.brewer.all()
ggplot(df,aes(broadAnn,perc,fill = GATA1s_status))+
  geom_col()+
  facet_grid(vars(donorID),vars(timePoint))+
  scale_fill_manual(values = c(grey(0.8),grey(0.4),col25[1],rev(brewer.pal(n=5,'OrRd')[-1]),rev(brewer.pal(n=5,'Greens')[-1])))+
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1)) + 
  ylab('Percentage of # Cells in the group')

## Generate the annotation group to extract BAM files
mdat = mlds@meta.data
write.csv(mdat,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jul23/MLDS_annot_jul23_noUnknowns_withGATA1s_mdat.csv')
mdat = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/1_annotation_scRNA_MLDS/jul23/MLDS_annot_jul23_noUnknowns_withGATA1s_mdat.csv')



##--------------------------------##
##    EZH2 mutation in L038     ####
##--------------------------------##

tumCells_annot = data.frame()
mlds$EZH2_cov = '?'
mlds$n_EZH2reads = '?'

unsure.EZH2status.all = data.frame()

for(donorID in c('L038')){
  print(donorID)
  if(donorID %in% c('L041','CC1')){next()}
  channels = unique(as.character(mlds$orig.ident[mlds$donorID == donorID]))
  if(donorID == 'CC2'){
    channels = gsub('\\.','_',channels)
    channels = paste0('cellranger700_count_48097_',channels)
  }
  
  for(channel in channels){
    print(channel)
    
    ##-- Which cells have EZH2 reads covering that region?
    reads = read.delim(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/',donorID,
                                 paste0(donorID,'_',channel,'_EZH2_reads.txt')))
    w = apply(reads, 2, FUN=function(x){sum(grepl('^CB\\:Z\\:',x))})
    reads$cellID = apply(reads,MARGIN = 1,FUN = function(x){paste0(channel,'_',gsub('^CB:Z:','',x[grepl('^CB\\:Z\\:',x)]))})
    table(reads$cellID %in% mlds$cellID)
    mlds$EZH2_cov[mlds$cellID %in% reads$cellID] = TRUE
    
    reads = reads %>% group_by(cellID) %>% summarise(n=n())
    
    ##-- EZH2 reads with mutation
    fp = file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads',donorID,
                   paste0(donorID,'_',channel,'_EZH2mut_reads.txt'))
    if(file.size(fp) == 0){
      mut = data.frame(readID = c(),
                       V25 = c(),
                       cellID = c(),
                       containsMUT=c())
    }else{
      mut = read.delim(fp,sep = '\t',header = F,strip.white = T,fill = T,col.names = paste0('V',seq(1:29)))
      w = apply(mut, 2, FUN=function(x){sum(grepl('^CB\\:Z\\:',x))})
      mut$cellID = apply(mut,MARGIN = 1,FUN = function(x){paste0(channel,'_',gsub('^CB:Z:','',x[grepl('^CB\\:Z\\:',x)]))})
      colnames(mut)[1] = c('readID')
      mut = mut[,c('readID','V25','cellID')]
      mut$containsMUT = T
      # check how many cellIDs are included in the seurat object
      table(mut$cellID %in% mlds$cellID)
      
      # There are some reads where CB is not in field 25, probably because these reads have some different tags (eg. AN for reads mapped to antisense strand)
      #mut$group = ifelse(grepl('^CB\\:Z',mut$V25),'norm','odd')
      # check if these Mut_odd CB are unique - not included in mut_norm
      #table(mut$cellID[mut$group == 'odd'] %in% mlds$cellID)
      #table(mut$cellID[mut$group == 'odd'] %in% mut$cellID[mut$group != 'odd'])
      
    }
    
    
    ##-- EZH2 reads withOUT mutation
    wt = read.delim(file.path('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads',donorID,
                              paste0(donorID,'_',channel,'_EZH2_WT_reads.txt')),sep = '\t',header = F)
    wt$cellID = apply(wt,MARGIN = 1,FUN = function(x){paste0(channel,'_',gsub('^CB:Z:','',x[grepl('^CB\\:Z\\:',x)]))})
    colnames(wt)[1] = c('readID')
    wt = wt[,c('readID','V25','cellID')]
    wt$containsMUT = F
    
    # check how many cellIDs are included in the seurat object
    table(wt$cellID %in% mlds$cellID)
    
    # There are some reads where CB is not in field 25, probably because these reads have some different tags (eg. AN for reads mapped to antisense strand)
    #wt$group = ifelse(grepl('^CB\\:Z',wt$V25),'norm','odd')
    # check if these wt_odd CB are unique - not included in wt_norm
    #table(wt$cellID[wt$group == 'odd'] %in% mlds$cellID)
    # if(sum(wt$cellID[wt$group == 'odd'] %in% wt$cellID[wt$group != 'odd'])>0){
    #   w = wt[wt$cellID %in% wt_odd$cellID[wt_odd$cellID %in% wt_norm$cellID],]
    #   print(w)
    # }
    
    ## check if there are duplicated readID being called between MUT vs WT
    if(sum(mut$readID %in% wt$readID) > 0){
      stop(sprintf('Duplicated readID for channel %s from patient %s. Please check!',channel,donorID))
      View(mut[mut$readID %in% wt$readID,])
      View(wt[wt$readID %in% mut$readID,])
    }
    
    
    ## summarise by cellID
    df = rbind(wt,mut)
    # df = df %>% group_by(cellID,containsMUT,group) %>% summarise(nReads = n_distinct(readID)) %>% 
    #   group_by(cellID) %>% mutate(mut_group = sum(containsMUT),
    #                               nGroup = n_distinct(paste0(containsMUT,':',group)),
    #                               annot = ifelse(mut_group == 1 & containsMUT,'Tum',
    #                                              ifelse((mut_group == 1 & !containsMUT)|(mut_group == 0 & containsMUT),'?',
    #                                                     ifelse(mut_group == 0 & !containsMUT,'Norm','?'))))
    df = df %>% group_by(cellID,containsMUT) %>% summarise(nReads = n_distinct(readID)) %>%
      group_by(cellID) %>% mutate(cat = n_distinct(containsMUT),
                                  #nGroup = n_distinct(paste0(containsMUT,':',group)),
                                  annot = ifelse(cat == 1 & containsMUT,'Tum',
                                                 ifelse((cat == 1 & !containsMUT),'Norm','?')))
    
    #df$id = seq(1:nrow(df))
    
    # Cells with '?' are most likely those with both normal and mutated reads
    df.unsure = df[df$cat == 2,]
    df.unsure$containsMUT = ifelse(df.unsure$containsMUT,'Mut','WT')
    df.unsure = pivot_wider(df.unsure,id_cols = c('cellID','cat','annot'),names_from = 'containsMUT',values_from = 'nReads')
    if(nrow(df.unsure) >= 1){
      df.unsure$annot = ifelse(df.unsure$WT <= 5 & df.unsure$Mut > 2*df.unsure$WT,'Tum',
                               ifelse(df.unsure$Mut <= 5 & df.unsure$WT > 2*df.unsure$Mut,'WT','unsure'))
      # df.unsure$annot = sapply(seq(1:nrow(df.unsure)),FUN = function(x){
      #   cellID = df.unsure$cellID[x]
      #   if(n_distinct(df.unsure$containsMUT[df.unsure$cellID == cellID]) == 2){
      #     nNorm_read = df.unsure$nReads[df.unsure$cellID == cellID & df.unsure$containsMUT == F]
      #     nMut_read = df.unsure$nReads[df.unsure$cellID == cellID & df.unsure$containsMUT == T]
      #     if(nMut_read > 2*nNorm_read){
      #       return('Tum')
      #     }else if (nNorm_read > 1.5*nMut_read){
      #       return('Norm')
      #     }else{
      #       return('unsure')
      #     }
      #   }else{
      #     return('unsure')
      #   }
      # 
      # })
      df.unsure$channelID = channel
      df.unsure$donorID = donorID
      unsure.EZH2status.all = rbind(unsure.EZH2status.all,df.unsure)
      df$annot[df$cellID %in% df.unsure$cellID] = df.unsure$annot[match(df$cellID[df$cellID %in% df.unsure$cellID],df.unsure$cellID)]
    }
    
    #df$annot[match(df.unsure$id,df$id)] = 'unsure' 
    df$channelID = channel
    df$donorID = donorID
    
    tumCells_annot = rbind(tumCells_annot,df)
  }
}

#tumCells_annot$cellID[tumCells_annot$donorID == 'CC2'] = gsub('^.*_MY_200531_1440628','MY.200531.1440628',tumCells_annot$cellID[tumCells_annot$donorID == 'CC2'])

# Save the output
write.csv(tumCells_annot, '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/EZH2mut_reads_summary_perCell.csv')



DimPlot(mlds,cells.highlight = tumCells_annot$cellID[tumCells_annot$annot == 'Tum'])

mlds$EZH2_status = ifelse(mlds$cellID %in% tumCells_annot$cellID[tumCells_annot$annot == 'Tum'],'Mut',
                          ifelse(mlds$cellID %in% tumCells_annot$cellID[tumCells_annot$annot == 'Norm'],'WT',
                                 ifelse(mlds$cellID %in% tumCells_annot$cellID[tumCells_annot$annot == 'unsure'],'unsure',
                                        ifelse(mlds$EZH2_cov==T,'cov_noInfo','noCov'))))
mlds$n_EZH2_supporting_reads = '?'
df = tumCells_annot[tumCells_annot$annot == 'Tum',]
mlds$n_EZH2_supporting_reads[mlds$EZH2_status =='Mut'] = df$nReads[match(mlds$cellID[mlds$EZH2_status =='Mut'],df$cellID)]
df = tumCells_annot[tumCells_annot$annot == 'Norm',]
mlds$n_EZH2_supporting_reads[mlds$EZH2_status =='WT'] = df$nReads[match(mlds$cellID[mlds$EZH2_status =='WT'],df$cellID)]

View(table(mlds$EZH2_status,mlds$EZH2_cov))
View(table(mlds$EZH2_status,mlds$sex,mlds$donorID))
a = as.data.frame(table(mlds$EZH2_status[mlds$donorID == 'L038'],mlds$GATA1s_status[mlds$donorID == 'L038'],mlds$timePoint[mlds$donorID == 'L038']))
a = a[a$Freq >0,]
View(a[a$Var1 %in% c('WT','Mut') | a$Var2%in% c('WT','Mut') ,])

DimPlot(mlds,cells.highlight = mlds$cellID[mlds$EZH2_status == 'Mut' & mlds$GATA1s_status == 'WT'],cells = mlds$cellID[mlds$donorID == 'L038'])+ggtitle('L038 - EZH2mut & GATA1wt')
DimPlot(mlds,cells.highlight = mlds$cellID[mlds$EZH2_status == 'WT' & mlds$GATA1s_status == 'Mut'],cells = mlds$cellID[mlds$donorID == 'L038'])+ggtitle('L038 - EZH2wt & GATA1s')

write.csv(mlds@meta.data,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/MLDS_231127_GATA1s.EZH2.status_mdat.csv')



