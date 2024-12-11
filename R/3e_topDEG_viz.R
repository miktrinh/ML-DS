##    Vizualization of 3e - DEGs in fLiver + fKidney + fAdrenal     ##

outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/apr24/without_published2n/geno'
setwd(outDir)



#############
# Libraries #
#############

library(Seurat)
library(GenomicFeatures)
library(DESeq2)
library(ComplexHeatmap)
library(reshape2)
library(zoo)
library(tidyverse)
library(RColorBrewer)
library(readxl)
source("~/lustre_mt22/generalScripts/utils/misc.R")
source("~/lustre_mt22/generalScripts/utils/sc_utils.R")
source("~/lustre_mt22/generalScripts/utils/pseudobulk.R")



#### Parameters ####
tgtChrs=c(1:22,'X')

#Define genomic coordinates
gtf = '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf'
txdb = makeTxDbFromGFF(gtf)
gns = genes(txdb)

geneMap = read.delim('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Data/Donor15680/liver/fLiver_MY_200531_10043298/filtered_feature_bc_matrix/features.tsv.gz',header = F)
colnames(geneMap) = c('ensID','geneSym','GEX')
geneMap$geneSym[duplicated(geneMap$geneSym)] = paste0(geneMap$geneSym[duplicated(geneMap$geneSym)],'.1')
geneMap$geneSym = gsub('_','-',geneMap$geneSym)
geneMap$chr = as.character(seqnames(gns[match(geneMap$ensID,gns$gene_id)]))

coords = gns[geneMap$ensID]

# Linearise genomic position
tgtChrs = c(1:22,'X')
w = which(as.character(seqnames(coords)) %in% tgtChrs)
coords = coords[w]
seqlevels(coords) = tgtChrs

# Get genomic coordinates
#Order aa and bb by genomic coords
#If it's positive stranded or unknown use start, else end
coords$TSS = ifelse(strand(coords)=='-',end(coords),start(coords))
o = order(seqnames(coords),coords$TSS)
coords = coords[o]

xpos = coords$TSS
chrLens = sapply(split(xpos,as.character(seqnames(coords))),max)
chrLens = chrLens[seqlevels(coords)]
offset = c(0,cumsum(as.numeric(chrLens))[-length(chrLens)])
names(offset) = names(chrLens)
xpos = xpos + offset[as.character(seqnames(coords))]+1e6 #1mb is correction for approximate chromosome length
coords$xpos_linearised = xpos

chrLens_pos = sapply(split(xpos,names(xpos)),max)


##-------------------------------##
##      Helper functions      ####
##-------------------------------##
# These functions have been moved to ~/lustre_mt22/generalScripts/utils/pseudobulk.R

# filter_sratObj = function(tissue,ageMatch=T,remove_cyclingCells=T){
#   
#   ## Import AK_REF merged seurat object for the tissue under consideration
#   srat_in_fp = ifelse(tissue == 'liver', '~/lustre_mt22/Aneuploidy/Results/2_annotation/liver/sept23/liver_liverREFmerged_clean_processed_annotated_noUnknowns_0923.RDS',
#                       file.path('~/lustre_mt22/Aneuploidy/Results/2_annotation',tissue,'sept23',paste0(tissue,'_',tissue,'REFmerged_clean_processed_annotated.RDS')))  
#   
#   
#   
#   if(!file.exists(srat_in_fp)){
#     warning(sprintf('Cannot find input directory for annotated AK-REF merged sratObj for tissue %s \nThe path is: \n%s',tissue, srat_in_fp))
#   }else{
#     big.srat = readRDS(srat_in_fp)
#     tissue_toKeep = tissue
#     big.srat$tissue = tissue_toKeep
#     
#     big.srat$Sex[big.srat$Sex == 'female'] = 'XX'
#     big.srat$Sex[big.srat$Sex == 'male'] = 'XY'
#     
#     if(!'finalAnn_broad' %in% colnames(big.srat@meta.data)){
#       big.srat$finalAnn = as.character(big.srat$cluster_ann)
#     }else{
#       big.srat$finalAnn = as.character(big.srat$finalAnn_broad)
#     }
#     
#     if(tissue != 'liver'){
#       ## Remove cells from donor 15806
#       big.srat = subset(big.srat,subset = donorID %in% unique(big.srat$donorID[big.srat$donorID != '15806']))
#       big.srat = subset(big.srat,subset = finalAnn %in% unique(big.srat$finalAnn[!big.srat$finalAnn %in% c('bad_cluster_15806','badcluster_15806')]))  
#     }
#     
#     
#     
#     ### check that only cells in G1 are retained
#     #if(tissue != 'liver'){
#     if(remove_cyclingCells){
#       if(n_distinct(big.srat@meta.data$Phase)>1){
#         message(sprintf('Removing cycling cells from big.srat for tissue %s',tissue))
#         big.srat = subset(big.srat, subset = Phase == 'G1')
#       }  
#     }
#     
#     # bin Gestational age into age group
#     if(tissue == 'kidney'){
#       
#       big.srat@meta.data$ageGroup = ifelse(big.srat@meta.data$gestationalAge %in% c('7+6','8+1','9+1','10pcw'),'7_10',
#                                            ifelse(big.srat@meta.data$gestationalAge %in% c('11pcw','12+0', '12pcw','13+6','14pcw'),'11_14','15_17'))
#       big.srat@meta.data$finalAnn_detailed = as.character(big.srat@meta.data$finalAnn)
#       
#       big.srat@meta.data$finalAnn[big.srat@meta.data$finalAnn_detailed %in% c('RVCSB','SSBm.d','SSBpr','SSBpod')] = 'RV_C_S_SB'
#       big.srat@meta.data$finalAnn[big.srat@meta.data$finalAnn_detailed %in% c('UBCD','CnT')] = 'UBCD_CnT'
#     }else if(tissue == 'liver'){
#       #big.srat@meta.data$gestationalAge[is.na(big.srat@meta.data$gestationalAge)] = '14pcw'
#       big.srat$finalAnn = as.character(big.srat$finalAnn_broad)
#       big.srat$finalAnn[big.srat$finalAnn %in% c('earlyMK')] = 'MK'
#       big.srat$finalAnn[big.srat$finalAnn %in% c('promyelocyte','myelocyte')] = 'myelocyte'
#       big.srat$finalAnn[big.srat$finalAnn %in% c('LMPP_ELP','pro.B.cell','pre.B.cell')] = 'B.cell.prog'
#       big.srat$finalAnn[big.srat$finalAnn %in% c('ILC.precursor','T.cell','NK_T')] = 'NK.T'
#       if(ageMatch){
#         big.srat@meta.data$ageGroup = ifelse(big.srat@meta.data$gestationalAge %in% c('7+6','8+1'),'8pcw',
#                                              ifelse(big.srat@meta.data$gestationalAge %in% c('9+1','9+5'),'10pcw',
#                                                     ifelse(big.srat@meta.data$gestationalAge %in% c('11+3', '11pcw'),'11pcw',
#                                                            ifelse(big.srat@meta.data$gestationalAge %in% c('12', '12pcw'),'12pcw',
#                                                                   ifelse(big.srat@meta.data$gestationalAge %in% c('13+6','14+3', '14pcw'),'14pcw',
#                                                                          ifelse(big.srat@meta.data$gestationalAge %in% c('15pcw'),'15pcw',
#                                                                                 ifelse(big.srat@meta.data$gestationalAge %in% c('16','16+2', '16pcw'),'16pcw','17pcw')))))))
#         #big.srat@meta.data$ageGroup[big.srat@meta.data$ageGroup %in% c('11pcw','12pcw')] = '11_12pcw'
#         big.srat@meta.data$ageGroup[big.srat@meta.data$ageGroup %in% c('15pcw','16pcw')] = '15_16pcw'
#       }else{
#         big.srat@meta.data$ageGroup = ifelse(big.srat@meta.data$gestationalAge %in% c('7+6','8+1','9+1','9+5'),'7_10',
#                                              ifelse(big.srat@meta.data$gestationalAge %in% c('11+3', '11pcw','12', '12pcw','13+6','14+3', '14pcw'),'11_14','15_18'))
#         
#       }
#       
#       #big.srat@meta.data$finalAnn_broad = as.character(big.srat@meta.data$finalAnn_broad)
#       #big.srat@meta.data$finalAnn_broad[big.srat@meta.data$finalAnn_broad %in% c('Pre.pro.B.cell','pre.B.cell','pro.B.cell','Immature.B.cell')] = 'B_cell_progenitor'
#       #big.srat@meta.data$finalAnn_broad[big.srat@meta.data$finalAnn_broad %in% c('Mono.Mac','Monocyte')] = 'Mono.Mac'
#       #big.srat@meta.data$finalAnn_detailed = as.character(big.srat@meta.data$finalAnn)
#       #big.srat@meta.data$finalAnn = as.character(big.srat@meta.data$finalAnn_broad)
#       
#       
#     }else if(tissue == 'adrenal'){
#       big.srat@meta.data$termination[is.na(big.srat@meta.data$termination)] = '??'
#       big.srat@meta.data$ageGroup = ifelse(big.srat@meta.data$gestationalAge %in% c('8','8+6','10pcw','10+5'),'7_10',
#                                            ifelse(big.srat@meta.data$gestationalAge %in% c('11','11pcw','12pcw'),'11_14',
#                                                   ifelse(big.srat@meta.data$gestationalAge %in% c('15pcw','16pcw', '17pcw'),'15_18','18+'))) 
#       
#       big.srat@meta.data$finalAnn_detailed = as.character(big.srat@meta.data$finalAnn)
#       big.srat@meta.data$finalAnn[big.srat@meta.data$finalAnn_detailed %in% c('Mesenchyme - early','Mesenchyme')] = 'Mesenchyme'
#       big.srat@meta.data$finalAnn[big.srat@meta.data$finalAnn_detailed %in% c('SCPs','Bridge')] = 'SCP_Bridge'
#       
#       big.srat$Genotype[is.na(big.srat$Genotype)] = 'T18'
#       big.srat$Sex[is.na(big.srat$Sex)] = 'XY'
#       big.srat$gestationalAge[is.na(big.srat$gestationalAge)] = '11pcw'
#       big.srat$assay[is.na(big.srat$assay)] = unique(big.srat$assay[big.srat$Genotype == 'T21'])
#     }
#     
#   }
#   
#   # # Change the annotation label as these are essentially the same cell population
#   # big.srat@meta.data$finalAnn[big.srat@meta.data$finalAnn == 'Mesenchyme - early'] = 'MSC'
#   
#   big.srat@meta.data$donorID = as.character(big.srat@meta.data$donorID)
#   
#   ### Number of samples per genotype
#   big.srat@meta.data$donorID = as.character(big.srat@meta.data$donorID)
#   sum_CT = big.srat@meta.data %>% 
#     dplyr::select(c(Genotype,donorID,finalAnn,ageGroup,Sex,Phase)) %>% 
#     group_by(Genotype) %>% mutate(nDonor_perGeno = n_distinct(donorID),
#                                   nFemale_perGeno = n_distinct(donorID[Sex == 'XX']),
#                                   nMale_perGeno = n_distinct(donorID[Sex == 'XY']),
#                                   # nAge_u10_perGeno = n_distinct(donorID[ageGroup == '8pcw']),
#                                   # nAge_u14_perGeno = n_distinct(donorID[ageGroup == '11_14']),
#                                   # nAge_u18_perGeno = n_distinct(donorID[ageGroup == '15_18']),
#                                   # nAge_u21_perGeno = n_distinct(donorID[ageGroup == '17+'])
#     ) %>% 
#     group_by(Genotype,finalAnn) %>% mutate(nCells_perGeno = n()) %>% 
#     group_by(Genotype,finalAnn,donorID) %>% summarise(nDonor_perGeno=unique(nDonor_perGeno),
#                                                       nFemale_perGeno = unique(nFemale_perGeno),
#                                                       nMale_perGeno = unique(nMale_perGeno),
#                                                       # nAge_u10_perGeno = unique(nAge_u10_perGeno),
#                                                       # nAge_u14_perGeno = unique(nAge_u14_perGeno),
#                                                       # nAge_u18_perGeno = unique(nAge_u18_perGeno),
#                                                       # nAge_u21_perGeno = unique(nAge_u21_perGeno),
#                                                       nCells_perGeno = unique(nCells_perGeno),
#                                                       nCells_perDonor = n()) %>% 
#     group_by(Genotype,finalAnn) %>% mutate(nDonor_with100cells=sum(nCells_perDonor>=100),
#                                            nDonor_with50cells=sum(nCells_perDonor>=50),
#                                            nDonor_lost_with100cells = nDonor_perGeno - nDonor_with100cells,
#                                            nDonor_lost_with50cells = nDonor_perGeno - nDonor_with50cells)
#   
#   
#   
#   sum_CT$tissue = tissue
#   #nCell_summary_unfiltered = rbind(nCell_summary_unfiltered,sum_CT)
#   
#   ### To get a decent control - we need at least 3 diploid-donors to have at least 3 values for the boxplots
#   # sum_CT = sum_CT[(sum_CT$finalAnn %in% sum_CT$finalAnn[sum_CT$Genotype=='diploid' & sum_CT$nDonor_with100cells >=3]) &
#   #                   !(sum_CT$finalAnn %in% c('Other','Leukocytes','Erythroblasts','?')),]
#   # 
#   # ct_toRemove = sum_CT[,c('finalAnn','Genotype','nDonor_with100cells')] %>% 
#   #   mutate(type = ifelse(Genotype == 'diploid','diploid','triploid')) %>% group_by(finalAnn, type) %>% 
#   #   summarise(nDonor = ifelse(type == 'diploid',(unique(nDonor_with100cells)>=3),sum(nDonor_with100cells>=1))) %>% 
#   #   distinct() %>% group_by(finalAnn) %>% summarise(toRemove = sum(nDonor ==0)>0)
#   # 
#   # sum_CT = sum_CT[sum_CT$finalAnn %in% ct_toRemove$finalAnn[ct_toRemove$toRemove==F],]
#   # 
#   # nCell_summary_filtered = rbind(nCell_summary_filtered,sum_CT)
#   # 
#   # ct_toKeep = unique(sum_CT$finalAnn)
#   # 
#   # 
#   # ### Subset to keep cells from relevant tissue only (both AK and REF)
#   # #big.srat = standard_clustering(big.srat)
#   # # DimPlot(big.srat,group.by = 'finalAnn',label = T)+NoLegend()
#   # #   xlim(5,15)+ylim(-2,6)
#   # # DimPlot(big.srat,cells.highlight = rownames(big.srat@meta.data[big.srat@meta.data$finalAnn == 'Endothelium',]))+NoLegend()
#   # #   
#   # 
#   # big.srat = subset(big.srat, subset = finalAnn %in% ct_toKeep)
#   
#   return(list(big.srat, sum_CT))
# }
# 
# 
# import_pbDEGresults = function(outDir,pCut = 0.05,max_pct_cells=10,allGenes=F,tissue=''){
#   ###### Import the results ######
#   out = list()
#   for(f in list.files(file.path(outDir,tissue),pattern = 'RDS$',full.names = T)){
#     geno = gsub('_.*$','',basename(f))
#     res = readRDS(f)
#     res = res[grepl('unmerged',names(res))]
#     
#     if(allGenes){
#       res_out = do.call(rbind,lapply(1:length(res),FUN = function(i){
#         ct = gsub('_unmerged$','',names(res)[i])
#         de = as.data.frame(res[[i]][['res']])  
#         de$ensID = rownames(de)
#         de$geno = geno
#         de$ct = ct
#         return(de)
#       }))
#       
#     }else{
#       res_out = do.call(rbind,lapply(1:length(res),FUN = function(i){
#         ct = gsub('_unmerged$','',names(res)[i])
#         de = res[[i]][['de']]
#         de$geno = geno
#         de$ct = ct
#         ## Removal rubbish genes
#         de = de[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',de$geneSym),]
#         #de$padj = p.adjust(de$pvalue)
#         de = de[de$padj < pCut,]
#         return(de)
#       }))
#     }
#     
#     out[[geno]] = res_out
#   }
#   
#   if(allGenes){
#     df = do.call(rbind,out)
#     table(df$geno,df$ct)
#     df$geneSym[is.na(df$geneSym)] = geneMap$geneSym[match(df$ensID[is.na(df$geneSym)],geneMap$ensID)]
#     return(df)
#   }else{
#     allDEGs = do.call(rbind,out)
#     table(allDEGs$geno,allDEGs$ct)
#     allDEGs$geneSym[is.na(allDEGs$geneSym)] = geneMap$geneSym[match(allDEGs$ensID[is.na(allDEGs$geneSym)],geneMap$ensID)]
#     allDEGs$direction = ifelse(allDEGs$log2FoldChange > 0, 'AK_up','AK_down')
#     allDEGs$ct_geno = paste0(allDEGs$ct,":",allDEGs$geno)
#     dim(allDEGs[is.na(allDEGs$geneSym),])
#     
#     ## Genes expressed in <=10% of cells of a Geno_CT is considered NOT expressed
#     allDEGs$max_pct = pmax(allDEGs$cellFrac_g1,allDEGs$cellFrac_g2)
#     allDEGs.sub = allDEGs[allDEGs$max_pct > max_pct_cells,]
#     #allDEGs.sub$ct_geno = gsub('diploid','2n',allDEGs.sub$ct_geno)
#     allDEGs.sub$ct_geno = gsub('complete$','3n',allDEGs.sub$ct_geno)
#     allDEGs.sub$geno[allDEGs.sub$geno == 'complete'] = '3n'
#     
#     
#     
#     # Looks like this step is too memory intesive - haven't been able to get it to work
#     # ## Or at least is epxressed in 10 cells
#     # # Calculate number of cells per Geno per CT epxressing each gene
#     # mtx = big.srat@assays$RNA@counts#[,big.srat$cellID[big.srat$Phase == 'G1']]
#     # big.srat$group = paste0(big.srat$finalAnn,':',big.srat$Genotype)
#     # mtx_byGenoCT = do.call(cbind,lapply(split(colnames(mtx),big.srat$group[match(colnames(mtx),big.srat$cellID)]),function(e) rowSums(as.matrix(mtx[,e] > 0.5))))
#     # mtx_byGenoCT = as.data.frame(mtx_byGenoCT)
#     # mtx_byGenoCT$gene = rownames(mtx_byGenoCT)
#     # mtx_byGenoCT = pivot_longer(mtx_byGenoCT, cols = 1:(ncol(mtx_byGenoCT) - 1),names_to = 'group',values_to = 'nCell_expr')
#     # mtx_byGenoCT$ct = gsub(':.*$','',mtx_byGenoCT$group)
#     # mtx_byGenoCT$group = gsub('complete_trisomy','3n',mtx_byGenoCT$group)
#     # mtx_byGenoCT_2n = mtx_byGenoCT[grepl('2n|diploid',mtx_byGenoCT$group),]
#     # colnames(mtx_byGenoCT_2n)[3] = 'nCell_expr_2n'
#     # 
#     # table(allDEGs.sub$ct_geno %in% mtx_byGenoCT$group)
#     # table(allDEGs.sub$ct %in% mtx_byGenoCT_2n$ct)
#     # #table(mtx_byGenoCT$group[!grepl('2n',mtx_byGenoCT$group) & !mtx_byGenoCT$group %in% allDEGs.sub$ct_geno])
#     # table(allDEGs.sub$geneSym %in% mtx_byGenoCT$gene)
#     # # Add this info to the previous data frame
#     # allDEGs.sub = merge(allDEGs.sub,mtx_byGenoCT,by.x = c('geneSym','ct_geno','ct'),by.y = c('gene','group','ct'),all.x = T)
#     # allDEGs.sub = merge(allDEGs.sub,mtx_byGenoCT_2n[,colnames(mtx_byGenoCT_2n) != 'group'],by.x = c('geneSym','ct'),by.y = c('gene','ct'),all.x = T)
#     # 
#     # allDEGs.sub$max_nCell = pmax(allDEGs.sub$nCell_expr,allDEGs.sub$nCell_expr_2n)
#     # allDEGs.sub = allDEGs.sub[allDEGs.sub$max_nCell > 10,]
#     return(allDEGs.sub)
#   }
#   
# }
# 

##-------------------------------------------##
##    1. Import dataset and DEG output     ####
##-------------------------------------------##
tissue = 'liver'

# Import dataset
big.srat = filter_sratObj(tissue,ageMatch = T,remove_cyclingCells = F)
big.srat = big.srat[[1]]

DimPlot(big.srat,group.by = 'finalAnn',label = T,label.box = T,cols = c(col25,pal34H))+theme(legend.position = 'none')
DimPlot(big.srat,group.by = 'Genotype',label = T,label.box = T,cols = c(col25,pal34H))+theme(legend.position = 'none')

resultDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/jan24/without_published2n/geno/'
# Import DEG output
out = import_pbDEGresults(outDir = resultDir,tissue = 'liver')
head(out)


##----------------------------------------------------------##
##    Plot log2FC of each gene along genomic position     ####
##----------------------------------------------------------##
# Import log2FC output

log2FC_out = import_pbDEGresults(outDir = resultDir, allGenes = T,tissue='liver')
log2FC_out$geneSym = geneMap$geneSym[match(log2FC_out$ensID,geneMap$ensID)]
log2FC_out$chr = geneMap$chr[match(log2FC_out$ensID,geneMap$ensID)]
# log2FC_out$xpos = coords$xpos_linearised[match(log2FC_out$ensID,coords$gene_id)]
# log2FC_out$xpos_byChr = coords$TSS[match(log2FC_out$ensID,coords$gene_id)]
# log2FC_out$isDE = ifelse(is.na(log2FC_out$padj) | log2FC_out$padj > 0.05,F,T) # To extract real DE, I also include criteria of being expressed in >= 10% cells in either group


# a = log2FC_out %>% filter(!is.na(log2FoldChange)) %>% group_by(ct,geno,chr) %>% summarise(med = median(log2FoldChange))
# ggplot(log2FC_out[log2FC_out$ct == 'Endo' & log2FC_out$chr %in% c('18','21','22','X'),],aes(y=log2FoldChange))+
#   geom_point(aes(x=xpos_byChr,col=isDE),size=0.1,alpha=0.5)+
#   scale_color_manual(values = c('grey','red'))+
#   facet_wrap(vars(geno),ncol = 1)+
#   facet_grid(vars(geno),vars(chr),scales = 'free_x',space='free_x')+
#   #geom_vline(xintercept = chrLens_pos[names(chrLens_pos) %in% c('18','21','22','X')],col='grey',lty=2,size=0.3,alpha=0.5) +
#   theme_classic()+ylim(-2,2)
# 
# df = log2FC_out[log2FC_out$ct == 'Hepatocyte',]
# df$geno[df$geno == 'complete'] = 'Triploid'
# df$geno = factor(df$geno,c('T18','T21','T22','MX','Triploid'))
# df$chr2 = ifelse(df$chr %in% c('18','21','22','X'),df$chr,'elsewhere')
# df$chr2 = factor(df$chr2,c('18','21','22','X','elsewhere'))
# 
# df$chr3 = ifelse(as.character(df$chr) %in% as.character(c(1:17)),'1-17',
#                  ifelse(as.character(df$chr) %in% as.character(c(19,20)),'19-20', 
#                         as.character(df$chr)))
# table(df$chr3)
# df$chr3 = factor(df$chr3,c('chr 1 - chr 17','18','chr 19 - chr 20','21','22','X'))
# df$chr3 = factor(df$chr3,c('1-17','18','19-20','21','22','X'))
# df$chr = factor(df$chr,levels = c(as.character(seq(1:22)),'X'))
# 
# 
# ggplot(df,aes(x=chr,y=log2FoldChange))+
#   geom_hline(yintercept = 0,col='grey',lty=2,linewidth=0.3,alpha=1) +
#   #geom_point(aes(col=isDE),size=0.1,alpha=0.5)+
#   geom_jitter(aes(col=isDE),size=0.01,alpha=0.5,width = 0.3)+
#   scale_color_manual(values = c('grey','red'))+
#   facet_wrap(vars(geno),ncol = 1)+
#   facet_grid(vars(geno),vars(chr3),scales = 'free_x',space = 'free_x')+
#   #ylim(-4,4) + 
#   #geom_vline(xintercept = chrLens_pos,col='black',lty=1,size=0.3,alpha=1) +
#   scale_y_continuous(breaks=seq(-4,4,by=4),limits = c(-3,3)) +
#   theme_bw(base_size = 13)+theme(panel.grid = element_blank(),
#                    #axis.text.x = element_blank(),
#                    #axis.ticks.x = element_blank(),
#                    axis.text.x = element_text(size=8),
#                    axis.text.y = element_text(size=8)) + #xlab('Genomic position') + 
#   xlab('Chromosome')+ 
#   ylab('log2 expression fold change (AK / diploid)')
# 


plotDir = file.path(outDir,paste0(tissue,'_plots'))
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}

##---- Heatmap of median log2FC per chromosome per genotype
df = log2FC_out#[log2FC_out$ct %in% c('HSC_MPP','MEMP_MEP'),]
df$geno[df$geno == 'complete'] = 'Triploid'
df$geno = factor(df$geno,c('T18','T21','T22','MX','Triploid'))
df$chr = factor(df$chr,levels = c(as.character(seq(1:22)),'X'))

write.csv(df,file.path(outDir,paste0(tissue,'_log2FC_allGenes.csv')))

data = df %>% filter(ct != 'doublets') %>% group_by(chr,geno,ct) %>% summarise(median_log2FC = mean(log2FoldChange,na.rm = T))
data = pivot_wider(data,id_cols = c('geno','ct'),names_from = 'chr',values_from = 'median_log2FC')

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.15), c('#255fa8','white','#b8393d'))
celltype_order = c('HSC_MPP','MEMP_MEP','EE','ME','LE','Mast.cell','MK',
                   'B.cell.prog','B.cell',
                   'CMP_GMP','proMono','MOP','Monocyte','Macrophage','Kupffer.cell',
                   'DC2','pDC','myelocyte','Myelocyte','NK.T',
                   'Endo','Fibroblast','Hepatocyte')



for(g in unique(df$geno)){
  d = data[data$geno == g,]
  d$ct = factor(d$ct,levels = celltype_order)
  mtx = d[,!colnames(d) %in% c('geno','ct')] %>% as.data.frame()
  rownames(mtx) = paste0(d$ct)
  
  hm = Heatmap(as.matrix(mtx),name = paste0('Average log2FC - ',g),border = F,
          column_split = factor(colnames(mtx),levels = c(1:22,'X')),
          row_split = d$ct,
          show_row_names = T,show_column_names = F,
          column_title_gp = gpar(fontsize=10),
          row_title_gp = gpar(fontsize=0),
          row_names_gp = gpar(fontsize=10),
          row_gap = unit(0.3,'cm'),
          col = col_fun,
          column_names_rot = 0,row_names_side = 'left',
          cluster_rows = F,cluster_columns = F)
  
  pdf(file.path(plotDir,paste0(tissue,'_',g,'_','meanLog2FC_hm.pdf')),width = 8,height = 5.2)
  draw(hm)
  dev.off()
  
  # ## Corrected log2FC for Triploid case
  # if(g %in% c('Triploid','3n')){
  #   mtx2 = mtx - mtx[,'X']
  #   
  #   hm = Heatmap(as.matrix(mtx2),name = paste0('Average log2FC - ',g),border = F,
  #                column_split = factor(colnames(mtx),levels = c(1:22,'X')),
  #                row_split = d$ct,
  #                show_row_names = T,show_column_names = F,
  #                column_title_gp = gpar(fontsize=10),
  #                row_title_gp = gpar(fontsize=0),
  #                row_names_gp = gpar(fontsize=10),
  #                row_gap = unit(0.3,'cm'),
  #                col = col_fun,
  #                column_names_rot = 0,row_names_side = 'left',
  #                cluster_rows = F,cluster_columns = F)
  #   
  #   pdf(file.path(plotDir,paste0(tissue,'_',g,'_','corrected_meanLog2FC_hm.pdf')),width = 8,height = 5.2)
  #   draw(hm)
  #   dev.off()
  # }
}


##Version 2: all Geno, subset for example celltypes
ct_toKeep = c('MEMP_MEP','MK','B.cell.prog','Hepatocyte')
data = data[data$ct %in% ct_toKeep,]

# # Correct log2fc for Triploid
# d = data[data$geno == g,]
# d$ct = factor(d$ct,levels = celltype_order)
# mtx_3n = d[,!colnames(d) %in% c('geno','ct')] %>% as.data.frame()
# rownames(mtx_3n) = paste0(d$ct)
# mtx_3n = mtx_3n - mtx_3n[,'X']
# 
# data2 = data[data$geno != g,]
# mtx_3n$geno = g
# mtx_3n$ct = d$ct[d$geno == g]
# mtx_3n = as.data.frame(mtx_3n)
# mtx_3n = mtx_3n[,c('geno','ct',colnames(mtx_3n)[!colnames(mtx_3n) %in% c('geno','ct')])]

# data2 = rbind(data2,mtx_3n)
# data = data2
data$ct = factor(data$ct,levels = ct_toKeep)
mtx = data[,!colnames(data) %in% c('geno','ct')] %>% as.data.frame()
rownames(mtx) = paste0(data$ct,':',data$geno)

hm=Heatmap(as.matrix(mtx),name = paste0('Average log2FC'),border = F,
             column_split = factor(colnames(mtx),levels = c(1:22,'X')),
             row_split = data$ct,
             show_row_names = T,show_column_names = F,
             column_title_gp = gpar(fontsize=11),
             row_title_gp = gpar(fontsize=0),
             row_names_gp = gpar(fontsize=10),
             row_gap = unit(0.4,'cm'),
             col = col_fun,
             column_names_rot = 0,row_names_side = 'left',
             cluster_rows = F,cluster_columns = F)

pdf(file.path(plotDir,paste0(tissue,'_allGeno_subsetCT_meanLog2FC_hm.pdf')),width = 9,height = 5.2)
draw(hm)
dev.off()

# ##----------------------------------------------##
# ##    Some other random investigations....    ####
# ##----------------------------------------------##
# 
# ggplot(log2FC_out[log2FC_out$ct == 'Endo' & log2FC_out$chr %in% c('18','21','22','X'),],aes(x=chr,y=log2FoldChange))+
#   geom_boxplot(outlier.size = 0.001,outlier.color = 'white',width=0.7)+
#   scale_color_manual(values = c('grey','red'))+
#   facet_wrap(vars(geno),ncol = 1)+
#   #facet_grid(vars(geno),vars(chr),scales = 'free_x',space='free_x')+
#   #geom_vline(xintercept = chrLens_pos[names(chrLens_pos) %in% c('18','21','22','X')],col='grey',lty=2,size=0.3,alpha=0.5) +
#   theme_classic()+ylim(-0.5,0.5)
# 
# 
# df = log2FC_out[log2FC_out$geno == 'T21' & 
#                   log2FC_out$ct %in% c('MEMP_MEP','Hepatocyte') & #log2FC_out$chr %in% c('1','2') &
#                   #log2FC_out$isDE ==F & 
#                   !is.na(log2FC_out$log2FoldChange),]
# ggplot(df,aes(x = chr,y=log2FoldChange))+
#   geom_boxplot(outlier.colour = 'white')+
#   geom_hline(yintercept = 0)+
#   facet_grid(vars(geno),vars(chr),scales = 'free_x',space='free_x') +
#   theme_classic()+ylim(-0.2,0.2)
# 
# 
# df2 = df[,c('log2FoldChange','ct','geneSym','isDE')]
# df2 = pivot_wider(df2,id_cols = c('geneSym'),names_from = 'ct',values_from = 'log2FoldChange')
# df2$isDE_MEMP = ifelse(df2$MEMP_MEP > 0,
#                        ifelse(df2$geneSym %in% out$geneSym[out$ct_geno == 'MEMP_MEP:T21' & out$direction == 'AK_up'],T,F),
#                        ifelse(df2$geneSym %in% out$geneSym[out$ct_geno == 'MEMP_MEP:T21' & out$direction == 'AK_down'],T,F))
# df2$isDE_Hep = ifelse(df2$Hepatocyte > 0,
#                       ifelse(df2$geneSym %in% out$geneSym[out$ct_geno == 'Hepatocyte:T21' & out$direction == 'AK_up'],T,F),
#                       ifelse(df2$geneSym %in% out$geneSym[out$ct_geno == 'Hepatocyte:T21' & out$direction == 'AK_down'],T,F))
# ggplot(df2,aes(Hepatocyte,MEMP_MEP,col=isDE))+
#   geom_point(size=0.01,alpha=0.5)+
#   geom_abline(col='grey',size=0.4)
# 
# Idents(big.srat) = big.srat$finalAnn
# big.srat$tmp = paste0(big.srat$Genotype,':',big.srat$ageGroup,':',big.srat$donorID)
# MEMP_MEP
# DotPlot(big.srat,idents= c('Hepatocyte'),group.by = 'tmp',
#         features = c('FAM89A'
#         )) + RotatedAxis() +
#   theme(axis.text.x = element_text(size=8))
# 
# 
# 
# 
# ## Top 50 DEGs in each genotype --> score in different lineages
# 
# topGenes_up = df[df$inPB == T & df$pct.diff > 0.1 & df$avg_log2FC > 0.5,]
# topGenes_down = df[df$inPB == T & df$pct.diff > 0.1 & df$avg_log2FC < -0.5,]
# topDEG = df[df$pct.diff > 0.1 & abs(df$avg_log2FC) > 0.5 & !df$gene %in% nCT_toRemove$gene,]
# 
# topDEG = rbind(topGenes_up,topGenes_down)
# topDEG$ensID = geneMap$ensID[match(topDEG$gene,geneMap$geneSym)]
# topDEG = annotateGenes(topDEG,geneMap = geneMap)
# topDEG = topDEG[!topDEG$gene %in% nonsenseGenes,]
# View(topDEG[!topDEG$gene %in% nonsenseGenes,])
# table(topGenes_up$geno)
# DotPlot(big.srat,idents = c('HSC_MPP','MEMP_MEP','CMP_GMP','LMPP_ELP','B.cell.prog','Hepatocyte','B.cell','MK'),
#         group.by = 'ann',features = unique(c(#topDEG$gene[topDEG$nGeno == 5 & topDEG$direction == 'AK_up'],
#           #topDEG$gene[topDEG$nGeno == 4 & topDEG$direction == 'AK_up'],
#           #topDEG$gene[topDEG$nGeno == 3 & topDEG$direction == 'AK_up'],
#           topDEG$gene[topDEG$nGeno == 2 & topDEG$direction == 'AK_down']
#           #topDEG$gene[topDEG$nGeno == 4 & topDEG$direction == 'AK_down']
#         ))) + RotatedAxis() + 
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))
# 
# pdf('~/lustre_mt22/Aneuploidy/Results/14_primitive_haematopoiesis/HSC_DEGs_DotPlot.pdf',width = 12,height = 5)
# 
# big.srat$ann2 = big.srat$ann
# big.srat$ann2 = gsub('complete_trisomy','triploid',big.srat$ann2)
# big.srat$ann2 = factor(big.srat$ann2,c('HSC_MPP/diploid','HSC_MPP/MX','HSC_MPP/T21','HSC_MPP/T18','HSC_MPP/T22','HSC_MPP/triploid',
#                                        unique(big.srat$ann2[!big.srat$ann2 %in% c('HSC_MPP/diploid','HSC_MPP/MX','HSC_MPP/T21','HSC_MPP/T18','HSC_MPP/T22','HSC_MPP/triploid')])))
# DotPlot(big.srat,idents = c('HSC_MPP'),cols = c(grey(0.8),'#2D4372'),
#         group.by = 'ann2',features = c(unique(topGenes_up$gene[topGenes_up$nGeno == 5]),
#                                        unique(topGenes_up$gene[topGenes_up$nGeno == 4]),
#                                        unique(topGenes_up$gene[topGenes_up$nGeno == 3]),
#                                        unique(topGenes_up$gene[topGenes_up$nGeno == 2]),
#                                        unique(topGenes_up$gene[topGenes_up$nGeno == 1]))) + RotatedAxis() + 
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7),legend.position = 'top')+xlab('')+ylab('')
# 
# dev.off()
# genes_toExclude = c('APOB','ALB','GNAS','TF','IGF2','JUND')
# 
# big.srat$lineage = big.srat$finalAnn
# big.srat$lineage[big.srat$finalAnn %in% c('CMP_GMP','MOP','Macrophage','Kupffer.cell','Monocyte','proMono','DC1','DC2','pDC','myelocyte')] = 'Myeloid lineage'
# big.srat$lineage[big.srat$finalAnn %in% c('B.cell.prog','B.cell')] = 'B lineage'
# big.srat$lineage[big.srat$finalAnn %in% c('ILC.precursor','NK.T')] = 'NK/T lineage'
# big.srat$lineage[big.srat$finalAnn %in% c('EE','ME','LE','MEMP_MEP','MK','Mast.cell')] = 'EE/MK/Mast lineage'
# big.srat$lineage[big.srat$finalAnn %in% c('Endo','Fibroblast','Hepatocyte','Mesenchyme','NPC')] = 'Stromal'
# 
# big.srat$lineage = factor(big.srat$lineage,c('HSC_MPP','EE/MK/Mast lineage','Myeloid lineage','B lineage','NK/T lineage','Stromal'))
# 
# 
# avgExp = AverageExpression(big.srat,group.by = 'ann',features = unique(topDEG$gene))
# avgExp = avgExp$RNA
# library(ComplexHeatmap)
# mtx = avgExp[unique(topDEG$gene[topDEG$direction == 'AK_down']),grepl('diploid',colnames(avgExp))]
# colnames(mtx) = gsub('/diploid','',colnames(mtx))
# library(circlize)
# q = quantile(mtx,seq(0,1,0.25))
# genes = c('RUNX1','LMO2','CEBPA','CEBPB','ETV6','CEBPA','EP300','LYN','LMO4','IGF1R','IGF2')
# rowAnno = rowAnnotation(mark = anno_mark(at = which(rownames(mtx) %in% genes),
#                                          labels = rownames(mtx)[which(rownames(mtx) %in% 
#                                                                         c('RUNX1','LMO2','CEBPA','CEBPB','ETV6','CEBPA','EP300','LYN','LMO4','IGF1R','IGF2'))],
#                                          which = 'row',labels_gp=gpar(fontsize=5)))
# 
# 
# 
# pdf('~/lustre_mt22/Aneuploidy/Results/14_primitive_haematopoiesis/HSC_DEGs_expressionIn2n.pdf',width = 7,height = 8)
# hm = Heatmap(t(scale(t(mtx))),name = 'Average expression',
#              show_row_dend = F,show_column_dend = F,show_row_names = T,show_column_names = T,
#              cluster_rows = T,cluster_columns = T,cluster_column_slices = F,km=4,
#              #col=colorRamp2(q, rev(c('#162668','#1A2D7A','#3A4B8D','#7C87B3','#9CA5C6','#BDC3D9','#DEE1EC'))),
#              col = colorRamp2(c(-4, 0, 4), c('white','#FED976','#E31A1C')),
#              column_split = big.srat$lineage[match(colnames(mtx),big.srat$finalAnn)],
#              column_title_rot = 90,column_title_gp = gpar(fontsize=17),
#              #right_annotation = rowAnno,
#              row_names_gp = gpar(fontsize=6),column_names_gp = gpar(fontsize=12),column_names_rot = 90)
# ht = draw(hm)
# nonsenseGenes = rownames(mtx)[row_order(ht)[['4']]]
# dev.off()
# 
# 
# geno_cols = c('2n' = grey(0.6),
#               'T21' = '#68389A',#A487C2
#               'T18' = '#EAB38A',
#               'T22' = '#AC872D',
#               'T13' = '#4E73BE',
#               'MX' = '#5E803F',
#               '3n' = '#B02318')
# df$geno[df$geno == '3n'] = 'Triploid'
# df$geno = factor(df$geno,levels = c('MX','T21','T18','T22','Triploid'))
# ggplot(df[abs(df$avg_log2FC) > 0.5,],aes(avg_log2FC,-log10(p_val_adj),col=geno))+
#   geom_point(size=0.3,alpha=0.4,col=grey(0.8)) +
#   geom_point(data = df[df$gene %in% genes[1],],aes(avg_log2FC,-log10(p_val_adj),col=geno),size=0.8) +
#   geom_vline(xintercept = c(0.5,-0.5),lty=2,lwd=0.3)+
#   #xlim(-0.5,max(topGenes_up$avg_log2FC)) +
#   scale_color_manual(values =  geno_cols) +
#   theme_bw() + facet_wrap(vars(geno)) + theme(panel.grid = element_blank())
# 
# pdf('~/lustre_mt22/Aneuploidy/Results/14_primitive_haematopoiesis/HSC_DEGs_pct.pdf',width = 9,height = 6)
# p = ggplot(df[abs(df$avg_log2FC) > 0.5,],aes(pct.2,pct.1,col=avg_log2FC))+
#   geom_point(size=0.6,alpha=0.8,col=grey(0.8)) +
#   geom_point(data = df[df$gene %in% genes,],aes(pct.2,pct.1,col=avg_log2FC),size=1) +
#   geom_text(data = df[df$gene %in% genes,],
#             label=df[df$gene %in% genes,]$gene, 
#             nudge_x = 0.1, nudge_y = 0.1, size=4,
#             check_overlap = T)+
#   geom_vline(xintercept = c(0.5),lty=2,lwd=0.3,alpha=0.5)+
#   geom_hline(yintercept = c(0.5),lty=2,lwd=0.3,alpha=0.5)+
#   xlab('Fraction of diploid cells') +
#   ylab('Fraction of AK cells') +
#   #xlim(-0.5,max(topGenes_up$avg_log2FC)) +
#   #scale_color_manual(values =  geno_cols) +
#   theme_bw(base_size = 18) + facet_wrap(vars(geno)) + theme(panel.grid = element_blank(),axis.text = element_text(size=11))
# 
# print(p)
# 
# dev.off()
# 
# 
# ## Score for the expression of these genes
# 
# 
# ##------------------------
# 
# ##      ME / LE   ##
# library(ComplexHeatmap)
# 
# celltype = 'ME'
# df = allDEGs.sub[allDEGs.sub$ct %in% c('ME','LE'),]
# df$group = paste0(df$geno,':',df$direction)
# 
# # Top genes up-regulated: expressed in >=50% AK but  <= 50% 2n
# df$fracDiff = abs(df$cellFrac_g1 - df$cellFrac_g2)
# df = df %>% group_by(geneSym,direction) %>% mutate(nGeno = n_distinct(geno))
# topDEG = df[df$fracDiff > 10,]
# topDEG$gene = topDEG$geneSym
# 
# 
# 
# 
# 
# avgExp = AverageExpression(big.srat,group.by = 'ann',features = unique(topDEG$geneSym))
# avgExp = avgExp$RNA
# 
# mtx = avgExp[nonsenseGenes,grepl('diploid',colnames(avgExp))]
# colnames(mtx) = gsub('/diploid','',colnames(mtx))
# library(circlize)
# q = quantile(mtx,seq(0,1,0.25))
# 
# 
# 
# pdf('~/lustre_mt22/Aneuploidy/Results/14_primitive_haematopoiesis/HSC_DEGs_expressionIn2n.pdf',width = 7,height = 8)
# 
# 
# hm = Heatmap(t(scale(t(mtx))),name = 'Average expression',
#              show_row_dend = F,show_column_dend = F,show_row_names = T,show_column_names = T,
#              cluster_rows = T,cluster_columns = T,cluster_column_slices = F,km=4,
#              #col=colorRamp2(q, rev(c('#162668','#1A2D7A','#3A4B8D','#7C87B3','#9CA5C6','#BDC3D9','#DEE1EC'))),
#              col = colorRamp2(c(-4, 0, 4), c('white','#FED976','#E31A1C')),
#              column_split = big.srat$lineage[match(colnames(mtx),big.srat$finalAnn)],
#              column_title_rot = 90,column_title_gp = gpar(fontsize=17),
#              #right_annotation = rowAnno,
#              row_names_gp = gpar(fontsize=6),column_names_gp = gpar(fontsize=12),column_names_rot = 90)
# ht = draw(hm)
# nonsenseGenes = rownames(mtx)[row_order(ht)[['1']]]
# 
# 
# ## EnrichR
# up <- enrichr(unique(topDEG$gene[topDEG$direction == 'AK_up' & topDEG$nGeno >= 3 & topDEG$avg_log2FC >= 0.5 & topDEG$pct.2<0.5]), dbs)
# plotEnrich(up[[5]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
# 
# down <- enrichr(unique(topDEG$gene[topDEG$direction == 'AK_down' & topDEG$nGeno >= 1 & topDEG$log2FoldChange <= -0.1]), dbs)
# plotEnrich(down[[3]], showTerms = 50, numChar = 40, y = "Count", orderBy = "P.value")
# View(down[[6]])
# 
# 
# 
# 
# DotPlot(big.srat,idents = c('HSC_MPP','CMP_GMP','EE','ME','LE','Hepatocyte','B.cell','MK'),
#         group.by = 'ann',features = unique(c(#topDEG$gene[topDEG$nGeno == 5 & topDEG$direction == 'AK_up'],
#           #topDEG$gene[topDEG$nGeno == 4 & topDEG$direction == 'AK_up'],
#           #topDEG$gene[topDEG$nGeno == 3 & topDEG$direction == 'AK_up'],
#           #topDEG$gene[topDEG$direction == 'AK_up' & topDEG$nGeno >= 3 & topDEG$log2FoldChange >= 0.5]
#           'SLC25A37','SLC6A9',
#           strsplit('SLC22A4;GYPA;ALAS2;EPB42;GYPB;SLC2A1;CPOX;FOXO3;CLCN3;NUDT4;TRAK2;SELENBP1;RNF19A;MXI1;TSPAN5;
#                    BNIP3L;HBZ;LMO2;NCOA4;XPO7;BPGM;HMBS;TRIM58',split=';')[[1]]
#           #topDEG$gene[topDEG$nGeno == 4 & topDEG$direction == 'AK_down']
#         ))) + RotatedAxis() + 
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))
# 
# Idents(big.srat) = big.srat$annot_sept23
# DotPlot(big.srat,idents = c('HSC_MPP','CMP_GMP','EE','ME','LE','Hepatocyte','B.cell','MK'),
#         group.by = 'ann',features = unique(#c(topDEG$gene[topDEG$nGeno == 5 & topDEG$direction == 'AK_up' & topDEG$isTF == T]
#           genes
#           #topDEG$gene[topDEG$nGeno == 4 & topDEG$direction == 'AK_up'],
#           #topDEG$gene[topDEG$nGeno == 3 & topDEG$direction == 'AK_up'],
#           #topDEG$gene[topDEG$direction == 'AK_up' & topDEG$nGeno >= 3 & topDEG$log2FoldChange >= 0.5]
#           
#           #topDEG$gene[topDEG$nGeno == 4 & topDEG$direction == 'AK_down']
#         )) + RotatedAxis() + 
#   theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1,size=7))
# 
# 
# 
# 
# df = df[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',df$gene),]
# topDEG = topDEG[!grepl('^RPL|^RPS|^MT-|^MALAT1|^NEAT1|^GAPDH|^GABARAP|^MRPL|^MRPS',topDEG$gene),]
# 
# geno_cols = c('2n' = grey(0.6),
#               'T21' = '#68389A',#A487C2
#               'T18' = '#EAB38A',
#               'T22' = '#AC872D',
#               'T13' = '#4E73BE',
#               'MX' = '#5E803F',
#               '3n' = '#B02318')
# df$geno[df$geno == '3n'] = 'Triploid'
# df$geno = factor(df$geno,levels = c('MX','T21','T18','T22','Triploid'))
# genes = c('KCNH2','EP300',
#           'GPC3','SLC2A1','SLC14A1','KLF13','MXD4',
#           'FOXO3','TSPAN5','HBZ','BPGM','BMP2K')
# df$geno = factor(df$geno,c('MX','T21','T18','T22','Triploid'))
# df$log2FoldChange = df$avg_log2FC
# df$cellFrac_g1 = df$pct.1*10
# df$cellFrac_g2 = df$pct.2*10
# df$geneSym = df$gene
# genes = genes[genes %in% df[(df$log2FoldChange) > 0.5 & df$pct.diff > 0.1 & df$pct.1 > df$pct.2,]$gene]
# pdf('~/lustre_mt22/Aneuploidy/Results/14_primitive_haematopoiesis/HSC_DEGs_pct.pdf',width = 9,height = 6)
# p = ggplot(df[(df$log2FoldChange) > 0.5 & df$pct.diff > 0.1 & df$pct.1 > df$pct.2,],aes(cellFrac_g2/10,cellFrac_g1/10,col=log2FoldChange))+
#   geom_point(size=0.6,alpha=0.8,col=grey(0.8)) +
#   geom_point(data = df[df$geneSym %in% genes,],aes(cellFrac_g2/10,cellFrac_g1/10,col=log2FoldChange),size=1) +
#   geom_point(data = df[df$geneSym %in% c('HBZ'),],aes(cellFrac_g2/10,cellFrac_g1/10,col=log2FoldChange),size=3,col='red') +
#   geom_text(data = df[df$geneSym %in% c('HBZ'),],
#             label=df[df$geneSym %in% c('HBZ'),]$geneSym,
#             nudge_x = 0.05, nudge_y = 0.05, size=5,col='darkred',
#             check_overlap = T)+
#   # geom_point(data = df[df$geneSym %in% c(genes,'HBZ'),],aes(cellFrac_g2/10,cellFrac_g1/10,col=log2FoldChange),size=1) +
#   # geom_text(data = df[df$geneSym %in% c(genes,'HBZ'),],
#   #           label=df[df$geneSym %in% c(genes,'HBZ'),]$geneSym,
#   #           nudge_x = 0.05, nudge_y = 0.05, size=4,
#   #           check_overlap = T)+
#   geom_vline(xintercept = c(0.5),lty=2,lwd=0.3,alpha=0.5)+
#   geom_hline(yintercept = c(0.5),lty=2,lwd=0.3,alpha=0.5)+
#   xlab('Fraction of diploid cells') +
#   ylab('Fraction of AK cells') +
#   #xlim(-0.5,max(topGenes_up$avg_log2FC)) +
#   #scale_color_manual(values =  geno_cols) +
#   theme_bw(base_size = 18) + facet_wrap(vars(geno)) + theme(panel.grid = element_blank(),axis.text = element_text(size=11))
# 
# print(p)
# 
# dev.off()
# 
# #####
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##----------------------------------------
# ### DGE analysis by celltype by genotype ###
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##===============       Plot Figure 1 Heatmap of mean_log2FC of genes per Chr per Cell type     =================##
# ##---------------------------##
# ##    Log2FC - fLiver      ####
# ##---------------------------##
# # Import log2FC output
# log2FC_out = import_pbDEGresults(outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/oct23_ageMatched',allGenes = T,tissue = 'liver')
# log2FC_out$geneSym = geneMap$geneSym[match(log2FC_out$ensID,geneMap$ensID)]
# log2FC_out$chr = geneMap$chr[match(log2FC_out$ensID,geneMap$ensID)]
# log2FC_out$xpos = coords$xpos_linearised[match(log2FC_out$ensID,coords$gene_id)]
# log2FC_out$xpos_byChr = coords$TSS[match(log2FC_out$ensID,coords$gene_id)]
# log2FC_out$isDE = ifelse(is.na(log2FC_out$padj) | log2FC_out$padj > 0.05,F,T)
# 
# 
# 
# ##----------------------------##
# ##    Log2FC - fKidney      ####
# ##----------------------------##
# kid_log2FC_out = import_pbDEGresults(outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/oct23_ageMatched',allGenes = T,tissue = 'kidney')
# kid_log2FC_out$geneSym = geneMap$geneSym[match(kid_log2FC_out$ensID,geneMap$ensID)]
# kid_log2FC_out$chr = geneMap$chr[match(kid_log2FC_out$ensID,geneMap$ensID)]
# kid_log2FC_out$xpos = coords$xpos_linearised[match(kid_log2FC_out$ensID,coords$gene_id)]
# kid_log2FC_out$xpos_byChr = coords$TSS[match(kid_log2FC_out$ensID,coords$gene_id)]
# kid_log2FC_out$isDE = ifelse(is.na(kid_log2FC_out$padj) | kid_log2FC_out$padj > 0.05,F,T)
# 
# 
# table(kid_log2FC_out$ct,kid_log2FC_out$geno)
# 
# ##-----------------------------##
# ##    Log2FC - fAdrenal      ####
# ##-----------------------------##
# adr_log2FC_out = import_pbDEGresults(outDir = '~/lustre_mt22/Aneuploidy/Results/3e_pseudoBulk_DESeq2_withCyclCells/oct23_ageMatched',allGenes = T,tissue = 'adrenal')
# adr_log2FC_out$geneSym = geneMap$geneSym[match(adr_log2FC_out$ensID,geneMap$ensID)]
# adr_log2FC_out$chr = geneMap$chr[match(adr_log2FC_out$ensID,geneMap$ensID)]
# adr_log2FC_out$xpos = coords$xpos_linearised[match(adr_log2FC_out$ensID,coords$gene_id)]
# adr_log2FC_out$xpos_byChr = coords$TSS[match(adr_log2FC_out$ensID,coords$gene_id)]
# adr_log2FC_out$isDE = ifelse(is.na(adr_log2FC_out$padj) | adr_log2FC_out$padj > 0.05,F,T)
# 
# table()
# 
# 
# ##-----------------------------##
# ##    Log2FC - Heatmap       ####
# ##-----------------------------##
# View(table(adr_log2FC_out$ct,adr_log2FC_out$geno))
# example_celltype = c('MK','ICb','med_Sympathoblast')
# 
# ##---- Heatmap of mean log2FC per chromosome per genotype
# df = rbind(log2FC_out[log2FC_out$ct %in% c('MK'),],
#            rbind(kid_log2FC_out[kid_log2FC_out$ct %in% c('NPC'),],
#                  adr_log2FC_out[adr_log2FC_out$ct %in% c('med_Sympathoblast'),]))
# 
# df$geno[df$geno == 'complete'] = 'Triploid'
# df = df[df$geno != 'diploidTest',]
# #df$geno = factor(df$geno,c('T18','T21','T22','MX','Triploid'))
# df$chr = factor(df$chr,levels = c(as.character(seq(1:22)),'X'))
# 
# data = df %>% group_by(chr,geno,ct) %>% summarise(mean_log2FC = mean(log2FoldChange,na.rm = T))
# data = pivot_wider(data,id_cols = c('geno','ct'),names_from = 'chr',values_from = 'mean_log2FC')
# data$ct = factor(data$ct,c('MK','NPC','med_Sympathoblast'))
# 
# mtx = data[,!colnames(data) %in% c('geno','ct')]
# rownames(mtx) = paste0(data$geno,':',data$ct)
# library(ComplexHeatmap)
# library(circlize)
# col_fun = colorRamp2(c(-0.1, 0, 0.06,0.15), c('#255fa8','white','#e8c5c6','#a4282c'))
# 
# Heatmap(as.matrix(mtx),name = 'Average log2FC',border = F,
#         column_split = factor(colnames(mtx),levels = c(1:22,'X')),
#         row_split = data$ct,
#         show_row_names = T,show_column_names = F,
#         column_title_gp = gpar(fontsize=10),
#         row_names_gp = gpar(fontsize=10),
#         row_gap = unit(0.3,'cm'),
#         col = col_fun, 
#         column_names_rot = 0,row_names_side = 'left',
#         cluster_rows = F,cluster_columns = F)
# 
# 
# 
# ##---- Heatmap of median log2FC per chromosome per genotype for ALL cells
# log2FC_out$tissue = 'fLiver'
# kid_log2FC_out$tissue = 'fKidney'
# adr_log2FC_out$tissue = 'fAdrenal'
# df = rbind(log2FC_out,
#            rbind(kid_log2FC_out,
#                  adr_log2FC_out))
# 
# 
# 
# df$geno[df$geno == 'complete'] = 'Triploid'
# df = df[df$geno != 'diploidTest',]
# #df$geno = factor(df$geno,c('T18','T21','T22','MX','Triploid'))
# df$chr = factor(df$chr,levels = c(as.character(seq(1:22)),'X'))
# 
# data = df %>% group_by(chr,geno,ct) %>% summarise(mean_log2FC = mean(log2FoldChange,na.rm = T))
# 
# ## Reset Triploid log2FC
# triploid_x = data[data$geno == 'Triploid' & data$chr == 'X',]
# data$chrX_3n = triploid_x$mean_log2FC[match(data$ct,triploid_x$ct)]
# data$mean_log2FC_v2 = ifelse(data$geno == 'Triploid',data$mean_log2FC - data$chrX_3n,data$mean_log2FC )
# data$mean_log2FC = data$mean_log2FC_v2
# 
# data = pivot_wider(data,id_cols = c('geno','ct'),names_from = 'chr',values_from = 'mean_log2FC')
# #data$ct = factor(data$ct,c('MK','NPC','med_Sympathoblast'))
# 
# 
# 
# 
# mtx = data[,!colnames(data) %in% c('geno','ct')]
# rownames(mtx) = paste0(data$geno,':',data$ct)
# library(ComplexHeatmap)
# library(circlize)
# col_fun = colorRamp2(c(-0.1, 0, 0.06,0.15), c('#255fa8','white','#e8c5c6','#a4282c'))
# 
# Heatmap(as.matrix(mtx[data$geno == 'Triploid',]),name = 'Average log2FC',border = F,
#         column_split = factor(colnames(mtx),levels = c(1:22,'X')),
#         row_split = data$ct[data$geno == 'Triploid'],
#         show_row_names = T,show_column_names = F,
#         column_title_gp = gpar(fontsize=10),
#         row_names_gp = gpar(fontsize=10),
#         row_title_gp = gpar(fontsize=10),row_title_rot = 0,
#         row_gap = unit(0.3,'cm'),
#         col = col_fun, 
#         column_names_rot = 0,row_names_side = 'left',
#         cluster_rows = F,cluster_columns = F)
# 
# 
