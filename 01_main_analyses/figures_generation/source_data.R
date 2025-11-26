# Generate Source data

library(openxlsx)

figure_generation_location = c(
  fig1A = 'figures_final.R',
  fig1B = 'figures_final.R',
  fig1C = 'figures_final.R',
  fig1D = 'figures_final.R',
  fig1E = '01_main_analyses/01.2_fetal_aneuploidy_analyses/3.2_DEanalysis_plots.R',
  fig1F = '01_main_analyses/01.2_fetal_aneuploidy_analyses/3.2_DEanalysis_plots.R',
  
  fig2A = 'figures_final.R',
  fig2B = 'figures_final.R',
  fig2C = 'figures_final.R',
  fig2D = '~/ML-DS/01_main_analyses/04_derive_transcriptional_modules/x3_MLDS_imprint_in_fLiver_T21_v2.R',
  fig2E = '~/ML-DS/01_main_analyses/04_derive_transcriptional_modules/x4_GATA1s_module.R',
  fig2F = '~/ML-DS/01_main_analyses/04_derive_transcriptional_modules/x6_TAM_vs_MLDS.R',
  
  fig3A = '~/ML-DS/01_main_analyses/figures_generation/figures_final_filtered.R',
  fig3B = '~/ML-DS/01_main_analyses/05_transcriptional_modules_specificity/xx01_moduleScore_terminal.R',
  fig3C = '~/ML-DS/01_main_analyses/05_transcriptional_modules_specificity/xx01_moduleScore_terminal.R',
  fig3D = '~/ML-DS/01_main_analyses/05_transcriptional_modules_specificity/xx01_moduleScore_terminal.R',
  
  fig4A = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l076_relapse.R',
  fig4B = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l076_deepdive.R',
  #fig4C = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l076_relapse.R',
  fig4D = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l038_relapse.R',
  fig4E = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l038_deepdive.R',
  
  figS1 = 'figures_final.R',
  figS2 = '~/ML-DS/01_main_analyses/01.2_fetal_aneuploidy_analyses/3.2_DEanalysis_plots.R',
  figS3 = 'figures_final.R',
  figS4 = 'figures_final.R',
  figS5 = '~/ML-DS/01_main_analyses/04_derive_transcriptional_modules/x4_GATA1s_module.R',
  figS6 = 'figures_final.R',
  figS7 = c('~/ML-DS/01_main_analyses/05_transcriptional_modules_specificity/xx01_moduleScore_terminal.R',
            '~/ML-DS/01_main_analyses/05_transcriptional_modules_specificity/xx01_moduleScore_bulkRNA_helperFunctions.R'),
  figS8 = '',
  figS9 = c('figures_final.R',
            '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l076_deepdive.R')
  )

# Load in raw data for each figure --------
## Figure 1 ------
fig1D <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure1/Fig1_2nAKLiv_mean.celltypeFractionByDonorByLineage_barPlot_rawData.tsv',sep = '\t')
fig1D <- fig1D %>% 
  dplyr::rename('Lineage'='broadLineage3',
                'cell_proportion'='ctFrac',
                Donor_ID = donorID,
                n_cell = nCell,
                total_cell = totalCell) %>% 
  dplyr::filter(Lineage %in% c('B lineage', 'Ery/MegK/Mast lineage')) %>% 
  dplyr::mutate(Karyotype = factor(Genotype,c('Diploid','T21','T18','T22','Triploid','MX')),
                Lineage = factor(Lineage, c('Ery/MegK/Mast lineage', 'B lineage'))) %>% 
  dplyr::select(Lineage, Karyotype, Donor_ID, n_cell, total_cell, cell_proportion) %>% 
  dplyr::arrange(Lineage, Karyotype, cell_proportion)

fig1E <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure1/Fig1_2nAKLiver_GlobalImpact_DEGs_hor_rawData.tsv',sep = '\t')
fig1E <- fig1E %>% 
  dplyr::rename(Karyotype = Genotype,
                cell_type = Celltype,
                global_impact = globalImpact) %>% 
  dplyr::select(-c(nCell_2n,nCell_AK))
  

fig1F <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure1/FigSup2X_nSharedDEGs_acrossGeno_rawData.tsv',sep = '\t')
fig1F <- fig1F %>% 
  dplyr::select(geno, direction, ct, nGeno, nGene) %>% 
  dplyr::filter(geno == 'T21') %>% 
  dplyr::mutate(direction = ifelse(direction == 'AK_up','AK up-regulated','AK down-regulated'),
                direction = factor(direction, c('AK up-regulated','AK down-regulated'))) %>% 
  dplyr::rename('Direction' = 'direction',
                Karyotype = 'geno',
                'Cell type' = 'ct',
                'Number of other AK-vs-diploid comparisons' = nGeno,
                'Number of DEG' = 'nGene') %>% 
  dplyr::arrange(Direction, `Number of other AK-vs-diploid comparisons`, desc(`Number of DEG`))

## Figure 2 ------
fig2B <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots/Fig2_MLDS_TumNorm_UMAP_rawData.tsv',sep = '\t')
fig2C <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots/Fig2_MLDS_GATA1s_UMAP_v2_rawData.tsv',sep = '\t')
shared_columns = intersect(colnames(fig2B), colnames(fig2C[,!colnames(fig2C) %in% c('broadLineage')]))
fig2bc <- dplyr::full_join(fig2B, fig2C[,!colnames(fig2C) %in% c('broadLineage')], by = shared_columns)
checkmate::assert_true(nrow(fig2bc) == nrow(fig2B))

fig2bc <- fig2bc %>% dplyr::select(cellID, orig.ident, donorID, broadLineage, annot, GATA1s_status2,    
                                   disease, timePoint, tissue, clinicalOutcome, UMAP_1, UMAP_2) %>% 
  dplyr::arrange(broadLineage, annot, GATA1s_status2, donorID, tissue, timePoint, orig.ident,clinicalOutcome) %>% 
  dplyr::rename(cell_type_annotation = annot,
                '10X_sample_id'=orig.ident,
                GATA1_mutation_status = GATA1s_status2)

fig2d <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure2/Fig3A_MLDS.T21.DEGs_barplot_simplified_sameScaleGATA1smodule_rawData.tsv',sep = '\t')
fig2e <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure2/Fig3B_MLDS.goodTAM.DEGs_barplot_simplified_rawData.tsv',sep = '\t')
fig2f <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure2/Fig3C_gTAM.allMLDS_barplot_rawData.tsv',sep = '\t')
fig2def <- do.call(rbind,list(fig2d,fig2e,fig2f))
fig2def <- fig2def %>% dplyr::mutate(module = dplyr::case_when(group == 'T21 vs 2n_MEMP_fLiver' ~ 'Trisomy 21-leukaemia',
                                                               group == 'TAM vs T21_MEMP' ~ 'GATA1-leukaemia',
                                                               group == 'MLDS_vs_canonicalTAM' ~ 'ML-DS',
                                                               .default = '?'),
                                     n_DEG = abs(nGene),
                                     comparison = group) %>% 
  dplyr::select(module,comparison,direction,n_DEG) 


## Figure 3 ------
fig3C <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure3/FigS7x_medianScorePerGroup_MLDS_topGenes_2410_cf2_rawData.tsv',sep = '\t')
fig3C <- fig3C %>% dplyr::mutate(karyotype = Genotype) %>% 
  dplyr::select(group,dataset_ID,karyotype,medianScore) %>% 
  dplyr::mutate(group = factor(group,c("MEMP_MEP","TAM / MLDS","Other leukaemia")),
                dataset_ID = dplyr::case_when(dataset_ID == 'FL' ~ 'Fetal liver',
                                              dataset_ID == 'FBM' ~ 'Fetal BM',
                                              dataset_ID == 'TAM' ~ 'Conventional TAM',
                                              dataset_ID == 'MLDS' ~ 'ML-DS',
                                              dataset_ID == 'pAML' ~ 'AML',
                                              dataset_ID == 'pBALL' ~ 'B-ALL',
                                              .default = dataset_ID),
                dataset_ID = factor(dataset_ID,c('Fetal liver',"Fetal BM", 
                                                 "Conventional TAM","Recurrent TAM","ML-DS",
                                                 "MDS", "AMKL", "AEL","AML","B-ALL", "LPD"))) %>% 
  dplyr::arrange(group,dataset_ID,desc(medianScore),karyotype)

fig3D <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots/Sup.FigSxx_MLDS_topGenes_moduleScore_bulkSamples_rawData.tsv',sep = '\t')
fig3D <- fig3D %>% dplyr::filter(moduleType=='MLDS_topGenes_all') %>% 
  dplyr::select(category,sampleGroup, Sample,Source,TotalScore) %>% 
  dplyr::mutate(category = factor(category,c('Normal', 'TAM / MLDS', 'Other leukaemia')),
                sampleGroup = ifelse(sampleGroup == 'MLDS','ML-DS',sampleGroup),
                sampleGroup = factor(sampleGroup,c("HSC","MPP", "MEP","CMP","GMP",
                                                   "Conventional TAM","Progressive TAM","ML-DS",
                                                   "AMKL","MLL"))) %>% 
  dplyr::rename(tissue = Source,
                enrichment_score = TotalScore,
                sample_id = Sample,
                sample_group=sampleGroup) %>% 
  dplyr::arrange(category,sample_group, sample_id,enrichment_score)


## Figure 4 ------
fig4A <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure4/Fig4x_mSignAct_L076_rawData.tsv',sep = '\t')
fig4A <- fig4A %>% dplyr::select(PDID,sig,prop) %>% 
  dplyr::arrange(PDID,desc(prop)) %>% 
  dplyr::rename('WGS_ID (cosine similarity)' = PDID,
                'SBS signature' = sig,
                'contribution' = prop) 

fig4B <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure4/Fig4B_L076_GATA1s_UMAP_rawData.tsv') %>% 
  dplyr::select(-c('UMAP_1','UMAP_2','AIres'))


fig4D <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure4/Fig4x_mSignAct_L038_rawData.tsv',sep = '\t')
fig4D <- fig4D %>% dplyr::select(PDID,sig,prop) %>% 
  dplyr::arrange(PDID,desc(prop)) %>% 
  dplyr::rename('WGS_ID (cosine similarity)' = PDID,
                'SBS signature' = sig,
                'contribution' = prop) 

fig4E_umap <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure4/Fig4E_L038_GATA1s_UMAP_rawData.tsv') %>% 
  dplyr::select(-c('UMAP_1','UMAP_2')) %>% 
  dplyr::rename('sampleID' = orig.ident)

fig4E_baf_chr5 <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure4/Fig4E_L038Tum_rawBAF_chr5_rawData.tsv')
fig4E_baf_chr17 <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure4/Fig4E_L038Tum_rawBAF_chr17_rawData.tsv')
fig4E_baf <- rbind(fig4E_baf_chr5,fig4E_baf_chr17) %>% 
  dplyr::select(c(clusterID,chr,pos,altFreq, totCount, clusterID_nCell, phasingAssign)) %>% 
  dplyr::rename(clone = clusterID,clone_nCell = clusterID_nCell)


# Supplementary Figures --------------------------------------------------------

## SuppFig S1 ------
figS1 <- data.table::fread('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/FigureS1/Supp.Fig2B_2nAKLiv_celltypeFraction_byDonorID_barPlot_rawData.tsv') %>% 
  dplyr::rename(Lineage = broadLineage3,
                "Proportion of cells" = ctFrac,
                Karyotype = Genotype) %>% 
  dplyr::mutate(Lineage = factor(Lineage,c("HSC_MPP","Ery/MegK/Mast lineage","Myeloid lineage","B lineage","T/NK lineage",
                                           "Stromal","others")),
                Karyotype = factor(Karyotype,c('Diploid','T21','T18','T22','MX','Triploid'))) %>% 
  dplyr::arrange(Karyotype,donorID,Lineage)

## SuppFig S2 ------
figS2 <- do.call(rbind,lapply(list.files('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/FigureS2/',pattern = '_meanLog2FC_withX_hm_rawData.tsv',recursive = T,full.names = T),data.table::fread)) %>% 
  dplyr::rename(Karyotype = geno,cell_type = ct) %>% 
  dplyr::select(c('Karyotype','cell_type',as.character(c(1:22,'X')))) %>% 
  dplyr::mutate(Karyotype = factor(Karyotype,c('T21','T18','T22','MX','Triploid'))) %>% 
  dplyr::arrange(Karyotype,cell_type)

## SuppFig S4 ------
figS4 <- data.table::fread('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/FigureS4/Figxx_MLDS_GATA1s_status_barplot_hor_rawData.tsv') %>% 
  dplyr::filter(donorID !='CC3') %>% 
  dplyr::rename(annotation_category = broadAnno,
                GATA1s_status = GATA1s_status2,
                nCell = n,
                cell_fraction = fracCells) %>% 
  dplyr::mutate(timePoint = factor(timePoint,c('Diagnostic','TP1','TP2','TP4',
                                               'D.Relapse', 'D.Relapse2')),
                GATA1s_status = factor(GATA1s_status,c('No GATA1 expression',
                                                       'Uninformative','GATA1 wild type','GATA1s mutation')),
                disease = factor(disease,c('TAM','MLDS')),
                annotation_category = ifelse(annotation_category == 'Tumour','Blast',annotation_category),
                annotation_category = factor(annotation_category,c("Blast","MEP","Megakaryocytes","Erythroblasts","others"))) %>% 
  dplyr::arrange(disease,donorID,annotation_category,GATA1s_status)

## SuppFig S5 ------
figS5B <- data.table::fread('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/FigureS5/Fig3C_MLDS.DEG_byDonor_barplot_rawData.tsv') %>% 
  dplyr::mutate(nGene=abs(nGene)) %>% 
  dplyr::rename(n_DEG = nGene)

## SuppFig S6 ------
figS6C <- data.table::fread('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/FigureS6/SFigxx_L067.MDS_Tum_rawBAF_sameSNPs_red_rawData.tsv') %>% 
  dplyr::select(c(celltype,chr,pos,altFreq, totCount,cov)) %>% 
  dplyr::mutate(celltype = ifelse(celltype == 'Normal','Normal cells', 'Blasts'),
                celltype = factor(celltype,c('Normal cells', 'Blasts')),
                chr = factor(chr,c('chr3','chr21','chr2'))) %>% 
  dplyr::arrange(chr) %>% 
  dplyr::rename(cell_type = celltype)

## SuppFig S7 ------
figS7B <- data.table::fread('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots/Sup.FigS6x_MLDSimprints.T21FL_topGenes_moduleScore_bulkSamples_newBulkMdat_1124_rawData.tsv') %>% 
  dplyr::filter(!Category %in% c('unclassified'),
                moduleType == 'MLDS_imprints_T21FL_all') %>% 
  dplyr::select(c(Sample, sampleGroup, TotalScore, Subgroup,Group,Category,Source)) %>% 
  dplyr::rename(enrichment_score = TotalScore) %>% 
  dplyr::mutate(Subgroup = ifelse(Subgroup == 'TMD','TAM',Subgroup),
                Group = ifelse(Group == 'Nomal_HSPCs','Normal_HSPCs',Group),
                Group = factor(Group,c('Normal_HSPCs','TAM_MLDS','Other Leukaemia')),
                sampleGroup = factor(sampleGroup,c('HSC','MPP','MEP','CMP','GMP','TAM','MLDS','AMKL','MLL'))) %>% 
  dplyr::arrange(Group,sampleGroup)
head(figS7B)
checkmate::assert_true(nrow(figS7B) == 127)


figS7C <- data.table::fread('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots/Sup.FigS6x_GATA1s_topGenes_moduleScore_bulkSamples_newBulkMdat_1124_cf20perc_rawData.tsv') %>% 
  dplyr::filter(!Category %in% c('unclassified'),
                moduleType == 'GATA1s_topGenes_all') %>% 
  dplyr::select(c(Sample, sampleGroup, TotalScore, Subgroup,Group,Category,Source)) %>% 
  dplyr::rename(enrichment_score = TotalScore) %>% 
  dplyr::mutate(Subgroup = ifelse(Subgroup == 'TMD','TAM',Subgroup),
                Group = ifelse(Group == 'Nomal_HSPCs','Normal_HSPCs',Group),
                Group = factor(Group,c('Normal_HSPCs','TAM_MLDS','Other Leukaemia')),
                sampleGroup = factor(sampleGroup,c('HSC','MPP','MEP','CMP','GMP','TAM','MLDS','AMKL','MLL'))) %>% 
  dplyr::arrange(Group,sampleGroup)
head(figS7C)
checkmate::assert_true(nrow(figS7C) == 127)



## SuppFig S9 ------
figS9A <- do.call(rbind,lapply(list.files('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/FigureS9/',pattern = 'L076_PD.*tsv$',recursive = T,full.names = T),
                               function(x){
                                 d = data.table::fread(x)
                                 d$sample = gsub('Sup.Figxx_L076_|_purpleCNprofile_rawData.tsv','',basename(x))
                                 d
                               })) %>% 
  dplyr::select(c(sample,chromosome,start, end,copyNumber,minorAlleleCopyNumber,x_start,x_end,chr_odd)) %>% 
  dplyr::rename(majorAlleleCopyNumber = copyNumber,
                x_start_for_plot = x_start,
                x_end_for_plot = x_end)

figS9B <- data.table::fread('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/FigureS9/FigS9B_L076_rawBAF_rawData.tsv') %>% 
  dplyr::mutate(clusterID = as.character(clusterID),
    clusterID = dplyr::case_when(clusterID == "D.Relapse2:BM (n=12256)" ~ "Relapse2:BM (n=12252)",
                                 clusterID == "D.Relapse:BM (n=695)" ~ "Relapse1:BM (n=693)",
                                 clusterID == "Normal (n=9469)" ~ "Normal (n=9463)",
                                             .default = clusterID),
                clusterID = factor(clusterID,c("Diagnostic:Blood (n=1040)",
                                               "Diagnostic:BM (n=2008)",
                                               "Relapse1:BM (n=693)",
                                               "Relapse2:BM (n=12252)",
                                               "Normal (n=9463)")),
    chr = factor(chr,c('chr5','chr8','chr13','chr21'))) %>% 
  dplyr::select(c('clusterID','chr','pos','altFreq','totCount','phasingAssign')) %>% 
  dplyr::arrange(clusterID,chr,pos)

checkmate::assert_true(all(!is.na(figS9B$clusterID))) 


figS9C <- do.call(rbind,lapply(list.files('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/FigureS9/',pattern = 'L038_PD.*tsv$',recursive = T,full.names = T),
                               function(x){
                                 d = data.table::fread(x)
                                 d$sample = gsub('Sup.Figxx_L038_|_purpleCNprofile_rawData.tsv','',basename(x))
                                 d
                               })) %>% 
  dplyr::select(c(sample,chromosome,start, end,copyNumber,minorAlleleCopyNumber,x_start,x_end,chr_odd)) %>% 
  dplyr::rename(majorAlleleCopyNumber = copyNumber,
                x_start_for_plot = x_start,
                x_end_for_plot = x_end)



# Create workbook ------
wb <- createWorkbook()
# Add sheets and write data
addWorksheet(wb, "Fig1D")
writeData(wb, "Fig1D", fig1D)
addWorksheet(wb, "Fig1E")
writeData(wb, "Fig1E", fig1E)
addWorksheet(wb, "Fig1F")
writeData(wb, "Fig1F", fig1F)

addWorksheet(wb, "Fig2B-C")
writeData(wb, "Fig2B-C", fig2bc)
addWorksheet(wb, "Fig2D-E-F")
writeData(wb, "Fig2D-E-F", fig2def)

addWorksheet(wb, "Fig3C")
writeData(wb, "Fig3C", fig3C)
addWorksheet(wb, "Fig3D")
writeData(wb, "Fig3D", fig3D)

addWorksheet(wb, "Fig4A")
writeData(wb, "Fig4A", fig4A)
addWorksheet(wb, "fig4B")
writeData(wb, "fig4B", fig4B)
addWorksheet(wb, "fig4D")
writeData(wb, "fig4D", fig4D)
addWorksheet(wb, "fig4E_umap")
writeData(wb, "fig4E_umap", fig4E_umap)
addWorksheet(wb, "fig4E_baf")
writeData(wb, "fig4E_baf", fig4E_baf)

addWorksheet(wb, "figS1B")
writeData(wb, "figS1B", figS1)
addWorksheet(wb, "figS2")
writeData(wb, "figS2", figS2)
addWorksheet(wb, "figS4")
writeData(wb, "figS4", figS4)
addWorksheet(wb, "figS5B")
writeData(wb, "figS5B", figS5B)
addWorksheet(wb, "figS6C")
writeData(wb, "figS6C", figS6C)
addWorksheet(wb, "figS7B")
writeData(wb, "figS7B", figS7B)
addWorksheet(wb, "figS7C")
writeData(wb, "figS7C", figS7C)
addWorksheet(wb, "figS9A")
writeData(wb, "figS9A", figS9A)
addWorksheet(wb, "figS9B")
writeData(wb, "figS9B", figS9B)
addWorksheet(wb, "figS9C")
writeData(wb, "figS9C", figS9C)


# Save file
saveWorkbook(wb, "~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Source_Data.xlsx", overwrite = TRUE)


# Prepare TAM/ML-DS/other_Leuk bulk RNA-seq count table ------------------------

# samples used in the study
bulk_mdat = readxl::read_excel('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Supplemental table 1.xlsx',skip = 2) %>% 
  dplyr::filter(`single-cell or bulk RNA-seq` == "bulk RNA-seq")
# Import count table
rawCnt = readxl::read_excel('~/lustre_mt22/MLDS_scRNAseq/bulkRNA_Data/Down_Leukemia_Klusmann_Apr23.xlsx',sheet = 'StarSalmon_Counts') %>% 
  dplyr::select(c("gene_id", "gene_name","EffectiveLength", gsub('-','.',bulk_mdat$`Sample ID`)))

checkmate::assert_true(all(c("gene_id", "gene_name","EffectiveLength", gsub('-','.',bulk_mdat$`Sample ID`)) %in% colnames(rawCnt)))
checkmate::assert_true(all(colnames(rawCnt) %in% c("gene_id", "gene_name","EffectiveLength", gsub('-','.',bulk_mdat$`Sample ID`))))
checkmate::assert_true(ncol(rawCnt) == 130)

write.csv(rawCnt,'~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/bulkRNAseq_count_table_for_Zenodo.csv',row.names = F)


# Determine number of individuals in each dataset ------------------------------
fliver_published= read.csv('~/Aneuploidy/Results/0_reprocessing_fetalREF/muzz_fetal_liver_obs.csv.gz')
table(fliver_published$fetal.ids[fliver_published$cell.labels=='MEMP'])
dplyr::n_distinct(fliver_published$fetal.ids)
dplyr::n_distinct(fliver_published$fetal.ids[fliver_published$cell.labels=='MEMP'])==dplyr::n_distinct(fliver_published$fetal.ids)

fBM = read.csv('~/Aneuploidy/Results/0_reprocessing_fetalREF/laura_fBM_obs.csv.gz')
dplyr::n_distinct(fBM$orig.ident)
table(fBM$orig.ident[fBM$cell.labels %in% c('MEMP','MEP')])
dplyr::n_distinct(fBM$orig.ident) == dplyr::n_distinct(fBM$orig.ident[fBM$cell.labels %in% c('MEMP','MEP')])
fBM_t21 = read.csv('~/Aneuploidy/Results/0_reprocessing_fetalREF/laura_fBM_T21_obs.csv.gz')
table(fBM_t21$cell.labels)
table(fBM_t21$orig.ident[fBM_t21$cell.labels %in% c('MEMP','MEP')])
dplyr::n_distinct(fBM_t21$orig.ident) == dplyr::n_distinct(fBM_t21$orig.ident[fBM_t21$cell.labels %in% c('MEMP','MEP')])

# Fetal liver MEMP/MEP -------
# Internal data: 
# Popescu 2019: 14 individuals

# Fetal BM MEMP/MEP ------
# Jardine 2021: 8 diploid individuals + 4 T21 individuals

# TAM ------
# sc data:
# bulk data

# ML-DS ------


