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
  fig4B = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l076_relapse.R',
  fig4C = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l076_relapse.R',
  fig4D = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l076_relapse.R',
  fig4E = '~/ML-DS/01_main_analyses/06_MLDS_refractory_relapse/l076_relapse.R'
  )

# Load in raw data for each figure --------
## Figure 1 ------
fig1D <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure1/Fig1_2nAKLiv_mean.celltypeFractionByDonorByLineage_barPlot_rawData.tsv',sep = '\t')
fig1D <- fig1D %>% 
  dplyr::rename('Lineage'='broadLineage3',
                'cellProportion'='ctFrac') %>% 
  dplyr::filter(Lineage %in% c('B lineage', 'Ery/MegK/Mast lineage')) %>% 
  dplyr::mutate(Genotype = factor(Genotype,c('Diploid','T21','T18','T22','Triploid','MX')),
                Lineage = factor(Lineage, c('Ery/MegK/Mast lineage', 'B lineage'))) %>% 
  dplyr::select(Lineage, Genotype, donorID, nCell, totalCell, cellProportion) %>% 
  dplyr::arrange(Lineage, Genotype, cellProportion)

fig1E <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure1/Fig1_2nAKLiver_GlobalImpact_DEGs_hor_rawData.tsv',sep = '\t')
fig1E <- fig1E %>% 
  dplyr::select(-nCell_2n, -nCell_AK)
  

fig1F <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure1/FigSup2X_nSharedDEGs_acrossGeno_rawData.tsv',sep = '\t')
fig1F <- fig1F %>% 
  dplyr::select(geno, direction, ct, nGeno, nGene) %>% 
  dplyr::filter(geno == 'T21') %>% 
  dplyr::mutate(direction = ifelse(direction == 'AK_up','AK up-regulated','AK down-regulated'),
                direction = factor(direction, c('AK up-regulated','AK down-regulated'))) %>% 
  dplyr::rename('Direction' = 'direction',
                Genotype = 'geno',
                'Cell type' = 'ct',
                'Number of other AK-vs-diploid comparisons' = nGeno,
                'Number of DEG' = 'nGene') %>% 
  arrange(Direction, `Number of other AK-vs-diploid comparisons`, desc(`Number of DEG`))

## Figure 2 ------
fig2B <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots/Fig2_MLDS_TumNorm_UMAP_rawData.tsv',sep = '\t')
fig2C <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2507/Plots/Fig2_MLDS_GATA1s_UMAP_v2_rawData.tsv',sep = '\t')
shared_columns = intersect(colnames(fig2B), colnames(fig2C[,!colnames(fig2C) %in% c('broadLineage')]))
fig2bc <- dplyr::full_join(fig2B, fig2C[,!colnames(fig2C) %in% c('broadLineage')], by = shared_columns)
checkmate::assert_true(nrow(fig2bc) == nrow(fig2B))

fig2bc <- fig2bc %>% dplyr::select(broadLineage, cellID, UMAP_1, UMAP_2,  GATA1s_status2, annot,  orig.ident, donorID, disease, timePoint, tissue, clinicalOutcome) %>% 
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
fig3C <- fig3C %>% dplyr::select(group,dataset_ID,Genotype,medianScore) %>% 
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
  dplyr::arrange(group,dataset_ID,desc(medianScore),Genotype)

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

fig4D <- read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Plots/Figure4/Fig4x_mSignAct_L038_rawData.tsv',sep = '\t')
fig4D <- fig4D %>% dplyr::select(PDID,sig,prop) %>% 
  dplyr::arrange(PDID,desc(prop)) %>% 
  dplyr::rename('WGS_ID (cosine similarity)' = PDID,
                'SBS signature' = sig,
                'contribution' = prop) 

  

# Create workbook
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
addWorksheet(wb, "fig4D")
writeData(wb, "fig4D", fig4D)

# Save file
saveWorkbook(wb, "~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Source_Data_wip.xlsx", overwrite = TRUE)
