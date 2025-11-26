## Pindel analysis ##

outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/Indels'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

library(VariantAnnotation)
source('~/lustre_mt22/generalScripts/caveman_readvcf.R')
source('~/lustre_mt22/generalScripts/binom_mix_model_Tim.R')
source('~/lustre_mt22/generalScripts/utils/wgs_analysis_helperFunctions.R')
