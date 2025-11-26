## Mutations in L156 (WGS and t-NS)

outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}



##----------------##
##   Libraries  ####
##----------------##

library(VariantAnnotation)
library(tidyverse)
source('~/lustre_mt22/generalScripts/caveman_readvcf.R')
source('~/lustre_mt22/generalScripts/binom_mix_model_Tim.R')
source('~/lustre_mt22/generalScripts/utils/wgs_analysis_helperFunctions.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')


##------------------------------------------------##
##    Import Caveman output - Diagnostic        ####
##------------------------------------------------##
d_caveman = caveman2(sample='PD60301a',projectid = 3030)
d_caveman$timepoint = 'Diagnostic'
d_caveman = d_caveman[d_caveman$Flag == T,]
d_caveman$snv_ID = paste0(d_caveman$Chr,':',d_caveman$Pos,'_',d_caveman$Ref,'/',d_caveman$Alt)



##-----------------------------------------------##
##      Import Caveman output - TP1            ####
##-----------------------------------------------##
tp1_caveman = caveman2(sample='PD61846a',projectid = 3030)
tp1_caveman$timepoint = 'TP1'
tp1_caveman = tp1_caveman[tp1_caveman$Flag == T,]
tp1_caveman$snv_ID = paste0(tp1_caveman$Chr,':',tp1_caveman$Pos,'_',tp1_caveman$Ref,'/',tp1_caveman$Alt)



