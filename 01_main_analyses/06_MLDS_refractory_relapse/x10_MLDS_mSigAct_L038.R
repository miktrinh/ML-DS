## Run SBS signature extraction on L076 and L038

##------------------------------------##
##    SBS signature extraction      ####
##------------------------------------##
outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/signatures'
if(!dir.exists(outDir)){
  dir.create(outDir)
}
setwd(outDir)

source('~/lustre_mt22/generalScripts/utils/wgs_analysis_helperFunctions.R')


##------ Generate trinucleotide matrix    ##

# Read in mutations after filtering
mutations = read.csv('~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/L038_D.vs.TP1_finalSomaticVariants.csv')
mutations$PD60301a_true_mut = ifelse(mutations$somaticVar_type != 'unique_TP1:shearwater_Failed',1,0)
mutations$PD61846a_true_mut = ifelse(grepl('unique_D:',mutations$somaticVar_type),0,1)
rownames(mutations) = gsub(':|/','_',mutations$snv_ID)

generate_trinucleotide_matrix(mutations,outDir=file.path(outDir,'01_Input'),patientID='L038')


##------ Prep for mSigAct farm job    ##
inputDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/signatures/01_Input'
dir.create(inputDir)
outputDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L038/signatures/02_Output'
dir.create(outputDir)



system(sprintf("cp /lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/signatures/mSigAct/01_Input/COSMIC_v3.4_SBS_GRCh38.txt %s",inputDir))

