#!/usr/bin/env bash

CORES=20
RAM="200G"
QUEUE="hugemem"
GROUP="team274"
#SCRIPT="/software/R-4.1.0/bin/Rscript /lustre/scratch125/casm/team274sb/mt22/farmTest/genotyping_farmTest.R"
#SCRIPT="Rscript /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/finalScripts/x1.2_LRv2_on_inhouse_MLDS_data.R"
#SCRIPT="bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/finalScripts/GATA1_mut_reads_v2.sh"
SCRIPT="Rscript /lustre/scratch126/casm/team274sb/aw35/misc_scripts/mSigAct_wrapper_v2.R GRCh38 ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L076/signatures/01_Input/ ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L076/signatures/02_Output/ SBS1_SBS5_SBS13_SBS15_SBS18_SBS31_SBS35_SBS44_SBS84_SBS86 FALSE FALSE"

# L038: SBS1_SBS5_SBS13_SBS15_SBS18_SBS31_SBS35_SBS39_SBS44_SBS84_SBS86

#module load R/4.3.1-jupyter-v0

bsub \
-G "${GROUP}" -q "${QUEUE}"  -n ${CORES} \
-M ${RAM} -R "select[mem>${RAM}] rusage[mem=${RAM}]" \
-o "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/mSigAct/%J.output" -e "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/mSigAct/%J.error" \
"${SCRIPT}"
