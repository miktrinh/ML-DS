#!/usr/bin/env bash

CORES=5
RAM="30G"
QUEUE="normal"
GROUP="team274-grp"
SCRIPT="bash ~/lustre_mt22/Aneuploidy/scripts/activeScritps/WGS_analysis/x10_MLDS_cgpVaf.sh L076 PD62331c_PD62331a_PD64665a_PD66167a ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L076/cgpVaf_input ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L076/cgpVaf_output ~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/L076/L076_unmatched_cavemanVar.bed"

bsub \
-G "${GROUP}" -q "${QUEUE}"  -n ${CORES} \
-M ${RAM} -R "select[mem>${RAM}] rusage[mem=${RAM}]" \
-o "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/farmJob_logs/%J.output" -e "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/farmJob_logs/%J.error" \
"${SCRIPT}"
