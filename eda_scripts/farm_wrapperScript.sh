#!/usr/bin/env bash

CORES=25
RAM="200G"
QUEUE="normal"
GROUP="team274-grp"
#SCRIPT="/software/R-4.1.0/bin/Rscript /lustre/scratch125/casm/team274sb/mt22/farmTest/genotyping_farmTest.R"
#SCRIPT="Rscript /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/finalScripts/x1.2_LRv2_on_inhouse_MLDS_data.R"
#SCRIPT="bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/finalScripts/GATA1_mut_reads_v2.sh"
#SCRIPT="Rscript /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/finalScripts/trajectoryProjection.R"
SCRIPT="Rscript /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/xx19_combine_leukObj.R"

#conda activate alleleIntegrator

bsub \
-G "${GROUP}" -q "${QUEUE}"  -n ${CORES} \
-M ${RAM} -R "select[mem>${RAM}] rusage[mem=${RAM}]" \
-o "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/%J.output" -e "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/%J.error" \
"${SCRIPT}"
