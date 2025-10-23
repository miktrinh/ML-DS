#!/usr/bin/env bash

CORES=80
RAM="200G"
QUEUE="hugemem"
GROUP="team274-grp"
#SCRIPT="/software/R-4.1.0/bin/Rscript /lustre/scratch125/casm/team274sb/mt22/farmTest/genotyping_farmTest.R"
SCRIPT="Rscript /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/finalScripts/0_genotyping_fetalSamples.R"

#conda activate alleleIntegrator

bsub \
-G "${GROUP}" -q "${QUEUE}"  -n ${CORES} \
-M ${RAM} -R "select[mem>${RAM}] rusage[mem=${RAM}]" \
-o "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/finalScripts/%J.output" -e "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/finalScripts/%J.error" \
"${SCRIPT}"
