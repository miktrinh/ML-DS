#!/usr/bin/env bash

CORES=5
RAM="10G"
QUEUE="normal"
GROUP="team274-grp"
SCRIPT="bash /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/16_tNS_foetalLivers/gata1_variants_GATK.sh"

#module load cellgen/gatk
bsub \
-G "${GROUP}" -q "${QUEUE}"  -n ${CORES} \
-M ${RAM} -R "select[mem>${RAM}] rusage[mem=${RAM}]" \
-o "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/farmJob_logs/%J.output" -e "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/farmJob_logs/%J.error" \
"${SCRIPT}"
