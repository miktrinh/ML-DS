#!/bin/bash

#bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R'select[mem>3000] rusage[mem=3000]' -M3000 -J cgpvaf'[1-24]' ./cgpvaf.sh $patient snp
#bash /lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/MD_cgpvaf.sh /lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/02_DNA_Processing/Variant_calls/Subs/02_LCM snp

#cgpVaf.pl -d /lustre/scratch119/casm/team274sb/to3/fetal_mapping/02_DNA_processing/cgpvaf/subs -o test -a snp -mq 0 -bq 0 -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -be .sample.dupmarked.bam -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz -bo 1 -b /lustre/scratch119/casm/team274sb/to3/fetal_mapping/02_DNA_processing/cgpvaf/subs/PD53943_subs.bed -nn PDv38is_wgs -tn PD53943ac_lo0008 -chr chr1

#/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_input/L076
#/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076
#/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L076_unmatched_cavemanVar.bed
inputDir=$1
outputDir=$2
bedFile=$3
sampleFile=$4
REF_fasta='/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa'
hdr_file='/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz'

module use --append /software/CASM/modules/modulefiles
module load vafcorrect

#REF=/lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa

#SAMPLES=$(cat samples_all.txt | grep $patient | tr '\n' ' ')

#SAMPLES=$(cat /lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/02_DNA_Processing/Variant_calls/Subs/02_LCM/samples_all.txt | tr '\n' ' ')

#cd /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/

SAMPLES=$(cat $sampleFile | cut -f 2 |  grep "PD" | tr '\n' ' ')

cgpVaf.pl \
  -d $inputDir \
  -o $outputDir  \
  -a snp \
  -mq 30 \
  -bq 25 \
  -g $REF_fasta \
  -be .sample.dupmarked.bam \
  -hdr $hdr_file \
  -bo 1 \
  -b $bedFile \
  -nn PDv38is_wgs -tn $SAMPLES #\
  #-chr ${arr[$LSB_JOBINDEX-1]}
  
  
# for s in $SAMPLES
# do
#   echo $s
#   
#   #declare -a arr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
#   
#   
# done

#cgpVaf.pl -d /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_input/L076 
# -o /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076 -a snp -mq 30 -bq 25 
# -g /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -be .sample.dupmarked.bam 
# -hdr /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz 
# -bo 1 -b /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L076_somaticVariants.bed -nn PDv38is_wgs -tn PD60301a

