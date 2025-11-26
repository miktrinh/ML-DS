#!/bin/bash

#bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R'select[mem>3000] rusage[mem=3000]' -M3000 -J cgpvaf'[1-24]' ./cgpvaf.sh $patient snp
#bash /lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/MD_cgpvaf.sh /lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/02_DNA_Processing/Variant_calls/Subs/02_LCM snp

#cgpVaf.pl -d /lustre/scratch119/casm/team274sb/to3/fetal_mapping/02_DNA_processing/cgpvaf/subs -o test -a snp -mq 0 -bq 0 -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -be .sample.dupmarked.bam -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz -bo 1 -b /lustre/scratch119/casm/team274sb/to3/fetal_mapping/02_DNA_processing/cgpvaf/subs/PD53943_subs.bed -nn PDv38is_wgs -tn PD53943ac_lo0008 -chr chr1

patient=$1
SAMPLES=$2
inDir=$3
outDir=$4
bed_file=$5
REF_fasta='/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa'
hdr_file='/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz'

module use --append /software/CASM/modules/modulefiles
module load vafcorrect


#REF=/lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa
#SAMPLES=$(cat samples_all.txt | grep $patient | tr '\n' ' ')

#SAMPLES=$(cat /lustre/scratch126/casm/team274sb/md39/Tim_Fetal_Analysis/02_DNA_Processing/Variant_calls/Subs/02_LCM/samples_all.txt | tr '\n' ' ')

cd $outDir

#SAMPLES=$(cat /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L076_samples.csv |grep "PD"| tr '\n' ' ')

echo $REF_fasta

echo $patient

if echo $SAMPLES | grep -q '_' ; then
    SAMPLES=$(echo $SAMPLES | tr '_' ' ')
fi
echo $SAMPLES

#declare -a arr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

if [[ ! -d $outDir ]]; then
    mkdir $outDir
    echo $outDir 
fi


# cgpVaf.pl \
# -d $inDir \
# -o $outDir \
# -a snp \
# -mq 30 \
# -bq 25 \
# -g $REF_fasta \
# -be .sample.dupmarked.bam \
# -hdr $hdr_file \
# -bo 1 \
# -b $bed_file \
# -nn PDv38is_wgs -tn $SAMPLES #\
# #-chr ${arr[$LSB_JOBINDEX-1]}

#cgpVaf.pl -d /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_input/L076 
# -o /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/cgpVaf_output/L076 -a snp -mq 30 -bq 25 
# -g /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -be .sample.dupmarked.bam 
# -hdr /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz 
# -bo 1 -b /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/L076_somaticVariants.bed -nn PDv38is_wgs -tn PD60301a



## cgpVaf with no filter mq = 0 & bq = 0
noFilter_outDir=$outDir"_nofilter"

if [[ ! -d $noFilter_outDir ]]; then
    mkdir $noFilter_outDir
    echo $noFilter_outDir 
fi

cgpVaf.pl \
-d $inDir \
-o $noFilter_outDir \
-a snp \
-mq 0 \
-bq 0 \
-g $REF_fasta \
-be .sample.dupmarked.bam \
-hdr $hdr_file \
-bo 1 \
-b $bed_file \
-nn PDv38is_wgs -tn $SAMPLES #\
# #cgpVaf.pl
# #'-g' genome must be defined
# #Usage:
# # cgpVaf.pl [-h] -d -a -g -tn -nn -e -o [ -tb -nb -b -t -c -r -m -ao -mq
# # -pid -bo -vcf -v]
# 
# #  Required Options (inputDir and variant_type must be defined):
# 
# #   --variant_type         (-a)     variant type (snp or indel) [default snp]
# #  --inputDir             (-d)     input directory path containing bam and vcf files
# # --genome               (-g)     genome fasta file name (default genome.fa)
# #--tumour_name          (-tn)    Tumour sample name(s), space separated [ co-located bam/cram, index and bas files expected, unless -tb specified]
# #--normal_name          (-nn)    Normal sample name [ co-located bam/cram, index and bas files expected, unless -nb specified ]
# #--outDir               (-o)     Output folder
# #
# #     Optional
# #     --vcfExtension         (-e)     vcf file extension string after the sample name - INCLUDE preceding dot (default: .caveman_c.annot.vcf.gz)
# #                                    - optional if -bo 1
# #                                   - ignored  when -vcf defined
# #  --tumour_bam           (-tb)    tumour bam/cram file(s) space separated list [ optional if -tn is specified]
# #                                 - if not defined will be deduced from tumour_name
# #                                - should be specified in same order as tumour sample names
# # --normal_bam           (-nb)    normal bam/cram file [optional if -nn is specified]
# #                                - if not defined will be deduced from --normal_name
# #--infoTags             (-t)     comma separated list of tags to be included in the tsv output, vcf file by default includes all data
# #                                 (default: VD,VW,VT,VC for Vagrent annotations)
# #       --bedIntervals         (-b)     tab separated file containing list of intervals in the form of <chr><pos> <ref><alt> (e.g 1  14000  A  C)
# #       --restrict_flag        (-r)     restrict analysis (possible values 1 : PASS or 0 : ALL) [default 1 ]
# #       --chromosome           (-chr)   restrict analysis chromosomes [space separated]
# #       --concat               (-ct)    concat per chromosome results to a single vcf file
# #       --augment              (-m)     Augment original vcf file (valid for indels)
# #                                       - this will add additional fields[ MTR, WTR, AMB] to FORMAT column of NORMAL and TUMOUR samples ] (default FALSE)
# #       --output_vcfExtension  (-oe)    Extension for augmented VCF, see --augment [vaf.vcf]
# #       --map_quality          (-mq)    read mapping quality threshold
# #       --base_quality         (-bq)    base quality threshold for snp
# #       --exonerate_pct        (-exp)   report alignment over a percentage of the maximum score attainable by each query (exonerate/indel parameter) [default 92]
# #       --exonerate_mb         (-emb)   Max memory to alow Exonerate to use (indel) [default: 6000]
# #       --depth                (-dp)    comma separated list of field(s) as specified in FORMAT field representing total depth at given location
# #       --high_depth_bed       (-hdr)   High Depth Region(HDR) bed file (tabix indexed) to mask high depth regions in the genome
# #       --bed_only             (-bo)    Only analyse bed intervals in the file (default 0: analyse vcf and bed interval)
# #       --vcf                  (-vcf)   user defined input vcf file path(s) [ optional if -tn is specified ]
# #                                       - if not defined will be deduced from --tumour_name and --vcfExtension
# #                                       - should be specified in same order as tumour sample names
# #       --filter_exc           (-fexc)  Sam flag values to exclude when checking reads for read length
# #       --filter_inc           (-finc)  Sam flag values to include when checking reads for read length
# #       --id_int_project       (-pid)   Internal project id [WTSI only]
# 
# --help                 (-h)     Display this help message
# --version              (-v)     provide version information for vaf
# 
# Examples:
#   #Merge vcf files to create single vcf containing union of all the variant sites and provides pileup output for each location
#   perl cgpVaf.pl -d tmpvcfdir -o testout -a snp -g genome.fa -e .caveman_c.annot.vcf.gz -nn PD21369b -tn PD26296a PD26296c2
# #Merge vcf files to create single vcf containing union of all the variant sites and provides allele count for underlying indel location
# perl cgpVaf.pl -d tmpvcfdir -o testout -a indel -g genome.fa -e .caveman_c.annot.vcf.gz -nn PD21369b -tn PD26296a PD26296c2
# # Run per chromosome analysis
# perl cgpVaf.pl -d tmpvcfdir -o testout -a indel -g genome.fa -e .caveman_c.annot.vcf.gz -nn sampleb -tn samplea samplec -chr 1
# perl cgpVaf.pl -d tmpvcfdir -o testout -a indel -g genome.fa -e .caveman_c.annot.vcf.gz -nn sampleb -tn samplea samplec -chr 2
# # concatenate per chromosome output to single vcf
# perl cgpVaf.pl -d tmpvcfdir -o testout -a indel -g genome.fa -e .caveman_c.annot.vcf.gz -nn sampleb -tn samplea samplec -ct 1
# 
# 
