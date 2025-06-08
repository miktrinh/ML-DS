############### Download all AK fAdrenal (22 chanels) + fKidney (22 chanels) scRNAseq files ###############
samples=(
# '/seq/illumina/runs/37/37295/cellranger/cellranger302_count_37295_MY_200531_10043294_GRCh38-1_2_0'
# '/seq/illumina/runs/37/37295/cellranger/cellranger302_count_37295_MY_200531_10043295_GRCh38-1_2_0'
# '/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10192865_GRCh38-1_2_0'
# '/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10192866_GRCh38-1_2_0'
# '/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10220136_GRCh38-1_2_0'
# '/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10220137_GRCh38-1_2_0'
# '/seq/illumina/runs/39/39041/cellranger/cellranger302_count_39041_MY_200531_10260914_GRCh38-1_2_0'
# '/seq/illumina/runs/39/39041/cellranger/cellranger302_count_39041_MY_200531_10260915_GRCh38-1_2_0'
# '/seq/illumina/runs/39/39041/cellranger/cellranger302_count_39041_MY_200531_10260916_GRCh38-1_2_0'
# '/seq/illumina/runs/39/39041/cellranger/cellranger302_count_39041_MY_200531_10260917_GRCh38-1_2_0'
# '/seq/illumina/runs/40/40665/cellranger/cellranger302_count_40665_MY_200531_10621687_GRCh38-1_2_0'
# '/seq/illumina/runs/40/40665/cellranger/cellranger302_count_40665_MY_200531_10621688_GRCh38-1_2_0'
# '/seq/illumina/runs/40/40665/cellranger/cellranger302_count_40665_MY_200531_10683132_GRCh38-1_2_0'
# '/seq/illumina/runs/40/40665/cellranger/cellranger302_count_40665_MY_200531_10683133_GRCh38-1_2_0'
# '/seq/illumina/runs/41/41456/cellranger/cellranger302_count_41456_MY_200531_10864083_GRCh38-1_2_0'
# '/seq/illumina/runs/41/41456/cellranger/cellranger302_count_41456_MY_200531_10864084_GRCh38-1_2_0'
# '/seq/illumina/runs/42/42835/cellranger/cellranger302_count_42835_MY_200531_11528565_GRCh38-1_2_0'
# '/seq/illumina/runs/42/42835/cellranger/cellranger302_count_42835_MY_200531_11528566_GRCh38-1_2_0'
# '/seq/illumina/runs/44/44847/cellranger/cellranger302_count_44847_MY_200531_12750606_GRCh38-1_2_0'
# '/seq/illumina/runs/44/44847/cellranger/cellranger302_count_44847_MY_200531_12750607_GRCh38-1_2_0'
# '/seq/illumina/runs/44/44258/cellranger/cellranger302_count_44258_MY_200531_12408350_GRCh38-1_2_0'
# '/seq/illumina/runs/44/44258/cellranger/cellranger302_count_44258_MY_200531_12408351_GRCh38-1_2_0'

# Kidney
'/seq/illumina/runs/37/37295/cellranger/cellranger302_count_37295_MY_200531_10043296_GRCh38-1_2_0'
'/seq/illumina/runs/37/37295/cellranger/cellranger302_count_37295_MY_200531_10043297_GRCh38-1_2_0'
'/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10220040_GRCh38-1_2_0'
'/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10220041_GRCh38-1_2_0'
'/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10192867_GRCh38-1_2_0'
'/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10192868_GRCh38-1_2_0'
'/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10220138_GRCh38-1_2_0'
'/seq/illumina/runs/37/37929/cellranger/cellranger302_count_37929_MY_200531_10220139_GRCh38-1_2_0'
'/seq/illumina/runs/39/39041/cellranger/cellranger302_count_39041_MY_200531_10260912_GRCh38-1_2_0'
'/seq/illumina/runs/39/39041/cellranger/cellranger302_count_39041_MY_200531_10260913_GRCh38-1_2_0'
'/seq/illumina/runs/39/39041/cellranger/cellranger302_count_39041_MY_200531_10260918_GRCh38-1_2_0'
'/seq/illumina/runs/39/39041/cellranger/cellranger302_count_39041_MY_200531_10260919_GRCh38-1_2_0'
'/seq/illumina/runs/40/40665/cellranger/cellranger302_count_40665_MY_200531_10621685_GRCh38-1_2_0'
'/seq/illumina/runs/40/40665/cellranger/cellranger302_count_40665_MY_200531_10621686_GRCh38-1_2_0'

'/seq/illumina/runs/40/40665/cellranger/cellranger302_count_40665_MY_200531_10683130_GRCh38-1_2_0'
'/seq/illumina/runs/40/40665/cellranger/cellranger302_count_40665_MY_200531_10683131_GRCh38-1_2_0'
'/seq/illumina/runs/41/41456/cellranger/cellranger302_count_41456_MY_200531_10864081_GRCh38-1_2_0'
'/seq/illumina/runs/41/41456/cellranger/cellranger302_count_41456_MY_200531_10864082_GRCh38-1_2_0'
'/seq/illumina/runs/42/42835/cellranger/cellranger302_count_42835_MY_200531_11528563_GRCh38-1_2_0'
'/seq/illumina/runs/42/42835/cellranger/cellranger302_count_42835_MY_200531_11528564_GRCh38-1_2_0'
'/seq/illumina/runs/44/44258/cellranger/cellranger302_count_44258_MY_200531_12408352_GRCh38-1_2_0'
'/seq/illumina/runs/44/44258/cellranger/cellranger302_count_44258_MY_200531_12408353_GRCh38-1_2_0'
)



#outDir=Data/fAdrenals_scRNAseq
outDir=Data/fKidneys_scRNAseq

for irod_path in ${samples[@]}; do
   sampleName=$(basename $irod_path)
   mkdir -p $outDir/$sampleName
   
   irods iget -r $irod_path/filtered_feature_bc_matrix $outDir/$sampleName/filtered_feature_bc_matrix
   irods iget -r $irod_path/raw_feature_bc_matrix $outDir/$sampleName/raw_feature_bc_matrix
    
#    irods iget -r $irod_path/possorted_genome_bam.bam $outDir/$sampleName/possorted_genome_bam.bam
#    irods iget -r $irod_path/possorted_genome_bam.bam.bai $outDir/$sampleName/possorted_genome_bam.bam.bai



done