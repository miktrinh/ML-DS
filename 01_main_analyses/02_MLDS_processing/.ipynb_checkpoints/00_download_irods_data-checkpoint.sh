############### TAM / MLDS - scRNAseq - cellranger700 ###############
samples=(
# # L039_D
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234194_GRCh38-2020-A'
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234195_GRCh38-2020-A'
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234196_GRCh38-2020-A'
# # L039_TP1
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415925_GRCh38-2020-A'
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415926_GRCh38-2020-A'
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415927_GRCh38-2020-A'
# # L039_TP4
# '/seq/illumina/runs/45/45842/cellranger/cellranger700_count_45842_SB_Leuk13104278_GRCh38-2020-A'

# L038_D
'/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234191_GRCh38-2020-A'
'/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234192_GRCh38-2020-A'
'/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234193_GRCh38-2020-A'
# L038_TP1 (Refractory)
'/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415922_GRCh38-2020-A'
'/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415923_GRCh38-2020-A'
'/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415924_GRCh38-2020-A'

# # L019_D
# '/seq/illumina/runs/45/45842/cellranger/cellranger700_count_45842_SB_Leuk13104279_GRCh38-2020-A'
# # L019_TP2
# '/seq/illumina/runs/47/47089/cellranger/cellranger700_count_47089_SB_Leuk13645528_GRCh38-2020-A'
# # L019_TP4
# '/seq/illumina/runs/47/47089/cellranger/cellranger700_count_47089_SB_Leuk13645529_GRCh38-2020-A'

# # L040_D
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234197_GRCh38-2020-A'
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234198_GRCh38-2020-A'
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234199_GRCh38-2020-A'
# # L040_TP1
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415928_GRCh38-2020-A'
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415929_GRCh38-2020-A'
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415930_GRCh38-2020-A'
# # L040_TP4
# '/seq/illumina/runs/47/47089/cellranger/cellranger700_count_47089_SB_Leuk13645524_GRCh38-2020-A'

# # L041_D (not good quality)
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234200_GRCh38-2020-A'
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234201_GRCh38-2020-A'
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234202_GRCh38-2020-A'
# # L041_TP1
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415931_GRCh38-2020-A'
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415932_GRCh38-2020-A'
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415933_GRCh38-2020-A'

# # L042_D
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234203_GRCh38-2020-A'
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234204_GRCh38-2020-A'
# '/seq/illumina/runs/46/46351/cellranger/cellranger700_count_46351_SB_Leuk13234205_GRCh38-2020-A'
# # L042_TP1
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415934_GRCh38-2020-A'
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415935_GRCh38-2020-A'
# '/seq/illumina/runs/46/46682/cellranger/cellranger700_count_46682_SB_Leuk13415936_GRCh38-2020-A'

# # L075_D
# '/seq/illumina/runs/47/47260/cellranger/cellranger700_count_47260_SB_Leuk13697518_GRCh38-2020-A'

# # L076_D_BM
# '/seq/illumina/runs/47/47580/cellranger/cellranger700_count_47580_SB_Leuk13760338_GRCh38-2020-A'
# L076_D_PB
# '/seq/illumina/runs/47/47260/cellranger/cellranger700_count_47260_SB_Leuk13697519_GRCh38-2020-A'
# # L076_Relapse1_D
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14635833_GRCh38-2020-A'
# # L076_Relapse2_D
# '/seq/illumina/runs/49/49031/cellranger/cellranger700_count_49031_ALeuk_RNA14832000_GRCh38-2020-A'
# '/seq/illumina/runs/49/49031/cellranger/cellranger700_count_49031_ALeuk_RNA14832001_GRCh38-2020-A'

# # L091_D
# '/seq/illumina/runs/47/47510/cellranger/cellranger700_count_47510_SB_Leuk13778112_GRCh38-2020-A'

# # L156_D
# '/seq/illumina/runs/48/48209/cellranger/cellranger700_count_48209_ALeuk_RNA14564480_GRCh38-2020-A'
# '/seq/illumina/runs/48/48209/cellranger/cellranger700_count_48209_ALeuk_RNA14564481_GRCh38-2020-A'
# # L156_Relapse
# '/seq/illumina/runs/48/48983/cellranger/cellranger700_count_48983_MY_200531_14784986_GRCh38-2020-A'
# '/seq/illumina/runs/48/48983/cellranger/cellranger700_count_48983_MY_200531_14784987_GRCh38-2020-A'

# # L178_D
# '/seq/illumina/runs/49/49031/cellranger/cellranger700_count_49031_ALeuk_RNA14831998_GRCh38-2020-A'
# '/seq/illumina/runs/49/49031/cellranger/cellranger700_count_49031_ALeuk_RNA14831999_GRCh38-2020-A'
# # L178_TP1
# '/seq/illumina/runs/49/49226/cellranger/cellranger700_count_49226_ALeuk_RNA14910210_GRCh38-2020-A'
# '/seq/illumina/runs/49/49226/cellranger/cellranger700_count_49226_ALeuk_RNA14910211_GRCh38-2020-A'

# # L182_D
# '/seq/illumina/runs/49/49087/cellranger/cellranger700_count_49087_ALeuk_RNA14872269_GRCh38-2020-A'
# '/seq/illumina/runs/49/49087/cellranger/cellranger700_count_49087_ALeuk_RNA14872270_GRCh38-2020-A'

# # L114_D
# '/seq/illumina/runs/49/49087/cellranger/cellranger700_count_49087_ALeuk_RNA14872267_GRCh38-2020-A'
# '/seq/illumina/runs/49/49087/cellranger/cellranger700_count_49087_ALeuk_RNA14872268_GRCh38-2020-A'

# # CC1
# '/seq/illumina/runs/48/48097/cellranger/cellranger700_count_48097_MY_200531_14406280_GRCh38-2020-A'
# '/seq/illumina/runs/48/48097/cellranger/cellranger700_count_48097_MY_200531_14406281_GRCh38-2020-A'

# # CC2
# '/seq/illumina/runs/48/48097/cellranger/cellranger700_count_48097_MY_200531_14406282_GRCh38-2020-A'
# '/seq/illumina/runs/48/48097/cellranger/cellranger700_count_48097_MY_200531_14406283_GRCh38-2020-A'

# # CC3
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14635834_GRCh38-2020-A'

# # CC4
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14635835_GRCh38-2020-A'

# # CC5
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14635836_GRCh38-2020-A'
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14635837_GRCh38-2020-A'

# # CC6
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14637036_GRCh38-2020-A'
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14637037_GRCh38-2020-A'

# # CC7
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14637038_GRCh38-2020-A'
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14637039_GRCh38-2020-A'

# # CC8
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14637040_GRCh38-2020-A'
# '/seq/illumina/runs/48/48665/cellranger/cellranger700_count_48665_MY_200531_14637041_GRCh38-2020-A'

)


outDir=Data/MLDS_scRNAseq
mkdir -p $outDir
for irod_path in ${samples[@]}; do
   sampleName=$(basename $irod_path)
   mkdir -p $outDir/$sampleName
   
#    irods iget -r $irod_path/filtered_feature_bc_matrix $outDir/$sampleName/filtered_feature_bc_matrix
#    irods iget -r $irod_path/raw_feature_bc_matrix $outDir/$sampleName/raw_feature_bc_matrix
    
   irods iget -r $irod_path/possorted_genome_bam.bam $outDir/$sampleName/possorted_genome_bam.bam
   irods iget -r $irod_path/possorted_genome_bam.bam.bai $outDir/$sampleName/possorted_genome_bam.bam.bai



done




############### Other Leukaemia - scRNAseq - cellranger700 ###############
# samples=(#
# # DS-brain lymphoma
# '/seq/illumina/runs/48/48097/cellranger/cellranger700_count_48097_CG_SB_NB14406184_GRCh38-2020-A'
# '/seq/illumina/runs/48/48097/cellranger/cellranger700_count_48097_CG_SB_NB14406185_GRCh38-2020-A'


# '/seq/illumina/runs/49/49031/cellranger/cellranger700_count_49031_ALeuk_RNA14832000_GRCh38-2020-A'
# '/seq/illumina/runs/49/49031/cellranger/cellranger700_count_49031_ALeuk_RNA14832001_GRCh38-2020-A'
# )



