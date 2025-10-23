############### Fetal Adrenal Reference cellranger ###############
samples=(
# Adrenal
# '/seq/27340/cellranger/cellranger302_count_27340_5388STDY7717452_GRCh38-1_2_0'
# '/seq/27340/cellranger/cellranger302_count_27340_5388STDY7717453_GRCh38-1_2_0'
# '/seq/27340/cellranger/cellranger302_count_27340_5388STDY7717454_GRCh38-1_2_0'
# '/seq/27340/cellranger/cellranger302_count_27340_5388STDY7717455_GRCh38-1_2_0'
# '/seq/27340/cellranger/cellranger302_count_27340_5388STDY7717456_GRCh38-1_2_0'
# '/seq/27340/cellranger/cellranger302_count_27340_5388STDY7717457_GRCh38-1_2_0'
# '/seq/27340/cellranger/cellranger302_count_27340_5388STDY7717458_GRCh38-1_2_0'
# '/seq/27340/cellranger/cellranger302_count_27340_5388STDY7717459_GRCh38-1_2_0'
# '/seq/illumina/runs/28/28495/cellranger/cellranger302_count_28495_5698STDY7839907_GRCh38-1_2_0'
# '/seq/illumina/cellranger/cellranger302_count_5698STDY7839909_GRCh38-1_2_0'
# '/seq/illumina/runs/28/28495/cellranger/cellranger302_count_28495_5698STDY7839917_GRCh38-1_2_0'
# '/seq/30070/cellranger/cellranger302_count_30070_WSSS8012017_GRCh38-1_2_0'
# '/seq/30714/cellranger/cellranger302_count_30714_WSSS8011223_GRCh38-1_2_0'
# '/seq/illumina/runs/32/32644/cellranger/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-1_2_0'
# '/seq/illumina/runs/32/32644/cellranger/cellranger302_count_32644_WSSS_F_Adr8710633_GRCh38-1_2_0'
# '/seq/illumina/runs/32/32644/cellranger/cellranger302_count_32644_WSSS_F_Adr8710634_GRCh38-1_2_0'
# '/seq/illumina/runs/32/32644/cellranger/cellranger302_count_32644_WSSS_F_Adr8710635_GRCh38-1_2_0'
# '/seq/illumina/runs/33/33701/cellranger/cellranger302_count_33701_WSSS_F_Adr8768489_GRCh38-1_2_0'
# '/seq/illumina/runs/33/33701/cellranger/cellranger302_count_33701_WSSS_F_Adr8768490_GRCh38-1_2_0'


# Kidney
'/seq/23208/cellranger/cellranger302_count_23208_4834STDY7002875_GRCh38-1_2_0'
'/seq/23208/cellranger/cellranger302_count_23208_4834STDY7002876_GRCh38-1_2_0'
'/seq/illumina/cellranger/cellranger302_count_4834STDY7002881_GRCh38-1_2_0'
'/seq/illumina/cellranger/cellranger302_count_4834STDY7002885_GRCh38-1_2_0'
'/seq/illumina/cellranger/cellranger302_count_4834STDY7002886_GRCh38-1_2_0'
'/seq/25871/cellranger/cellranger302_count_25871_FCAImmP7462242_GRCh38-1_2_0'
'/seq/25871/cellranger/cellranger302_count_25871_FCAImmP7462243_GRCh38-1_2_0'
'/seq/26087/cellranger/cellranger302_count_26087_FCAImmP7528292_GRCh38-1_2_0'
'/seq/26087/cellranger/cellranger302_count_26087_FCAImmP7528293_GRCh38-1_2_0'
'/seq/26280/cellranger/cellranger302_count_26280_FCAImmP7555849_GRCh38-1_2_0'
'/seq/26322/cellranger/cellranger302_count_26322_FCAImmP7555850_GRCh38-1_2_0'
'/seq/26321/cellranger/cellranger302_count_26321_FCAImmP7555859_GRCh38-1_2_0'
'/seq/26421/cellranger/cellranger302_count_26421_FCAImmP7579214_GRCh38-1_2_0'
'/seq/26421/cellranger/cellranger302_count_26421_FCAImmP7579215_GRCh38-1_2_0'

# Liver
# '/seq/23208/cellranger/cellranger302_count_23208_4834STDY7002877_GRCh38-1_2_0'
# '/seq/23208/cellranger/cellranger302_count_23208_4834STDY7002878_GRCh38-1_2_0'
# '/seq/illumina/cellranger/cellranger302_count_4834STDY7002882_GRCh38-1_2_0'
# '/seq/23395/cellranger/cellranger302_count_23395_4834STDY7038750_GRCh38-1_2_0'
# '/seq/23395/cellranger/cellranger302_count_23395_4834STDY7038751_GRCh38-1_2_0'
# '/seq/24479/cellranger/cellranger302_count_24479_FCAImmP7179363_GRCh38-1_2_0'
# '/seq/24479/cellranger/cellranger302_count_24479_FCAImmP7179364_GRCh38-1_2_0'
# '/seq/24543/cellranger/cellranger302_count_24543_FCAImmP7198434_GRCh38-1_2_0'
# '/seq/24543/cellranger/cellranger302_count_24543_FCAImmP7198628_GRCh38-1_2_0'
# '/seq/24543/cellranger/cellranger302_count_24543_FCAImmP7198629_GRCh38-1_2_0'
# '/seq/24543/cellranger/cellranger302_count_24543_FCAImmP7198630_GRCh38-1_2_0'
# '/seq/24543/cellranger/cellranger302_count_24543_FCAImmP7198631_GRCh38-1_2_0'
# '/seq/24966/cellranger/cellranger302_count_24966_FCAImmP7277552_GRCh38-1_2_0'
# '/seq/24966/cellranger/cellranger302_count_24966_FCAImmP7277553_GRCh38-1_2_0'
# '/seq/24967/cellranger/cellranger302_count_24967_FCAImmP7277560_GRCh38-1_2_0'
# '/seq/24967/cellranger/cellranger302_count_24967_FCAImmP7277561_GRCh38-1_2_0'
# '/seq/25102/cellranger/cellranger302_count_25102_FCAImmP7316894_GRCh38-1_2_0'
# '/seq/25102/cellranger/cellranger302_count_25102_FCAImmP7316895_GRCh38-1_2_0'
# '/seq/25102/cellranger/cellranger302_count_25102_FCAImmP7316889_GRCh38-1_2_0'
# '/seq/25102/cellranger/cellranger302_count_25102_FCAImmP7316890_GRCh38-1_2_0'
# '/seq/25102/cellranger/cellranger302_count_25102_FCAImmP7316891_GRCh38-1_2_0'
# '/seq/25102/cellranger/cellranger302_count_25102_FCAImmP7316892_GRCh38-1_2_0'
# '/seq/25102/cellranger/cellranger302_count_25102_FCAImmP7316893_GRCh38-1_2_0'
# '/seq/25507/cellranger/cellranger302_count_25507_FCAImmP7352192_GRCh38-1_2_0'
# '/seq/25507/cellranger/cellranger302_count_25507_FCAImmP7352193_GRCh38-1_2_0'
# '/seq/25507/cellranger/cellranger302_count_25507_FCAImmP7352194_GRCh38-1_2_0'
# '/seq/25507/cellranger/cellranger302_count_25507_FCAImmP7352195_GRCh38-1_2_0'
# '/seq/25507/cellranger/cellranger302_count_25507_FCAImmP7352196_GRCh38-1_2_0'
# '/seq/25871/cellranger/cellranger302_count_25871_FCAImmP7462237_GRCh38-1_2_0'
# '/seq/25871/cellranger/cellranger302_count_25871_FCAImmP7462238_GRCh38-1_2_0'
# '/seq/25871/cellranger/cellranger302_count_25871_FCAImmP7462239_GRCh38-1_2_0'
# '/seq/26109/cellranger/cellranger302_count_26109_FCAImmP7528286_GRCh38-1_2_0'
# '/seq/26109/cellranger/cellranger302_count_26109_FCAImmP7528287_GRCh38-1_2_0'
# '/seq/26109/cellranger/cellranger302_count_26109_FCAImmP7528288_GRCh38-1_2_0'
# '/seq/26087/cellranger/cellranger302_count_26087_FCAImmP7528289_GRCh38-1_2_0'
# '/seq/26087/cellranger/cellranger302_count_26087_FCAImmP7528294_GRCh38-1_2_0'
# '/seq/26322/cellranger/cellranger302_count_26322_FCAImmP7555846_GRCh38-1_2_0'
# '/seq/26280/cellranger/cellranger302_count_26280_FCAImmP7555847_GRCh38-1_2_0'
# '/seq/26321/cellranger/cellranger302_count_26321_FCAImmP7555856_GRCh38-1_2_0'
# '/seq/26321/cellranger/cellranger302_count_26321_FCAImmP7555857_GRCh38-1_2_0'
# '/seq/26421/cellranger/cellranger302_count_26421_FCAImmP7579210_GRCh38-1_2_0'
# '/seq/26421/cellranger/cellranger302_count_26421_FCAImmP7579211_GRCh38-1_2_0'
# '/seq/26441/cellranger/cellranger302_count_26441_FCAImmP7579222_GRCh38-1_2_0'
# '/seq/26431/cellranger/cellranger302_count_26431_FCAImmP7579223_GRCh38-1_2_0'
# '/seq/26431/cellranger/cellranger302_count_26431_FCAImmP7579226_GRCh38-1_2_0'
# '/seq/26431/cellranger/cellranger302_count_26431_FCAImmP7579227_GRCh38-1_2_0'
)



#outDir=Data/fAdrenals_scRNAseq
outDir=Data/fetal_references/fKidneys_REF

for irod_path in ${samples[@]}; do
   sampleName=$(basename $irod_path)
   mkdir -p $outDir/$sampleName
   
   irods iget -r $irod_path/filtered_feature_bc_matrix $outDir/$sampleName/filtered_feature_bc_matrix
   irods iget -r $irod_path/raw_feature_bc_matrix $outDir/$sampleName/raw_feature_bc_matrix
    
#    irods iget -r $irod_path/possorted_genome_bam.bam $outDir/$sampleName/possorted_genome_bam.bam
#    irods iget -r $irod_path/possorted_genome_bam.bam.bai $outDir/$sampleName/possorted_genome_bam.bam.bai



done