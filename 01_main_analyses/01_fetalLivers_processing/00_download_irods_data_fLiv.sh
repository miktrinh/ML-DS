############### TAM / MLDS - scRNAseq - cellranger700 ###############
samples=(
# # Hsb32	Diploid	Male (CD45-/+)
# MY_200531_13679026
# MY_200531_13679027
# MY_200531_13679028
# MY_200531_13679029

# # Hsb31	Diploid	Female (CD45-/+)
# MY_200531_13651140
# MY_200531_13651141
# MY_200531_13651142
# MY_200531_13651143

# # Hsb35	Diploid	Female (CD45-/+)
# MY_200531_14552808
# MY_200531_14552809
# MY_200531_14552810
# MY_200531_14552811

# # Hsb40	Diploid	Female (CD45-/+)
# MY_200531_14816326
# MY_200531_14816327
# MY_200531_14816328
# MY_200531_14816329

# # Hsb33	Trisomy 21	Male (CD45-/+)
# MY_200531_13777913
# MY_200531_13777914
# MY_200531_13777915
# MY_200531_13777916

# #Hsb37	Trisomy 21	Female (CD45-/+)
# MY_200531_14665735
# MY_200531_14665736
# MY_200531_14665737
# MY_200531_14665738

# # Hsb38	Trisomy 21	Female (CD45-/+)
# MY_200531_14784978
# MY_200531_14784979
# MY_200531_14784982
# MY_200531_14784983

# # 15724	Trisomy 21	Male
# MY_200531_10220042
# MY_200531_10220043

# # Hsb36	Trisomy 21	Female (CD45-/+)
# MY_200531_14635737
# MY_200531_14635738
# MY_200531_14635739
# MY_200531_14635740

# # Hsb34	Trisomy 21	Female (CD45-/+)
# MY_200531_13777917
# MY_200531_13777918
# MY_200531_13777919
# MY_200531_13777920

# # 15877	Trisomy 21	Male
# MY_200531_10621689
# MY_200531_10621690

# # Hsb39	Trisomy 21	Female (CD45-/+)
# MY_200531_14784980
# MY_200531_14784981
# MY_200531_14784984
# MY_200531_14784985

# # 16049	Trisomy 18	Male
# MY_200531_11528567
# MY_200531_11528568

# # 15756	Trisomy 18	Female
# MY_200531_10220140
# MY_200531_10220141

# 15806	Trisomy 18	Female
MY_200531_10260920
MY_200531_10260921

# Hsb22	Trisomy 18	Female
MY_200531_12750608
MY_200531_12750609

# 15733	Trisomy 22	Male
MY_200531_10192869
MY_200531_10192870

# 15905	Mono X	Female
MY_200531_10683134
MY_200531_10683135

# Hsb21	Mono X	Female
MY_200531_12408354
MY_200531_12408355

# 15680	Triploid	Female
MY_200531_10043298
MY_200531_10043299
)


outDir=Data/fLivers_scRNAseq
mkdir -p $outDir
for sample in ${samples[@]}; do
   irod_path=$(irods imeta qu -z seq -C sample = "$sample" | sed 's/^collection: //' | head -n 1)
   sampleName=$(basename $irod_path)
   mkdir -p $outDir/$sampleName
   
   irods iget -r $irod_path/filtered_feature_bc_matrix $outDir/$sampleName/filtered_feature_bc_matrix
   irods iget -r $irod_path/raw_feature_bc_matrix $outDir/$sampleName/raw_feature_bc_matrix
    
#    irods iget -r $irod_path/possorted_genome_bam.bam $outDir/$sampleName/possorted_genome_bam.bam
#    irods iget -r $irod_path/possorted_genome_bam.bam.bai $outDir/$sampleName/possorted_genome_bam.bam.bai



done



