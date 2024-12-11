# Mi's attempt to detect GATA1 mutation in scRNAseq data

#samtools view -b /lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_45842_SB_Leuk13104279_GRCh38-2020-A/possorted_genome_bam.bam chrX:48791295-48791295 | samtools fillmd -e - /nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa | grep -v "^@"| awk -v pos="48791295" 'BEGIN {OFS = FS = "\t" } ; {n=split($10,a,"") ; if(a[(pos-$4)+1] != "=" ) print pos,(pos-$4)+1, a[(pos-$4)+1], $1, $4, $10 }' > /lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/L019/L019_Leuk13104279_GATA1mut_readID_method2.txt

extract_reads(){
    outDir=$1
    channels=$2
    donorID=$3
    range=$4
    mut_pat=$5
    wt_pat=$6
    
    echo -e "\n$donorID"
    
    mkdir -p $outDir
    for sample in ${channels[@]}; do
        if echo $sample | grep '_MY_' ; then
            sampleName1=$(basename $sample | awk '{gsub(/_GRCh38-2020-A.*$/,"")} 1')
            sampleName2=$(echo $sampleName1 | awk '{gsub(/^.*_MY_/,"MY_")} 1')
            sampleName=$(echo $sampleName2 | awk '{gsub(/_/,".")} 1')
        else
            sampleName=$(basename $sample | awk '{gsub(/^.*SB_|^.*ALeuk_|_GRCh38-2020-A.*$/,"")} 1')
        fi

        echo $sampleName

        # Grep for GATA1 reads
        samtools view $sample/possorted_genome_bam.bam $range | grep "GATA1" | grep "CB\:Z\:" > $outDir/$donorID"_"$sampleName"_GATA1_reads.txt"

        # Grep for mutant reads
        samtools view $sample/possorted_genome_bam.bam $range | cut -f 1,10 | grep $mut_pat | cut -f 1 > $outDir/$donorID"_"$sampleName"_GATA1mut_readID.txt"    
        samtools view $sample/possorted_genome_bam.bam $range | grep -f $outDir/$donorID"_"$sampleName"_GATA1mut_readID.txt" | grep "CB\:Z\:" > $outDir/$donorID"_"$sampleName"_GATA1mut_reads.txt"   


        # Grep for WT reads
        samtools view $sample/possorted_genome_bam.bam $range | cut -f 1,10 | grep $wt_pat | cut -f 1 > $outDir/$donorID"_"$sampleName"_GATA1_WT_readID.txt"
        samtools view $sample/possorted_genome_bam.bam $range | grep -f $outDir/$donorID"_"$sampleName"_GATA1_WT_readID.txt" | grep "CB\:Z\:" > $outDir/$donorID"_"$sampleName"_GATA1_WT_reads.txt"


    done
}




extract_EZH2_reads(){
    outDir=$1
    channels=$2
    donorID=$3
    range=$4
    mut_pat=$5
    wt_pat=$6
    
    echo -e "\n$donorID"
    
    mkdir -p $outDir
    for sample in ${channels[@]}; do

        sampleName=$(basename $sample | awk '{gsub(/^.*SB_|_GRCh38-2020-A.*$/,"")} 1')
        echo $sampleName
        echo $sample
        # Grep for GATA1 reads
        samtools view $sample/possorted_genome_bam.bam $range | grep "EZH2" | grep "CB\:Z\:" > $outDir/$donorID"_"$sampleName"_EZH2_reads.txt"

        # Grep for mutant reads
        samtools view $sample/possorted_genome_bam.bam $range | cut -f 1,10 | grep $mut_pat | cut -f 1 > $outDir/$donorID"_"$sampleName"_EZH2mut_readID.txt"    
        samtools view $sample/possorted_genome_bam.bam $range | grep -f $outDir/$donorID"_"$sampleName"_EZH2mut_readID.txt" | grep "CB\:Z\:" > $outDir/$donorID"_"$sampleName"_EZH2mut_reads.txt"   


        # Grep for WT reads
        samtools view $sample/possorted_genome_bam.bam $range | cut -f 1,10 | grep $wt_pat | cut -f 1 > $outDir/$donorID"_"$sampleName"_EZH2_WT_readID.txt"
        samtools view $sample/possorted_genome_bam.bam $range | grep -f $outDir/$donorID"_"$sampleName"_EZH2_WT_readID.txt" | grep "CB\:Z\:" > $outDir/$donorID"_"$sampleName"_EZH2_WT_reads.txt"


    done
}







# allowing position being -/+ 70 bases either end


# L039 - with 11bp insertion T.GGCCTACTACA.GG after chrX:48791288

donorID='L039'
mut_pat='TGGCCTACTACAGGC\|GGCCTACTACAGGC\|GCCTACTACAGGC\|CCTACTACAGGC\|CTACTACAGGCC\|TACTACAGGCCT\|ACTACAGGCCT\|CTACAGGCCT\|TACAGGCCTACTACAGGG\|ACAGGCCTACTACAGGG\|CAGGCCTACTACAGGG\|AGGCCTACTACAGGG'
wt_pat='TGGCCTACTACAGGGA\|TGGCCTACTACAGGGACG\|TGGCCTACTACAGGG'
#range="chrX:48791178-48791398"
range="chrX:48791218-48791358"

channels=('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234194_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234195_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234196_GRCh38-2020-A'

'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415925_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415926_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415927_GRCh38-2020-A'

'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_45842_SB_Leuk13104278_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat





#+/- 70 bp either end
## L019 with point mutation at 48791295
donorID='L019'
mut_pat='CCTAGTAC\|GCACTGGCCTAGTACAGGGACGCT\|GCACTGGCCTAGTAC\|GCACTGGCCTAGTA\|GCACTGGCCTAGT\|GCACTGGCCTAG\|CTGGCCTAGTACAGGGACGCT\|GGCCTAGTACAGGGACGCT\|CTAGTACAGGGACGCT\|TAGTACAGGGACGCT\|AGTACAGGGACGCT\|GTACAGGGACGCT\|GTACAGGGACGCT'
wt_pat='CCTACTAC\|GCACTGGCCTAC\|CTACAGGGACGCT'
range="chrX:48791225-48791365"

channels=('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_45842_SB_Leuk13104279_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47089_SB_Leuk13645528_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47089_SB_Leuk13645529_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat






## L038 with deletion at chrX: 48791248..48791259
donorID='L038'
mut_pat='TCCTAGCACAG\|GCTTCCTAGCACA\|GCTTCCTAGCAC\|GCTTCCTAGCA\|GCTTCCTAGC\|GCAGCTTCCTAGCA\|TCCTAGCACAGCC'
wt_pat='CCTCCACTGCCCCGAGCACAG\|CACTGCCCCGAGCAC\|CTGCCCCGAGCAC\|CCACTGCCCCG\|GCTTCCTCCA\|GCAGCTTCCTCCA\|GCAGCTTCCTCC\|GCAGCTTCCTC'
#range="chrX:48791138-48791369"
range="chrX:48791178-48791329"

echo -e "\n$donorID"

channels=('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234191_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234192_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234193_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415922_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415923_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415924_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID



# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat
#TGCATCCA C CACAAAATC
EZH2_mut_pat='TGCATCCATCACAAAATC\|TGCATCCATCACAAAA\|TGCATCCATCACAA\|CATCCATCACAAAA\|ATCCATCACAAAATC\|ATCCATCACAAAA'
EZH2_wt_pat='TGCATCCACCACAAAATC\|TGCATCCACCACAAAA\|TGCATCCACCACAA\|CATCCACCACAAAA\|ATCCACCACAAAATC\|ATCCACCACAAAA'
EZH2_range="chr7:148809312-148809452"
#extract_EZH2_reads $outDir $channels $donorID $EZH2_range $EZH2_mut_pat $EZH2_wt_pat








## L042 with point insertion after 48791179
donorID='L042'
mut_pat='CTCTGGTTGT\|CCTGCTCTGGTT\|TTGTCCTCCACA'
wt_pat='CTCTGGTGT\|GTGTCCTCCACACCA\|TCCTGCTCTGGTG'
#range="chrX:48791069-48791289"
range="chrX:48791109-48791249"


channels=('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234203_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234204_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234205_GRCh38-2020-A'

'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415934_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415935_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415936_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID

# extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat 
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat






## L040 with insertion mutation at 48791296
donorID='L040'
mut_pat='CCTACTTAC\|ACTTACAGGGACGCT\|CTTACAGGGACGCT\|TTACAGGGACGCTGAG\|GCGGCACTGGCCTACTTA\|GCGGCACTGGCCTACTT'
wt_pat='CCTACTAC\|GCGGCACTGGCCTACTA\|CTACAGGGACGCTGAG'
#range="chrX:48791186-48791406"
range="chrX:48791226-48791366"

channels=('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234197_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234198_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46351_SB_Leuk13234199_GRCh38-2020-A'

'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415928_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415929_GRCh38-2020-A'
'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_46682_SB_Leuk13415930_GRCh38-2020-A'

'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47089_SB_Leuk13645524_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat






## L076 with insertion mutation at  c.160_161insGCTCA after 48791269
donorID='L076'
mut_pat='AGCCAGCTCACCGCT\|GCCAGCTCACCGCT\|AGCCAGCTCACCGC\|CCGAGCACAGCCAG\|TCACCGCTGCAGCTGCG'
wt_pat='AGCCACCGCT\|CCGAGCACAGCCACCG\|GCCACCGCTGCAGCTGCG'
#range="chrX:48791159-48791379"
range="chrX:48791199-48791339"

channels=(#'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47260_SB_Leuk13697519_GRCh38-2020-A'
#'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47580_SB_Leuk13760338_GRCh38-2020-A'

#'/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47580_SB_Leuk13760338_GRCh38-2020-A'
#'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14635833_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49031_ALeuk_RNA14832000_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49031_ALeuk_RNA14832001_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat








## L075 with deletion/insertion mutations at c.31_45delinsT 48791142 - 48791158
donorID='L075'
mut_pat='TCCCTGGGGTCCCCAG\|CCTGGGGTCCCCAGTTT\|TGGGGTCCCCAGTTT\|GGGTCCCCAGTTT\|GTCCCCAGTTTGTG\|CTGGGTCCCTGGGGT\|TCCCCAGTTTGTG'
wt_pat='TCCCTGGGGACCTCAGAGCCCCTCCCCCAG\|CCTGGGGACCTCAGAGCCCCTCCCCCAG\|GGACCTCAGAGCCCCTCCCCCAG\|GGACCTCAGAGCCCCTCCCC\|CTGGGTCCCTGGGGA\|TCCCCCAGTTTGTG'
#range="chrX:48791032-48791268"
range="chrX:48791072-48791228"

channels=('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47260_SB_Leuk13697518_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat







## L091 with GATA1 c.113_128del chrX:48791222-48791236
donorID='L091'
mut_pat='TCTGGGCAGCAGCT\|CTGGGCAGCAGCT\|TCTGGGCAGCA\|TCTGGGCAG\|TGGGCAGCAGCTTCT\|GGCAGCAGCTTCT\|GCAGCAGCTTCT\|CAGCAGCTTCTTCCACT\|CCCTCTGGGCAGCAG\|CCCTCTGGGCAGCAGC\|CCCTCTGGGCAGCA\|CCCTCTGGGCAGC\|CCCTCTGGGCAG\|TTCCCCTCTGGGCA\|GGCAGCAGCTTCCTCC'
wt_pat='GGCCTGAGGGCTTGGATGCAGC\|CCTGAGGGCTTGGATGCAGC\|GGCCTGAGGGCTTGGATGCA\|CCTGAGGGCTTGGATGCAGCAGCT\|TGCAGCAGCTTCCTCC'
#range="chrX:48791112-48791346"
range="chrX:48791152-48791306"

channels=('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/MLDS_GEX/cellranger700_count_47510_SB_Leuk13778112_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat



## L041






## CC2 with insertion mutation of C after 48791222
donorID='CC2'
mut_pat='CCCTCTGGGCCCTGAGGGCTTGGA\|CCCTCTGGGCCCTGAGGGCTTG\|CCCTCTGGGCCCTGAGGGCT\|CCCTCTGGGCCCTGAGGG\|TTCTTCCCCTCTGGGCCCTGA\|TTCTTCCCCTCTGGGCCCTG\|TTCTTCCCCTCTGGGCCCT\|TTCTTCCCCTCTGGGCCC'
wt_pat='CCCTCTGGGCCTGAGGGCTTGGA\|CCCTCTGGGCCTGAGGGCTTG\|CCCTCTGGGCCTGAGGGCT\|CCCTCTGGGCCTGAGGG\|TTCTTCCCCTCTGGGCCTGA\|TTCTTCCCCTCTGGGCCTG\|TTCTTCCCCTCTGGGCCT'

range="chrX:48791152-48791292"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48097_MY_200531_14406282_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48097_MY_200531_14406283_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat






## L156 with 16bp duplication [c158_173 duplication p.(Ala59HisfsTer14)], starting at chrX:48791282

donorID='L156'
mut_pat='ACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCCACCGCTGCAGCTGCGGCACTGGCCTACTAC\|ATGCAGCAGCTTCCTCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCCACCGCTGCAGCTGC\|ATGCAGCAGCTTCCTCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCCACCGCTGCAGC\|ATGCAGCAGCTTCCTCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCCACCGCTGC\|ACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCCACCGCT\|ACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCCACCG\|TCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCCAC\|TCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCCAC\|TGCCCACCGCTGCAGCTGCGGCACTGGCCTACTACAGG\|AGCTGCCCACCGCTGCAGCTGCGGCACTGGCCTACTACAGG\|TGCAGCTGCCCACCGCTGCAGCTGCGGCACTGGCCTACTAC\|CGCTGCAGCTGCCCACCGCTGCAGCTGCGGCACTGGCC\|CCACCGCTGCAGCTGCCCACCGCTGCAGCTGCGGCACTGGCC\|GCACAGCCACCGCTGCAGCTGCCCACCGCTGCAGCTGCGGCACTGGCC\|CAGCCACCGCTGCAGCTGCCCACCGCTGCAGCTGCGGCACTGGCC\|CAGCCACCGCTGCAGCTGCCCACCGCTGCAGCTGCGGCA\|ATGCAGCAGCTTCCTCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCCC'

wt_pat='ACTGCCCCGAGCACAGCCACCGCTGCAGCTGCGGCACTGGCCTACTACAGG\|GCCCCGAGCACAGCCACCGCTGCAGCTGCGGCACTGGCCTACTAC\|CCGAGCACAGCCACCGCTGCAGCTGCGGCACTGGCC\|AGCACAGCCACCGCTGCAGCTGCGGCACTG\|CAGCCACCGCTGCAGCTGCGGCACTGGCCTACTACAGG\|ACTGCCCCGAGCACAGCCACCGCTGCAGCTGCGGCA\|ATGCAGCAGCTTCCTCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCG'

range="chrX:48791232-48791350"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48209_ALeuk_RNA14564480_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48209_ALeuk_RNA14564481_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48983_MY_200531_14784986_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48983_MY_200531_14784987_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat







## CC3 with 18bp deletion (chrX: 48,791,301-48,791,318) p.aa.60-64 
## Note: [c.227_230dup accroding to Henning's team. However, no coverage here in scRNAseq data]

donorID='CC3'
mut_pat='GCCTACTACAGACACTCCCCAGTCTTT\|GCCTACTACAGACACTCC'

wt_pat='GCCTACTACAGGGACGCTGAGGCCTACAGA\|TACAGACACTCCCCAGTCTTTCAGGTG\|AGGGACGCTGAGGCCTACAGACACTCCCCA\|GGGACGCTGAGGCCTACAGACACTCCCCAGTC'

range="chrX:48791231-48791388"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14635834_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat





## CC4 with 4bp insertion (after chrX: 48,791,240) c.128_131dup
donorID='CC4'
mut_pat='GCTTGGATGCAGCCAGCAGCTTCCTCCACTGCC\|CCAGCAGCTTCCTCCACTGCC\|TCTGGGCCTGAGGGCTTGGATGCAGCC'

wt_pat='GCCTGAGGGCTTGGATGCAGCAGCTTCCTCCACTGCC\|GCAGCAGCTTCCTCCACTGCC'

range="chrX:48791170-48791320"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14635835_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat




## CC5 with 11bp duplication (after chrX: 48,791,283) c.164_174dup
donorID='CC5'
mut_pat='ACCGCTGCAGCTGCGCTGCAGCTGCGGCACTGGCCTACTACAGG\|GCGCTGCAGCTGCGGCACTGGCCTACTACAGGGACGCT\|CCGAGCACAGCCACCGCTGCAGCTGCGC'

wt_pat='CCGAGCACAGCCACCGCTGCAGCTGCGG\|CCGCTGCAGCTGCGGCACTGGCCTACTACAGGGACGCT\|GCCACCGCTGCAGCTGCGGCACTGGCCTACTACAGGGAC'

range="chrX:48791213-48791353"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14635836_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14635837_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat





## CC6 with 31bp deletion (chrX: 48,791,228 – 48,791,258) c125_150del
donorID='CC6'
mut_pat='TCAGGGGTTTTCTTCCCCTCTGGGCCTGAGGGA\|GGAGCACAGCCACCGCTGCAGCTGCG\|GGGCCTGAGGGAGCACAGCCA'

wt_pat='TCAGGGGTTTTCTTCCCCTCTGGGCCTGAGGGC\|GAGGGCTTGGATGCAGCAGCTTCCTCCACTGCCCCGAGCACAGC\|CGAGCACAGCCACCGCTGCAGCTGCGGCACTGGCCT'

range="chrX:48791158-48791328"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14637036_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14637037_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat



## CC7 with 10bp insertion (after chrX: 48,791,298) c.188_189insGTGCCTACT
donorID='CC7'
mut_pat='GCAGCTGCGGCACTGGCCTACTAGT\|CTGGCCTACTAGTGCCTACTACAGGGACGCT\|TAGTGCCTACTACAGGGACGCTGAGGCCTACAGACACTCCCCA'

wt_pat='GCAGCTGCGGCACTGGCCTACTACA\|GCTGCGGCACTGGCCTACTACAGGGACGCTGAGGCC\|CTGGCCTACTACAGGGACGCTGAGGCCTACAGACACTCCCCA'

range="chrX:48791228-48791368"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14637038_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14637039_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat




## CC8 with 26bp deletion (after chrX: 48,791,287 - 48,791,313) c.180_205del
donorID='CC8'
mut_pat='ACAGCCACCGCTGCAGCTGCGGCACTAC\|CGGCACTACAGACACTCCCCAG\|CACTACAGACACTCCCCAGTCTTTCA'

wt_pat='GCACTGGCCTACTAGTGCCTACTACAGGGACGCT\|TGCGGCACTGGCCTACTACAGGGACGCTGAGGCCTACAGA\|GCCTACAGACACTCCCCAGTCTTTCA'

range="chrX:48791217-48791383"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14637040_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_48665_MY_200531_14637041_GRCh38-2020-A'
)

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat






## L067 - MDS 11bp insertion c.153_154insGCAGCTGCAGC (after chrX: 48,791,262) 
donorID='L067'
mut_pat='TTCCTCCACTGCCCCGAGCGCAGCTGCAGCACAGCCACCGCTGCAGCTGCG\|TTCCTCCACTGCCCCGAGCG\|CAGCACAGCCACCGCTGCAGCTGCG'

wt_pat='TTCCTCCACTGCCCCGAGCACAGCCACCGCTGCAGCTGCG\|TTCCTCCACTGCCCCGAGCA\|GAGCACAGCCACCGCTGCAGCTGCG'

range="chrX:48791192-48791332"

channels=('/lustre/scratch126/casm/team274sb/bl10/B-ALL/Data/GEX/cellranger700_count_47089_SB_Leuk13645530_GRCh38-2020-A')

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat



## L178 - MLDS 1bp insertion insG (after chrX: 48,791,112) 
donorID='L178'
mut_pat='AATCCCCAGAGGCTCCATGGGAGTTCCCTGGCCTGGGGTCC\|GGGAGTTCCCTGGCCTGGGGTCCCTGGGGA\|AATCCCCAGAGGCTCCATGGG'

wt_pat='AATCCCCAGAGGCTCCATGGAGTTCCCTGGCCTGGGGTCC\|TGGAGTTCCCTGGCCTGGGGTCCCTGGGGA\|AATCCCCAGAGGCTCCATGGA'

range="chrX:48791012-48791212"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49031_ALeuk_RNA14831998_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49031_ALeuk_RNA14831999_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49226_ALeuk_RNA14910210_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49226_ALeuk_RNA14910211_GRCh38-2020-A')

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat



## L182 - MLDS 2bp deletion AG (chrX: 48,791,199 - 48,791,200) 
donorID='L182'
mut_pat='TGTCCTCCACACCAGAATCGGGTTTTCTTCCCCTCTGG\|TGGATCCTGCTCTGGTGTCCTCCACACCAGAATCG\|CGGGTTTTCTTCCCCTCTGGGCCTGAGGGCTTG'

wt_pat='TGTCCTCCACACCAGAATCAGGGGTTTTCTTCCCCTCTGG\|TGGATCCTGCTCTGGTGTCCTCCACACCAGAATCA\|GGGGTTTTCTTCCCCTCTGGGCCTGAGGGCTTG'

range="chrX:48791099-48791300"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49087_ALeuk_RNA14872269_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49087_ALeuk_RNA14872270_GRCh38-2020-A')

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat



## L114 - TAM
donorID='L114'
mut_pat=''

wt_pat=''

range="chrX:48791099-48791300"

channels=('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49087_ALeuk_RNA14872269_GRCh38-2020-A'
'/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/MLDS_scRNAseq_Data/cellranger700_count_49087_ALeuk_RNA14872270_GRCh38-2020-A')

outDir='/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/GATA1_reads/scRNAseq_GATA1mut_reads/'$donorID


# Run extract reads
#extract_reads $outDir $channels $donorID $range $mut_pat $wt_pat
