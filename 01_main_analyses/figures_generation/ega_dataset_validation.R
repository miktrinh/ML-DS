# Check EGA dataset

# Read in the list of samples provided by Data-Release team
sc = read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/EGAD00001015452_sample_list',sep = '')
wgs = read.delim('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/EGAD00001015453_sample_list',sep = '')

# My current Supplementary Table S1

# Check scRNA-seq dataset ---------
s1 = read_excel('~/lustre_mt22/MLDS_scRNAseq/manuscript_revision_2510/Supplemental table 1.xlsx',skip = 2)
table(s1$`Sample ID` %in% sc$sample_name,s1$`single-cell or bulk RNA-seq`) # All scRNAseq channels are included, including 2 MDS samples L061 and L067
checkmate::assert_true(all(sc$sample_name %in% s1$`Sample ID`))
table(table(sc$sample_name)) # some chanels have >4 fastq files  --> Ask Nathan how to find out what these are... Shall we ask Data-Release team to remove these?
table(sc$sample_name)[table(sc$sample_name) >=6]

# Check WGS dataset ---------------
table(s1$`WGS ID` %in% wgs$supplier_name,s1$`WGS ID` == '-')
s1$`WGS ID`[!s1$`WGS ID` %in% wgs$supplier_name & s1$`WGS ID` != '-'] # PD61857a (L067) is missing from EGA WGS list
checkmate::assert_true(all(wgs$supplier_name %in% s1$`WGS ID`))
table(table(wgs$supplier_name))
table(wgs$supplier_name)[table(wgs$supplier_name) >=4] # Same here... Ask Nathan what these are!

wgs_extra = read.csv('~/lustre_mt22/Aneuploidy/manuscriptDraft_1124/noCC3/WGS_data_info_for_EGA.csv') %>% dplyr::filter(Data_type == 'WGS')
# To ask Data-Release Team:
# Add PD61857a to WGS dataset EGAD00001015453
# Remove 4602STDY7920965 from the sc dataset (as data for P9 should have already been published with Ellie's infant B-ALL paper)
# ?.Remove potential TCR/BCR files (sc dataset) and what are the WGS extra files?
