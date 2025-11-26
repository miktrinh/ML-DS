## TAM / ML-DS Mutation analysis

##--------------##
##      WGS   ####
##--------------##

outDir = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}



library(VariantAnnotation)
library(tidyverse)
source('~/lustre_mt22/generalScripts/caveman_readvcf.R')
source('~/lustre_mt22/generalScripts/binom_mix_model_Tim.R')
source('~/lustre_mt22/generalScripts/utils/wgs_analysis_helperFunctions.R')
source('~/lustre_mt22/generalScripts/utils/misc.R')


##------------------------------##
##      Global paraneters     ####
##------------------------------##

variantImpact_toKeep = c( "missense","ess_splice","inframe","splice_region", "cds_disrupted", "nonsense","frameshift","start loss","start lost","start_lost","start_loss")
TSG_impact = variantImpact_toKeep[!variantImpact_toKeep %in% c("missense", "inframe","splice_region")]

##--- Import list of ML-DS driver genes -------##
potential_driver_genes = c("GATA1", "STAG2", "RAD21", "CTCF", "SMC1A", "NIPBL", "EZH2", "GSE1", "KANSL1", "KMT2C",
                           "SUZ12",'EED', "EP300", "KDM6A", "BCOR",
                           "JAK2", "JAK3", "JAK1", "EPOR","SH2B3", "MPL","STAT5B","GNB1",
                           "NRAS", "KRAS","PTPN11","BRAF",
                           'DCAF7','CSF2RB','TP53','CDKN2A',
                           "DAZAP1", "MBNL1", "SRSF2",
                           "ZBTB7A", "RUNX1",'IRX1','IRF2',"NFE2","NFIA","IKZF1", "WT1","MYC","GATA2")
labuhn_drivers = data.frame(Gene = c('SMC3','KMT2A','KMT2C','NAT6','TET2','DNMT3A','PTPRD','KIT','PTEN','NF1','CREBBP','ATRX','SF3B1','ITGB1','DLEC1','IL3'),
                            Method = 'Labuhn_etal_2019')
sato_2024_MLDS_driver_genes =readxl::read_excel('~/lustre_mt22/generalResources/Sato_etal_2024_MLDS_driverLandscape.xlsx',sheet = 'Table 9',skip = 2)
sato_2024_MLDS_driver_genes$Method = paste0(sato_2024_MLDS_driver_genes$Method,' (Sato_etal_2024)')
driverGenes = rbind(labuhn_drivers,sato_2024_MLDS_driver_genes)

write.csv(driverGenes,'~/lustre_mt22/generalResources/MLDS_drivers.csv',row.names = F)
## Combine both resouces
driverGenes = unique(c(potential_driver_genes,sato_2024_MLDS_driver_genes$Gene,labuhn_drivers$Gene))

cosmicGenes = read.delim('~/lustre_mt22/generalResources/COSMIC_v100_202408/Cosmic_CancerGeneCensus_v100_GRCh38.tsv.gz',sep = '\t')
cosmic = read.delim('~/lustre_mt22/generalResources/COSMIC_v100_202408/Cosmic_MutantCensus_v100_GRCh38.tsv.gz',sep = '\t')

## Match analysis
sampleMani = data.frame(PDID=c('PD60302a','PD61847a','PD58851a','PD61854a','PD60303a','PD60304a'),
                        projectid=c(3030,3030,3030,3030,3030,3030),
                        timePoint=c('Diagnostic','TP1','Diagnostic','TP2','Diagnostic','Diagnostic'),
                        donorID=c('L039','L039','L019','L019','L040','L042'))



##--------------------------##
##      Subs - Indels     ####
##--------------------------##

# 1. Subs 
sub = import_multiple_caveman(sampleMani,sampleID_name='PDID',ASMD_CLPM_filter=F,filter=T)
sub$varClass = 'sub'

# 2. Indels
indels = import_multiple_pindel(sampleMani,sampleID_name='PDID',intv=10,filter_lowQual = T)
indels$varClass = 'indel'
table(indels$Impact)

# Combine subs and indels
shared_columns = intersect(colnames(sub),colnames(mcols(indels)))
sub_indels = rbind(sub[,shared_columns],as.data.frame(mcols(indels))[,shared_columns])
table(sub_indels$PDID,sub_indels$varClass)

# assign groups, with preference given to TSG then Wilms tumour genes
sub_indels$group = ifelse(sub_indels$Gene %in% driverGenes,'knownDriverGene',
                          ifelse(sub_indels$Gene %in% cosmicGenes$GENE_SYMBOL[cosmicGenes$ROLE_IN_CANCER !=''],paste0('cosmic_',cosmicGenes$ROLE_IN_CANCER[match(sub_indels$Gene,cosmicGenes$GENE_SYMBOL)]),
                                 ifelse(sub_indels$Gene %in% cosmicGenes$GENE_SYMBOL[cosmicGenes$ROLE_IN_CANCER ==''],'cosmicGene_unknownRole','unknown')))
                                 

# keep TSG mutations if functional and remove splice region mutation in OG (WHY?? --> ask Taryn)
sub_indels$isDriver = F
sub_indels$isDriver[sub_indels$group == 'knownDriverGene' & sub_indels$Impact %in% variantImpact_toKeep] = T
sub_indels$isDriver[sub_indels$group %in% c('cosmic_TSG','cosmic_TSG, fusion') & sub_indels$Impact %in% TSG_impact] = T

sub_indels$comment = 'non-functional'
sub_indels$comment[sub_indels$isDriver] = sub_indels$Impact[sub_indels$isDriver]



##------------------------------##
##      Hotspot mutations     ####
##------------------------------##
# 3. Hotspots 
# check if hotspot in either missense in TSG or oncogene or unclassified
subIndels_toCheck = rbind(sub_indels[grepl('TSG',sub_indels$group) & !sub_indels$Impact %in% TSG_impact & sub_indels$isDriver==F & sub_indels$Gene %in% cosmic$GENE_SYMBOL & sub_indels$AAchange != '-',],
                          sub_indels[grepl('oncogene',sub_indels$group) & sub_indels$Impact != "splice_region" & sub_indels$isDriver==F & sub_indels$Gene %in% cosmic$GENE_SYMBOL & sub_indels$AAchange != '-',],
                          sub_indels[grepl('unclassified|unknown|cosmicGene_unknownRole',sub_indels$group) & sub_indels$isDriver==F & sub_indels$Gene %in% cosmic$GENE_SYMBOL & sub_indels$AAchange != '-',])  

if(nrow(subIndels_toCheck) > 0){
  ## Check those with aa change
  
  subIndels_toCheck$aa_pos = paste0(subIndels_toCheck$Gene,'_',gsub('^p.|\\S$','',subIndels_toCheck$AAchange))
  cosmic_toCheck = cosmic[cosmic$GENE_SYMBOL %in% subIndels_toCheck$Gene[!subIndels_toCheck$AAchange %in% c('-','p.?')],]  
  cosmic_toCheck$aa_pos = paste0(cosmic_toCheck$GENE_SYMBOL,'_',gsub('^p.|\\S$','',cosmic_toCheck$MUTATION_AA))
  cosmic_toCheck = cosmic_toCheck[cosmic_toCheck$aa_pos %in% subIndels_toCheck$aa_pos,]
  if(nrow(cosmic_toCheck) > 0 ){
    message('Please continue working on this')
    table(cosmic_toCheck$GENE_SYMBOL)
  }else{
    message('No matching mutations found in COSMIC')
  }
  
  
}else{
  message('No qualified sub/indel to check')
}



##------------------------------##
##    Ess_splice mutations    ####
##------------------------------##
# 4.Ess_splice 
# check on varsome, can also us MTBP if not obviously benign on varsome
varsome = sub_indels[sub_indels$Impact %in% c("ess_splice", "splice_region") & sub_indels$isDriver == F,]
varsome = print(gsub("_", " ", varsome))

# remove all ess_ and splice mutations once have checked benign in varsome
# BE CAREFUL TO CHECK BEFORE THIS LINE
sub_indels = sub_indels[!c(sub_indels$Impact == "splice_region" | sub_indels$Impact == "ess_splice"),]



##-----------------------##
##    TERT mutations   ####
##-----------------------##

# 5. TERT ####
#check for TERT as per TERT script
# hotspot mutations 
promoter = c(1295228,1295250,1295242,1295243,1295228,1295229)

# read Caveman outputs, not limited to PASS variants
sub_unfiltered = import_multiple_caveman(sampleMani,sampleID_name='PDID',ASMD_CLPM_filter=F,filter=F)
table(sub_unfiltered$Pos %in% promoter & sub_unfiltered$Chr == 'chr5')






## Save FULL output
sub_indels = sub_indels[,c('donorID','sampleID','PDID','timePoint','varClass','varID','Gene','Impact','AAchange','Type','group','isDriver','comment','Chr','Pos','Ref','Alt','ID','Filter','projectid')]
write.csv(sub_indels,'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/wgs_mutations.driversCalling_MatchedAnalysis.csv')

driver = sub_indels[sub_indels$isDriver == T,]
table(driver$PDID,driver$varClass)

##-------------##
## Check the list of non-driver mutations to see if there might be some important genes in there ####
# View(table(sub_indels$Gene[sub_indels$isDriver ==F]))
# PTPRN2 - leukaemia genes, something to do with methylation
# KCNIP4 - associated with ML-DS - hypermethylated
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3790517/


mlds_srat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns.RDS'
mlds_mdat_fp = '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/1_scProcessing_leukaemia/MLDS/MLDS_clean_annotated_2408_noUnknowns_mdat.csv'
library(Seurat)
sub_indels[!is.na(sub_indels$Gene) & sub_indels$Gene == 'ROBO2',]

FeaturePlot(mlds,c('PTPRN2','ROBO2','KCNIP4','ACSL1'))
DimPlot(mlds,group.by='donorID',label=T,cols = col25)


##---------------------------------------------------##
##                targetted Nanoseq                ####
##---------------------------------------------------##

##--- 1. import all variants called by tNS (VAF > 0.1)
library(VariantAnnotation)
## 1. Proportion of SNPs being dbSNPs
sampleMani_tNS = data.frame(PDID=c('PD60302a','PD61847a','PD60302c','PD58851a','PD61854a','PD61856a','PD60303a','PD61851a','PD60304a','PD61850a',
                                   'PD60301a','PD61846a','PD62336a','PD62331c','PD62331a'),
                        projectid=c(3327,3327,3327,3327,3327,3327,3327,3327,3327,3327,3327,3327,3327,3327,3327),
                        timePoint=c('Diagnostic','TP1','TP4','Diagnostic','TP2','TP4','Diagnostic','TP4','Diagnostic','TP1','Diagnostic','TP1','Diagnostic','Diagnostic','Diagnostic'),
                        donorID=c('L039','L039','L039','L019','L019','L019','L040','L040','L042','L042','L038','L038','L075','L076','L076'),
                        normal_reference = c('PD60302c','PD60302c','PD60302c','PD61856a','PD61856a','PD61856a','PD61851a','PD61851a','PD61850a','PD61850a',
                                             'PD60301a','PD61846a','PD62336a','PD62331c','PD62331a'))
sampleMani_tNS$canapp_projectID = sampleMani_tNS$projectid
sampleMani_tNS = cbind(sampleMani_tNS,mlds@meta.data[match(sampleMani_tNS$donorID,mlds$donorID),c('sex','age_yrs')])
sampleMani_tNS$PDID = paste0(sampleMani_tNS$PDID,'_tds0001')
sampleMani_tNS$normal_reference = paste0(sampleMani_tNS$normal_reference,'_tds0001')
write.csv(sampleMani_tNS[,c('donorID','sex','age_yrs','timePoint','PDID','normal_reference','projectid')],'~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/tNS_MLDS_manifest.csv',row.names = F)

cmd = ('perl ~/lustre_mt22/generalScripts/annovar/table_annovar.pl %s ~/lustre_mt22/generalScripts/annovar/humandb -buildver hg19 -out %s -remove -protocol refGene -operation g -nastring . -polish -vcfinput')


# 
# '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/tNS_allVariantsmutations.bed',
# '~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/Results/x10.1_WGS_variantAnalysis/tNS/tNS_allVariantsmutations_annot'
outDir_tNS = file.path(outDir,'tNS')
if(!dir.exists(outDir_tNS)){
  dir.create(outDir_tNS,recursive = T)
}

all_var = data.frame()
for(sample in sampleMani_tNS$PDID){
  projectid = sampleMani_tNS$projectid[sampleMani_tNS$PDID == sample]
  sample = paste0(sample,'_tds0001')
  vcf_fp = paste0(path,projectid,"/",sample,"/",sample,".tnanoseq.vcf.gz")
  txt_annot_fp = file.path(outDir,'tNS',paste0(sample,'_allVariants_annotVCF.hg19_multianno.txt'))
  if(!file.exists(txt_annot_fp)){
    vcf_annot_fp = file.path(outDir,'tNS',paste0(sample,'_allVariants_annotVCF'))
    cmd_toRun = sprintf(cmd,vcf_fp,vcf_annot_fp)
    system(cmd_toRun)  
  }
  
  # Import annotation
  annot = read.delim(txt_annot_fp,header = T,sep = '\t')
  annot$varID = paste0(annot$Otherinfo4,':',annot$Otherinfo5,'_',annot$Otherinfo7,'/',annot$Otherinfo8)
  
  data = readVcf(vcf_fp)
  var = rowRanges(data)
  info = info(data)
  mcols(var) = cbind(mcols(var),info[match(names(var),rownames(info)),])
  var$varID = as.character(names(var))
  var = as.data.frame(mcols(var))
  var$sample = sample
  var$ALT = sapply(var$ALT,as.character)
  ## Check that all variants had annotation
  print(nrow(var) == nrow(annot))
  table(annot$varID %in% var$varID)
  var = cbind(var,annot[match(var$varID,annot$varID),grepl('refGene',colnames(annot))])
  all_var = rbind(all_var,var)
}

all_var = cbind(all_var,sampleMani_tNS[match(all_var$sample,paste0(sampleMani_tNS$PDID,'_tds0001')),!colnames(sampleMani_tNS) %in% colnames(all_var)])

write.csv(all_var,file.path(outDir_tNS,'tNS_allVariants_allSamples_annot.csv'))





















output_unclassified_check = sub_indels[sub_indels$group == "unclassified" & sub_indels$Gene %in% cosmic$Gene.Symbol,]
output_OG = sub_indels[sub_indels$Gene %in% OG & !sub_indels$Impact =="splice_region",,]
output_OG = rbind(output_TSG_check, output_OG, output_unclassified_check)

# outputs hotspot plot to review and dataframe
OG_check = c()
check=c()
pdf("hotspot.pdf")
for (i in 1:nrow(subIndels_toCheck)){
  gene = subIndels_toCheck[i,"Gene"]
  AA =  gsub('^p.','',subIndels_toCheck[i,"AAchange"])
  type = subIndels_toCheck[i,"varClass"]
  if(type == "sub"){
    number = as.numeric(gsub("[^0-9.-]", "", AA))
    first = substr(AA, 1,1)
    spot = paste0(first,number)
  }else{
    spot = gsub("*fs.*", "", gsub("*del.*", "", gsub("*_.*", "", AA)))
  }
  cosmic3 = cosmic[cosmic$GENE_SYMBOL == gene,]
  cosmic3 = data.frame(table(cosmic3$MUTATION_AA))
  cosmic3$AA = as.numeric(gsub("[^[:digit:], ]", "", cosmic3$Var))
  temp = cosmic3[grepl(spot, cosmic3$Var1),] 
  cosmic3 = na.omit(cosmic3[!grepl("p.[*]", cosmic3$Var1) & !grepl("fs", cosmic3$Var1)  & !grepl("_", cosmic3$Var1),])
  cosmic5=c()
  for(j in unique(cosmic3$AA)){
    cosmic4 = na.omit(cosmic3[cosmic3$AA == j,])
    if(nrow(cosmic4)>1){
      sum = sum(cosmic4$Freq)
    }else{sum = cosmic4$Freq}
    cosmic4 = data.frame("NA", sum, j)
    cosmic5 = rbind(cosmic5,cosmic4)
  }
  df = data.frame("AA", 1, number)
  colnames(df) = colnames(cosmic5) = c("Var1", "Freq", "AA")
  df = rbind(df, cosmic5)
  df$AA = as.numeric(df$AA)
  df$col = ifelse(df$Var1 == "AA", "withcolor", NA)
  # review hotspot plot, remove any that are not the hotspot 
  p<-ggplot(data=df, aes(x=AA, y=Freq, fill=col)) +
    geom_bar(stat="identity", width=4, show.legend=FALSE) +
    ggtitle(paste0(gene, " ", AA) ) + theme_bw() +
    geom_segment(aes(y=max(Freq), x=number, yend=max(Freq)-5, xend=number), arrow = arrow(length=unit(0.5, 'cm')),      color='red')
  print(p)
  if (nrow(temp)>0){
    if(sum(temp$Freq) >2){
      check = output_OG[output_OG$Gene == gene & output_OG$AAchange == AA ,]
      check$Freq = sum(temp$Freq)
    }
  }
  print(gene)
  OG_check = rbind(OG_check, check)
}
dev.off()









##------------------------------------------------##
##    Import Caveman output - Diagnostic        ####
##------------------------------------------------##
d_caveman = caveman2(sample='PD62331c',projectid = 3030)
d_caveman$timepoint = 'Diagnostic'
d_caveman = d_caveman[d_caveman$Flag == T,]
d_caveman$snv_ID = paste0(d_caveman$Chr,':',d_caveman$Pos,'_',d_caveman$Ref,'/',d_caveman$Alt)
View(d_caveman[d_caveman$Gene %in% potential_driver_genes,])
table(d_caveman[d_caveman$Gene %in% potential_driver_genes,]$Gene,
      d_caveman[d_caveman$Gene %in% potential_driver_genes,]$Impact
      )
## Prepare "mutations" table for dndscv: 5 columns in the following order: sampleID, chr, pos, ref, alt
wgs_mutations = d_caveman[,c('Chr','Pos','Ref','Alt')]
wgs_mutations$sampleID = 'PD62331c'
colnames(wgs_mutations) = c('chr', 'pos', 'ref', 'alt','sampleID')
wgs_mutations = wgs_mutations[,c('sampleID', 'chr', 'pos', 'ref', 'alt')]
wgs_mutations$chr = gsub('chr','',wgs_mutations$chr)
wgs_mutations = dndscv(mutations = wgs_mutations,cv=NULL,max_muts_per_gene_per_sample = 1e9,
                    refdb = '~/lustre_mt22/Aneuploidy/Results/16_GATA1_variants/RefCDS_human_GRCh38_GencodeV18_recommended.rda')


# Including FAILED Pindel variants too

intv = 10

# pindel filter
library(data.table)
indels = pindel(sample = 'PD62331c',projectid=3030,filter_lowQual=F)

# Read in t-NS result
tNS_var = readVcf('/nfs/cancer_ref01/nst_links/live/3327/PD60302c_tds0001/PD60302c_tds0001.tnanoseq.vcf.gz')
tNS_var = as.data.frame(info(tNS_var))

path="/nfs/cancer_ref01/nst_links/live/"
projectid = 3327 #foetal livers

samples = c('PD62336a_tds0001','PD60301a_tds0001','PD61846a_tds0001','PD61856a_tds0001','PD61851a_tds0001','PD61850a_tds0001','PD62331c_tds0001','PD62331a_tds0001')
all_var = data.frame()
for(sample in samples){
  data = readVcf(paste0(path,projectid,"/",sample,"/",sample,".tnanoseq.vcf.gz"))
  var = rowRanges(data)
  info = info(data)
  mcols(var) = cbind(mcols(var),info[match(names(var),rownames(info)),])
  var$varID = as.character(names(var))
  var = as.data.frame(mcols(var))
  var$sample = sample
  all_var = rbind(all_var,var)
}



## Prepare "mutations" table for dndscv: 5 columns in the following order: sampleID, chr, pos, ref, alt
all_var$chr = gsub(':.*$','',all_var$varID)
all_var$pos = gsub('^.*:|_.*$','',all_var$varID)
mutation = all_var[,c('sample','chr','pos','REF','ALT')]
#mutation$chr = paste0('chr',mutation$chr)
mutation = dndscv::dndscv(mutations = mutation,max_muts_per_gene_per_sample = 1e9,refdb = "hg19")
aggCnt.sub = dndscv(mutations = aggCnt.sub,cv=NULL,max_muts_per_gene_per_sample = 1e9,
                    refdb = '~/lustre_mt22/Aneuploidy/Results/16_GATA1_variants/RefCDS_human_GRCh38_GencodeV18_recommended.rda')



mutations = aggCnt.sub[['annotmuts']]
mutations$aachange2 = paste0('p.',mutations$aachange)
table(mutations$aachange2 %in% cosmic_gata1$AA.Mutation)
mutations$inCosmic = (mutations$aachange2 %in% cosmic_gata1$AA.Mutation)
mutations$varID = rownames(mutations)
mutations = cbind(mutations,aggCnt[match(mutations$varID,aggCnt$varID),!colnames(aggCnt) %in% colnames(mutations)])

