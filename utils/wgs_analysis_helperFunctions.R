## WGS analysis - helper functions ##


# Written by Taryn, with Mi's edits
# function to read Caveman file
# output df with Chr_Pos_Ref_Alt, NV, NR, VAF=NV/NR, Gene
# filter for PASS, ASMD, CLMP
caveman2 = function(sample, projectid, path="/nfs/cancer_ref01/nst_links/live/",filter=TRUE,ASMD_value=140,CLMP_value=0){
  require(VariantAnnotation)
  caveman_fp = paste0(path,projectid,"/",sample,"/",sample,".caveman_c.annot.vcf.gz")
  if(!file.exists(caveman_fp)){
    warning(sprintf('Caveman output for %s in project %s does not exist.',sample,projectid))
    return()
  }
  data = readVcf(caveman_fp)
  # Number of reads 
  reads = cbind((geno(data)$FAZ)[,2]+(geno(data)$RAZ)[,2],(geno(data)$FCZ)[,2]+(geno(data)$RCZ)[,2],
                (geno(data)$FGZ)[,2]+(geno(data)$RGZ)[,2],(geno(data)$FTZ)[,2]+(geno(data)$RTZ)[,2])
  Var = data.frame(Chr=as.character(seqnames(rowRanges(data))),
                   Pos=start(ranges(rowRanges(data))),
                   Ref=as.character(ref(data)))
  Alt_tmp = CharacterList(alt(data))
  Delete = which(sum(alt(data)%in%c("A","C","T","G"))!=1) # check all Variants have an alt allele
  Alt_tmp[Delete]="N"
  Var$Alt = as.character(unlist(Alt_tmp))
  Var$NR=rowSums(reads) #Total number of reads
  Var$NV=NA
  Var$Gene=info(data)$VD
  Var$Depth=info(data)$DP
  Var$Gene=gsub("\\|.*","",Var$Gene)
  Var$Impact=as.character(info(data)$VC)
  Var$Snp=info(data)$SNP
  Var$AAchange = sapply(strsplit(info(data)$VW,'\\|'),'[',5)
  Var$Type = sapply(strsplit(info(data)$VW,'\\|'),'[',6)
  Var$Filter = rowRanges(data)$FILTER
  Var$QUAL = rowRanges(data)$QUAL
  colnames(reads)=c("A","C","G","T")
  for (k in c("A","C","G","T")){
    Var$NV[Var$Alt==k] = reads[Var$Alt==k,k]  #number of reads with the alternate base
  }
  if(filter){
    select = rowRanges(data)$FILTER=="PASS"
    Var$Flag = info(data)$ASMD>=ASMD_value & info(data)$CLPM==CLMP_value
    Var = Var[select,]
  }
  Var$ID = paste(Var$Chr,Var$Pos,Var$Ref,Var$Alt,sep = "_")
  Var$VAF=Var$NV/Var$NR
  Var$sampleID = sample
  Var$varID = paste0(Var$Chr,':',Var$Pos,'_',Var$Ref,'/',Var$Alt)
  return(Var)
}






snvBinomTest=function(gender,data,genotype='diploid',pCut=0.05,useChr=NULL){
  # Function to filter out germline variants based on unmatched variant calls 
  # data is a data.frame where: rownames = snv_ID (chr1:pos_REF/ALT), required columns are refCnt, altCnt, DP (total depth)
  # Returns a logical vector, 
  # TRUE if mutation is likely to be germline.
  
  if(is.null(useChr)){
    if(all(grepl('^chr',data$snv_ID))){
      useChr = T
    }else if(all(grepl('^Chr',data$snv_ID))){
      useChr = T
      data$snv_ID = gsub('^Chr','chr',data$snv_ID)
    }else if(all(!grepl('^chr|^Chr',data$snv_ID))){
      useChr = F
    }
  }
  
  ## Check if there's cn region
  if('CN_tot2major' %in% colnames(data)){
    cnRegions = (data$CN_tot2major != '2:1')
  }else{
    cnRegions = F
  }
  
  
  if(grepl('T',genotype)){
    if(useChr){
      trisomyChr = gsub('T','chr',genotype)
    }else{
      trisomyChr = gsub('T','',genotype)
    }
    #print(unique(trisomyChr))
    trisomy = grepl(paste0('^',unique(trisomyChr)),data$snv_ID)
  }else{
    trisomy = F
  }
  
  if(gender %in% c("male",'Male','M')){
    XY_chromosomal = grepl("X|Y",data$snv_ID)
  }else{
    XY_chromosomal = F
  }
  
  
  autosomal = !(XY_chromosomal | cnRegions | trisomy)
  
  ## Test for autosomal 
  data$pval = NA
  data$pval[autosomal] = sapply(seq(1:nrow(data[autosomal,])),function(s){
    #print(s)
    binom.test(x=data$altCnt[autosomal][s],n=data$DP[autosomal][s], p=0.5,alt='less')$p.value
  })
  
  if(sum(XY_chromosomal) > 0){
    data$pval[XY_chromosomal] = sapply(seq(1:nrow(data[XY_chromosomal,])),function(s){
      binom.test(x=data$altCnt[XY_chromosomal][s],n=data$DP[XY_chromosomal][s], p=0.95,alt='less')$p.value
    })  
  }
  
  if(grepl('T',genotype)){
    data$pval[trisomy] = sapply(seq(1:nrow(data[trisomy,])),function(s){
      if(data$altCnt[trisomy][s] > 0.5* data$DP[trisomy][s]){
        binom.test(x=data$altCnt[trisomy][s],n=data$DP[trisomy][s], p=2/3,alt='two.sided')$p.value
      }else{
        binom.test(x=data$altCnt[trisomy][s],n=data$DP[trisomy][s], p=1/3,alt='two.sided')$p.value
      }
    })
  }
  
  # if(gender %in% c("female",'Female','F')){
  #   data$pval = NA
  #   if(genotype == 'diploid'){
  #     data$pval = sapply(seq(1:nrow(data)),function(s){
  #       binom.test(x=data$altCnt[s],n=data$DP[s], p=0.5,alt='less')$p.value
  #     })  
  #   }else if(grepl('T',genotype)){
  #     not_trisomy = !trisomy
  #     data$pval[not_trisomy] = sapply(seq(1:nrow(data[not_trisomy,])),function(s){
  #       binom.test(x=data$altCnt[not_trisomy][s],n=data$DP[not_trisomy][s], p=0.5,alt='less')$p.value
  #     })
  #     
  #     data$pval[trisomy] = sapply(seq(1:nrow(data[trisomy,])),function(s){
  #       if(data$altCnt[trisomy][s] > 0.5* data$DP[trisomy][s]){
  #         binom.test(x=data$altCnt[trisomy][s],n=data$DP[trisomy][s], p=2/3,alt='two.sided')$p.value  
  #       }else{
  #         binom.test(x=data$altCnt[trisomy][s],n=data$DP[trisomy][s], p=1/3,alt='two.sided')$p.value  
  #       }
  #       
  #     })
  #   }
  #   
  # }
  # # For male, split test in autosomal and XY chromosomal part
  # if(gender %in% c("male",'Male','M')){
  #   data$pval = NA
  #   
  #   if(genotype == 'diploid'){
  #     data$pval[autosomal] = sapply(seq(1:nrow(data[autosomal,])),function(s){
  #       binom.test(x=data$altCnt[autosomal][s],n=data$DP[autosomal][s], p=0.5,alt='less')$p.value
  #     })  
  #   }else if(grepl('T',genotype)){
  #     not_trisomy = !(trisomy|XY_chromosomal)
  #     data$pval[not_trisomy] = sapply(seq(1:nrow(data[not_trisomy,])),function(s){
  #       binom.test(x=data$altCnt[not_trisomy][s],n=data$DP[not_trisomy][s], p=0.5,alt='less')$p.value
  #     })
  #     
  #     data$pval[trisomy] = sapply(seq(1:nrow(data[trisomy,])),function(s){
  #       if(data$altCnt[trisomy][s] > 0.5* data$DP[trisomy][s]){
  #         binom.test(x=data$altCnt[trisomy][s],n=data$DP[trisomy][s], p=2/3,alt='two.sided')$p.value  
  #       }else{
  #         binom.test(x=data$altCnt[trisomy][s],n=data$DP[trisomy][s], p=1/3,alt='two.sided')$p.value  
  #       }
  #       
  #     })
  #   }
  # 
  #   
  #   if(sum(XY_chromosomal) > 0){
  #     data$pval[XY_chromosomal] = sapply(seq(1:nrow(data[XY_chromosomal,])),function(s){
  #       binom.test(x=data$altCnt[XY_chromosomal][s],n=data$DP[XY_chromosomal][s], p=0.95,alt='less')$p.value
  #     })  
  #   }
  # }
  # 
  
  data$qval = p.adjust(data$pval,method="BH")
  data$isGermline = NA
  data$isGermline = data$qval > pCut
  # data$isGermline[autosomal] = (data$qval[autosomal]>pCut)  
  # data$isGermline[!autosomal] = (data$qval[!autosomal]<pCut)  
    
  # if(genotype == 'diploid'){
  #   data$isGermline = data$qval>pCut  
  # }else{
  #   data$isGermline[!trisomy] = data$qval[!trisomy]>pCut  
  #   data$isGermline[trisomy] = data$qval[trisomy]>pCut*  
  # }
  
  return(data)
}

plot_binomMixModel = function(data,res,sampleID=NULL,mode='Full'){
  # data is a dataframe used as input for binom_mix function. Data must have NV and NR columns
  # res is the output of binom_mix function
  NV = data$NV
  NR = data$NR
  p=hist(NV/NR,breaks=20,xlim=c(0,1),col='gray',freq=F,xlab="Variant Allele Frequency",
         main=ifelse(is.null(sampleID),paste0("n=",length(NV)),paste0(sampleID,' (n=',length(NV),')')))
  cols=c("red","blue","green","magenta","cyan")
  
  y_coord=max(p$density)-0.5
  y_intv=y_coord/5
  text(y=y_coord,x=0.9,label='Data')
  segments(lwd=2,lty='dashed',col='black',y0=y_coord+0.25,x0=0.85,x1=0.95)
  
  for (i in 1:res$n){
    depth=rpois(n=5000,lambda=median(NR))
    sim_NV=unlist(lapply(depth,rbinom,n=1,prob=res$p[i]))
    sim_VAF=sim_NV/depth
    if (mode=="Truncated") sim_VAF=sim_VAF[sim_NV>3]
    dens=density(sim_VAF)
    lines(x=dens$x,y=res$prop[i]*dens$y,lwd=2,lty='dashed',col=cols[i])
    y_coord=y_coord-y_intv/2
    text(y=y_coord,x=0.9,label=paste0("p1: ",round(res$p[i],digits=2),', ',100*round(res$prop[i],digits=2),'%'))
    segments(lwd=2,lty='dashed',col=cols[i],y0=y_coord+y_intv/4,x0=0.85,x1=0.95)
  }
  
}





# function to read Pindel file
# filter for PASS
pindel = function(sample, projectid, path="/nfs/cancer_ref01/nst_links/live/",filter_lowQual=T){
  pindel_fp = paste0(path,projectid,"/",sample,"/",sample,".pindel.annot.vcf.gz")
  if(!file.exists(pindel_fp)){
    warning(sprintf('Pindel output not found for sample %s in project %d',sample,projectid))
    return()
  }
  data = VariantAnnotation::readVcf(pindel_fp)
  if(filter_lowQual){
    data = data[mcols(data)$FILTER == 'PASS' & mcols(data)[['QUAL']] > 300]
  }
  
  Var = data.frame(Chr=as.character(seqnames(rowRanges(data))),
                   Pos=start(ranges(rowRanges(data))),
                   Ref=as.character(ref(data)),
                   Alt=as.character(unlist(alt(data))),
                   Filter = mcols(data)[['FILTER']],
                   Qual = mcols(data)[['QUAL']],
                   Gene = gsub("\\|.*","",info(data)$VD),
                   AAchange = sapply(strsplit(info(data)$VW,'\\|'),'[',5),
                   Type = sapply(strsplit(info(data)$VW,'\\|'),'[',6))
  Var$Impact = sapply(info(data)$VC,function(i){paste(i,collapse = ':')})
  Var$ID = paste(Var$Chr,Var$Pos,Var$Ref,Var$Alt,sep = "_")
  Var$varID = paste0(Var$Chr,':',Var$Pos,'_',Var$Ref,'/',Var$Alt)
  Var$End=Var$Pos+nchar(Var$Ref)-1
  Var$sampleID = sample

  return(Var)
}



pindel_all = function(sample, projectid, path="/nfs/cancer_ref01/nst_links/live/"){
  data2 = fread(paste0(path,projectid,"/",sample,"/",sample,".pindel.annot.vcf.gz"), skip="CHROM", select=c(1,2,4,5,7,6,8))
  #select = data$FILTER=="PASS" & data$QUAL >=300
  #data = data[select,]
  data$ID = paste(data$`#CHROM`,data$POS,data$REF,data$ALT,sep = "_")
  data$Gene = gsub(".*VD=","",data$INFO)
  data$Gene=gsub("\\|.*","",data$Gene)
  data$Impact = gsub(".*VC=","",data$INFO)
  data$Impact=gsub("\\;.*","",data$Impact)
  data$Snp=data[,grepl("SNP", data$INFO)]
  data$AAchange = gsub(".*p\\.","",data$INFO)
  data$AAchange=gsub("\\|.*","",data$AAchange)
  data$Impact= as.character(data$Impact)
  return(data)
}




depthFilter = function(data,gender,genotype,depth_cutoff = NULL,useChr=NULL,upperCutOff=F){
  
  if(is.null(useChr)){
    if(all(grepl('^chr',data$snv_ID))){
      useChr = T
    }else if(all(grepl('^Chr',data$snv_ID))){
      useChr = T
      data$snv_ID = gsub('^Chr','chr',data$snv_ID)
    }else if(all(!grepl('^chr|^Chr',data$snv_ID))){
      useChr = F
    }
  }
  
  
  
  # If female, only need to consider the Y-chr, as female would also have 2 copies of chrX
  # If male, having 1 copy of ChrX and 1 copy of ChrY --> need to exclude both from depth distribution
  if (gender %in% c("male",'Male','M','XY')) {
    XY_chr = c('X','Y')
  }else{
    XY_chr = c('Y')
  }
  
  if(!genotype %in% c('2n','diploid') & grepl('^T|trisomy',genotype)){
    trisomyChr = gsub('^T|Trisomy|trisomy','',genotype)
  }else{
    trisomyChr = c()
  }
  
  chr_toExclude = c(XY_chr,trisomyChr)
  if(useChr){
    chr_toExclude = paste0('chr',chr_toExclude)
  }
  
  ## Set cutoff threshold without chr to exclude (i.e. regions with CN changes)
  df = data[!grepl(paste(chr_toExclude,collapse = '|'),data$snv_ID),]
  
  if(is.null(depth_cutoff)){
    depth_cutoff = c(mean(df$sample_DEP) - 3*sd(df$sample_DEP),
                     mean(df$sample_DEP) + 3*sd(df$sample_DEP))
  }
  
  
  # Plot histogram of depth at each variant in both samples
  p = ggplot(data,aes(sample_DEP))+
    geom_histogram(binwidth = 1)+
    geom_vline(xintercept = c(depth_cutoff),lty=2,col='red')+
    theme_classic() + ggtitle('cgpVAF DEPTH distribution from all samples',
                              subtitle = paste0(n_distinct(data$snv_ID),' variants'))
  print(p)
  
  # Plot boxplot/hist of depth distribution by chromosome
  if(useChr){
    l = paste0('chr',c(1:22,'X','Y'))
  }else{
    l = c(1:22,'X','Y')
  }
  
  
  data$Chrom = factor(data$Chrom,levels = l)
  p2 = ggplot(data,aes(Chrom,sample_DEP))+
    geom_boxplot(outlier.colour = 'white')+
    #geom_histogram(binwidth = 1)+
    geom_hline(yintercept = c(depth_cutoff),lty=2,col='red')+
    theme_classic() + ggtitle('cgpVAF DEPTH distribution from all samples',
                              subtitle = paste0(n_distinct(data$snv_ID),' variants')) + xlab('')
  
  print(p2)
  
  data$depthFilter = 'NA'
  data$depthFilter[data$snv_ID %in% df$snv_ID & data$sample_DEP < min(depth_cutoff)] = 'tooLow_DEP'
  if(upperCutOff == F){
    depth_cutoff = min(depth_cutoff)
  }else if(upperCutOff == T & length(depth_cutoff) > 1){
    data$depthFilter[data$snv_ID %in% df$snv_ID & data$sample_DEP > max(depth_cutoff)] = 'tooHigh_DEP'
  }
  data$depthFilter[data$snv_ID %in% df$snv_ID & data$depthFilter == 'NA'] = 'PASS'
  
  
  
  return(data)
}


df2granges = function(data, chr, start, end,id,sampleID=NULL){
  require(GenomicRanges)
  data_gr = GRanges(data[[chr]],IRanges(data[[start]],data[[end]]))
  mcols(data_gr)[[id]] = data[[id]]
  if(!is.null(sampleID)){
    mcols(data_gr)[[sampleID]] = data[[sampleID]]
  }
  
  mcols(data_gr) = cbind(mcols(data_gr),data[match(paste0(mcols(data_gr)[[id]],':',mcols(data_gr)[[sampleID]]),
                                                   paste0(data[[id]],':',data[[sampleID]])),])
  if(length(data_gr) == nrow(data)){
    return(data_gr)
  }else{
    stop('something has gone wrong...Please check!')
  }
  
}






##--------------------------------------------------------------------------------##
##      HELPER FUNCTIONS for running Tim's version of Shearwater filter         ####
##--------------------------------------------------------------------------------##

cgpVaf_onNormPanel = function(outDir,somaticVar,normPanel_fp = "/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt",
                              script = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/scripts/activeScritps/WGS_analysis/x10_MLDS_cgpVaf_forShearwater.sh'){
  
  ##--------------------------------------------------------------------##
  ##    STEP1: cgpVAF of variants of interst in panel of normal samples ##
  ##--------------------------------------------------------------------##
  print('STEP1: cgpVAF of variants of interst in panel of normal samples')
  
  ## Create cgpVaf output folder
  outDir_list = file.path(outDir,c('1_cgpVaf/cgpVaf_input','1_cgpVaf/cgpVaf_output','2_allelecount'))
  for(d in outDir_list){
    if(!dir.exists(d)){
      dir.create(d,recursive = T)
    }  
  }
  
  
  ## Run cgpVaf on the normal samples
  ## 1. Create bed files
  # requried format:
  #   chr pos ref_allele alt_allele
  #   filename ends with .bed, delim='\t'
  bed = data.frame(snv_ID = unique(somaticVar$snv_ID))
  bed = merge(bed,somaticVar[,c('Chrom','Pos','Ref','Alt','snv_ID')],by='snv_ID',all.x=T)
  bed = bed[!duplicated(bed),]
  dim(bed)
  n_distinct(somaticVar$snv_ID)
  bed = bed[order(bed$Chrom),]
  write_delim(bed[,colnames(bed) != 'snv_ID'],file.path(outDir,'1_potential_somaticVar_forShearwater.bed'),col_names = F,delim = '\t')
  
  
  ## Import list of normal samples
  inDir = file.path(outDir,'1_cgpVaf/cgpVaf_input')
  
  normal_samples=read.table(normPanel_fp, header=T) 
  for(i in 1:nrow(normal_samples)){
    sampleID = normal_samples$Sample_ID[i]
    projectID = normal_samples$Project_ID[i]
    path = file.path('/nfs/cancer_ref01/nst_links/live',projectID,sampleID,paste0(sampleID,'.sample.dupmarked.bam'))
    if(!file.exists(path)){
      stop(sprintf('[%s] File does not exist: %s',sampleID,path))
    }
    
    system(sprintf('cp -s %s %s',path,inDir))
    system(sprintf('cp -s %s %s',paste0(path,'.bai'),inDir))
    
  }
  
  system(sprintf('cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam %s/PDv38is_wgs.sample.dupmarked.bam',inDir))
  system(sprintf('cp -s /nfs/cancer_ref01/nst_links/live/2480/PDv38is_wgs/PDv38is_wgs.sample.dupmarked.bam.bai %s/PDv38is_wgs.sample.dupmarked.bam.bai',inDir))
  
  
  cgpVaf_outDir = file.path(outDir,'1_cgpVaf/cgpVaf_output/')
  bedFile = file.path(outDir,'1_potential_somaticVar_forShearwater.bed')
  sampleFile = normPanel_fp#"/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt"
  #sampleFile = "/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/L038/test_sample.txt"
  print(sprintf('bash %s %s %s %s %s',script,inDir,cgpVaf_outDir,bedFile,sampleFile))
  
}






alleleCount_onNormPanel = function(outDir,allele_count,
                                   normPanel_fp = "/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt"){
  
  ##-------------------------------------------------------------------------------##
  ##    STEP2: alleleCount at each variant position across panel of normal samples ##
  ##-------------------------------------------------------------------------------##
  
  print('STEP2: alleleCount at each variant position across panel of normal samples')
  
  allelecount_dir = file.path(outDir,'2_allelecount')
  if(!dir.exists(allelecount_dir)){
    stop(sprintf('Please check: allelecount_dir does not exist: %s',allelecount_dir))
  }  
  
  
  # This will be used for the shearwater analysis
  normal_samples=read.table(normPanel_fp, header=T) 
  f = list.files(file.path(outDir,'1_cgpVaf/cgpVaf_output'),pattern = 'PDv38is_wgs_.*_snp_vaf.tsv$',full.names = T)
  n = system(sprintf("cat %s | grep '^#' | wc -l",f),intern = T) %>% as.integer()
  normSamples_cgpVaf = read.table(f, header = TRUE, sep = '\t', stringsAsFactors = F, comment.char = "",skip = (n-1))
  colnames(normSamples_cgpVaf)[which(colnames(normSamples_cgpVaf)=="VariantID")]="ID"
  normSamples_cgpVaf=normSamples_cgpVaf[!duplicated(normSamples_cgpVaf[,c("Chrom", "Pos")]),]
  row.names(normSamples_cgpVaf)=paste0(sub('.*chr', '', normSamples_cgpVaf$Chrom), "_", normSamples_cgpVaf$Pos)
  #Reorder normSamples_cgpVaf to have the same order of rows as allele count
  table(row.names(allele_count) %in% row.names(normSamples_cgpVaf))
  normSamples_cgpVaf=normSamples_cgpVaf[match(row.names(allele_count), row.names(normSamples_cgpVaf)), ]
  
  for (i in 1:nrow(normal_samples)) {
    case=normal_samples$Case[i]
    print(case)
    num_samples=normal_samples$Sample_ID[startsWith(normal_samples$Sample_ID, case)]
    
    for (j in 1:length(num_samples)) {
      sample=num_samples[j]
      print(sample)
      mut_temp=normSamples_cgpVaf[,grep(paste0(sample,"_"), colnames(normSamples_cgpVaf))]
      #reset the allele counts to NA before starting a new sample. Not strictly necessary but doing it for my sanity/piece of mind
      allele_count$Count_A=NA
      allele_count$Count_C=NA
      allele_count$Count_G=NA
      allele_count$Count_T=NA
      
      allele_count$Count_A = mut_temp[,grepl('FAZ',colnames(mut_temp))] + mut_temp[,grepl('RAZ',colnames(mut_temp))]
      allele_count$Count_C = mut_temp[,grepl('FCZ',colnames(mut_temp))] + mut_temp[,grepl('RCZ',colnames(mut_temp))]
      allele_count$Count_G = mut_temp[,grepl('FGZ',colnames(mut_temp))] + mut_temp[,grepl('RGZ',colnames(mut_temp))]
      allele_count$Count_T = mut_temp[,grepl('FTZ',colnames(mut_temp))] + mut_temp[,grepl('RTZ',colnames(mut_temp))]
      
      write.table(allele_count, file = file.path(allelecount_dir,paste0(sample, "_allelecounts.txt")),
                  sep = "\t", row.names = F, quote = F, col.names = T)
    }
  }
  
  print('STEP2: Completed')
}


runShearWaterLike = function(samples_ID,outDir,normPanel_fp = "/lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/mutations/normal_samples.txt",
                             FDR=0.001,patient=NULL){
  
  ##----------------------------------------------------------------------##
  ##    STEP3: Run Tim's shearwater-like filter to remove false mutations ##
  ##----------------------------------------------------------------------##
  
  print('STEP3: Run Tims shearwater-like filter to remove false mutations')
  
  ##---- Inputs
  # Vector of all sample names
  samples_ID=samples_ID
  # Vector of normal panel (only blood from other project [3110])
  normal_samples=read.table(normPanel_fp, header=T) 
  
  # "Bed" file of all mutations to be considered (across all patients)
  # Format: Chr Pos Ref Alt
  bedFile = file.path(outDir,'1_potential_somaticVar_forShearwater.bed')
  variants = read.delim(bedFile,header = F)
  variants=variants[!duplicated(variants[,1:2]),]
  coords=paste(variants$V1,variants$V2,sep="_")
  
  
  # Read in data from AlleleCounter/cgpVAF
  allelecount_dir = file.path(outDir,'2_allelecount')
  
  all_counts = array(0,dim=c(length(samples_ID),length(coords),4),
                     dimnames=list(samples_ID,coords,c("A","C","G","T")))
  print(length(samples_ID))
  for (k in 1:length(samples_ID)){
    #Read in allele counts per sample
    if(file.exists(file.path(allelecount_dir, paste0(samples_ID[k],"_allelecounts.txt")))){
      print(samples_ID[k])
      data=read.table(file.path(allelecount_dir, paste0(samples_ID[k],"_allelecounts.txt")),comment.char = '',header=T)
      rownames(data)=paste(data$Chrom,data$Pos,sep="_")
      all_counts[k,,]=as.matrix(data[coords,3:6])
    }
  }
  
  # Read in data from AlleleCounter/cgpVAF for all samples in the normal panel
  norm_all_counts = array(0,dim=c(nrow(normal_samples),length(coords),4),
                          dimnames=list(normal_samples$Sample_ID,coords,c("A","C","G","T")))
  
  for (k in 1:nrow(normal_samples)){
    #Read in allele counts per sample
    if(file.exists(file.path(allelecount_dir, paste0(normal_samples$Sample_ID[k],"_allelecounts.txt")))){
      print(normal_samples$Sample_ID[k])
      data=read.table(file.path(allelecount_dir, paste0(normal_samples$Sample_ID[k],"_allelecounts.txt")),comment.char = '',header=T)
      rownames(data)=paste(data$Chrom,data$Pos,sep="_")
      norm_all_counts[k,,]=as.matrix(data[coords,3:6])
    }
  }
  
  
  
  
  
  
  #------- Run shearwater filter
  #setwd("~/lustre_mt22/Aneuploidy/MLDS_scRNAseq/x10.1_WGS_variantAnalysis/shearwater/")
  source('~/lustre_mt22/generalScripts/shearwater/shearwater_flt_2019_WT_mt22.R')
  colnames(variants) = c('Chr','Pos','Ref','Alt')
  if(!is.null(patient)){
    outFile = file.path(outDir,paste0(patient,"_shearwater_pval_mat.txt"))
  }else{
    outFile = file.path(outDir,"shearwater_pval_mat.txt")
  }
  pval_mat = shearwater_probability(patient, save=outFile,allVar=variants,
                                    all_counts=all_counts,norm_all_counts=norm_all_counts,case_samples=samples_ID, rho=10^-3)
  
  
  qval_mat=apply(pval_mat,2,function(x) p.adjust(x,method="BH",n = length(pval_mat)))
  # Note that if MTR=DEP, the p-value will be so low that it is listed as NA(n). Change these values to 0
  qval_mat[which(qval_mat=="NaN")]=0
  row.names(qval_mat)=sub('.*chr', '', row.names(qval_mat))
  
  # Add the qval info to the mutations object (make sure rows are in the same order first)
  mutations = variants
  row.names(mutations)=paste0(sub('.*chr', '', mutations$Chr), "_", mutations$Pos)#, "_", mutations$Ref,"_", mutations$Alt)
  table(row.names(mutations) %in% row.names(qval_mat))
  reorder_idx=match(row.names(mutations), row.names(qval_mat))
  qval_mat=qval_mat[reorder_idx, ]
  mutations=cbind(mutations, qval_mat)
  
  return(mutations)
}









import_multiple_caveman=function(caveman_samples,sampleID_name='PDID',ASMD_CLPM_filter=T,filter=T){
  out_list = list()
  for(s in unique(caveman_samples[[sampleID_name]])){
    caveman = caveman2(sample=s,projectid = unique(caveman_samples$projectid[caveman_samples[[sampleID_name]] == s]),filter = filter)
    caveman$PDID = s
    caveman = cbind(caveman,caveman_samples[match(caveman$PDID,caveman_samples[[sampleID_name]]),!colnames(caveman_samples) %in% colnames(caveman)])
    if(ASMD_CLPM_filter){
      caveman = caveman[caveman$Flag == T,]  
    }
    caveman$snv_ID = paste0(caveman$Chr,':',caveman$Pos,'_',caveman$Ref,'/',caveman$Alt)
    
    out_list[[s]] = caveman
  }
  #------ Combine all variants across all samples
  out = do.call(rbind,out_list)
}




import_cgpVAF_output = function(output_wFilter_file,caveman_samples,sampleID_name='sampleID',output_noFilter_file=NULL,
                                columns = c('snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect')){
  
  if(length(columns) == 1 & columns[1] == 'all'){
    grep_pattern = 'F*Z|R*Z|WTR|MTR|DEP|VAF' 
    columns = c('snv_ID','Chrom','Pos','Ref','Alt','Gene','Type','Effect')
  }else{
    grep_pattern = 'WTR|MTR|DEP|VAF' 
  }
  ## all samples are included in this cgpVAF output
  
  n = system(sprintf("cat %s | grep '^#' | wc -l",output_wFilter_file),intern = T) %>% as.integer()
  
  ##---- With quality parameters
  cgpVaf_withReadFilters_output = read.delim(output_wFilter_file,sep = '\t',skip = (n-1))
  cgpVaf_withReadFilters_output$snv_ID = paste0(cgpVaf_withReadFilters_output$Chrom,':',cgpVaf_withReadFilters_output$Pos,'_',cgpVaf_withReadFilters_output$Ref,'/',cgpVaf_withReadFilters_output$Alt)
  samples = gsub('_MTR$','',colnames(cgpVaf_withReadFilters_output)[grepl('MTR',colnames(cgpVaf_withReadFilters_output))])
  samples = samples[samples != 'PDv38is_wgs']
  ## Convert to a long table
  cgpVaf_withReadFilters_output = do.call(rbind,lapply(samples,function(s){
    #print(s)
    tmp = cgpVaf_withReadFilters_output[,unique(c('snv_ID',columns,colnames(cgpVaf_withReadFilters_output)[grepl(grep_pattern,colnames(cgpVaf_withReadFilters_output)) & grepl(s,colnames(cgpVaf_withReadFilters_output))]))]
    colnames(tmp) = gsub(s,'wFilter_sample',colnames(tmp))
    tmp$PDID = s
    tmp[[sampleID_name]] = caveman_samples[[sampleID_name]][match(tmp$PDID,caveman_samples$PDID)]
    return(tmp)
  }))
  
  if(!is.null(output_noFilter_file)){
    ##---- WithOUT quality parameters
    cgpVaf_noReadFilters_output = read.delim(output_noFilter_file,sep = '\t',skip = (n-1))
    cgpVaf_noReadFilters_output$snv_ID = paste0(cgpVaf_noReadFilters_output$Chrom,':',cgpVaf_noReadFilters_output$Pos,'_',cgpVaf_noReadFilters_output$Ref,'/',cgpVaf_noReadFilters_output$Alt)
    
    ## Convert to a long table
    cgpVaf_noReadFilters_output = do.call(rbind,lapply(samples,function(s){
      #print(s)
      tmp = cgpVaf_noReadFilters_output[,unique(c('snv_ID',columns,colnames(cgpVaf_noReadFilters_output)[grepl(grep_pattern,colnames(cgpVaf_noReadFilters_output)) & grepl(s,colnames(cgpVaf_noReadFilters_output))]))]
      colnames(tmp) = gsub(s,'noFilter_sample',colnames(tmp))
      tmp$PDID = s
      tmp[[sampleID_name]] = caveman_samples[[sampleID_name]][match(tmp$PDID,caveman_samples$PDID)]
      return(tmp)
    }))
    
    ## Combine with and without read quality filter
    out = cbind(cgpVaf_withReadFilters_output[,c(columns,'PDID',sampleID_name,colnames(cgpVaf_withReadFilters_output)[!colnames(cgpVaf_withReadFilters_output) %in% c(columns,'PDID',sampleID_name)])],
                cgpVaf_noReadFilters_output[match(paste0(cgpVaf_withReadFilters_output$snv_ID,'_',cgpVaf_withReadFilters_output$PDID),
                                                  paste0(cgpVaf_noReadFilters_output$snv_ID,'_',cgpVaf_noReadFilters_output$PDID)),colnames(cgpVaf_noReadFilters_output)[!colnames(cgpVaf_noReadFilters_output) %in% c(columns,'PDID',sampleID_name)]]
    )
  }else{
    out = cgpVaf_withReadFilters_output[,c(columns,'PDID',sampleID_name,colnames(cgpVaf_withReadFilters_output)[!colnames(cgpVaf_withReadFilters_output) %in% c(columns,'PDID',sampleID_name)])]
  }
  
  
  
  
  return(out)
}






import_multiple_pindel = function(caveman_samples,sampleID_name='PDID',intv=10,filter_lowQual=F){
  library(data.table)
  out_list = list()
  for(s in unique(caveman_samples[[sampleID_name]])){
    indels = pindel(sample = s,projectid=unique(caveman_samples$projectid[caveman_samples[[sampleID_name]] == s]),filter_lowQual=filter_lowQual)
    dim(indels)
    n_distinct(indels$varID)
    indels$End=indels$Pos+nchar(indels$Ref)-1
    indels$PDID = s
    indels = cbind(indels,caveman_samples[match(indels$PDID,caveman_samples[[sampleID_name]]),!colnames(caveman_samples) %in% colnames(indels)])
    
    out_list[[s]] = indels
  }
  
  #------ Combine all variants across all samples
  indels = do.call(rbind,out_list)
  
  indels_gr = GRanges(indels$Chr,IRanges(indels$Pos - intv,indels$End + intv),
                      varID = indels$varID,sampleID=indels$sampleID)
  mcols(indels_gr) = cbind(mcols(indels_gr),indels[match(paste0(indels_gr$varID,':',indels_gr$sampleID),
                                                         paste0(indels$varID,':',indels$sampleID)),])
  length(indels_gr) == nrow(indels)
  
  indels=indels_gr
  
  return(indels)
}









generate_trinucleotide_matrix = function(mutations,outDir,patientID,tool='mSigAct',genomeVersion = c('hg38','hg37')){
  #Annotate the trinuc ref
  #NOTE! Need to change the genomeFile if you have hg37 data!
  
  if(genomeVersion == 'hg38'){
    genomeFile="/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa"  
  }else if(genomeVersion == 'hg37'){
    genomeFile="/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/genome.fa"  
  }
  
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$Chrom, IRanges(as.numeric(mutations$Pos)-1, 
                                                                                       as.numeric(mutations$Pos)+1)))) #Note naming of the columns
  # Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$Ref,mutations$Alt,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$Ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$Ref[j]],ntcomp[mutations$Alt[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # Save trinucleotide matrix for the patient by looping through each sample 
  samples_pat = gsub('_true_mut','',colnames(mutations)[grepl('^PD\\d+.*_true_mut',colnames(mutations))])
  freqs_all=data.frame(matrix(NA, nrow = 96, ncol=length(samples_pat)))
  for(i in 1:length(samples_pat)) {
    sample=samples_pat[i]
    print(sample)
    idx=which(colnames(mutations)==paste0(sample, "_true_mut")) #I have a column called sample_true_mut for each sample to indicate whether they have passed all filters. Your identifier might be called something else!
    mut_temp=mutations[mutations[,idx]==1,] #Only pick the variants that are true_mut for this sample
    
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    if(tool!='mSigAct'){
      freqs = table(paste(mut_temp$sub,paste(substr(mut_temp$trinuc_ref_py,1,1),substr(mut_temp$trinuc_ref_py,3,3),sep="-"),sep=","))
      ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
      full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")  
    }else{
      freqs = table(paste0(substr(mut_temp$trinuc_ref_py,1,1),'[',mut_temp$sub,']',substr(mut_temp$trinuc_ref_py,3,3)))
      full_vec = paste0(rep(rep(c("A","C","G","T"),each=4),times=6),'[',rep(sub_vec,each=16),']',rep(rep(c("A","C","G","T"),times=4),times=6))
    }
    
    freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
    df=as.data.frame(freqs_full)
    freqs_all[,i]=df$Freq
    colnames(freqs_all)[i]=sample
  }
  
  
  
  row.names(freqs_all)=full_vec
  write.table(freqs_all, 
              file = file.path(outDir,paste0(patientID,'_tri_Matrix.txt')),
              sep = "\t", row.names = T, quote = F, col.names = T)
  
}









plot_spectrum = function(mutations,save,add_to_title="",style=1,genomeVersion=c('hg38','hg37')){
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  
  if(genomeVersion == 'hg38'){
    genomeFile="/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa"  
  }else if(genomeVersion == 'hg37'){
    genomeFile="/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/genome.fa"  
  }
  #genomeFile="/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa"
  #mutations=as.data.frame(bed)
  
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$Chrom, IRanges(as.numeric(mutations$Pos)-1, 
                                                                                       as.numeric(mutations$Pos)+1))))
  
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$Ref,mutations$Alt,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$Ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$Ref[j]],ntcomp[mutations$Alt[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
  
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  if(is.null(save)){dev.new(width=12,height=4)}
  #if(is.null(save)){dev.new(width=16,height=6)} #for tiff
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full; maxy = max(y)
  t = ifelse(add_to_title =='',paste0(add_to_title," - Total number of mutations: ",nrow(mutations)),paste0("Total number of mutations: ",nrow(mutations)))
  if(style==1){
    
    h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.lab=1.5, cex.axis = 1.3, cex.names=0.9, names.arg=xstr, ylab="Number of mutations",main=t)
    #h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number of mutations", main=paste0("Total number of mutations: ",nrow(mutations)))
    for (j in 1:length(sub_vec)) {
      xpos = h[c((j-1)*16+1,j*16)]
      rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
      text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
    }  
  }else if(style==2){
    df = as.data.frame(freqs_full)
    colnames(df) = c('tri_nt_context','nMutation')
    df$mutationGroup = gsub(',.*$','',df$tri_nt_context)
    df$mutationGroup = factor(df$mutationGroup,sub_vec)
    df$tri_nt_context = xstr
    
    colvec = unique(colvec)
    names(colvec) = sub_vec
    
    library(ggh4x)
    strip <- strip_themed(background_x = elem_list_rect(fill = unique(colvec)))
    
    p = ggplot(df,aes(tri_nt_context,nMutation,fill=mutationGroup))+
      geom_col(width = 0.4)+
      #facet_wrap(vars(mutationGroup),nrow=1,scales = 'free_x')+
      facet_wrap2(vars(mutationGroup),nrow=1,scales = 'free_x',strip = strip)+
      scale_fill_manual(values = colvec)+
      theme_classic()+
      theme(axis.text.x = element_text(angle=90,size=6,vjust = 0.5,hjust = 1,colour = 'black'),
            axis.text.y = element_text(size=11,colour = 'black'),
            axis.title.y = element_text(size=12.5,colour = 'black'),
            panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
            strip.background = element_blank(),
            legend.position = 'none',
            panel.spacing.x = unit(0.1,'cm'),
            #strip.background.x = element_rect(fill=colvec),
            strip.text = element_text(colour='white',size=12),
            axis.ticks = element_line(colour = 'black')) + xlab('') + ylab('Number of mutations')+ggtitle(t)
    print(p)
  }
   
  
  
  if(!is.null(save)){ 
    #dev.copy(tiff,save,width=16,height=6, units="in", compression="lzw", res=300); dev.off()} #for tiff image
    dev.copy(pdf,save,width=7,height=3); dev.off()}
}
