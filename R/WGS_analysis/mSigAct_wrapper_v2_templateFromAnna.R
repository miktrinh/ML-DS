# Wrapper for assigning COSMIC signatures to samples based on a mutational matrix as input
# Mutational matrix and cosmic ref signatures imported below
# Run as:
# module load R/4.3.1-jupyter-v0
# export R_HOME=$(R RHOME)
# prepend_path R_LIBS_USER $HOME/R-modules-4.3.1
# cd /lustre/scratch126/casm/team274sb/aw35/hepatoblastoma/WGS/signatures/mSigAct/
# bsub -o $PWD/log.%J -e $PWD/err.%J -G team274-grp -q normal -R'select[mem>10000] rusage[mem=10000]' -M10000 -n 20 Rscript /lustre/scratch126/casm/team274sb/aw35/misc_scripts/mSigAct_wrapper_v2.R GRCh37 01_Input/ 02_Output/ SBS1_SBS5_SBS18_SBS31_SBS35_SBS40a_SBS40b_SBS40c FALSE FALSE

# Arguments:
# hg_version : GRCh37 or GRCh38
# input_prefix: folder containing your trinculeotide matrix and SBS COSMIC REF [use 01_Input/]
# output_prefix: output folder [use 02_Output/]
# sigs: SBS1_SBS5_SBS31_SBS35 etc (must use underscore between sigs!)
# use_sig_presence_to_assign_sig: TRUE/FALSE (I usually use FALSE as this parameter generally doesn't make a difference)
# drop_low_mut_samples: TRUE/FALSE (if TRUE, will remove any samples with less than 100 mutations)

#NOTE. The trinucleotide matrix must be named something including Matrix. See line 52

#-------------------------------------------------------------------------------------

# Get arguments

#-------------------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 6) {
  cat("Error: improper arguments. See header of the script for details.\n")
  exit();
}

hg_version = args[1]
input_prefix=args[2]
output_prefix = args[3]
sigs=args[4]
use_sig_presence_to_assign_sig=args[5]
drop_low_mut_samples=args[6]

cat("Arguments loaded \n")

#-------------------------------------------------------------------------------------

library(mSigAct)

cat("Library loaded successfully \n")

#-------------------------------------------------------------------------------------

#work_dir=getwd()

#Load the Mutational matrix of your samples
mut_mat_file=list.files(input_prefix)[grep("Matrix",list.files(input_prefix))]

mut_mat=as.matrix(read.table(paste0(input_prefix,mut_mat_file), header = T, row.names = 1, sep = '\t'))

cat("Mutational matrix:", mut_mat_file, " has been loaded \n")


#Load the COSMIC SBS REF for your genome version
cosmic_file=list.files(input_prefix)[grep("COSMIC",list.files(input_prefix))]
cosmic_hg=cosmic_file[grep(hg_version, cosmic_file)] #Grep only of correct hg version
cosmic_hg_newest=sort(cosmic_hg, decreasing = TRUE)[1]
cosmic_ref=as.matrix(read.table(paste0(input_prefix,cosmic_hg_newest), header=TRUE, row.names = 1))


# Put the mutation types into the same order in the cosmic_ref and the mutational matrix
reorder_idx=match(row.names(mut_mat), row.names(cosmic_ref))
cosmic_ref=cosmic_ref[reorder_idx, ]

cat("COSMIC ref ", cosmic_hg_newest, " has been loaded  \n")


# Need to change the row names (mutation names) to the format that mSigAct accepts (ACAA, ACCA etc)

row.names(mut_mat)=paste0(substr(row.names(mut_mat), 1,1),
                     substr(row.names(mut_mat), 3,3),
                     substr(row.names(mut_mat), 7,7),
                     substr(row.names(mut_mat), 5,5))

row.names(cosmic_ref)=paste0(substr(row.names(cosmic_ref), 1,1),
                          substr(row.names(cosmic_ref), 3,3),
                          substr(row.names(cosmic_ref), 7,7),
                          substr(row.names(cosmic_ref), 5,5))


#Pick out sigs to use

sigs_to_use=unlist(strsplit(sigs, "_"))

cat("Signatures that mSigAct will try to assign to the samples:", sigs_to_use, "\n")

cat("Use signature presence test for assigning the signatures:", use_sig_presence_to_assign_sig, "\n")

cat("Drop samples with <100 mutations:", drop_low_mut_samples, "\n")


#Run mSigAct signature assignment

sparse_out=SparseAssignActivity(spectra=mut_mat, 
                                sigs=cosmic_ref[, sigs_to_use], 
                                output.dir = output_prefix,
                                max.level=length(sigs_to_use)-1,
                                p.thresh=0.05/length(sigs_to_use),
                                max.subsets=1000,
                                num.parallel.samples = 2,
                                mc.cores.per.sample=20,
                                seed=8717,
                                use.sig.presence.test=use_sig_presence_to_assign_sig, #default setting is FALSE
                                drop.low.mut.samples = drop_low_mut_samples) #removes samples with a total of <100 mutations IF = TRUE

cat("Preparing results \n")

results=sparse_out$proposed.assignment
results=rbind(results,colSums(results)) # Add a row with the total number of mutations
row.names(results)[nrow(results)]="total_mutations"

# State the result in percentage instead of number of mutations attributed to each signature
for (i in 1:ncol(results)) {
  results[1:(nrow(results)-1),i]=results[1:(nrow(results)-1),i]/results[nrow(results),i]  #State the result in percentage
  
}


# Add a row to the result file with the cosine similarity
recon_dist=sparse_out$reconstruction.distances
prop_assignment=recon_dist$proposed.assignment
results=rbind(results, prop_assignment$cosine)
row.names(results)[nrow(results)]="cosine"


#Run a presence test of each of the signatures and append to the results file
new_row=matrix(data=NA, ncol=ncol(results), nrow=length(sigs_to_use))
row.names(new_row)=paste0(sigs_to_use, "_p-value")


cat("Running signature presence test for all included sigs \n")
for (i in 1:length(sigs_to_use)) { #Loop through and do a presence test for each of the included sigs
  presence=SignaturePresenceTest(spectra=mut_mat, sigs=cosmic_ref[, sigs_to_use], target.sig.index = sigs_to_use[i], seed=8717, mc.cores=20)
  for (j in 1:length(presence)) { #Loop through each sample
    new_row[i,j]=unlist(presence[[j]][4])
  }
}

results=rbind(results, new_row)

cat("Saving results \n")

write.table(results, file=paste0(output_prefix,"mSigAct_assignment.txt"), row.names=TRUE, col.names = TRUE, sep="\t")

save(sparse_out, file=paste0(output_prefix,"sparse_out.Rdata"))


cat("Done \n")