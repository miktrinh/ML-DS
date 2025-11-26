module use --append /software/CASM/modules/modulefiles
module load jbrowse_rasterize
shopt -s expand_aliases

alias jbrowse-rasterize='jbrowse-rasterize-exec node /app/jbrowse_rasterize.js'

alias jbrowse-rasterize-exec='singularity exec --cleanenv --home /nfs/users/nfs_m/mt22:/nfs/users/nfs_m/mt22 --bind /nfs:/nfs --bind /lustre:/lustre --bind /nfs/users/nfs_m/mt22:/nfs/users/nfs_m/mt22 /software/CASM/singularity/jbrowse_rasterize/jbrowse_rasterize_2.4.0.sif'

alias jbrowse-rasterize-shell='singularity shell --cleanenv --home /nfs/users/nfs_m/mt22:/nfs/users/nfs_m/mt22 --bind /nfs:/nfs --bind /lustre:/lustre --bind /nfs/users/nfs_m/mt22:/nfs/users/nfs_m/mt22 /software/CASM/singularity/jbrowse_rasterize/jbrowse_rasterize_2.4.0.sif'

bedFile=$1
outDir=$2

jbrowse-rasterize --locs=$bedFile --outdir=$outDir --width=1200 --imgType=png --passwdFile=/lustre/scratch125/casm/team274sb/mt22/generalScripts/passwdFile.txt --timeout 30 --highlight
