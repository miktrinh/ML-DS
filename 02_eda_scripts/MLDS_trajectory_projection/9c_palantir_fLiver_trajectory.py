# Import modules
import palantir
import scanpy as sc
import numpy as np
import pandas as pd
import os

# Plotting 
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

import scipy.stats
import anndata
import matplotlib as mpl
from matplotlib.axes._axes import _log as matplotlib_axes_logger
from scipy import sparse
matplotlib_axes_logger.setLevel('ERROR')

# Inline plotting
#%matplotlib inline

sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['image.cmap'] = 'Spectral_r'
warnings.filterwarnings(action="ignore", module="matplotlib", message="findfont")

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 800)


sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
# Set up the plot config for viewing the annotation clearly.
sc.settings.set_figure_params(dpi=120, dpi_save=1000)
#sc.logging.print_versions()


# Reset random seed
np.random.seed(5)


##--------------------------##
##    Helper functions    ####
##--------------------------##

# ad = adata object
# ann_key = column name containing cell type annotation
# start_cell = user's defined cellID as the start of the trajectory
# marker_genes = list of key marker genes for plotting
# result_out_fp = path to output the results
# num_waypoints = number of waypoints (palantir parameter)
# n_top_genes = number of top genes for clustering
# n_components = number of PCs
# runHarmony = logical, should harmony be run for integration
# harmonyVar = if runHarmony == T, what are the variables for batch correction

## output: palantir_results, adata_object


def build_palantir_traj_v2(ad,ann_key,start_cell,marker_genes=['CD34', 'MPO', 'GATA1', 'IRF8','HBB','PF4','KLF1','FLI1'],
                           result_out_fp=None,num_waypoints=1200, n_top_genes=1500,n_comps=50,n_components=5,terminal_states=None,
                           runHarmony=False,harmonyVar=None,
                           use_early_cell_as_start=False,determine_termState=False,terminal_celltypes=[]):
    
    ### Data processing - assuming input anndata contains raw count (not normalized / scaled)
    print('Data processing...')
    # Normalisation
    sc.pp.normalize_per_cell(ad)
    # log1p transformation, palantir function uses a pseudocount of 0.1 instead of 1.
    palantir.preprocess.log_transform(ad)
    
    # Select highly variable genes
    sc.pp.highly_variable_genes(ad, n_top_genes=n_top_genes)
    
    
    ### PCA
    print('Computing PCA...')
    sc.pp.pca(ad,n_comps=n_comps)
    sc.pl.pca_scatter(ad, color='CD34')
    
    if runHarmony:
        ### Harmony
        print('Running Harmony')
        sc.external.pp.harmony_integrate(ad,harmonyVar)
        
    
    ### Compute Diffusion maps
    print('Computing Diffusion maps...') 
    # Run diffusion maps
    if runHarmony:
        dm_res = palantir.utils.run_diffusion_maps(ad, n_components=n_components,pca_key="X_pca_harmony")
        sc.pp.neighbors(ad,use_rep='X_pca_harmony')     
    else:
        dm_res = palantir.utils.run_diffusion_maps(ad, n_components=n_components)
        sc.pp.neighbors(ad)  
    
    # The low dimensional embeddeing of the data is estimated based on the eigen gap using the following function
    ms_data = palantir.utils.determine_multiscale_space(ad)
    
    ### Visualization
    # Recommend using FDG, but UMAP is a good alternative for visualization
    
    sc.tl.umap(ad)
    # Use scanpy functions to visualize umaps or FDL
    sc.pl.embedding(ad, basis='umap',color=ann_key,palette='tab10')
    
    # Diffusion map visualisation
    palantir.plot.plot_diffusion_components(ad)
    plt.show()
    
    ### MAGIC imputation
    print('MAGIC imputation...')
    #imputed_X = palantir.utils.run_magic_imputation(ad)
    #sc.pl.embedding(ad, basis='umap', layer='MAGIC_imputed_data',
    #               color=marker_genes)
    
    
    ### Plot start cells
    #palantir.plot.highlight_cells_on_umap(ad, {'start_cell':start_cell})
    ad.obs['isStartCell'] = [x == start_cell for x in ad.obs.index.tolist()]
    palantir.plot.highlight_cells_on_umap(ad,'isStartCell')
    
    plt.show()
    
    if determine_termState and (len(terminal_celltypes) > 0) and not terminal_states:
        terminal_states = []
        for ct in terminal_celltypes:
            termCell = palantir.utils.early_cell(ad, ct, ann_key)
            terminal_states.append(termCell)
   
    #ad.obs['isStartCell'] = [x == start_cell for x in ad.obs.index.tolist()]
    #palantir.plot.highlight_cells_on_umap(ad,'isStartCell')
    
    ### Running Palantir
    print('Computing trajectory with Palantir...')
    pr_res = palantir.core.run_palantir(ad, start_cell, num_waypoints=num_waypoints,
                                        terminal_states=terminal_states,
                                       use_early_cell_as_start=use_early_cell_as_start)
    
    ### Branch selection
    masks = palantir.presults.select_branch_cells(ad, eps=0)
    palantir.plot.plot_branch_selection(ad)
    plt.show()
    
    ### Extract Palantir results
    mask_output = pd.DataFrame.from_records(ad.obsm['branch_masks'],index=ad.obsm['branch_masks'].index)
    results = pr_res.branch_probs
    
    mdat = pd.merge(ad.obs,mask_output,left_index=True, right_index=True)
    mdat = pd.merge(mdat,results,left_index=True, right_index=True)
    
    ## Save output
    if result_out_fp:
        mdat.to_csv(result_out_fp)
            
    palantir.plot.plot_palantir_results(ad, s=3)
    plt.show()
    
    return(mdat,ad)



def standard_clustering(adata,ann_key=[],n_neighbors=10, n_pcs=40):
    import scanpy as sc
    print('Data processing...')
    # Normalisation
    #sc.pp.filter_genes(adata, min_cells=100)|
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata)#, min_mean=0.0125, max_mean=3, min_disp=0.5)
    #sc.pl.highly_variable_genes(adata)
    sc.pp.scale(adata)#, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    if len(ann_key) > 0:
        sc.pl.umap(adata,color=ann_key)
    else:
        sc.pl.umap(adata)
        
    return adata


##---------------------##
##    Import data    ####
##---------------------##

adata = sc.read_h5ad('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/Results/2_annotation/liver/jun24/liver_clean_processed_annotated_noUnknowns_0724.h5ad')
print(adata.obs['Phase'].value_counts())
print(adata.obs.donorID[adata.obs.Genotype == 'T21'].value_counts())

# Plot UMAP
sc.pl.embedding(adata, basis='umap',color=['donorID'],palette='tab10')

cells_toKeep = pd.read_csv('/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/jul24/wholeTranscriptome/1_palantir_output/fLiver_cells_downsampled.csv')
print(adata.obs['cellID'].isin(cells_toKeep.cellID[cells_toKeep.cells_toKeep==True]).value_counts())

##--------------------------------------------------##
##    Build Trajectory with Palantir - 2n only    ####
##--------------------------------------------------##

# # Subset to keep revelant cell types only

# celltypes_toKeep = ['HSC_MPP','MEMP_MEP','EE','ME','LE',
#                     'Mast.cell','earlyMK','MK',
#                    'LMPP_ELP','pro.B.cell','pre.B.cell','B.cell',
#                     'CMP_GMP','proMono','Monocyte']
#                     #'Macrophage','myelocyte','promyelocyte','Kupffer.cell']


# adata_2n = adata[(adata.obs['Genotype'] == 'diploid') & 
#                  (adata.obs['finalAnn_broad'].isin(celltypes_toKeep)) &
#                  (adata.obs['cellID'].isin(cells_toKeep.cellID[cells_toKeep.cells_toKeep==True]))
#                 ]
# print(adata_2n.obs.donorID.value_counts())
# print(adata_2n.obs.assay.value_counts())




# num_waypoints = 1200
# runHarmony = True
# harmonyVar = ['donorID','Phase','assay']
# terminal_states = None
# annotation_key = 'annot_jun24'
# use_early_cell_as_start=False
# n_components=10

# ad = adata_2n.copy()
# #result_out_fp = '/home/jovyan/sc126/tmp_2312/1_palantir_output/diploid_haemTraj.csv' # main output is a csv file
# result_out_fp = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/jul24/wholeTranscriptome/1_palantir_output/inhouse2n_haemTraj_2407.csv'

# # define start cell
# # eg. I defined start cells as HSC with highest expression of CD34
# cd34_expr = ad.X[ad.obs['finalAnn_broad'] == 'HSC_MPP',np.where(ad.var.index == 'CD34')[0].tolist()[0]].toarray()
# # Find cell ID with highest expression level for CD34
# start_cell = ad.obs.loc[ad.obs['finalAnn_broad'] == 'HSC_MPP',:].index[cd34_expr.argmax()]

# ## Define terminal states
# terminal_celltypes = ['Monocyte','LE','MK','Mast.cell','B.cell']
# determine_termState = True
# palantir_results, ad = build_palantir_traj_v2(ad,ann_key=annotation_key,start_cell=start_cell,
#                                            result_out_fp=result_out_fp,num_waypoints=num_waypoints,n_components=n_components,
#                                            runHarmony=runHarmony,harmonyVar=harmonyVar,
#                                               terminal_states=terminal_states,
#                                               determine_termState = determine_termState,
#                                               terminal_celltypes=terminal_celltypes)




##---------------------------------------------------------------------##
##    Build Trajectory with Palantir - EE/MK/Mast lineage 2n only    ####
##---------------------------------------------------------------------##

# Subset to keep revelant cell types only

celltypes_toKeep = ['HSC_MPP','MEMP_MEP','EE','ME','LE',
                    'Mast.cell','earlyMK','MK']

adata_2n = adata[(adata.obs['Genotype'] == 'diploid') & 
                 (adata.obs['finalAnn_broad'].isin(celltypes_toKeep)) &
                 (adata.obs['cellID'].isin(cells_toKeep.cellID[cells_toKeep.cells_toKeep==True]))
                ]
print(adata_2n.obs.donorID.value_counts())
print(adata_2n.obs.assay.value_counts())




num_waypoints = 1200
runHarmony = True
harmonyVar = ['donorID','Phase','assay']
terminal_states = None
annotation_key = 'annot_jun24'
use_early_cell_as_start=False
n_components=10

ad = adata_2n.copy()
#result_out_fp = '/home/jovyan/sc126/tmp_2312/1_palantir_output/diploid_haemTraj.csv' # main output is a csv file
result_out_fp = '/lustre/scratch125/casm/team274sb/mt22/Aneuploidy/MLDS_scRNAseq/Results/trajectoryProjection/fLiver2n_allcells/jul24/wholeTranscriptome/1_palantir_output/inhouse2n_MK.LE.Mast.Traj_2407.csv'

# define start cell
# eg. I defined start cells as HSC with highest expression of CD34
cd34_expr = ad.X[ad.obs['finalAnn_broad'] == 'HSC_MPP',np.where(ad.var.index == 'CD34')[0].tolist()[0]].toarray()
# Find cell ID with highest expression level for CD34
start_cell = ad.obs.loc[ad.obs['finalAnn_broad'] == 'HSC_MPP',:].index[cd34_expr.argmax()]

## Define terminal states
terminal_celltypes = ['LE','MK','Mast.cell']
determine_termState = True
palantir_results, ad = build_palantir_traj_v2(ad,ann_key=annotation_key,start_cell=start_cell,
                                           result_out_fp=result_out_fp,num_waypoints=num_waypoints,n_components=n_components,
                                           runHarmony=runHarmony,harmonyVar=harmonyVar,
                                              terminal_states=terminal_states,
                                              determine_termState = determine_termState,
                                              terminal_celltypes=terminal_celltypes)
