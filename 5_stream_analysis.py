# singularity shell -B $PWD/../../../ ../../../../shared/software/containers/stream

# ### DO NOT RUN
# #export the processed object from R to an anndata object
# #run in a ../software/r4.3.2-seuratv5_v2 R session
# seu.obj <- readRDS("../output/s3/allCells_clean_highRes_integrated.harmony_res1.4_dims45_dist0.15_neigh20_S3.rds")
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "majorID")
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "celltype")
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "clusterID2_integrated.harmony", metaAdd = "majCol")
# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/allCells_ID_disconected_highRes.csv", groupBy = "majCol", metaAdd = "labCol")
# 
# seu.obj[["RNA"]] <- as(object = seu.obj[["RNA"]], Class = "Assay")
# SaveH5Seurat(seu.obj, filename = "../output/s3/allCells_clean_highRes_integrated.harmony_res1.4_dims45_dist0.15_neigh20_S3.h5Seurat", verbose = TRUE, overwrite = T)
# Convert("../output/s3/allCells_clean_highRes_integrated.harmony_res1.4_dims45_dist0.15_neigh20_S3.h5Seurat", dest = "h5ad", overwrite = T)
# write.table(seu.obj$clusterID2_integrated.harmony, file = "./metaData/stream_lab.tsv", sep="\t", col.names = F)
# write.table(seu.obj$majCol, file = "./metaData/stream_col.tsv", sep="\t", col.names = F)


#load modules
import scanpy as sc
sc.__version__
import stream as st
st.__version__
from matplotlib import pyplot as plt
import pandas as pd


#set params and read in data
st.set_figure_params(dpi=80,style='white',figsize=[5.4,4.8],
                     rc={'image.cmap': 'viridis'})

adata = sc.read_h5ad("../output/s3/allCells_clean_highRes_integrated.harmony_res1.4_dims45_dist0.15_neigh20_S3.h5ad")
st.add_cell_labels(adata, file_name='./metaData/stream_col.tsv')
adata.obs[['label_color']] = adata.obs[['label']]
st.add_cell_labels(adata, file_name='./metaData/stream_lab.tsv')

colz_df = pd.read_csv('./metaData/allCells_ID_disconected_highRes.csv')
colz_df[['clusterID2_integrated']] = colz_df[['clusterID2_integrated.harmony']].astype(str)
colz_dict = dict(zip(colz_df.clusterID2_integrated, colz_df.majCol.astype(str)))

# adata.obs[['label']] = adata.obs[['celltype']]
# adata.obs[['label_color']] = adata.obs[['majCol']]
st.set_workdir(adata,'../output/stream_results')

#not working
# st.add_cell_colors(adata,file_name='./cell_label_color.tsv.gz')
# st.add_cell_labels(adata,file_path='./metaData/',file_name='stream_col.tsv')
# st.add_cell_labels(adata,file_name='./metaData/stream_lab.tsv')


### All skipped b/c already preprocessed data
# #clean anndata
# adata.var_names_make_unique()
# adata.obs_names_make_unique()

# #QC filter the dataset
# st.cal_qc(adata,assay='rna')
# st.filter_features(adata,min_n_cells = 5)
# st.filter_cells(adata,min_n_features= 100)
# st.normalize(adata,method='lib_size')
# st.log_transform(adata)
# st.remove_mt_genes(adata)

# #dim reduction
# st.select_variable_genes(adata,loess_frac=0.01,n_genes=2000)
# st.select_top_principal_components(adata,feature='var_genes',first_pc=True,n_pc=40)
# st.dimension_reduction(adata,method='se',feature='top_pcs',n_components=40,n_neighbors=15,n_jobs=8)

# st.plot_dimension_reduction(adata,color=['label'],
#                             show_graph=False,show_text=False)

# plt.savefig('../output/stream_results/foo.png')

# #trajectory analysis
# st.seed_elastic_principal_graph(adata,n_clusters=20)
# st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=False)
# st.plot_branches(adata,show_text=True)
# st.elastic_principal_graph(adata,epg_alpha=0.02,epg_mu=0.05,epg_lambda=0.01)


#use the Seurat embeddings
adata.obsm['top_pcs'] = adata.obsm['X_pca']
adata.obsm['X_dr'] = adata.obsm['X_umap.integrated.harmony']
adata.obsm['X_vis_umap'] = adata.obsm['X_umap.integrated.harmony'][:,:2]

#learn trajectories directly on UMAP from seurat
st.plot_visualization_2D(adata,method='umap',color=['majorID'],use_precomputed=True)
st.seed_elastic_principal_graph(adata,n_clusters=50,use_vis=True)
st.elastic_principal_graph(adata,epg_alpha=0.0001,epg_mu=0.025,epg_lambda=0.03,epg_n_processes=8)
st.plot_dimension_reduction(adata,color=['majorID'],n_components=2,show_graph=True,show_text=False)
plt.savefig('../output/stream_results/foo.png')

#extend leaf branches
st.extend_elastic_principal_graph(adata, epg_ext_mode='WeigthedCentroid',epg_ext_par=0.8)
st.plot_dimension_reduction(adata,color=['majorID'],n_components=2,show_graph=True,show_text=True)
plt.savefig('../output/stream_results/foo.png')

adata.uns['label_color'] = colz_dict
st.plot_stream(adata,root='S9',color=['label'],log_scale=True,fig_size=(8,6),fig_legend_ncol=1)
plt.savefig('../output/stream_results/foo.png')

st.plot_stream_sc(adata,root='S9',color=['label'],dist_scale=0.5,show_graph=True,show_text=False,fig_legend_ncol=3)
plt.savefig('../output/stream_results/foo.png')

