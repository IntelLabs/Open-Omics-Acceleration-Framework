import time
import daal4py
#from sklearnex import patch_sklearn
#patch_sklearn()

import os
n_cores = os.cpu_count()
print("Number of cores: ", n_cores)

import numpy as np
import scanpy as sc
from sklearn.cluster import KMeans

import mkl
import math

import fastpp
import sc_nbrs

import warnings
warnings.filterwarnings('ignore', 'Expected ')
warnings.simplefilter('ignore')

os.environ["NUMEXPR_MAX_THREADS"] = str(n_cores)
os.environ["NUMBA_NUM_THREADS"] = str(n_cores)

# Add path to input file here.
input_file = "/data/1M_brain_cells_10X.sparse.h5ad"
# USE_FIRST_N_CELLS = 10000
USE_FIRST_N_CELLS = -1           # -1 indicates use whole file

# marker genes
MITO_GENE_PREFIX = "mt-"              # Mitochondrial gene prefix
markers = ["Stmn2", "Hes1", "Olig1"]  # Marker genes for visualization

# filtering cells (Filter out genes outside [Min, Max])
min_genes_per_cell = 200 
max_genes_per_cell = 6000

# filtering genes
n_top_genes = 4000                    # Number of highly variable genes to retain

# PCA
n_components = 50

# Batched PCA
pca_train_ratio = 0.35               # percentage of cells to use for PCA training in batched PCA
n_pca_batches = 1

# t-SNE
tsne_n_pcs = 20                     # Number of principal components for t-SNE

# k-means
k = 35                              # Number of clusters for k-means

# KNN
n_neighbors = 15                    # Number of nearest neighbors for KNN graph
knn_n_pcs = 50                      # Number of principal components for nearest neighbors

# UMAP
umap_min_dist = 0.3 
umap_spread = 1.0

# Number jobs 
# sc.settings.n_jobs=112               # Set this value to number of hyperthread cores
sc.settings.n_jobs=n_cores               # Set this value to number of hyperthread cores

start = time.time()

pp_start = time.time()

adata = fastpp.loadpreprocess(input_file, USE_FIRST_N_CELLS, min_genes_per_cell, max_genes_per_cell, 1, 1e4, markers, n_top_genes)

mito_genes = adata.var_names.str.startswith(MITO_GENE_PREFIX)
n_counts = fastpp.sum(adata.X, axis=1)
adata.obs['percent_mito'] = np.array(np.sum(adata[:, mito_genes].X, axis=1)) / n_counts
adata.obs['n_counts'] = n_counts

# regression
fastpp.regress_out(adata, ['n_counts', 'percent_mito'])
fastpp.scale(adata, 10)

print("Total Preprocess time : %s" % (time.time() - pp_start))

pca_start = time.time()

from sklearn.decomposition import PCA

train_size = math.ceil(adata.X.shape[0] * pca_train_ratio)
pca = PCA(n_components=n_components).fit(adata.X[:train_size])

embeddings = np.zeros((adata.X.shape[0], n_components))
batch_size = int(embeddings.shape[0] / n_pca_batches)

for batch in range(n_pca_batches):
    start_idx = batch * batch_size
    end_idx = start_idx + batch_size

    if(adata.X.shape[0] - end_idx < batch_size):
        end_idx = adata.X.shape[0]

    embeddings[start_idx:end_idx,:] = np.asarray(pca.transform(adata.X[start_idx:end_idx]))

adata.obsm["X_pca"] = embeddings

pca_time = time.time()
print("Total pca time : %s" % (pca_time-pca_start))

tsne_time = time.time()

daal4py.daalinit(sc.settings.n_jobs)
# from scanpy.tools._utils import _choose_representation
from sklearn.manifold import TSNE
sc.tl.tsne(adata, n_pcs=tsne_n_pcs)
# X = _choose_representation(adata, n_pcs=tsne_n_pcs)
# X_tsne = TSNE().fit_transform(X.astype(np.float32))
# adata.obsm['X_tsne'] = X_tsne

#adata.obsm['X_tsne'] = TSNE().fit_transform(adata.obsm["X_pca"][:,:tsne_n_pcs].astype(np.float32))

print("TSNE time : %s" % (time.time() - tsne_time))

kmeans_time = time.time()
# print("Num of jobs for KMeans:", sc.settings.n_jobs)
kmeans = KMeans(n_clusters=k, random_state=0).fit(adata.obsm['X_pca'])
adata.obs['kmeans'] = kmeans.labels_.astype(str)
print("KMeans time : %s" % (time.time() - kmeans_time))

sc.pl.tsne(adata, color=["kmeans"], save="_kmeans.png")
sc.pl.tsne(adata, color=["Stmn2_raw"], color_map="Blues", vmax=1, vmin=-0.05, save="_Stmn2_raw.png")
sc.pl.tsne(adata, color=["Hes1_raw"], color_map="Blues", vmax=1, vmin=-0.05, save="_Hes1_raw.png")


neighbor_time = time.time()
sc_nbrs.neighbors(adata, n_neighbors=n_neighbors, n_pcs=knn_n_pcs, method='sklearn')
print("neighbors time : %s" % (time.time() - neighbor_time))

umap_time = time.time()
sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread)
print("UMAP time : %s" % (time.time() - umap_time))

########################## Adjacency matrix ################################

from scanpy._utils import _choose_graph

adjacency = _choose_graph(adata, obsp=None, neighbors_key=None)

sources, targets = adjacency.nonzero()
weights = adjacency[sources, targets]
if isinstance(weights, np.matrix):
    weights = weights.A1


##################### Katana ###############################################
louvain_time = time.time()

import pandas as pd
from natsort import natsorted

from katana.local import Graph
from katana.local.analytics import louvain_clustering, LouvainClusteringStatistics, LouvainClusteringPlan
from katana.local.import_data import from_edge_list_arrays

import katana.local
katana.local.initialize()
katana.set_active_threads(sc.settings.n_jobs)

property_dict = {"value": weights}
graph = from_edge_list_arrays(sources, targets, property_dict)

enable_vf = False
modularity_threshold_per_round = 0.0001
modularity_threshold_total = 0.0001
max_iterations = 100000
min_graph_size = 0

louvain_plan = LouvainClusteringPlan.do_all(enable_vf, modularity_threshold_per_round, modularity_threshold_total, max_iterations, min_graph_size)
# louvain_plan = LouvainClusteringPlan.deterministic(enable_vf, modularity_threshold_per_round, modularity_threshold_total, max_iterations, min_graph_size)

louvain_clustering(graph, "value", "output", plan=louvain_plan)
stats = LouvainClusteringStatistics(graph, "value", "output")
print(stats)
groups = graph.get_node_property("output").to_numpy().astype('int')

adata.obs['louvain'] = pd.Categorical(values=groups.astype('U'), categories=natsorted(map(str, np.unique(groups))),)
print("louvain time (Katana python): %s" % (time.time() - louvain_time))


leiden_time = time.time()
##################################### Katana ##########################################
from katana.local.analytics import leiden_clustering, LeidenClusteringStatistics, LeidenClusteringPlan
graph = from_edge_list_arrays(sources, targets, property_dict)

enable_vf = True
modularity_threshold_per_round = 0.0001
modularity_threshold_total = 0.0001  
max_iterations = 1000000
min_graph_size = 0
resolution = 1.0
randomness = 0
 
leiden_plan = LeidenClusteringPlan.deterministic(enable_vf, modularity_threshold_per_round, modularity_threshold_total, max_iterations, min_graph_size, resolution, randomness)
# leiden_plan = LeidenClusteringPlan.do_all(enable_vf, modularity_threshold_per_round, modularity_threshold_total, max_iterations, min_graph_size, resolution, randomness)
leiden_clustering(graph, "value", "leiden_output", plan=leiden_plan)
stats = LeidenClusteringStatistics(graph, "value", "leiden_output")
print(stats)
groups = graph.get_node_property("leiden_output").to_numpy().astype('int')

adata.obs['leiden'] = pd.Categorical(values=groups.astype('U'), categories=natsorted(map(str, np.unique(groups))),)
print("leiden time (Katana python): %s" % (time.time() - leiden_time))

sc.pl.umap(adata, color=["louvain"], save="_louvain.png")
sc.pl.umap(adata, color=["leiden"], save="_leiden.png")
sc.pl.umap(adata, color=["Stmn2_raw"], color_map="Blues", vmax=1, vmin=-0.05, save="_Stmn2_raw.png")
sc.pl.umap(adata, color=["Hes1_raw"], color_map="Blues", vmax=1, vmin=-0.05, save="_Hes1_raw.png")


re_start = time.time()

adata = adata[adata.obs["Hes1_raw"] > 0.0, :]
print(adata.X.shape)

sc.tl.pca(adata, n_comps=n_components)
sc_nbrs.neighbors(adata, n_neighbors=n_neighbors, n_pcs=knn_n_pcs, method='sklearn')
sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread)

#################### adjacency #############
from scanpy._utils import _choose_graph
adjacency = _choose_graph(adata, obsp=None, neighbors_key=None)

sources, targets = adjacency.nonzero()
weights = adjacency[sources, targets]
if isinstance(weights, np.matrix):
    weights = weights.A1


##################### Katana ########################################

property_dict = {"value": weights}
graph = from_edge_list_arrays(sources, targets, property_dict)

enable_vf = False
modularity_threshold_per_round = 0.0001
modularity_threshold_total = 0.0001   
max_iterations = 100000
min_graph_size = 0
resolution = 1.0
randomness = 0
 
leiden_plan = LeidenClusteringPlan.deterministic(enable_vf, modularity_threshold_per_round, modularity_threshold_total, max_iterations, min_graph_size, resolution, randomness)
leiden_clustering(graph, "value", "leiden_output", plan=leiden_plan)
groups = graph.get_node_property("leiden_output").to_numpy().astype('int')
adata.obs['leiden'] = pd.Categorical(values=groups.astype('U'), categories=natsorted(map(str, np.unique(groups))),)

sc.pl.umap(adata, color=["leiden"], save="_re_leiden.png")
sc.pl.umap(adata, color=["Olig1_raw"], color_map="Blues", vmax=1, vmin=-0.05, save="_re_Olig1_raw.png")

print("Total reanalysis time : %s" % (time.time() - re_start))

print("Total notebook time: %s" % (time.time() - start))

