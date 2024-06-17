import scanpy as sc

def prepare_adata(gene_mtx, df_counts, centroids, min_counts=10, min_cells=3):
    
    adata = sc.AnnData(X=gene_mtx.copy())
    df_counts = df_counts[df_counts.index.isin(gene_mtx.index)]

    df_counts = df_counts.reindex(gene_mtx.index)
    adata.obs = df_counts.copy()
    centroids = centroids.reindex(adata.obs.index)
    ### sort cell_id in centroid to match the order in adata.obs

    adata.obs[['x_location', 'y_location']] = centroids[['centroid_x', 'centroid_y']].copy()
    adata.var['genes'] = gene_mtx.columns.astype(str)
    
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)

    return adata


def process_adata(adata, target_sum=1e4, n_comps=21, n_neighbors=10, resolution=0.3):
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=21)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added='leiden'+str(resolution), resolution=0.3)
    return adata


