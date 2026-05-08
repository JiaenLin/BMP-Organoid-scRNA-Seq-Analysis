# ==============================================================================
# Figure 8 — Step 2: Map Query to scANVI Reference (Label Transfer)
# Input:  models/scVI_ref/               trained scVI model
#         models/scANVI_ref/             trained scANVI model
#         BMP_project_organoids.h5ad     raw query counts
# Output: output/query_mapped.h5ad       query with X_scANVI, predictions,
#                                        prediction_score, UMAP
#         output/scanvi_soft_probs.csv   per-cell soft probability matrix
# Author: Jiaen Lin
# Date: 08-May-2026
# ==============================================================================

import warnings; warnings.simplefilter("ignore")
import os
import numpy as np
import scanpy as sc
import scvi
import scipy.sparse as sp


# ==============================================================================
# CONFIG
# ==============================================================================

BASE             = "reference_mapping"
SCVI_MODEL_DIR   = f"{BASE}/models/scVI_ref"
SCANVI_MODEL_DIR = f"{BASE}/models/scANVI_ref"
QUERY_H5AD       = f"{BASE}/BMP_Project/data/BMP_project_organoids.h5ad"
REF_H5AD         = f"{BASE}/models/snMultiome_atlas_ref.h5ad"
OUTPUT_DIR       = f"{BASE}/BMP_Project/output"

MAX_EPOCHS_SCVI_QUERY   = 200
MAX_EPOCHS_SCANVI_QUERY = 100

os.makedirs(OUTPUT_DIR, exist_ok=True)
scvi.settings.seed = 94705


# ==============================================================================
# 1. LOAD REFERENCE + TRAINED MODELS
# ==============================================================================

adata_ref = sc.read(REF_H5AD)
adata_ref.var_names_make_unique()
print(f"Reference: {adata_ref.shape[0]:,} cells × {adata_ref.shape[1]:,} genes")

vae_ref      = scvi.model.SCVI.load(SCVI_MODEL_DIR,   adata_ref)
vae_ref_scan = scvi.model.SCANVI.load(SCANVI_MODEL_DIR, adata_ref)

# Retrieve keys used during reference training
scvi_batch_key    = vae_ref.registry["setup_args"].get("batch_key",  None)
scanvi_labels_key = vae_ref_scan.registry["setup_args"].get("labels_key", None)
print(f"  scVI  batch_key   = '{scvi_batch_key}'")
print(f"  scANVI labels_key = '{scanvi_labels_key}'")


# ==============================================================================
# 2. LOAD + PREPROCESS QUERY
# ==============================================================================

adata_query = sc.read(QUERY_H5AD)
adata_query.var_names_make_unique()

# Ensure raw integer counts are in layers["counts"]
if "counts" not in adata_query.layers:
    mat = adata_query.X
    arr = mat.toarray() if sp.issparse(mat) else np.array(mat)
    assert np.allclose(arr, arr.round(), atol=1e-2), \
        "adata.X does not appear to be raw integer counts. Assign counts layer manually."
    adata_query.layers["counts"] = adata_query.X.copy()

# Save raw before normalisation
if adata_query.raw is None:
    adata_query.raw = adata_query

sc.pp.normalize_total(adata_query, target_sum=1e4, exclude_highly_expressed=True)
sc.pp.log1p(adata_query)
sc.pp.highly_variable_genes(adata_query, n_top_genes=2500, subset=False)

# Cast all obs/var columns to str (required by scvi-tools)
for col in adata_query.obs.columns:
    adata_query.obs[col] = adata_query.obs[col].astype(str)
for col in adata_query.var.columns:
    adata_query.var[col] = adata_query.var[col].astype(str)

# Add required batch/label columns if absent
if scvi_batch_key and scvi_batch_key not in adata_query.obs.columns:
    adata_query.obs[scvi_batch_key] = "query"
if scanvi_labels_key and scanvi_labels_key not in adata_query.obs.columns:
    adata_query.obs[scanvi_labels_key] = "Unknown"


# ==============================================================================
# 3. MAP QUERY WITH scVI  (latent embedding)
# ==============================================================================

scvi.model.SCVI.prepare_query_anndata(adata_query, SCVI_MODEL_DIR)
vae_q = scvi.model.SCVI.load_query_data(adata_query, vae_ref)
vae_q.train(max_epochs=MAX_EPOCHS_SCVI_QUERY, plan_kwargs=dict(weight_decay=0.0))
adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()


# ==============================================================================
# 4. MAP QUERY WITH scANVI  (label transfer)
# ==============================================================================

scvi.model.SCANVI.prepare_query_anndata(adata_query, SCANVI_MODEL_DIR)
vae_q_scan = scvi.model.SCANVI.load_query_data(adata_query, vae_ref_scan)
vae_q_scan.train(
    max_epochs             = MAX_EPOCHS_SCANVI_QUERY,
    plan_kwargs            = dict(weight_decay=0.0),
    check_val_every_n_epoch = 10,
)

# Predicted labels + confidence scores
adata_query.obsm["X_scANVI"]        = vae_q_scan.get_latent_representation()
adata_query.obs["predictions"]      = vae_q_scan.predict()
soft                                 = vae_q_scan.predict(soft=True)
adata_query.obs["prediction_score"] = soft.max(axis=1).values
soft.index                           = adata_query.obs_names
soft.to_csv(f"{OUTPUT_DIR}/scanvi_soft_probs.csv")


# ==============================================================================
# 5. UMAP + SAVE
# ==============================================================================

sc.pp.neighbors(adata_query, use_rep="X_scANVI", key_added="scANVI_nn")
sc.tl.umap(adata_query, neighbors_key="scANVI_nn")

adata_query.write(f"{OUTPUT_DIR}/query_mapped.h5ad", compression="gzip")
print(f"Saved → {OUTPUT_DIR}/query_mapped.h5ad")
