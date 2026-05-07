# Figure 8 — Step 1: Train scVI/scANVI Reference Model
#
# Input:  snMultiome_atlas_reanalyzed.h5ad   (fetal brain reference)
# Output: models/scVI_ref/                   (scVI model)
#         models/scANVI_ref/                 (scANVI model)
#         models/adata_ref_trained.h5ad      (reference with UMAP embedding)
# ==============================================================================

import warnings; warnings.simplefilter("ignore")
import os
import scanpy as sc
import scvi


# ==============================================================================
# CONFIG
# ==============================================================================

BASE             = "/reference_mapping"
REF_H5AD         = f"{BASE}/snMultiome_atlas_reanalyzed.h5ad"
SCVI_MODEL_DIR   = f"{BASE}/models/scVI_ref"
SCANVI_MODEL_DIR = f"{BASE}/models/scANVI_ref"
REF_OUT_H5AD     = f"{BASE}/models/adata_ref_trained.h5ad"

CELLTYPE_KEY = "labels"   # obs column: cell type labels
BATCH_KEY    = "dataset"  # obs column: batch correction

N_TOP_GENES  = 2500
N_LATENT     = 30
N_LAYERS     = 2
MAX_EPOCHS_SCVI  = 400
MAX_EPOCHS_SCANVI = 20

os.makedirs(SCVI_MODEL_DIR,   exist_ok=True)
os.makedirs(SCANVI_MODEL_DIR, exist_ok=True)
scvi.settings.seed = 94705


# ==============================================================================
# 1. LOAD REFERENCE
# ==============================================================================

adata_ref = sc.read(REF_H5AD)
if adata_ref.raw is not None:
    adata_ref = adata_ref.raw.to_adata()

adata_ref.var_names_make_unique()
adata_ref.layers["counts"] = adata_ref.X.copy()
print(f"Reference: {adata_ref.shape[0]:,} cells × {adata_ref.shape[1]:,} genes")


# ==============================================================================
# 2. NORMALISE + SELECT HVGs
# ==============================================================================

sc.pp.normalize_total(adata_ref, exclude_highly_expressed=True)
sc.pp.log1p(adata_ref)
sc.pp.highly_variable_genes(adata_ref, n_top_genes=N_TOP_GENES,
                             batch_key=BATCH_KEY, subset=True)
print(f"HVGs selected: {adata_ref.n_vars}")


# ==============================================================================
# 3. TRAIN scVI
# ==============================================================================

scvi.model.SCVI.setup_anndata(adata_ref, batch_key=BATCH_KEY, layer="counts")

vae_ref = scvi.model.SCVI(
    adata_ref,
    use_layer_norm    = "both",
    use_batch_norm    = "none",
    encode_covariates = True,
    dropout_rate      = 0.2,
    n_layers          = N_LAYERS,
    n_latent          = N_LATENT,
)
vae_ref.train(max_epochs=MAX_EPOCHS_SCVI, early_stopping=True)
vae_ref.save(SCVI_MODEL_DIR, overwrite=True)
print(f"scVI model saved → {SCVI_MODEL_DIR}")


# ==============================================================================
# 4. TRAIN scANVI
# ==============================================================================

adata_ref.obs["labels_scanvi"] = adata_ref.obs[CELLTYPE_KEY].values

vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
    vae_ref,
    unlabeled_category = "Unknown",
    labels_key         = "labels_scanvi",
)
vae_ref_scan.train(max_epochs=MAX_EPOCHS_SCANVI, n_samples_per_label=100)
vae_ref_scan.save(SCANVI_MODEL_DIR, overwrite=True)
print(f"scANVI model saved → {SCANVI_MODEL_DIR}")


# ==============================================================================
# 5. EMBED + SAVE REFERENCE
# ==============================================================================

adata_ref.obsm["X_scANVI"] = vae_ref_scan.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scANVI")
sc.tl.umap(adata_ref)
adata_ref.write_h5ad(REF_OUT_H5AD)
print(f"Reference AnnData saved → {REF_OUT_H5AD}")