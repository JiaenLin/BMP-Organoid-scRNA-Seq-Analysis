# ==============================================================================
# Figure 8: Label Transfer Analysis — Control vs. Treatment
#           (Transcriptional Similarity to Fetal Brain Reference)
#
# Corresponds to: Figure 8A–E
#   8A  — scANVI prediction score distribution
#   8B  — Fraction-overlap heatmaps (Control | Treatment, side-by-side)
#   8C  — Per-condition heatmaps (Control only, Treatment only)
#   8D  — Treatment − Control difference heatmap
#   8E  — Pseudobulk Pearson-R heatmap + paired Wilcoxon boxplot
#
# Input:
#   query_mapped.h5ad           scANVI-mapped query (HVG space)
#   snMultiome_atlas_reanalyzed.h5ad   Fetal reference
#   markers_hvg_all_celltypes.csv      Pre-computed top markers per cell type
#
# Output:
#   figures/prediction_score_dist.pdf/.png
#   figures/label_transfer_overlap_ctrl_vs_treat.pdf/.png
#   figures/label_transfer_overlap_ctrl.pdf
#   figures/label_transfer_overlap_treat.pdf
#   figures/label_transfer_diff_heatmap.pdf/.png
#   figures/pearson_r_heatmap.pdf/.png
#   figures/pearson_r_boxplot.pdf/.png
#   stats/pearson_r_per_celltype.csv
# Author: Jiaen Lin
# Date: 08-May-2026
# ==============================================================================

import warnings; warnings.simplefilter("ignore")
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr, wilcoxon
import scanpy as sc
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns


# ==============================================================================
# 0. CONFIG
# ==============================================================================

BASE       = "/reference_mapping"
OUTPUT_DIR = f"{BASE}//output"

# ── Input files ───────────────────────────────────────────────────────────────
QUERY_H5AD     = f"{OUTPUT_DIR}/query_mapped.h5ad"
RAW_QUERY_H5AD = f"{BASE}//data/.h5ad"
REF_H5AD       = f"{BASE}/snMultiome_atlas_reanalyzed.h5ad"
MARKER_CSV     = f"{BASE}/markers/markers_hvg_all_celltypes.csv"

# ── Column keys ───────────────────────────────────────────────────────────────
REF_CELLTYPE_KEY   = "type"
GROUP_KEY          = "Treatment"
QUERY_CELLTYPE_KEY = "class"
PRED_KEY           = "predictions"
COUNTS_LAYER       = "counts"

CTRL_LABEL  = "Control"
TREAT_LABEL = "BMP4"

# ── Analysis parameters ───────────────────────────────────────────────────────
N_TOP_MARKERS        = 50
CONFIDENCE_THRESHOLD = 0.5
MIN_CELLS_PSEUDOBULK = 1

# ── Derived column names ──────────────────────────────────────────────────────
CTRL_COL  = f"r_{CTRL_LABEL}"
TREAT_COL = f"r_{TREAT_LABEL}"

# ── Cell types excluded from all analyses ─────────────────────────────────────
EXCLUDE_TYPES = [
    "Endothelial", "Pericyte", "Smooth muscle cell", "VSMC", "Vascular",
    "Microglia", "Macrophage", "T cell", "B cell",
    "Unknown", "unknown", "Doublet", "Fibroblast", "Mesenchymal",
]

# ── Reference cell-type column ordering for heatmaps ─────────────────────────
REF_ORDER = [
    "Astrocyte-Fibrous", "Astrocyte-Protoplasmic",
    "RG-tRG", "RG-oRG", "RG-vRG", "IPC-EN",
    "EN-Newborn", "EN-IT-Immature", "EN-Non-IT-Immature",
    "EN-L2_3-IT", "EN-L4-IT", "EN-L5-IT", "EN-L5-ET",
    "EN-L5_6-NP", "EN-L6b", "EN-L6-IT", "EN-L6-CT",
    "Cajal-Retzius cell", "IPC-Glia",
    "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-SNCG",
    "IN-CGE-VIP", "IN-MGE-Immature", "IN-MGE-SST",
    "OPC", "Oligodendrocyte-Immature", "Oligodendrocyte",
]

os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/stats",   exist_ok=True)


# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

def fraction_matrix(adata_subset, query_key, pred_key, threshold=0.02):
    """Cross-tabulate query vs. predicted labels; return raw and threshold-masked fractions."""
    obs      = adata_subset.obs[[query_key, pred_key]].copy()
    ct       = pd.crosstab(obs[query_key], obs[pred_key])
    frac     = ct.div(ct.sum(axis=1), axis=0)
    frac_plot = frac.where(frac >= threshold)
    return frac, frac_plot

def make_row_labels(adata_subset, query_key, row_order):
    counts = adata_subset.obs[query_key].value_counts()
    return [f"{r}  (n={counts.get(r, 0):,})" for r in row_order]

def make_col_labels(adata_subset, pred_key, col_order):
    counts = adata_subset.obs[pred_key].value_counts()
    return [f"{c}  (n={counts.get(c, 0):,})" for c in col_order]

def make_annot(df):
    """Signed-float annotation; blank for NaN or zero."""
    return df.applymap(lambda v: f"{v:+.2f}" if (pd.notna(v) and v != 0) else "")

def pseudobulk_lognorm_fast(adata, groupby, layer="counts"):
    """Per-group mean log-normalised expression via sparse matrix ops."""
    mat = adata.layers[layer] if layer in adata.layers else adata.X
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)
    elif not isinstance(mat, sp.csr_matrix):
        mat = mat.tocsr()
    mat = mat.astype(np.float32)

    groups = np.asarray(adata.obs[groupby].values, dtype=str)
    uniq_grps, inverse = np.unique(groups, return_inverse=True)

    indicator = sp.csr_matrix(
        (np.ones(mat.shape[0], np.float32), (inverse, np.arange(mat.shape[0]))),
        shape=(len(uniq_grps), mat.shape[0]),
    )
    group_sums   = (indicator @ mat).toarray().astype(np.float64)
    row_totals   = group_sums.sum(axis=1, keepdims=True)
    mean_lognorm = np.log1p(group_sums / (row_totals + 1e-9) * 1e4)

    gene_idx = {g: i for i, g in enumerate(adata.var_names)}
    pb_dict  = {grp: mean_lognorm[i] for i, grp in enumerate(uniq_grps)}
    return pb_dict, gene_idx

def transfer_obs_to_raw(adata_mapped, adata_raw, cols):
    """Transfer obs columns from mapped (HVG) query to raw query by cell barcode."""
    common = adata_mapped.obs_names.intersection(adata_raw.obs_names)
    if len(common) == 0:
        raise ValueError("No common barcodes between mapped and raw query.")
    sub = adata_raw[common].copy()
    for col in cols:
        if col in adata_mapped.obs.columns:
            sub.obs[col] = adata_mapped[common].obs[col].values
        else:
            print(f"  ⚠️  Column '{col}' not in mapped query — skipping")
    print(f"  Common barcodes: {len(common):,}  |  transferred: {cols}")
    return sub

def plot_single_condition(frac_plot, adata_sub, col_order, row_order,
                           title, savepath=None):
    """Single-condition fraction-overlap heatmap."""
    fig, ax = plt.subplots(
        figsize=(max(8, frac_plot.shape[1] * 0.40), max(4, frac_plot.shape[0] * 0.6)),
        constrained_layout=True,
    )
    annot_arr = frac_plot.applymap(lambda v: f"{v:.2f}" if (pd.notna(v) and v > 0) else "")
    sns.heatmap(frac_plot, ax=ax, cmap="Greys", vmin=0, vmax=1,
                linewidths=0.3, linecolor="white",
                annot=annot_arr, fmt="", annot_kws={"size": 6},
                mask=frac_plot.isna(), cbar=True,
                cbar_kws={"shrink": 0.4, "ticks": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]})
    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)
    ax.set_xticks(np.arange(len(col_order)) + 0.5)
    ax.set_xticklabels(make_col_labels(adata_sub, PRED_KEY, col_order),
                       rotation=45, ha="right", fontsize=7)
    ax.set_xlabel("Predicted cell type  (n = cells predicted to that type)", fontsize=8, labelpad=6)
    ax.set_yticks(np.arange(len(row_order)) + 0.5)
    ax.set_yticklabels(make_row_labels(adata_sub, QUERY_CELLTYPE_KEY, row_order),
                       rotation=0, fontsize=8)
    ax.set_ylabel("Original cell type  (n = cells in query)", fontsize=8, labelpad=6)
    ax.collections[0].colorbar.set_label("Fraction", fontsize=8)
    ax.collections[0].colorbar.ax.tick_params(labelsize=8)
    fig.suptitle("Fraction overlap of labels transferred cell types to the cell types defined in this study",
                 fontsize=11)
    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
        print(f"Saved → {savepath}")
    plt.close()


# ==============================================================================
# 2. LOAD & PREPROCESS
# ==============================================================================

adata_q_hvg = sc.read(QUERY_H5AD)
adata_raw   = sc.read(RAW_QUERY_H5AD)
adata_q     = transfer_obs_to_raw(adata_q_hvg, adata_raw, cols=[PRED_KEY, "prediction_score"])
sc.pp.filter_cells(adata_q, min_genes=200)

adata_ref = sc.read(REF_H5AD)
if COUNTS_LAYER not in adata_ref.layers:
    print(f"  ⚠️  '{COUNTS_LAYER}' layer missing — saving .X as raw counts before normalising")
    adata_ref.layers[COUNTS_LAYER] = adata_ref.X.copy()
    sc.pp.normalize_total(adata_ref, target_sum=1e4, exclude_highly_expressed=True)
    sc.pp.log1p(adata_ref)

# ── Filter: excluded types + low-confidence predictions ──────────────────────
n_before  = adata_q.shape[0]
type_mask = ~adata_q.obs[PRED_KEY].isin(EXCLUDE_TYPES)
cell_mask = ~adata_q.obs[QUERY_CELLTYPE_KEY].isin(EXCLUDE_TYPES)
conf_mask = adata_q.obs["prediction_score"] >= CONFIDENCE_THRESHOLD
adata_q   = adata_q[type_mask & conf_mask & cell_mask].copy()
print(f"  Filtered {n_before:,} → {adata_q.shape[0]:,} cells")

obs_ctrl  = adata_q[adata_q.obs[GROUP_KEY] == CTRL_LABEL]
obs_treat = adata_q[adata_q.obs[GROUP_KEY] == TREAT_LABEL]


# 3. Fig S9A — PREDICTION SCORE DISTRIBUTION
fig, ax = plt.subplots(figsize=(5, 4))
for grp, color in [(CTRL_LABEL, "#4393C3"), (TREAT_LABEL, "#D6604D")]:
    scores = adata_q.obs.loc[adata_q.obs[GROUP_KEY] == grp, "prediction_score"]
    ax.hist(scores, bins=40, alpha=0.6, color=color,
            label=f"{grp} (mean={scores.mean():.3f})", density=True)
ax.set_xlabel("scANVI Prediction Score")
ax.set_ylabel("Density")
ax.set_title("Label Transfer Confidence")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(f"{OUTPUT_DIR}/figures/Fig8A_prediction_score_dist.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print("Saved → Fig8A_prediction_score_dist.pdf/.png")

# 4. Fig 8 and Fig S9 C-D — FRACTION-OVERLAP HEATMAPS

frac_ctrl,  frac_ctrl_plot  = fraction_matrix(obs_ctrl,  QUERY_CELLTYPE_KEY, PRED_KEY)
frac_treat, frac_treat_plot = fraction_matrix(obs_treat, QUERY_CELLTYPE_KEY, PRED_KEY)

# Unified row / column ordering
row_order = sorted(set(frac_ctrl.index) | set(frac_treat.index))
all_cols  = set(frac_ctrl_plot.columns) | set(frac_treat_plot.columns)
col_order = [c for c in REF_ORDER if c in all_cols]
col_order += [c for c in sorted(all_cols) if c not in col_order]

frac_ctrl_plot  = frac_ctrl_plot.reindex(index=row_order,  columns=col_order)
frac_treat_plot = frac_treat_plot.reindex(index=row_order, columns=col_order)

# ── Fig S9 C-D: Per-condition heatmaps ────────────────────────────────────────────
plot_single_condition(frac_ctrl_plot,  obs_ctrl,  col_order, row_order,
                      f"This study — {CTRL_LABEL}",
                      f"{OUTPUT_DIR}/figures/Fig8C_label_transfer_overlap_ctrl.pdf")
plot_single_condition(frac_treat_plot, obs_treat, col_order, row_order,
                      f"This study — {TREAT_LABEL}",
                      f"{OUTPUT_DIR}/figures/Fig8C_label_transfer_overlap_treat.pdf")

# ── Fig 8E: Treatment − Control difference heatmap ───────────────────────────
frac_ctrl_full  = frac_ctrl.reindex(index=row_order,  columns=col_order).fillna(0)
frac_treat_full = frac_treat.reindex(index=row_order, columns=col_order).fillna(0)
frac_diff       = frac_treat_full - frac_ctrl_full
frac_diff_plot  = frac_diff.where((frac_ctrl_full != 0) | (frac_treat_full != 0))

abs_max = np.nanmax(np.abs(frac_diff_plot.values))
fig, ax = plt.subplots(figsize=(max(8, len(col_order) * 0.45), max(4, len(row_order) * 0.6)))
sns.heatmap(frac_diff_plot, ax=ax,
            cmap="RdBu_r", vmin=-abs_max, vmax=abs_max,
            linewidths=0.3, linecolor="white",
            annot=make_annot(frac_diff_plot), fmt="", annot_kws={"size": 6},
            mask=frac_diff_plot.isna(), cbar=True,
            cbar_kws={"shrink": 0.5, "ticks": np.round(np.linspace(-abs_max, abs_max, 9), 2)})
ax.set_title(f"Label transfer fraction difference ({TREAT_LABEL} − {CTRL_LABEL})",
             fontsize=11, fontweight="bold", pad=10)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
ax.set_xlabel("Reference cell type (fetal atlas)", fontsize=9)
ax.set_ylabel("Cell type (this study)", fontsize=9)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=8)
cbar.ax.set_ylabel(f"Δ fraction  ({TREAT_LABEL} − {CTRL_LABEL})", fontsize=8, labelpad=8)
cbar.ax.axhline(0, color="black", linewidth=0.8)
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(f"{OUTPUT_DIR}/figures/Fig8D_label_transfer_diff_heatmap.{ext}",
                dpi=300, bbox_inches="tight")
plt.close()
print("Saved → Fig8D_label_transfer_diff_heatmap.pdf/.png")


# 5. Fig 8G — PSEUDOBULK PEARSON-R VS. FETAL REFERENCE

# ── Load top markers from pre-computed CSV ────────────────────────────────────
print(f"Loading markers from: {MARKER_CSV}")
marker_df   = pd.read_csv(MARKER_CSV)
top_markers = {
    ct: [g for g in grp["names"].tolist()[:N_TOP_MARKERS] if g in adata_q.var_names]
    for ct, grp in marker_df.groupby("cell_type")
}
print(f"  Loaded markers for {len(top_markers)} cell types")

if "logfoldchanges" in marker_df.columns:
    down = marker_df[marker_df["logfoldchanges"] < 0]
    if len(down):
        print(f"  ⚠️  {len(down)} downregulated genes in marker CSV:")
        for _, row in down.iterrows():
            print(f"    {row['cell_type']:30s} | {row['names']:20s}  logFC={row['logfoldchanges']:+.3f}")
    else:
        print("  ✅ All loaded markers are upregulated (logFC > 0)")

# ── Pseudobulk expression ─────────────────────────────────────────────────────
pb_fetal, ref_gene_idx = pseudobulk_lognorm_fast(adata_ref,   REF_CELLTYPE_KEY, COUNTS_LAYER)
pb_ctrl,  q_gene_idx   = pseudobulk_lognorm_fast(obs_ctrl,    PRED_KEY,         COUNTS_LAYER)
pb_treat, _            = pseudobulk_lognorm_fast(obs_treat,   PRED_KEY,         COUNTS_LAYER)

counts_ctrl  = obs_ctrl.obs[PRED_KEY].value_counts()
counts_treat = obs_treat.obs[PRED_KEY].value_counts()

# ── Pearson-R per shared cell type ────────────────────────────────────────────
shared_types = sorted(set(pb_fetal) & set(pb_ctrl) & set(pb_treat))
cor_rows = []

for ct in shared_types:
    n_ctrl_ct  = counts_ctrl.get(ct, 0)
    n_treat_ct = counts_treat.get(ct, 0)
    if min(n_ctrl_ct, n_treat_ct) < MIN_CELLS_PSEUDOBULK:
        print(f"  Skipping {ct}: too few cells (ctrl={n_ctrl_ct}, treat={n_treat_ct})")
        continue
    markers = top_markers.get(ct, [])
    if len(markers) < 5:
        print(f"  Skipping {ct}: only {len(markers)} shared markers")
        continue
    ref_idx    = [ref_gene_idx[g] for g in markers]
    q_idx      = [q_gene_idx[g]   for g in markers]
    fetal_expr = pb_fetal[ct][ref_idx]
    r_ctrl,  _ = pearsonr(fetal_expr, pb_ctrl[ct][q_idx])
    r_treat, _ = pearsonr(fetal_expr, pb_treat[ct][q_idx])
    cor_rows.append({
        "celltype": ct, CTRL_COL: r_ctrl, TREAT_COL: r_treat,
        "n_markers": len(markers),
        f"n_{CTRL_LABEL}": n_ctrl_ct, f"n_{TREAT_LABEL}": n_treat_ct,
    })
    print(f"  {ct:35s}  {CTRL_LABEL}={r_ctrl:.3f}  {TREAT_LABEL}={r_treat:.3f}"
          f"  (markers={len(markers)}, ctrl_n={n_ctrl_ct}, treat_n={n_treat_ct})")

cor_df = pd.DataFrame(cor_rows)
cor_df.to_csv(f"{OUTPUT_DIR}/stats/pearson_r_per_celltype.csv", index=False)
print(f"Saved → stats/pearson_r_per_celltype.csv  ({len(cor_df)} cell types)")

# ── Fig 8E (i): Pearson-R heatmap ─────────────────────────────────────────────
heat_data = cor_df.set_index("celltype")[[CTRL_COL, TREAT_COL]].rename(
    columns={CTRL_COL: CTRL_LABEL, TREAT_COL: TREAT_LABEL}
)
y_labels = [
    f"{ct}  ({CTRL_LABEL}: n={counts_ctrl.get(ct,0):,}  |  {TREAT_LABEL}: n={counts_treat.get(ct,0):,})"
    for ct in heat_data.index
]

fig, ax = plt.subplots(figsize=(6, max(4, len(cor_df) * 0.35)))
sns.heatmap(heat_data, annot=True, fmt=".2f", cmap="RdBu_r",
            center=0, vmin=-1, vmax=1, linewidths=0.5, ax=ax,
            annot_kws={"size": 8}, cbar_kws={"shrink": 0.5})
ax.set_title("Pearson R  (pseudobulk vs. fetal reference)", fontsize=10)
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha="center", fontsize=10)
ax.set_yticks(np.arange(len(heat_data)) + 0.5)
ax.set_yticklabels(y_labels, rotation=0, fontsize=8)
ax.collections[0].colorbar.set_label("Pearson r", fontsize=8)
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(f"{OUTPUT_DIR}/figures/Fig8E_pearson_r_heatmap.{ext}", dpi=300, bbox_inches="tight")
plt.close()
print("Saved → Fig8G_pearson_r_heatmap.pdf/.png")

# 6. FIG 8F — BOXPLOT + PAIRED WILCOXON TEST
ctrl_vals  = cor_df[CTRL_COL].values
treat_vals = cor_df[TREAT_COL].values
_, p_val   = wilcoxon(treat_vals, ctrl_vals, alternative="greater")

BLUE = "#4393C3"; RED = "#D6604D"
np.random.seed(42)
jitter = np.random.uniform(-0.08, 0.08, len(ctrl_vals))

fig, ax = plt.subplots(figsize=(4.5, 5.5))
for xpos, vals, color in [(0, ctrl_vals, BLUE), (1, treat_vals, RED)]:
    ax.boxplot(vals, positions=[xpos], widths=0.35, patch_artist=True,
               boxprops=dict(facecolor=color+"44", color=color, linewidth=1.5),
               medianprops=dict(visible=False),
               whiskerprops=dict(color=color), capprops=dict(color=color),
               flierprops=dict(marker=""))
    ax.scatter(xpos + jitter, vals, color=color, s=55, zorder=3,
               alpha=0.85, edgecolors="white", linewidths=0.5)
    ax.plot([xpos - 0.175, xpos + 0.175], [vals.mean()] * 2,
            color="black", lw=2.5, zorder=5)
    ax.text(xpos, vals.mean() + 0.04, f"{vals.mean():.3f}",
            ha="center", va="bottom", fontsize=9, fontweight="bold")

for c, t, j in zip(ctrl_vals, treat_vals, jitter):
    ax.plot([0 + j, 1 + j], [c, t], color="grey", alpha=0.2, lw=0.8, zorder=1)

y_top = max(ctrl_vals.max(), treat_vals.max()) + 0.07
ax.plot([0, 0, 1, 1], [y_top - 0.02, y_top, y_top, y_top - 0.02], color="black", lw=1.2)
ax.text(0.5, y_top + 0.01,
        f"p = {p_val:.4f}" if p_val >= 1e-4 else "p < 0.0001",
        ha="center", va="bottom", fontsize=9)

ax.set_xticks([0, 1]); ax.set_xticklabels([CTRL_LABEL, TREAT_LABEL], fontsize=12)
ax.set_xlim(-0.5, 1.5); ax.set_ylabel("Pearson R", fontsize=11)
ax.set_title(
    f"Transcriptional Similarity to Fetal Brain\n"
    f"(n={len(cor_df)} cell types, top {N_TOP_MARKERS} markers each)", fontsize=11
)
ax.axhline(0, color="lightgrey", lw=0.8, ls="--", zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(handles=[
    mpatches.Patch(facecolor=BLUE + "44", edgecolor=BLUE, label=CTRL_LABEL),
    mpatches.Patch(facecolor=RED  + "44", edgecolor=RED,  label=TREAT_LABEL),
], frameon=False, fontsize=10, loc="lower right")
plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(f"{OUTPUT_DIR}/figures/Fig8F_pearson_r_boxplot.{ext}", dpi=300)
plt.close()
