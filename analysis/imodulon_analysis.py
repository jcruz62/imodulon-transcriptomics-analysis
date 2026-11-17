"""
iModulon analysis of Acinetobacter baumannii ATCC 19606
=======================================================

Reproduces and extends the Phase 4 analysis from the curated iModulon
matrices in phase2_metadata_curation/.

Run from the repository root:
    python analysis/imodulon_analysis.py

Outputs:
    figures/*.png                    plots used in the README
    results/summary.txt              key numbers printed by this script
    results/correlation_robustness.csv
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

# ---------------------------------------------------------------------------
# 0. Paths and settings
# ---------------------------------------------------------------------------
DATA_DIR = Path("phase2_metadata_curation")
FIGURE_DIR = Path("figures")
RESULTS_DIR = Path("results")
FIGURE_DIR.mkdir(exist_ok=True)   # exist_ok=True: no error if the folder already exists
RESULTS_DIR.mkdir(exist_ok=True)

RANDOM_SEED = 42
N_BOOTSTRAP = 2000

# Biological labels from the Phase 3 curation report.
CURATED_LABELS = {
    "IM_109": "Ribosome",
    "IM_139": "ATP synthase",
    "IM_104": "Arsenic resistance",
    "IM_45": "Type IV pili",
    "PmrA": "PmrA",
    "Fur": "Fur",
}

# The two relationships highlighted in the project.
KEY_PAIRS = [
    ("IM_109", "IM_139"),   # Ribosome vs ATP synthase (growth program)
    ("IM_104", "IM_45"),    # Arsenic resistance vs Type IV pili (trade-off)
]

plt.rcParams.update({
    "figure.dpi": 150,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
})

summary_lines = []


def report(line):
    """Print a line and also save it for results/summary.txt."""
    print(line)
    summary_lines.append(line)


def label(imodulon_name):
    """Return 'IM_109 (Ribosome)' style labels when a curated name exists."""
    if imodulon_name in CURATED_LABELS and CURATED_LABELS[imodulon_name] != imodulon_name:
        return f"{imodulon_name} ({CURATED_LABELS[imodulon_name]})"
    return imodulon_name


# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
# A matrix: rows = iModulons, columns = samples, values = activity level.
activity = pd.read_csv(DATA_DIR / "A.csv", index_col=0)

# iModulon table: one row per iModulon with its size and explained variance.
imodulon_table = pd.read_csv(DATA_DIR / "im_table.csv", index_col=0)
imodulon_table = imodulon_table.set_index(keys="name")

# Sample table: one row per RNA-seq run. Its index looks like
# "DRX489136_DRR505159"; the A matrix uses the first part ("DRX489136").
samples = pd.read_csv(DATA_DIR / "sample_table.csv", index_col=0)
samples = samples.set_index(keys="unique_identifier")
samples = samples.loc[activity.columns]   # put samples in the same order as A

report(f"iModulons: {activity.shape[0]} | samples: {activity.shape[1]} | "
       f"projects: {samples['project_name'].nunique()} | "
       f"conditions: {samples[['project_name', 'condition']].drop_duplicates().shape[0]}")

total_explained = imodulon_table["exp_var"].sum()
is_single_gene = imodulon_table["single_gene"].fillna(value=False).astype(bool)
report(f"Total explained variance: {total_explained:.1%}")
report(f"Single-gene iModulons: {is_single_gene.sum()} of {len(imodulon_table)} "
       f"({is_single_gene.mean():.0%})")

# ---------------------------------------------------------------------------
# 2. Figure 1: what the decomposition looks like
# ---------------------------------------------------------------------------
figure, (axis_left, axis_right) = plt.subplots(nrows=1, ncols=2, figsize=(11, 4.2))

# Left: how many genes each iModulon contains (log scale because sizes are skewed).
multi_gene_sizes = imodulon_table.loc[~is_single_gene, "n_genes"]
axis_left.hist(x=multi_gene_sizes, bins=np.logspace(start=np.log10(2), stop=np.log10(multi_gene_sizes.max()), num=20),
               color="#2E6FBA", edgecolor="white")
axis_left.set_xscale(value="log")
axis_left.set_xlabel("Genes per iModulon (log scale)")
axis_left.set_ylabel("Number of iModulons")
axis_left.set_title(f"{(~is_single_gene).sum()} multi-gene iModulons\n(+ {is_single_gene.sum()} single-gene, not shown)")

# Right: the 15 iModulons that explain the most variance.
top_by_variance = imodulon_table.sort_values(by="exp_var", ascending=False).head(n=15)
bar_colors = ["#D1495B" if name in CURATED_LABELS else "#9AA5B1" for name in top_by_variance.index]
axis_right.barh(y=[label(name) for name in top_by_variance.index][::-1],
                width=(top_by_variance["exp_var"] * 100)[::-1],
                color=bar_colors[::-1])
axis_right.set_xlabel("Explained variance (%)")
axis_right.set_title("Top 15 iModulons by explained variance\n(red = curated in this project)")
axis_right.tick_params(axis="y", labelsize=8)

figure.tight_layout()
figure.savefig(fname=FIGURE_DIR / "01_decomposition_overview.png", bbox_inches="tight")
plt.close(fig=figure)

# ---------------------------------------------------------------------------
# 3. Are the key correlations robust?
# ---------------------------------------------------------------------------
# Problem with the sample-level correlation: replicates of the same condition
# are not independent, and big projects contribute many near-identical points.
# Two fixes:
#   (a) average replicates first, so each condition counts once;
#   (b) bootstrap whole PROJECTS (cluster bootstrap) to get a confidence
#       interval that respects the study structure.
condition_key = samples["project_name"].astype(str) + " | " + samples["condition"].astype(str)
condition_activity = activity.T.groupby(by=condition_key.values).mean()   # rows = conditions
condition_project = (samples.assign(condition_key=condition_key.values)
                     .drop_duplicates(subset="condition_key")
                     .set_index(keys="condition_key")
                     .loc[condition_activity.index, "project_name"])

random_generator = np.random.default_rng(seed=RANDOM_SEED)
unique_projects = condition_project.unique()

robustness_rows = []
for first, second in KEY_PAIRS:
    sample_r = activity.loc[first].corr(other=activity.loc[second], method="pearson")
    condition_r = condition_activity[first].corr(other=condition_activity[second], method="pearson")
    condition_spearman = condition_activity[first].corr(other=condition_activity[second], method="spearman")

    bootstrap_values = []
    for _ in range(N_BOOTSTRAP):
        chosen_projects = random_generator.choice(a=unique_projects, size=len(unique_projects), replace=True)
        # Take every condition from each sampled project (a project drawn twice appears twice).
        row_positions = np.concatenate([np.where(condition_project.values == project)[0]
                                        for project in chosen_projects])
        resampled = condition_activity.iloc[row_positions]
        bootstrap_values.append(resampled[first].corr(other=resampled[second]))
    lower, upper = np.nanpercentile(a=bootstrap_values, q=[2.5, 97.5])

    # (c) Within-project correlations: does the relationship hold INSIDE
    # individual studies, or is it driven by differences between studies
    # (different strains, labs, media)? Only projects with >= 6 samples and
    # real variation in both iModulons are used.
    within_project_r = []
    for project_name, project_samples in samples.groupby(by="project_name"):
        sample_ids = project_samples.index
        first_values = activity.loc[first, sample_ids]
        second_values = activity.loc[second, sample_ids]
        if len(sample_ids) >= 6 and first_values.std() > 0.5 and second_values.std() > 0.5:
            within_project_r.append(first_values.corr(other=second_values))
    within_project_r = np.array(within_project_r)

    robustness_rows.append({
        "pair": f"{label(first)} vs {label(second)}",
        "sample_level_pearson_r": round(sample_r, 3),
        "condition_level_pearson_r": round(condition_r, 3),
        "condition_level_spearman_rho": round(condition_spearman, 3),
        "project_bootstrap_95ci_low": round(lower, 3),
        "project_bootstrap_95ci_high": round(upper, 3),
        "within_project_median_r": round(float(np.median(within_project_r)), 3),
        "within_project_fraction_same_sign": round(float(np.mean(np.sign(within_project_r) == np.sign(condition_r))), 2),
        "n_projects_within": len(within_project_r),
        "n_samples": activity.shape[1],
        "n_conditions": len(condition_activity),
    })
    report(f"{label(first)} vs {label(second)}: sample r = {sample_r:.3f}, "
           f"condition r = {condition_r:.3f} (95% project-bootstrap CI {lower:.2f} to {upper:.2f}), "
           f"Spearman = {condition_spearman:.3f}; within-project median r = {np.median(within_project_r):.2f} "
           f"(same sign in {np.mean(np.sign(within_project_r) == np.sign(condition_r)):.0%} of {len(within_project_r)} projects)")

robustness_table = pd.DataFrame(data=robustness_rows)
robustness_table.to_csv(path_or_buf=RESULTS_DIR / "correlation_robustness.csv", index=False)

# Figure 2: scatter plots at the condition level.
figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(11, 4.6))
for axis, (first, second), row in zip(axes, KEY_PAIRS, robustness_rows):
    axis.scatter(x=condition_activity[first], y=condition_activity[second],
                 s=10, alpha=0.5, color="#2E6FBA", edgecolors="none")
    # Least-squares line: np.polyfit returns [slope, intercept] for degree 1.
    slope, intercept = np.polyfit(x=condition_activity[first], y=condition_activity[second], deg=1)
    x_line = np.linspace(start=condition_activity[first].min(), stop=condition_activity[first].max(), num=50)
    axis.plot(x_line, slope * x_line + intercept, color="#D1495B", linewidth=1.5)
    axis.axhline(y=0, color="grey", linewidth=0.5)
    axis.axvline(x=0, color="grey", linewidth=0.5)
    axis.set_xlabel(f"{label(first)} activity")
    axis.set_ylabel(f"{label(second)} activity")
    axis.set_title(f"r = {row['condition_level_pearson_r']:.2f} "
                   f"[{row['project_bootstrap_95ci_low']:.2f}, {row['project_bootstrap_95ci_high']:.2f}]"
                   f"  ({row['n_conditions']} conditions)")
figure.suptitle("Key iModulon relationships (replicates averaged; CI = project-level bootstrap)", y=1.02)
figure.tight_layout()
figure.savefig(fname=FIGURE_DIR / "02_key_correlations.png", bbox_inches="tight")
plt.close(fig=figure)

# ---------------------------------------------------------------------------
# 4. Correlation structure among the strongest iModulons
# ---------------------------------------------------------------------------
# Use the 25 multi-gene iModulons with the largest explained variance.
top_multi_gene = (imodulon_table.loc[~is_single_gene]
                  .sort_values(by="exp_var", ascending=False)
                  .head(n=25).index.tolist())
for name in CURATED_LABELS:
    if name not in top_multi_gene and name in activity.index:
        top_multi_gene.append(name)

correlation_matrix = condition_activity[top_multi_gene].corr(method="pearson")

# Order rows/columns by hierarchical clustering so related modules sit together.
# Distance = 1 - |r|, so strongly positive OR negative pairs end up close.
# .to_numpy(copy=True) makes a writable copy we can edit in place.
distance_array = (1 - correlation_matrix.abs()).to_numpy(copy=True)
np.fill_diagonal(a=distance_array, val=0)
cluster_tree = linkage(y=squareform(X=distance_array, checks=False), method="average")
order = [correlation_matrix.index[i] for i in leaves_list(Z=cluster_tree)]
ordered = correlation_matrix.loc[order, order]

figure, axis = plt.subplots(figsize=(9, 8))
image = axis.imshow(X=ordered.values, cmap="RdBu_r", vmin=-1, vmax=1)
axis.set_xticks(ticks=range(len(order)), labels=[label(name) for name in order], rotation=90, fontsize=7)
axis.set_yticks(ticks=range(len(order)), labels=[label(name) for name in order], fontsize=7)
figure.colorbar(mappable=image, ax=axis, shrink=0.7, label="Pearson r (condition level)")
axis.set_title("Correlation of iModulon activities across conditions")
figure.tight_layout()
figure.savefig(fname=FIGURE_DIR / "03_correlation_heatmap.png", bbox_inches="tight")
plt.close(fig=figure)

# ---------------------------------------------------------------------------
# 5. PCA and k-means on condition-level activity profiles
# ---------------------------------------------------------------------------
pca = PCA(n_components=10, random_state=RANDOM_SEED)
pca_scores = pca.fit_transform(X=condition_activity.values)
report(f"PCA: PC1 {pca.explained_variance_ratio_[0]:.1%}, PC2 {pca.explained_variance_ratio_[1]:.1%}, "
       f"first 10 PCs {pca.explained_variance_ratio_.sum():.1%}")

# Choose k by silhouette score instead of fixing it by hand.
silhouette_by_k = {}
for number_of_clusters in range(2, 9):
    kmeans_trial = KMeans(n_clusters=number_of_clusters, n_init=20, random_state=RANDOM_SEED)
    trial_labels = kmeans_trial.fit_predict(X=pca_scores)
    silhouette_by_k[number_of_clusters] = silhouette_score(X=pca_scores, labels=trial_labels)
best_k = max(silhouette_by_k, key=silhouette_by_k.get)
report("Silhouette by k: " + ", ".join(f"k={k}: {v:.2f}" for k, v in silhouette_by_k.items())
       + f" -> chosen k = {best_k}")

kmeans = KMeans(n_clusters=best_k, n_init=20, random_state=RANDOM_SEED)
cluster_labels = kmeans.fit_predict(X=pca_scores)

# Which projects dominate each cluster?
cluster_summary = (pd.DataFrame(data={"cluster": cluster_labels, "project": condition_project.values})
                   .groupby(by="cluster")["project"]
                   .agg(func=lambda projects: ", ".join(projects.value_counts().head(n=3).index)))
cluster_sizes = pd.Series(data=cluster_labels).value_counts().sort_index()
for cluster_id in cluster_summary.index:
    report(f"  Cluster {cluster_id} ({cluster_sizes[cluster_id]} conditions): top projects = {cluster_summary[cluster_id]}")

figure, (axis_left, axis_right) = plt.subplots(nrows=1, ncols=2, figsize=(11, 4.6))
cluster_palette = ["#2E6FBA", "#D1495B", "#5BA45B", "#E8A33D", "#7B4FB0", "#3AAFA9", "#8C6D46", "#E377C2"]
for cluster_id in range(best_k):
    in_cluster = cluster_labels == cluster_id
    axis_left.scatter(x=pca_scores[in_cluster, 0], y=pca_scores[in_cluster, 1],
                      s=12, alpha=0.7, color=cluster_palette[cluster_id], label=f"Cluster {cluster_id}",
                      edgecolors="none")
axis_left.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
axis_left.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
axis_left.set_title(f"Conditions in iModulon-activity space (k-means, k = {best_k})")
axis_left.legend(frameon=False, fontsize=8)

# Right: which iModulons drive PC1 and PC2 (largest absolute loadings).
loadings = pd.DataFrame(data=pca.components_[:2].T, index=condition_activity.columns, columns=["PC1", "PC2"])
top_loading_names = loadings.abs().max(axis=1).sort_values(ascending=False).head(n=12).index
loadings.loc[top_loading_names].rename(index=label).plot.barh(ax=axis_right, color=["#2E6FBA", "#E8A33D"])
axis_right.axvline(x=0, color="grey", linewidth=0.5)
axis_right.set_title("iModulons with the largest PC1/PC2 loadings")
axis_right.set_xlabel("Loading")
axis_right.tick_params(axis="y", labelsize=8)
axis_right.legend(frameon=False)

figure.tight_layout()
figure.savefig(fname=FIGURE_DIR / "04_pca_kmeans.png", bbox_inches="tight")
plt.close(fig=figure)

# ---------------------------------------------------------------------------
# 6. Where is arsenic resistance active?
# ---------------------------------------------------------------------------
arsenic = condition_activity["IM_104"].sort_values()
extremes = pd.concat(objs=[arsenic.head(n=8), arsenic.tail(n=8)])
figure, axis = plt.subplots(figsize=(9, 5))
axis.barh(y=[name[:70] for name in extremes.index], width=extremes.values,
          color=["#2E6FBA"] * 8 + ["#D1495B"] * 8)
axis.axvline(x=0, color="grey", linewidth=0.5)
axis.set_xlabel("IM_104 (Arsenic resistance) activity")
axis.set_title("Conditions with the lowest and highest arsenic iModulon activity")
axis.tick_params(axis="y", labelsize=7)
figure.tight_layout()
figure.savefig(fname=FIGURE_DIR / "05_arsenic_conditions.png", bbox_inches="tight")
plt.close(fig=figure)

(RESULTS_DIR / "summary.txt").write_text(data="\n".join(summary_lines) + "\n")
print("\nDone. Figures in figures/, numbers in results/.")
