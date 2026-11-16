# ==============================================================================
# 07_celltype_composition_plot.py
# Plot of cell-type composition, both cohorts.
# Produces: Fig S9
# Run     : python3 scripts/06_dissociation_ap1/07_celltype_composition_plot.py   (from the repository root)
# ==============================================================================
"""Cell-type composition lollipop plot, discovery (GSE162798) vs replication (GSE244585).
One SHARED y axis of harmonised cell types, so cell types present in both cohorts vs only one
are visible at a glance; a cell type absent from a cohort's annotation is marked "not annotated". Sub-clusters that one
cohort splits (LASP-1/-2 in discovery; Endo-1/-2 and Adi-1/-2 in replication) are summed on the shared row, with the
breakdown printed under the count. PNG only.
Counts (post-QC cells used by the TF-GWAS):
  discovery   -- 05_motif_sets_upset_heatmap.R cell table (16,451 cells)
  replication -- replication/tfgwas_grouped/replication_atac_umap.csv (10,700 cells, CELLxGENE author labels)
"""
import sys, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = sys.argv[1]
plt.rcParams.update({"font.family": "Arial", "font.size": 9})

## shared rows, top to bottom: epithelial | stromal | immune ;  (row label, colour)
ROWS = [("LASP", "#E64B35"), ("LHS", "#F39B7F"), ("B-Myo", "#4DBBD5"),
        ("Fb", "#00A087"), ("Endo", "#3C5488"), ("Pc", "#8491B4"), ("Adipocyte", "#B09C85"),
        ("T", "#7E6148"), ("B", "#DC9A00"), ("Myeloid", "#91D1C2")]
GROUPS = [("Epithelial", 0, 2), ("Stromal", 3, 6), ("Immune", 7, 9)]
## cohort -> shared row -> list of (native sub-cluster label, n)
DISC = {"LASP": [("LASP-1", 4277), ("LASP-2", 1781)], "LHS": [("LHS", 2404)], "B-Myo": [("B-Myo", 1726)],
        "Fb": [("Fb", 2566)], "Endo": [("Endo", 1498)], "Pc": [("Pc", 794)], "T": [("T", 897)],
        "B": [("B", 216)], "Myeloid": [("Myl", 292)]}
REPL = {"LASP": [("LASP", 2964)], "LHS": [("LHS", 1350)], "B-Myo": [("BM", 3937)], "Fb": [("Fibroblasts", 906)],
        "Endo": [("Endo-1", 392), ("Endo-2", 392)], "Adipocyte": [("Adi-1", 190), ("Adi-2", 57)],
        "T": [("T-cells", 313)], "Myeloid": [("Macrophages", 199)]}
assert sum(n for v in DISC.values() for _, n in v) == 16451 and sum(n for v in REPL.values() for _, n in v) == 10700

ny = len(ROWS); ypos = {r: ny - 1 - i for i, (r, _) in enumerate(ROWS)}
fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.3), sharey=True)
for ax, (title, D) in zip(axes, [("Main (GSE162798)\n4 donors, 16,451 cells", DISC),
                                  ("Replication (GSE244585)\n15 donors, 10,700 cells", REPL)]):
    tot = sum(n for v in D.values() for _, n in v); xmax = max(sum(n for _, n in v) for v in D.values())
    for _, a, b in GROUPS[1:]:                                   # light separators between lineage groups
        ax.axhline(ypos[ROWS[a][0]] + 0.5, color="grey", lw=0.5, ls=":")
    for row, col in ROWS:
        y = ypos[row]
        if row not in D:
            ax.text(xmax * 0.02, y, "not annotated", va="center", fontsize=7.5, color="grey", style="italic")
            continue
        n = sum(k for _, k in D[row])
        ax.hlines(y, 0, n, color=col, lw=2); ax.plot(n, y, "o", ms=7, color=col, mec="white", mew=0.6)
        sub = D[row]
        lab = f"{n:,} ({100 * n / tot:.1f}%)"
        if len(sub) > 1:                                         # split sub-clusters: show the breakdown
            ax.text(n + xmax * 0.03, y + 0.13, lab, va="center", fontsize=7.5)
            ax.text(n + xmax * 0.03, y - 0.25, " + ".join(f"{s} {k:,}" for s, k in sub), va="center",
                    fontsize=6.3, color="dimgrey")
        else:
            ax.text(n + xmax * 0.03, y, lab, va="center", fontsize=7.5)
    ax.set_xlim(0, xmax * 1.55); ax.set_ylim(-0.6, ny - 0.4)
    ax.set_xlabel("Number of cells"); ax.set_title(title, fontsize=9.5, fontweight="bold")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: f"{int(x):,}"))
axes[0].set_yticks([ypos[r] for r, _ in ROWS]); axes[0].set_yticklabels([r for r, _ in ROWS])
for name, a, b in GROUPS:                                        # lineage-group brackets left of the y labels
    y0, y1 = ypos[ROWS[b][0]] - 0.3, ypos[ROWS[a][0]] + 0.3
    axes[0].annotate("", xy=(-0.30, y0), xytext=(-0.30, y1), xycoords=("axes fraction", "data"),
                     arrowprops=dict(arrowstyle="-", color="grey", lw=0.8), annotation_clip=False)
    axes[0].text(-0.34, (y0 + y1) / 2, name, transform=axes[0].get_yaxis_transform(), rotation=90,
                 ha="center", va="center", fontsize=8, color="dimgrey")
fig.tight_layout()
fig.savefig(f"{OUT}/CellType_distribution_lollipop.png", dpi=300, bbox_inches="tight")
print("wrote", f"{OUT}/CellType_distribution_lollipop.png")
