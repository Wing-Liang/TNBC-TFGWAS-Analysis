#!/usr/bin/env python3
# ==============================================================================
# 03b_extract_cellxgene_embeddings.py
# Extract embeddings and cell metadata from the CELLxGENE object.
# Produces: replication input
# Run     : python3 scripts/05_replication/03b_extract_cellxgene_embeddings.py   (from the repository root)
# ==============================================================================
"""
03b_extract_cellxgene_embeddings.py -- pull the AUTHORS' embeddings out of the CELLxGENE h5ad.

WHY
03a_extract_cellxgene_data.py took only the raw counts and obs$author_cell_type. The h5ad also
carries the embeddings the authors published:

    obsm: X_pca (n x 100)   X_lsi (n x 100)
          X_umap.rna (n x 2)  X_umap.atac (n x 2)  X_wnn.umap (n x 2)
    obs : ATAC.weight   (the WNN modality weight)

Importing them means 05b/05c/05d can plot the EXACT coordinates behind the source
publication's figures instead of re-deriving an embedding that will differ in orientation
and cluster shape. Same reasoning that led to using the authors' cell-type labels rather
than re-clustering.

TWO THINGS TO KNOW
1. The embeddings were computed on ALL 51,367 cells across 16 pools. Restricting to the
   three European pools yields valid coordinates, but the layout reflects the full atlas --
   cells sit where they sat among all donors. State this in Methods; it is not the same as
   embedding this subset independently.
2. MATCHING IS ON (barcode, Pool), NOT barcode alone. The h5ad index is
   'AAACAGCCAATTGAAG-1_1' -- barcode plus a per-pool integer suffix that is NOT our
   pool prefix. 88 barcodes recur across pools 9/11/12, so a barcode-only join would
   silently mis-assign them.

Verified before writing this: the h5ad barcode set per pool matches
data/cellxgene/pool{9,11,12}_cells.txt exactly (4417 / 3408 / 4644, zero discrepancy).

INPUT   data/cellxgene/breast_atlas.h5ad
OUTPUT  data/cellxgene/embeddings_replication.csv   cell, pool, umap.rna_1/2, umap.atac_1/2,
                                            wnn.umap_1/2, ATAC.weight
        data/cellxgene/obs_replication.csv          cell + all 69 obs metadata columns
        data/cellxgene/var_replication.csv          gene + feature_name/type/biotype/length
        data/cellxgene/pca_replication.csv          cell + 100 PCs
        data/cellxgene/lsi_replication.csv          cell + 100 LSI components
"""
import h5py
import numpy as np
import pandas as pd
import os

H5   = "data/cellxgene/breast_atlas.h5ad"
OUT  = "cellxgene"
POOLS = {"Pool9": "pool9", "Pool11": "pool11", "Pool12": "pool12"}

def decode(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])

def cat_column(obs, name):
    """Read an obs column that may be stored categorically or as a plain array."""
    col = obs[name]
    if isinstance(col, h5py.Group):                    # categorical
        cats = decode(col["categories"][:])
        return cats[col["codes"][:]]
    return decode(col[:]) if col.dtype.kind in "OS" else col[:]

def main():
    if not os.path.exists(H5):
        raise SystemExit(f"Not found: {H5}")
    f = h5py.File(H5, "r")
    obs = f["obs"]

    index = decode(obs["_index"][:])
    pool  = cat_column(obs, "Pool")
    print(f"h5ad: {len(index)} cells, {len(set(pool))} pools")

    keep = np.isin(pool, list(POOLS.keys()))
    print(f"European pools: {keep.sum()} cells")

    # barcode is everything before the LAST underscore; the suffix is a per-pool integer
    barcode = np.array([x.rsplit("_", 1)[0] for x in index])
    cell = np.array([f"{POOLS[p]}_{b}" for p, b in zip(pool[keep], barcode[keep])])
    if len(set(cell)) != len(cell):
        raise SystemExit("Pool-prefixed cell names are not unique -- refusing to write.")

    out = pd.DataFrame({"cell": cell, "pool": [POOLS[p] for p in pool[keep]]})

    # ---- 2D embeddings ----
    for key, stem in [("X_umap.rna", "umap.rna"),
                      ("X_umap.atac", "umap.atac"),
                      ("X_wnn.umap", "wnn.umap")]:
        if key in f["obsm"]:
            m = f["obsm"][key][:][keep]
            out[f"{stem}_1"] = m[:, 0]
            out[f"{stem}_2"] = m[:, 1]
            print(f"  {key}: {m.shape}")
        else:
            print(f"  {key}: ABSENT")

    # ---- WNN modality weight ----
    if "ATAC.weight" in obs:
        out["ATAC.weight"] = np.asarray(cat_column(obs, "ATAC.weight"), dtype=float)[keep]
        print("  ATAC.weight: ok")

    out.to_csv(f"{OUT}/embeddings_replication.csv", index=False)
    print(f"wrote {OUT}/embeddings_replication.csv  ({out.shape[0]} x {out.shape[1]})")

    # ---- FULL obs metadata (69 columns) ----
    # Everything the authors recorded per cell. Written separately from the embeddings so
    # the embedding file stays small and readable.
    #
    # WARNING carried into 05b: many of these names COLLIDE with metrics this pipeline
    # computes itself -- nCount_RNA, nCount_ATAC, nFeature_*, percent.mt, TSS.enrichment,
    # nucleosome_signal, high.tss, nucleosome_group. 05b prefixes every imported column
    # with `cxg_` so the authors' values sit ALONGSIDE ours rather than overwriting them.
    meta = pd.DataFrame({"cell": cell})
    skipped = []
    for k in sorted(obs.keys()):
        if k == "_index":
            continue
        try:
            v = cat_column(obs, k)
            v = np.asarray(v)
            if v.ndim != 1 or len(v) != len(keep):
                skipped.append(k); continue
            meta[k] = v[keep]
        except Exception as e:
            skipped.append(f"{k} ({e})")
    meta.to_csv(f"{OUT}/obs_replication.csv", index=False)
    print(f"wrote {OUT}/obs_replication.csv  ({meta.shape[0]} x {meta.shape[1]})")
    if skipped:
        print("  skipped columns:", ", ".join(skipped))

    # donor-level fields in this deposit are POOL-level strings, not per-cell assignments
    for k in ("donor_id", "donor_age", "self_reported_ethnicity", "Group"):
        if k in meta.columns:
            u = meta[k].astype(str).unique()
            print(f"  {k}: {len(u)} distinct value(s) across the 3 European pools")

    # ---- gene (var) annotation ----
    # feature_type carries the GENCODE biotype: protein_coding 19,267 / lncRNA 15,492 /
    # 718 others (mostly IG_* and TR_* segments, pseudogenes, 17 'artifact', 1 TEC).
    # 03a_extract_cellxgene_data.py did not carry this across, so the Seurat RNA assay has all
    # 35,477 genes with no way to filter. Written here so downstream scripts can restrict
    # to protein_coding + lncRNA.
    v = f["var"]
    vg = pd.DataFrame({"gene": decode(v["_index"][:])})
    for k in ("feature_name", "feature_type", "feature_biotype", "feature_length"):
        if k in v:
            try:
                vg[k] = cat_column(v, k)
            except Exception as e:
                print(f"  var/{k}: skipped ({e})")
    vg.to_csv(f"{OUT}/var_replication.csv", index=False)
    print(f"wrote {OUT}/var_replication.csv  ({vg.shape[0]} x {vg.shape[1]})")
    if "feature_type" in vg.columns:
        # NB: do NOT name this `keep` -- that is the boolean CELL mask used below to subset
        # obsm. Shadowing it turned the mask into an int and broke the PCA/LSI export.
        n_keep = int(vg["feature_type"].isin(["protein_coding", "lncRNA"]).sum())
        print(f"  protein_coding + lncRNA: {n_keep} of {len(vg)} "
              f"({100*n_keep/len(vg):.1f}%); {len(vg)-n_keep} other biotypes")

    # ---- high-dimensional reductions (used by 05d for the joint graph) ----
    for key, name in [("X_pca", "pca"), ("X_lsi", "lsi")]:
        if key not in f["obsm"]:
            print(f"  {key}: ABSENT")
            continue
        m = f["obsm"][key][:][keep]
        df = pd.DataFrame(m, columns=[f"{name}_{i+1}" for i in range(m.shape[1])])
        df.insert(0, "cell", cell)
        df.to_csv(f"{OUT}/{name}_replication.csv", index=False)
        print(f"wrote {OUT}/{name}_replication.csv  ({df.shape[0]} x {df.shape[1]})")

    f.close()
    print("=== done ===")

if __name__ == "__main__":
    main()
