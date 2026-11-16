#!/usr/bin/env python3
# 03a_extract_cellxgene_data.py
# Extract the CELLxGENE-derived inputs for the replication cohort (pools 9, 11, 12) from the
# source atlas object: per-pool cell barcodes, author cell types and raw RNA counts.
#
# SOURCE  : data/cellxgene/breast_atlas.h5ad
#   downloaded from CELLxGENE dataset 61af564d-e5ea-4d34-a0f3-2668a00db376
#   URL: https://datasets.cellxgene.cziscience.com/31798b24-255e-4726-8805-66507589ccb8.h5ad
#   (Bhat-Nakshatri et al., Nat Med 2024;30:3482-3494)
#
# PRODUCES:
#   data/cellxgene/pool{9,11,12}_cells.txt            per-pool base barcodes (12,469 QC'd cells; used by 03c_combine_fragments.sbatch)
#   celltypes/cellxgene_celltypes_replication.csv   cell,pool,author_cell_type (used by 04_build_combined_object.R)
#   data/cellxgene/rna_replication/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}
#                                                raw RNA counts /raw/X (used by 04_build_combined_object.R)
# Run from the repository root:  python3 scripts/05_replication/03a_extract_cellxgene_data.py
import os, gzip, csv, collections, h5py, numpy as np, scipy.sparse as sp, scipy.io
os.chdir(os.path.join(os.environ.get("TNBC_ROOT", "."), "replication"))
H5="data/cellxgene/breast_atlas.h5ad"
POOLMAP={"Pool9":"pool9","Pool11":"pool11","Pool12":"pool12"}   # European (non-Ashkenazi) pools
f=h5py.File(H5,"r")
def rd(node):
    if isinstance(node,h5py.Group) and "categories" in node:
        cats=[x.decode() if isinstance(x,bytes) else str(x) for x in node["categories"][:]]
        return [cats[c] if c>=0 else "NA" for c in node["codes"][:]]
    return [x.decode() if isinstance(x,bytes) else str(x) for x in node[:]]
Pool=np.array(rd(f["obs"]["Pool"])); act=np.array(rd(f["obs"]["author_cell_type"]))
idx=np.array(rd(f["obs"]["_index"]))
mask=np.isin(Pool,list(POOLMAP)); rows=np.where(mask)[0]
os.makedirs("cellxgene",exist_ok=True); os.makedirs("celltypes",exist_ok=True)
os.makedirs("data/cellxgene/rna_replication",exist_ok=True)

# per-pool base barcodes + cell types
ctrows=[]
per_pool=collections.defaultdict(list)
for i in rows:
    ours=POOLMAP[Pool[i]]; base="_".join(idx[i].split("_")[:-1])   # strip _N suffix
    per_pool[ours].append(base)
    ctrows.append((f"{ours}_{base}", ours, act[i]))
for p,bc in per_pool.items():
    open(f"data/cellxgene/{p}_cells.txt","w").write("\n".join(bc)+"\n")
with open("celltypes/cellxgene_celltypes_replication.csv","w",newline="") as o:
    w=csv.writer(o); w.writerow(["cell","pool","author_cell_type"]); w.writerows(ctrows)
print("cells:", {p:len(v) for p,v in per_pool.items()}, "total", len(ctrows))

# raw RNA counts (/raw/X: CSR cells x genes) -> genes x cells 10x mtx
X=f["raw"]["X"]; shape=tuple(X.attrs["shape"])
M=sp.csr_matrix((X["data"][:],X["indices"][:],X["indptr"][:]),shape=shape)
genes=rd(f["raw"]["var"]["_index"]); sym=rd(f["raw"]["var"]["feature_name"])
cellnames=[r[0] for r in ctrows]
sub=M[rows,:].T.tocsc()
scipy.io.mmwrite("data/cellxgene/rna_replication/matrix.mtx",sub,field="integer")
os.system("gzip -f data/cellxgene/rna_replication/matrix.mtx")
with gzip.open("data/cellxgene/rna_replication/barcodes.tsv.gz","wt") as o: o.write("\n".join(cellnames)+"\n")
with gzip.open("data/cellxgene/rna_replication/features.tsv.gz","wt") as o:
    o.write("\n".join(f"{g}\t{s}\tGene Expression" for g,s in zip(genes,sym))+"\n")
print("RNA genes x cells:", sub.shape, "nnz", sub.nnz)
