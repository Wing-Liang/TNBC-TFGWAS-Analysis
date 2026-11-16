# ==============================================================================
# 21_three_way_tf_overlap.py
# Three-way overlap of TFs across levels and cohorts.
# Produces: Table S8
# Run     : python3 scripts/05_replication/21_three_way_tf_overlap.py   (from the repository root)
# ==============================================================================
"""21_three_way_tf_overlap.py — TF-level overlap of the three shared evidence sets.

For one response, compares three sets at TF-name level:
  A  TFs of the motifs that are SUBTYPE-UNIQUE in both cohorts (motif model)
  B  TFs of every member motif of the JASPAR122 clusters significant in both cohorts
  C  TFs of every member motif of the TF families significant in both cohorts

B and C are full group memberships, not only the significant members.

  python3 21_three_way_tf_overlap.py [RESPONSE]     default TNBC
"""
import csv, io, os, sys, collections

REF = os.environ.get("GROUP_REF",   "reference")
V   = os.environ.get("OBS_RESULTS", "replication/tfgwas_grouped")
RESP = sys.argv[1] if len(sys.argv) > 1 else "TNBC"

def rd(fn):
    with io.open(os.path.join(REF, fn), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))
cl  = {r["motif_id"]: r["group"] for r in rd("group_JASPAR122.csv")}
fam = {r["motif_id"]: r["group"] for r in rd("group_TFfamily.csv")}
tf  = {r["motif_id"]: r["tf"] for r in rd("jaspar_tf_family.csv")}

SUB = ["TNBC","LuminalA","LuminalB","LuminalB_HER2Neg","HER2"]
def sig(c, lvl, k):
    p = os.path.join(V, f"obs_{c}_{lvl}_{k}.csv")
    with io.open(p, encoding="utf-8-sig") as fh:
        return {r["group"] for r in csv.DictReader(fh)
                if r["bonferroni_significant"].upper() == "TRUE"}

# A — subtype-unique motifs shared by both cohorts
def uniq(c):
    s = sig(c, "motif", RESP)
    return s - set().union(*[sig(c, "motif", x) for x in SUB if x != RESP])
A_motifs = uniq("discovery") & uniq("replication")
A = {tf[m] for m in A_motifs if m in tf}

# B — every motif in the clusters significant in both cohorts
Cshared = sig("discovery","JASPAR122",RESP) & sig("replication","JASPAR122",RESP)
B_motifs = {m for m, g in cl.items() if g in Cshared}
B = {tf[m] for m in B_motifs if m in tf}

# C — every motif in the families significant in both cohorts
Fshared = sig("discovery","TFfamily",RESP) & sig("replication","TFfamily",RESP)
C_motifs = {m for m, g in fam.items() if g in Fshared}
C = {tf[m] for m in C_motifs if m in tf}

show = lambda s: ", ".join(sorted(s)) if s else "(none)"
print(f"{RESP} — TF-level overlap of the three shared sets\n")
print(f"  A  subtype-unique motifs, both cohorts   : {len(A_motifs)} motifs, {len(A)} TFs")
print(f"  B  motifs in shared JASPAR122 clusters   : {len(B_motifs)} motifs, {len(B)} TFs "
      f"({len(Cshared)} clusters)")
print(f"  C  motifs in shared TF families          : {len(C_motifs)} motifs, {len(C)} TFs "
      f"({len(Fshared)} families)\n")
print(f"  A n B     ({len(A&B)}): {show(A&B)}")
print(f"  A n C     ({len(A&C)}): {show(A&C)}")
print(f"  B n C     ({len(B&C)}): {show(B&C)}")
print(f"\n  A n B n C ({len(A&B&C)}): {show(A&B&C)}")
print(f"\n  A only    ({len(A-B-C)}): {show(A-B-C)}")

# ---- same three sets, compared at TF-FAMILY level ----------------------------
import collections as _c
fA = {fam[m] for m in A_motifs if m in fam}
_cnt = {g: _c.Counter(fam[m] for m in cl if cl[m] == g and m in fam) for g in Cshared}
fB_dom = {c.most_common(1)[0][0] for c in _cnt.values() if c}
fB_all = {f for c in _cnt.values() for f in c}
fC = set(Fshared)

print("")
print("=" * 78)
print(f"{RESP} - same three sets compared at TF-FAMILY level")
print("=" * 78)
print(f"  A  families of the subtype-unique motifs   : {len(fA)}")
print(f"  B  dominant family of each shared cluster  : {len(fB_dom)}  (all families present: {len(fB_all)})")
print(f"  C  families significant in both cohorts    : {len(fC)}")
print("")
print(f"  A n B     ({len(fA&fB_dom)}): {show(fA&fB_dom)}")
print(f"  A n C     ({len(fA&fC)}): {show(fA&fC)}")
print(f"  B n C     ({len(fB_dom&fC)}): {show(fB_dom&fC)}")
print("")
print(f"  A n B n C ({len(fA&fB_dom&fC)}): {show(fA&fB_dom&fC)}")
print("")
print(f"  A only    ({len(fA-fB_dom-fC)}): {show(fA-fB_dom-fC)}")
print(f"  B only    ({len(fB_dom-fA-fC)}): {show(fB_dom-fA-fC)}")
print(f"  C only    ({len(fC-fA-fB_dom)}): {show(fC-fA-fB_dom)}")
