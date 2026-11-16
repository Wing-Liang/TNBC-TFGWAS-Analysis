# ==============================================================================
# 20_triple_convergence.py
# Motifs significant at motif, family and cluster level.
# Produces: Table S8
# Run     : python3 scripts/05_replication/20_triple_convergence.py   (from the repository root)
# ==============================================================================
"""20_triple_convergence.py — motifs supported at ALL THREE levels in BOTH cohorts.

A motif qualifies only if:
  (1) it is Bonferroni-significant at motif level in discovery AND replication;
  (2) its JASPAR122 cluster is significant in the cluster model in BOTH cohorts;
  (3) its TF family is significant in the family model in BOTH cohorts.

  python3 20_triple_convergence.py [RESPONSE]      default TNBC
"""
import csv, io, os, sys

REF = os.environ.get("GROUP_REF",   "reference")
V   = os.environ.get("OBS_RESULTS", "replication/tfgwas_grouped")
RESP = sys.argv[1] if len(sys.argv) > 1 else "TNBC"

def rd(fn):
    with io.open(os.path.join(REF, fn), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))
cl  = {r["motif_id"]: r["group"] for r in rd("group_JASPAR122.csv")}
fam = {r["motif_id"]: r["group"] for r in rd("group_TFfamily.csv")}
tf  = {r["motif_id"]: r["tf"] for r in rd("jaspar_tf_family.csv")} \
      if os.path.exists(os.path.join(REF, "jaspar_tf_family.csv")) else {}

def load(c, lvl):
    p = os.path.join(V, f"obs_{c}_{lvl}_{RESP}.csv")
    if not os.path.exists(p): sys.exit(f"missing: {p}")
    with io.open(p, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    return ({r["group"] for r in rows if r["bonferroni_significant"].upper()=="TRUE"},
            {r["group"]: float(r["t_observed"]) for r in rows if r["t_observed"] not in ("","NA")})

(md, td), (mr, tr) = load("discovery","motif"),     load("replication","motif")
(cd, _),  (cr, _)  = load("discovery","JASPAR122"), load("replication","JASPAR122")
(fd, _),  (fr, _)  = load("discovery","TFfamily"),  load("replication","TFfamily")

M, C, F = md & mr, cd & cr, fd & fr
step1 = sorted(M)
step2 = [m for m in step1 if cl.get(m) in C]
step3 = [m for m in step2 if fam.get(m) in F]

print(f"{RESP} — convergence across all three levels, both cohorts\n")
print(f"  motifs significant in both cohorts            : {len(M)}")
print(f"  ...whose cluster is significant in both       : {len(step2)}   (of {len(C)} shared clusters)")
print(f"  ...and whose family is significant in both    : {len(step3)}   (of {len(F)} shared families)\n")
step3.sort(key=lambda m: -td.get(m, 0))
print(f"{'motif':<11}{'TF':<15}{'t_disc':>8}{'t_repl':>8}  {'cluster':<12}{'TF family'}")
print("-"*95)
for m in step3:
    print(f"{m:<11}{tf.get(m,''):<15}{td.get(m,float('nan')):>+8.2f}{tr.get(m,float('nan')):>+8.2f}"
          f"  {cl.get(m,''):<12}{fam.get(m,'')[:40]}")
