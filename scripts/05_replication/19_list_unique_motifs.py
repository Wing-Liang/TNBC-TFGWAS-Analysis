# ==============================================================================
# 19_list_unique_motifs.py
# List subtype-unique motifs in both cohorts.
# Produces: Table S8
# Run     : python3 scripts/05_replication/19_list_unique_motifs.py   (from the repository root)
# ==============================================================================
"""19_list_unique_motifs.py — name the UNIQUE significant motifs shared by both cohorts.

Motif model (727 tests). Uniqueness is computed inside each cohort, then intersected:
a motif is listed only if it is unique in discovery AND unique in replication.

  python3 19_list_unique_motifs.py [RESPONSE] [MODE]
      RESPONSE  TNBC | LuminalA | LuminalB | LuminalB_HER2Neg | HER2   (default TNBC)
      MODE      sub     = unique vs the other four subtypes  (default)
                overall = unique vs overall BC

Locations from GROUP_REF / OBS_RESULTS; defaults assume the project layout under $TNBC_ROOT.
"""
import csv, io, os, sys

REF = os.environ.get("GROUP_REF",   "reference")
V   = os.environ.get("OBS_RESULTS", "replication/tfgwas_grouped")
RESP = sys.argv[1] if len(sys.argv) > 1 else "TNBC"
MODE = sys.argv[2] if len(sys.argv) > 2 else "sub"

def rd(fn):
    with io.open(os.path.join(REF, fn), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))
cl  = {r["motif_id"]: r["group"] for r in rd("group_JASPAR122.csv")}
fam = {r["motif_id"]: r["group"] for r in rd("group_TFfamily.csv")}
tf  = {r["motif_id"]: r["tf"] for r in rd("jaspar_tf_family.csv")} \
      if os.path.exists(os.path.join(REF, "jaspar_tf_family.csv")) else {}

SUB = ["TNBC","LuminalA","LuminalB","LuminalB_HER2Neg","HER2"]
def load(c, k):
    p = os.path.join(V, f"obs_{c}_motif_{k}.csv")
    if not os.path.exists(p): sys.exit(f"missing: {p}")
    with io.open(p, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    return ({r["group"] for r in rows if r["bonferroni_significant"].upper()=="TRUE"},
            {r["group"]: float(r["t_observed"]) for r in rows if r["t_observed"] not in ("","NA")})

S, T = {}, {}
for c in ("discovery","replication"):
    for k in SUB + ["OverallBC"]:
        S[(c,k)], T[(c,k)] = load(c, k)

def uniq(c):
    if MODE == "overall":
        return S[(c,RESP)] - S[(c,"OverallBC")]
    if MODE == "all":                      # every significant motif, no uniqueness filter
        return S[(c,RESP)]
    return S[(c,RESP)] - set().union(*[S[(c,x)] for x in SUB if x != RESP])

d, r = uniq("discovery"), uniq("replication")
both = sorted(d & r, key=lambda m: -T[("discovery",RESP)].get(m, 0))
label = {"sub":"unique vs other subtypes","overall":"unique vs overall BC","all":"all significant motifs"}[MODE]
print(f"{RESP} — motif model, {label}")
print(f"discovery {len(d)}, replication {len(r)}, shared {len(both)}\n")
print(f"{'motif':<11}{'TF':<15}{'t_disc':>8}{'t_repl':>8}  {'cluster':<12}{'TF family'}")
print("-"*95)
for m in both:
    print(f"{m:<11}{tf.get(m,''):<15}{T[('discovery',RESP)].get(m,float('nan')):>+8.2f}"
          f"{T[('replication',RESP)].get(m,float('nan')):>+8.2f}  {cl.get(m,''):<12}{fam.get(m,'')[:40]}")
