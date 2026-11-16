# ==============================================================================
# 18_list_shared_groups.py
# List TF families significant in both cohorts.
# Produces: Table S8
# Run     : python3 scripts/05_replication/18_list_shared_groups.py   (from the repository root)
# ==============================================================================
"""18_list_shared_groups.py — name the groups significant in BOTH cohorts.

Lists, for one response and one group level, the groups that are Bonferroni-
significant in discovery AND in replication under the group model tested
directly, with their member motifs, TF names and TFClass families.

  python3 18_list_shared_groups.py [RESPONSE] [LEVEL]
      RESPONSE  TNBC | LuminalA | LuminalB | LuminalB_HER2Neg | HER2 | OverallBC
      LEVEL     JASPAR122 | TFfamily        (default JASPAR122)

Locations from GROUP_REF / OBS_RESULTS; defaults assume the project layout under $TNBC_ROOT.
"""
import csv, io, os, sys, collections

REF = os.environ.get("GROUP_REF",   "reference")
V   = os.environ.get("OBS_RESULTS", "replication/tfgwas_grouped")
RESP  = sys.argv[1] if len(sys.argv) > 1 else "TNBC"
LEVEL = sys.argv[2] if len(sys.argv) > 2 else "JASPAR122"
MAPF  = {"JASPAR122": "group_JASPAR122.csv", "TFfamily": "group_TFfamily.csv"}[LEVEL]

def rd(fn, d=REF):
    with io.open(os.path.join(d, fn), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))

grp = {r["motif_id"]: r["group"] for r in rd(MAPF)}
fam = {r["motif_id"]: r["group"] for r in rd("group_TFfamily.csv")}
tf  = {r["motif_id"]: r["tf"] for r in rd("jaspar_tf_family.csv")} if \
      os.path.exists(os.path.join(REF, "jaspar_tf_family.csv")) else {}

def sig(cohort):
    p = os.path.join(V, f"obs_{cohort}_{LEVEL}_{RESP}.csv")
    if not os.path.exists(p): sys.exit(f"missing: {p}")
    with io.open(p, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    s = {r["group"] for r in rows if r["bonferroni_significant"].upper() == "TRUE"}
    t = {r["group"]: float(r["t_observed"]) for r in rows if r["t_observed"] not in ("", "NA")}
    return s, t

d, td = sig("discovery"); r, tr = sig("replication")
members = collections.defaultdict(list)
for m, g in grp.items(): members[g].append(m)

shared = sorted(d & r, key=lambda g: -abs(td.get(g, 0)))
print(f"{RESP} — {LEVEL}: discovery {len(d)}, replication {len(r)}, shared {len(shared)}\n")
for g in shared:
    ms = sorted(members[g])
    fams = collections.Counter(fam[m] for m in ms if m in fam)
    names = ", ".join(sorted({tf.get(m, m) for m in ms}))
    print(f"{g}  (n={len(ms)} motifs)   t_disc={td.get(g,float('nan')):+.2f}  t_repl={tr.get(g,float('nan')):+.2f}")
    print(f"   dominant family : {fams.most_common(1)[0][0] if fams else '-'}"
          f"   ({len(fams)} famil{'y' if len(fams)==1 else 'ies'})")
    print(f"   TFs             : {names[:300]}{'...' if len(names)>300 else ''}")
    print()
