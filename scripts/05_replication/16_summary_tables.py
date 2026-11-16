# ==============================================================================
# 16_summary_tables.py
# Summary counts of significant and subtype-unique units per level.
# Produces: Table S8
# Run     : python3 scripts/05_replication/16_summary_tables.py   (from the repository root)
# ==============================================================================
"""16_summary_tables.py — the per-response summary tables for one cohort.

Reproduces, for the motif (727) model:
  (1) Bonferroni-significant motifs
  (2) distinct JASPAR122 clusters those motifs belong to
  (3) distinct TF families those motifs belong to
  (4) motifs unique vs the other four subtypes   (+ clusters/families they span)
  (5) motifs unique vs overall BC                    (+ clusters/families they span)
and the same counts for the JASPAR122 and TF-family models tested directly.

"Unique vs other subtypes" = significant in this response and in none of the other
four subtypes. "Unique vs overall BC" = significant in this subtype but not in overall BC.
For the overall BC row, unique is against the union of all five subtypes and column (5)
is not defined.

Locations are configurable; defaults assume the project layout under $TNBC_ROOT. No hardcoded paths.
  GROUP_REF    group mappings          default reference
  OBS_RESULTS  model output            default replication/tfgwas_grouped
  COHORT       discovery | replication default discovery
"""
import csv, io, os, sys

REF    = os.environ.get("GROUP_REF",   "reference")
V      = os.environ.get("OBS_RESULTS", "replication/tfgwas_grouped")
COHORT = os.environ.get("COHORT", sys.argv[1] if len(sys.argv) > 1 else "discovery")

def load_map(fn):
    with io.open(os.path.join(REF, fn), encoding="utf-8-sig") as fh:
        return {r["motif_id"]: r["group"] for r in csv.DictReader(fh)}

cl  = load_map("group_JASPAR122.csv")
fam = load_map("group_TFfamily.csv")
NC, NF = len(set(cl.values())), len(set(fam.values()))

def sig(level, key):
    p = os.path.join(V, f"obs_{COHORT}_{level}_{key}.csv")
    if not os.path.exists(p):
        return None
    with io.open(p, encoding="utf-8-sig") as fh:
        return {r["group"] for r in csv.DictReader(fh)
                if r["bonferroni_significant"].upper() == "TRUE"}

RESP = [("OverallBC","overall BC"), ("LuminalA","Luminal A-like"), ("TNBC","TNBC"),
        ("LuminalB","Luminal B-like"), ("LuminalB_HER2Neg","Luminal B/HER2-neg"),
        ("HER2","HER2-enriched")]
SUB  = [k for k, _ in RESP if k != "OverallBC"]

M = {k: sig("motif", k) for k, _ in RESP}
if any(v is None for v in M.values()):
    sys.exit(f"missing motif-level results for {COHORT} in {V}")

span = lambda s, m: len({m[x] for x in s if x in m})

print(f"{COHORT.upper()} — motif level (727 tests, Bonferroni p < {0.05/727:.2e})\n")
hdr = f"{'Response':<21}{'Sig':>5}{'Clusters':>10}{'Families':>10}" \
      f"{'Unique vs subtypes':>21}{'Unique vs overall BC':>22}"
print(hdr); print("-" * len(hdr))
for k, lab in RESP:
    s = M[k]
    others = set().union(*[M[x] for x in SUB if x != k]) if k != "OverallBC" \
             else set().union(*[M[x] for x in SUB])
    u1 = s - others
    u2 = None if k == "OverallBC" else s - M["OverallBC"]
    f1 = f"{len(u1)} ({span(u1,cl)}c / {span(u1,fam)}f)"
    f2 = "—" if u2 is None else f"{len(u2)} ({span(u2,cl)}c / {span(u2,fam)}f)"
    print(f"{lab:<21}{len(s):>5}{span(s,cl):>10}{span(s,fam):>10}{f1:>21}{f2:>22}")

print(f"\n(clusters out of {NC}; families out of {NF})")

for lvl, N in (("JASPAR122", NC), ("TFfamily", NF)):
    G = {k: sig(lvl, k) for k, _ in RESP}
    if any(v is None for v in G.values()):
        continue
    print(f"\n{lvl} model tested directly ({N} tests, Bonferroni p < {0.05/N:.2e})\n")
    h2 = f"{'Response':<21}{'Sig':>5}{'% of tests':>12}{'Unique vs subtypes':>21}{'Unique vs overall BC':>22}"
    print(h2); print("-" * len(h2))
    for k, lab in RESP:
        s = G[k]
        others = set().union(*[G[x] for x in SUB if x != k]) if k != "OverallBC" \
                 else set().union(*[G[x] for x in SUB])
        u2 = "—" if k == "OverallBC" else str(len(s - G["OverallBC"]))
        print(f"{lab:<21}{len(s):>5}{100*len(s)/N:>11.0f}%{len(s-others):>21}{u2:>22}")
