# ==============================================================================
# 18_cross_cohort_overlap.py
# Overlap of significant units between cohorts.
# Produces: Table S8
# Run     : python3 scripts/05_replication/18_cross_cohort_overlap.py   (from the repository root)
# ==============================================================================
"""18_cross_cohort_overlap.py — discovery vs replication overlap, all three levels.

Motif model (727), per response:
  discovery / replication / overlap of Bonferroni-significant motifs, and the
  cohort-only remainders; the JASPAR122 clusters and TF families those motifs
  span in each cohort and the overlap of those group sets; and the overlap of the
  "unique vs other subtypes" and "unique vs overall BC" sets (i.e. motifs unique in
  discovery AND unique in replication).

Group models (JASPAR122, TF family) tested directly: discovery / replication /
overlap of significant groups.

Uniqueness is computed within each cohort separately, then intersected — a motif
counts only if it is unique in both.

Locations configurable; defaults assume the project layout under $TNBC_ROOT. No hardcoded paths.
  GROUP_REF    default reference
  OBS_RESULTS  default replication/tfgwas_grouped
"""
import csv, io, os, sys

REF = os.environ.get("GROUP_REF",   "reference")
V   = os.environ.get("OBS_RESULTS", "replication/tfgwas_grouped")

def load_map(fn):
    with io.open(os.path.join(REF, fn), encoding="utf-8-sig") as fh:
        return {r["motif_id"]: r["group"] for r in csv.DictReader(fh)}
cl, fam = load_map("group_JASPAR122.csv"), load_map("group_TFfamily.csv")

def sig(cohort, level, key):
    p = os.path.join(V, f"obs_{cohort}_{level}_{key}.csv")
    if not os.path.exists(p): return None
    with io.open(p, encoding="utf-8-sig") as fh:
        return {r["group"] for r in csv.DictReader(fh)
                if r["bonferroni_significant"].upper() == "TRUE"}

RESP = [("OverallBC","overall BC"), ("LuminalA","Luminal A-like"), ("TNBC","TNBC"),
        ("LuminalB","Luminal B-like"), ("LuminalB_HER2Neg","Luminal B/HER2-neg"),
        ("HER2","HER2-enriched")]
SUB  = [k for k, _ in RESP if k != "OverallBC"]
COH  = ("discovery", "replication")

M = {c: {k: sig(c, "motif", k) for k, _ in RESP} for c in COH}
for c in COH:
    if any(v is None for v in M[c].values()):
        sys.exit(f"missing motif-level results for {c} in {V}")

def uniq_sub(c, k):
    return M[c][k] - set().union(*[M[c][x] for x in SUB if x != k]) if k != "OverallBC" \
           else M[c][k] - set().union(*[M[c][x] for x in SUB])
def uniq_overall(c, k):
    return None if k == "OverallBC" else M[c][k] - M[c]["OverallBC"]
span = lambda s, m: {m[x] for x in s if x in m}

print(f"MOTIF MODEL (727 tests, Bonferroni p < {0.05/727:.2e}) — discovery vs replication\n")
h = (f"{'Response':<21}{'Disc':>6}{'Repl':>6}{'Both':>6}"
     f"{'J122_D':>8}{'J122_R':>8}{'J122_ovl':>9}"
     f"{'Fam_D':>7}{'Fam_R':>7}{'Fam_ovl':>8}")
print(h); print("-" * len(h))
for k, lab in RESP:
    d, r = M["discovery"][k], M["replication"][k]
    cd, cr = span(d, cl),  span(r, cl)
    fd, fr = span(d, fam), span(r, fam)
    print(f"{lab:<21}{len(d):>6}{len(r):>6}{len(d&r):>6}"
          f"{len(cd):>8}{len(cr):>8}{len(cd&cr):>9}"
          f"{len(fd):>7}{len(fr):>7}{len(fd&fr):>8}")

# ---- same layout, restricted to the UNIQUE significant motifs ----------------
# Uniqueness is computed inside each cohort, then the cohorts are compared:
#   Disc/Repl = size of that cohort's unique set
#   Both      = motifs unique in BOTH cohorts
#   J122/Fam  = clusters / families those unique motifs span, and their overlap
for tag, fn, rows in (
        ("UNIQUE vs OTHER SUBTYPES", uniq_sub, RESP),
        ("UNIQUE vs OVERALL BC",         uniq_overall, [(k, l) for k, l in RESP if k != "OverallBC"])):
    print("")
    print(f"MOTIF MODEL - {tag} (Bonferroni p < {0.05/727:.2e})")
    print("")
    h3 = (f"{'Response':<21}{'Disc':>6}{'Repl':>6}{'Both':>6}"
          f"{'J122_D':>8}{'J122_R':>8}{'J122_ovl':>9}"
          f"{'Fam_D':>7}{'Fam_R':>7}{'Fam_ovl':>8}")
    print(h3); print("-" * len(h3))
    for k, lab in rows:
        d, r = fn("discovery", k), fn("replication", k)
        cd, cr = span(d, cl),  span(r, cl)
        fd, fr = span(d, fam), span(r, fam)
        print(f"{lab:<21}{len(d):>6}{len(r):>6}{len(d&r):>6}"
              f"{len(cd):>8}{len(cr):>8}{len(cd&cr):>9}"
              f"{len(fd):>7}{len(fr):>7}{len(fd&fr):>8}")

for lvl in ("JASPAR122", "TFfamily"):
    G = {c: {k: sig(c, lvl, k) for k, _ in RESP} for c in COH}
    if any(v is None for c in COH for v in G[c].values()): continue
    n = len(set(cl.values())) if lvl == "JASPAR122" else len(set(fam.values()))
    print(f"\n{lvl} MODEL tested directly ({n} tests, Bonferroni p < {0.05/n:.2e})\n")
    h2 = f"{'Response':<21}{'Disc':>6}{'Repl':>6}{'Overlap':>9}{'DiscOnly':>10}{'ReplOnly':>10}"
    print(h2); print("-" * len(h2))
    for k, lab in RESP:
        d, r = G["discovery"][k], G["replication"][k]
        print(f"{lab:<21}{len(d):>6}{len(r):>6}{len(d&r):>9}{len(d-r):>10}{len(r-d):>10}")

    # uniqueness at the group level: the unit tested IS the group, so there are
    # no spanned-group sub-counts here. Computed inside each cohort, then intersected.
    def g_uniq_sub(c, k):
        return G[c][k] - set().union(*[G[c][x] for x in SUB if x != k]) if k != "OverallBC"                else G[c][k] - set().union(*[G[c][x] for x in SUB])
    def g_uniq_overall(c, k):
        return G[c][k] - G[c]["OverallBC"]
    for tag, fn, rows in (
            ("UNIQUE vs OTHER SUBTYPES", g_uniq_sub, RESP),
            ("UNIQUE vs OVERALL BC",         g_uniq_overall,
             [(a, b) for a, b in RESP if a != "OverallBC"])):
        print("")
        print(f"{lvl} MODEL - {tag}")
        print("")
        h4 = f"{'Response':<21}{'Disc':>6}{'Repl':>6}{'Both':>6}"
        print(h4); print("-" * len(h4))
        for k, lab in rows:
            d, r = fn("discovery", k), fn("replication", k)
            print(f"{lab:<21}{len(d):>6}{len(r):>6}{len(d&r):>6}")
