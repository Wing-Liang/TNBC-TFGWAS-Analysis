# ==============================================================================
# 14_motif_master_table.py
# Per-motif table: peak and variant counts, significance and subtype-uniqueness flags, both cohorts.
# Produces: Tables S7-S8
# Run     : python3 scripts/05_replication/14_motif_master_table.py   (from the repository root)
# ==============================================================================
"""Per-motif master table: peak/variant counts in both cohorts, plus Bonferroni
significance and subtype-uniqueness flags for every response in both cohorts."""
import csv, io, os

# Locations are configurable; defaults assume the project layout under $TNBC_ROOT so this runs
# server-side without editing. No hardcoded workstation paths.
REF = os.environ.get("GROUP_REF", "reference")
V   = os.environ.get("OBS_RESULTS", "replication/tfgwas_grouped")
SP  = os.environ.get("PEAK_TABLE_DIR", V)

RESP = [("TNBC","TNBC"), ("LuminalA","LumA"), ("LuminalB","LumB"),
        ("LuminalB_HER2Neg","LumB_H2n"), ("HER2","HER2"), ("OverallBC","OverallBC")]
SUBS = [k for k,_ in RESP if k != "OverallBC"]

def sigset(coh, key):
    p = os.path.join(V, f"obs_{coh}_motif_{key}.csv")
    return {r["group"] for r in csv.DictReader(io.open(p, encoding="utf-8"))
            if r["bonferroni_significant"].upper() == "TRUE"}

S = {c: {k: sigset(c,k) for k,_ in RESP} for c in ("discovery","replication")}

# uniqueness: significant in this subtype and in no other subtype; and vs overall BC
U_sub, U_overall = {}, {}
for c in ("discovery","replication"):
    U_sub[c], U_overall[c] = {}, {}
    for k in SUBS:
        U_sub[c][k] = S[c][k] - set().union(*[S[c][x] for x in SUBS if x != k])
        U_overall[c][k] = S[c][k] - S[c]["OverallBC"]

tf  = {r["motif_id"]: r["tf"] for r in csv.DictReader(io.open(os.path.join(REF,"jaspar_tf_family.csv"), encoding="utf-8-sig"))}
base = list(csv.DictReader(io.open(os.path.join(SP,"peaks_per_motif.csv"), encoding="utf-8-sig")))

cols = ["motif","motif_name","jaspar122","tf_family",
        "peaks_discovery","variants_discovery","peaks_replication","variants_replication",
        "peak_ratio","variant_ratio"]
for c,ct in (("discovery","disc"),("replication","repl")):
    cols += [f"sig_{ct}_{lab}" for _,lab in RESP]
for c,ct in (("discovery","disc"),("replication","repl")):
    cols += [f"uniqSub_{ct}_{lab}" for k,lab in RESP if k != "OverallBC"]
for c,ct in (("discovery","disc"),("replication","repl")):
    cols += [f"uniqOverall_{ct}_{lab}" for k,lab in RESP if k != "OverallBC"]

out = []
for r in base:
    m = r["motif"]
    row = {c: r.get(c,"") for c in cols[:10]}
    row["motif"] = m; row["motif_name"] = tf.get(m,"")
    for c,ct in (("discovery","disc"),("replication","repl")):
        for k,lab in RESP:
            row[f"sig_{ct}_{lab}"] = "Y" if m in S[c][k] else ""
        for k,lab in [(k,l) for k,l in RESP if k != "OverallBC"]:
            row[f"uniqSub_{ct}_{lab}"] = "Y" if m in U_sub[c][k] else ""
            row[f"uniqOverall_{ct}_{lab}"] = "Y" if m in U_overall[c][k] else ""
    out.append(row)

p = os.path.join(V, "motif_master_table.csv")
with io.open(p, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(out)
print("wrote", os.path.basename(p), "-", len(out), "rows,", len(cols), "columns")

print("\ncolumn check — flagged Y counts:")
print(f"{'response':<10}{'sig disc':>9}{'sig repl':>9}{'uniqSub disc':>14}{'uniqSub repl':>14}{'uniqOverall disc':>18}{'uniqOverall repl':>18}")
for k,lab in RESP:
    sd = sum(1 for r in out if r[f"sig_disc_{lab}"]=="Y")
    sr = sum(1 for r in out if r[f"sig_repl_{lab}"]=="Y")
    if k != "OverallBC":
        ud = sum(1 for r in out if r[f"uniqSub_disc_{lab}"]=="Y")
        ur = sum(1 for r in out if r[f"uniqSub_repl_{lab}"]=="Y")
        pd = sum(1 for r in out if r[f"uniqOverall_disc_{lab}"]=="Y")
        pr = sum(1 for r in out if r[f"uniqOverall_repl_{lab}"]=="Y")
    else:
        ud=ur=pd=pr=""
    print(f"{lab:<10}{sd:>9}{sr:>9}{str(ud):>14}{str(ur):>14}{str(pd):>18}{str(pr):>18}")
