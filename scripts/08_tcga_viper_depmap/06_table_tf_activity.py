# ==============================================================================
# 06_table_tf_activity.py
# Format the TF-activity table.
# Produces: Table S10
# Run     : python3 scripts/08_tcga_viper_depmap/06_table_tf_activity.py   (from the repository root)
# ==============================================================================
"""06_table_tf_activity.py -- inferred TF activity (msVIPER) of the 61 TNBC-specific TF genes in TCGA-BRCA,
in the two-row-header layout of Table S10.
Source (the final VIPER run used for the expression-vs-activity figure, 05_expression_activity_panel.R):
  results/tf_activity/viper_vst/f10n39_min5/VIPER_VST_msviper_61TF_4contrasts_ARACNe_BRCA.csv
  (DESeq2 VST of 836 female, White TCGA-BRCA samples; ARACNe-BRCA regulons; msVIPER minsize 5)
Gene-expression flag from results/tf_activity/TCGA_DE_vs_VIPER_61TF_concordance.csv (up_all4: DESeq2 FDR < 0.05 and
log2FC > 0 in all four comparisons, FDR within the 61 TF genes -- the 13 genes of Fig. 4A).
Per comparison: NES (positive = higher activity in basal-like), p, FDR (BH within the tested TFs), regulon targets tested,
and the call: Higher / Lower (FDR < 0.05) or blank. TFs not testable (< 5 regulon targets) are listed with "not tested".
Outputs: results/tf_activity/SuppTable_TCGA_TFactivity_VIPER.xlsx (+ .csv, flattened header)
Run from $TNBC_ROOT:  python3 scripts/08_tcga_viper_depmap/06_table_tf_activity.py
"""
import csv
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter
B = "results/tf_activity/"
V = {(r["TF"], r["contrast"]): r for r in csv.DictReader(open(B + "viper_vst/f10n39_min5/VIPER_VST_msviper_61TF_4contrasts_ARACNe_BRCA.csv"))}
C = list(csv.DictReader(open(B + "TCGA_DE_vs_VIPER_61TF_concordance.csv")))
up4 = {r["TF"]: r["up_all4"].upper() == "TRUE" for r in C}
TFS = sorted(up4)                                   # 61 TF genes
CON = [("LumA", "Basal-like vs luminal A"), ("LumB", "Basal-like vs luminal B"),
       ("Her2", "Basal-like vs HER2-enriched"), ("NAT", "Basal-like vs NAT")]
assert len(TFS) == 61
head1 = ["TF gene", "Upregulated gene expression in basal-like tumors (all four comparisons; Fig. 4A)"]
head2 = ["", ""]
for _, lab in CON:
    head1 += [lab, "", "", "", ""]; head2 += ["NES", "p", "FDR", "Regulon targets tested", "Inferred activity (FDR < 0.05)"]
head1 += ["No. of comparisons with higher activity", "Upregulated expression and higher activity"]
head2 += ["", ""]
rows = []
for tf in TFS:
    row = [tf, "Y" if up4[tf] else ""]; nh = 0; tested = False
    for k, _ in CON:
        r = V.get((tf, k))
        if r is None:
            row += ["not tested", "", "", "", ""]; continue
        tested = True
        nes, p, q = float(r["NES"]), float(r["p"]), float(r["FDR_61"])
        call = ("Higher" if nes > 0 else "Lower") if q < 0.05 else ""
        nh += call == "Higher"
        row += [nes, p, q, int(r["tested_size"]), call]
    row += [nh if tested else "not tested", "Y" if (up4[tf] and nh > 0) else ""]
    rows.append(row)
rows.sort(key=lambda x: (x[-2] == "not tested", -(x[-2] if isinstance(x[-2], int) else -1), x[1] != "Y", x[0]))
both = [r[0] for r in rows if r[-1] == "Y"]
print("TFs:", len(rows), "| tested:", sum(r[-2] != "not tested" for r in rows), "| up all4:", sum(r[1] == "Y" for r in rows),
      "| up + higher activity:", both)
OUT = B + "SuppTable_TCGA_TFactivity_VIPER"
with open(OUT + ".csv", "w", newline="") as f:
    w = csv.writer(f); flat, last = [], ""
    for a, b in zip(head1, head2):
        last = a or last; flat.append(f"{last}: {b}" if b else last)
    w.writerow(flat)
    for r in rows: w.writerow([("%.3g" % v) if isinstance(v, float) else v for v in r])
wb = Workbook(); ws = wb.active; ws.title = "TCGA-BRCA TF activity"
ws.append(head1); ws.append(head2)
c = 1
while c <= len(head1):
    span = 1
    while c + span <= len(head1) and head1[c + span - 1] == "" and head2[c + span - 1] != "": span += 1
    if head2[c - 1] == "": ws.merge_cells(start_row=1, start_column=c, end_row=2, end_column=c)
    elif span > 1: ws.merge_cells(start_row=1, start_column=c, end_row=1, end_column=c + span - 1)
    c += span
for row in ws.iter_rows(min_row=1, max_row=2):
    for cell in row: cell.font = Font(bold=True); cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
for r in rows: ws.append(r)
for r in ws.iter_rows(min_row=3):
    for cell in r:
        if isinstance(cell.value, float):
            cell.number_format = "0.00" if ws.cell(row=2, column=cell.column).value == "NES" else "0.00E+00"
for i in range(1, len(head1) + 1): ws.column_dimensions[get_column_letter(i)].width = 12
ws.column_dimensions["B"].width = 22; ws.freeze_panes = "C3"
wb.save(OUT + ".xlsx"); print("wrote", OUT + ".xlsx", "and .csv")
