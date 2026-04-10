"""kinship_qc.py — Compare KING kinship results against UKB reference (ukb_rel.dat).

Produces kinship and IBS0 comparison summaries + diagnostic plots.
Runs inside a Swiss Army Knife job on DNAnexus.

Outputs (written to working directory for SAK auto-upload):
  - kinship_comparison_summary.txt
  - kinship_comparison_plots.png
  - kinship_bland_altman.png
  - ibs0_comparison_summary.txt
  - ibs0_comparison_plots.png
  - ibs0_bland_altman.png
"""

import math
import os
import sys

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

print("Loading ukb_rel.dat ...")
ukb = {}
with open("/mnt/project/Bulk/Genotype Results/Genotype calls/ukb_rel.dat") as f:
    f.readline()  # header: ID1 ID2 HetHet IBS0 Kinship
    for line in f:
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        id1, id2 = parts[0], parts[1]
        key = (min(id1, id2), max(id1, id2))
        ukb[key] = {
            "HetHet": float(parts[2]),
            "IBS0": float(parts[3]),
            "Kinship": float(parts[4]),
        }
print(f"  UKB pairs: {len(ukb):,}")

DX_OUTPUT_DIR = os.environ.get("DX_OUTPUT_DIR", "/sbayesrc_genotypes")
_BASE = f"/mnt/project{DX_OUTPUT_DIR}"

print("Loading ukb_all_direct_rel.kin0 ...")
wgs = {}
with open(f"{_BASE}/kinship/ukb_all_direct_rel.kin0") as f:
    f.readline()  # header: #FID1 IID1 FID2 IID2 NSNP HETHET IBS0 KINSHIP
    for line in f:
        parts = line.strip().split("\t")
        iid1, iid2 = parts[1], parts[3]
        key = (min(iid1, iid2), max(iid1, iid2))
        wgs[key] = {
            "HetHet": float(parts[5]),
            "IBS0": float(parts[6]),
            "Kinship": float(parts[7]),
        }
print(f"  WGS pairs: {len(wgs):,}")

overlap_keys = sorted(set(ukb.keys()) & set(wgs.keys()))
n = len(overlap_keys)
print(f"  Overlapping pairs: {n:,}")
print(f"  WGS-only pairs: {len(wgs) - n:,}")
print(f"  UKB-only pairs: {len(ukb) - n:,}")

if n == 0:
    print("ERROR: no overlapping pairs")
    sys.exit(1)

# Build combined diff records
diffs = []
for key in overlap_keys:
    diffs.append({
        "id1": key[0], "id2": key[1],
        "wgs_kinship": wgs[key]["Kinship"], "ukb_kinship": ukb[key]["Kinship"],
        "kin_diff": wgs[key]["Kinship"] - ukb[key]["Kinship"],
        "kin_abs_diff": abs(wgs[key]["Kinship"] - ukb[key]["Kinship"]),
        "wgs_ibs0": wgs[key]["IBS0"], "ukb_ibs0": ukb[key]["IBS0"],
        "ibs0_diff": wgs[key]["IBS0"] - ukb[key]["IBS0"],
        "ibs0_abs_diff": abs(wgs[key]["IBS0"] - ukb[key]["IBS0"]),
        "wgs_hethet": wgs[key]["HetHet"], "ukb_hethet": ukb[key]["HetHet"],
    })

# Relationship bins (using UKB kinship)
bins = [
    ("MZ twins / duplicates", 0.354, 1.0),
    ("1st degree", 0.177, 0.354),
    ("2nd degree", 0.0884, 0.177),
    ("3rd degree", 0.0442, 0.0884),
    ("Below 3rd degree", 0.03, 0.0442),
]


def pearson(xs, ys):
    """Pearson correlation coefficient."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in xs) / n)
    sy = math.sqrt(sum((b - my) ** 2 for b in ys) / n)
    return cov / (sx * sy) if sx > 0 and sy > 0 else float("nan")


def percentile(sorted_vals, p):
    """Return the p-th percentile from a sorted list."""
    return sorted_vals[int(len(sorted_vals) * p)]


# ---------------------------------------------------------------------------
# Kinship comparison
# ---------------------------------------------------------------------------
print("\n=== Kinship Comparison ===")

kin_wgs = [d["wgs_kinship"] for d in diffs]
kin_ukb = [d["ukb_kinship"] for d in diffs]
kin_signed = [d["kin_diff"] for d in diffs]
kin_abs = sorted([d["kin_abs_diff"] for d in diffs])
kin_corr = pearson(kin_wgs, kin_ukb)
kin_mean_signed = sum(kin_signed) / n
kin_mean_abs = sum(kin_abs) / n

with open("kinship_comparison_summary.txt", "w") as out:
    out.write("=== Kinship Comparison: WGS vs UKB ===\n\n")
    out.write(f"WGS total pairs:   {len(wgs):,}\n")
    out.write(f"UKB total pairs:   {len(ukb):,}\n")
    out.write(f"Overlapping pairs: {n:,}\n")
    out.write(f"WGS-only pairs:    {len(wgs) - n:,}\n")
    out.write(f"UKB-only pairs:    {len(ukb) - n:,}\n\n")

    out.write("--- Kinship Difference (WGS - UKB) ---\n")
    out.write(f"Mean absolute difference:   {kin_mean_abs:.6f}\n")
    out.write(f"Median absolute difference: {kin_abs[n // 2]:.6f}\n")
    out.write(f"90th percentile abs diff:   {percentile(kin_abs, 0.9):.6f}\n")
    out.write(f"95th percentile abs diff:   {percentile(kin_abs, 0.95):.6f}\n")
    out.write(f"99th percentile abs diff:   {percentile(kin_abs, 0.99):.6f}\n")
    out.write(f"Max absolute difference:    {max(kin_abs):.6f}\n")
    out.write(f"Mean signed difference:     {kin_mean_signed:.6f}\n")
    out.write(f"Pearson correlation:         {kin_corr:.6f}\n\n")

    out.write("--- Pair Counts by Kinship Bin (using UKB kinship) ---\n")
    out.write(f"{'Relationship':<25} {'UKB kinship range':<20} {'Both':<10} {'WGS-only':<12} {'UKB-only':<10}\n")
    for label, lo, hi in bins:
        both = sum(1 for d in diffs if lo <= d["ukb_kinship"] < hi)
        wgs_only_bin = sum(1 for k, v in wgs.items() if k not in ukb and lo <= v["Kinship"] < hi)
        ukb_only_bin = sum(1 for k, v in ukb.items() if k not in wgs and lo <= v["Kinship"] < hi)
        out.write(f"{label:<25} [{lo:.4f}, {hi:.4f})    {both:<10} {wgs_only_bin:<12} {ukb_only_bin:<10}\n")
    out.write("\n")

    out.write("--- Top 50 Pairs by Absolute Kinship Discrepancy ---\n")
    out.write(f"{'ID1':<12} {'ID2':<12} {'WGS_Kin':>10} {'UKB_Kin':>10} {'Diff':>10} {'WGS_HetHet':>12} {'UKB_HetHet':>12} {'WGS_IBS0':>10} {'UKB_IBS0':>10}\n")
    top50 = sorted(diffs, key=lambda x: -x["kin_abs_diff"])[:50]
    for d in top50:
        out.write(
            f"{d['id1']:<12} {d['id2']:<12} {d['wgs_kinship']:>10.4f} {d['ukb_kinship']:>10.4f} "
            f"{d['kin_diff']:>10.4f} {d['wgs_hethet']:>12.4f} {d['ukb_hethet']:>12.4f} "
            f"{d['wgs_ibs0']:>10.4f} {d['ukb_ibs0']:>10.4f}\n"
        )

print("  kinship_comparison_summary.txt written.")

# ---------------------------------------------------------------------------
# IBS0 comparison
# ---------------------------------------------------------------------------
print("\n=== IBS0 Comparison ===")

ibs0_wgs = [d["wgs_ibs0"] for d in diffs]
ibs0_ukb = [d["ukb_ibs0"] for d in diffs]
ibs0_signed = [d["ibs0_diff"] for d in diffs]
ibs0_abs = sorted([d["ibs0_abs_diff"] for d in diffs])
ibs0_corr = pearson(ibs0_wgs, ibs0_ukb)
ibs0_mean_signed = sum(ibs0_signed) / n
ibs0_mean_abs = sum(ibs0_abs) / n

with open("ibs0_comparison_summary.txt", "w") as out:
    out.write("=== IBS0 Comparison: WGS vs UKB ===\n\n")
    out.write(f"Overlapping pairs: {n:,}\n\n")

    out.write("--- IBS0 Difference (WGS - UKB) ---\n")
    out.write(f"Mean absolute difference:   {ibs0_mean_abs:.6f}\n")
    out.write(f"Median absolute difference: {ibs0_abs[n // 2]:.6f}\n")
    out.write(f"90th percentile abs diff:   {percentile(ibs0_abs, 0.9):.6f}\n")
    out.write(f"95th percentile abs diff:   {percentile(ibs0_abs, 0.95):.6f}\n")
    out.write(f"99th percentile abs diff:   {percentile(ibs0_abs, 0.99):.6f}\n")
    out.write(f"Max absolute difference:    {max(ibs0_abs):.6f}\n")
    out.write(f"Mean signed difference:     {ibs0_mean_signed:.6f}\n")
    out.write(f"Pearson correlation:         {ibs0_corr:.6f}\n\n")

    out.write("--- Mean IBS0 Difference by Kinship Bin ---\n")
    out.write(f"{'Relationship':<25} {'Count':<10} {'Mean WGS IBS0':<16} {'Mean UKB IBS0':<16} {'Mean Diff':<12} {'Mean Abs Diff':<14}\n")
    for label, lo, hi in bins:
        bin_diffs = [d for d in diffs if lo <= d["ukb_kinship"] < hi]
        if bin_diffs:
            cnt = len(bin_diffs)
            mk = sum(d["wgs_ibs0"] for d in bin_diffs) / cnt
            mu = sum(d["ukb_ibs0"] for d in bin_diffs) / cnt
            md = sum(d["ibs0_diff"] for d in bin_diffs) / cnt
            mad = sum(d["ibs0_abs_diff"] for d in bin_diffs) / cnt
            out.write(f"{label:<25} {cnt:<10} {mk:<16.6f} {mu:<16.6f} {md:<12.6f} {mad:<14.6f}\n")
        else:
            out.write(f"{label:<25} 0\n")
    out.write("\n")

    out.write("--- Top 50 Pairs by Absolute IBS0 Discrepancy ---\n")
    out.write(f"{'ID1':<12} {'ID2':<12} {'WGS_IBS0':>10} {'UKB_IBS0':>10} {'Diff':>10} {'WGS_Kin':>10} {'UKB_Kin':>10}\n")
    top50 = sorted(diffs, key=lambda x: -x["ibs0_abs_diff"])[:50]
    for d in top50:
        out.write(
            f"{d['id1']:<12} {d['id2']:<12} {d['wgs_ibs0']:>10.4f} {d['ukb_ibs0']:>10.4f} "
            f"{d['ibs0_diff']:>10.4f} {d['wgs_kinship']:>10.4f} {d['ukb_kinship']:>10.4f}\n"
        )

print("  ibs0_comparison_summary.txt written.")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

print("\n=== Generating plots ===")

# --- Kinship plots ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.hist(kin_signed, bins=200, edgecolor="none", alpha=0.8, color="#4C72B0")
ax.set_xlabel("Kinship Difference (WGS - UKB)")
ax.set_ylabel("Number of Pairs")
ax.set_title(f"Distribution of Kinship Differences\n(n={n:,} overlapping pairs)")
ax.axvline(0, color="red", linestyle="--", linewidth=0.8)
ax.set_yscale("log")
p99_kin = percentile(sorted([abs(d) for d in kin_signed]), 0.99)
central = [d for d in kin_signed if abs(d) <= p99_kin]
ax_inset = ax.inset_axes([0.55, 0.55, 0.4, 0.4])
ax_inset.hist(central, bins=100, edgecolor="none", alpha=0.8, color="#55A868")
ax_inset.set_title("Central 99%", fontsize=8)
ax_inset.tick_params(labelsize=7)

ax = axes[1]
ax.scatter(kin_ukb, kin_wgs, s=0.3, alpha=0.3, color="#4C72B0", rasterized=True)
lims = [min(min(kin_ukb), min(kin_wgs)) - 0.01, max(max(kin_ukb), max(kin_wgs)) + 0.01]
ax.plot(lims, lims, "r--", linewidth=0.8, label="y = x")
ax.set_xlabel("UKB Kinship")
ax.set_ylabel("WGS Kinship")
ax.set_title(f"Kinship: WGS vs UKB\n(r = {kin_corr:.4f}, n = {n:,})")
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect("equal")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("kinship_comparison_plots.png", dpi=200)
print("  kinship_comparison_plots.png saved.")

# Kinship Bland-Altman
fig2, ax2 = plt.subplots(figsize=(8, 5))
means = [(k + u) / 2 for k, u in zip(kin_wgs, kin_ukb)]
sd_diff = math.sqrt(sum((d - kin_mean_signed) ** 2 for d in kin_signed) / n)
ax2.scatter(means, kin_signed, s=0.3, alpha=0.3, color="#4C72B0", rasterized=True)
ax2.axhline(kin_mean_signed, color="red", linestyle="-", linewidth=0.8,
            label=f"Mean diff = {kin_mean_signed:.4f}")
ax2.axhline(kin_mean_signed + 1.96 * sd_diff, color="gray", linestyle="--", linewidth=0.7,
            label=f"+1.96 SD = {kin_mean_signed + 1.96 * sd_diff:.4f}")
ax2.axhline(kin_mean_signed - 1.96 * sd_diff, color="gray", linestyle="--", linewidth=0.7,
            label=f"-1.96 SD = {kin_mean_signed - 1.96 * sd_diff:.4f}")
ax2.set_xlabel("Mean Kinship [(WGS + UKB) / 2]")
ax2.set_ylabel("Difference (WGS - UKB)")
ax2.set_title("Bland-Altman Plot: Kinship Agreement")
ax2.legend(fontsize=8)
plt.tight_layout()
plt.savefig("kinship_bland_altman.png", dpi=200)
print("  kinship_bland_altman.png saved.")

# --- IBS0 plots ---
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5))

ax = axes3[0]
ax.scatter(ibs0_ukb, ibs0_wgs, s=0.3, alpha=0.3, color="#DD8452", rasterized=True)
lims = [min(min(ibs0_ukb), min(ibs0_wgs)) - 0.001, max(max(ibs0_ukb), max(ibs0_wgs)) + 0.001]
ax.plot(lims, lims, "r--", linewidth=0.8, label="y = x")
ax.set_xlabel("UKB IBS0")
ax.set_ylabel("WGS IBS0")
ax.set_title(f"IBS0: WGS vs UKB\n(r = {ibs0_corr:.4f}, n = {n:,})")
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect("equal")
ax.legend(fontsize=8)

ax = axes3[1]
ax.hist(ibs0_signed, bins=200, edgecolor="none", alpha=0.8, color="#DD8452")
ax.set_xlabel("IBS0 Difference (WGS - UKB)")
ax.set_ylabel("Number of Pairs")
ax.set_title(f"Distribution of IBS0 Differences\n(n={n:,} overlapping pairs)")
ax.axvline(0, color="red", linestyle="--", linewidth=0.8)
ax.set_yscale("log")

plt.tight_layout()
plt.savefig("ibs0_comparison_plots.png", dpi=200)
print("  ibs0_comparison_plots.png saved.")

# IBS0 Bland-Altman
fig4, ax4 = plt.subplots(figsize=(8, 5))
ibs0_means = [(k + u) / 2 for k, u in zip(ibs0_wgs, ibs0_ukb)]
ibs0_sd = math.sqrt(sum((d - ibs0_mean_signed) ** 2 for d in ibs0_signed) / n)
ax4.scatter(ibs0_means, ibs0_signed, s=0.3, alpha=0.3, color="#DD8452", rasterized=True)
ax4.axhline(ibs0_mean_signed, color="red", linestyle="-", linewidth=0.8,
            label=f"Mean diff = {ibs0_mean_signed:.4f}")
ax4.axhline(ibs0_mean_signed + 1.96 * ibs0_sd, color="gray", linestyle="--", linewidth=0.7,
            label=f"+1.96 SD = {ibs0_mean_signed + 1.96 * ibs0_sd:.4f}")
ax4.axhline(ibs0_mean_signed - 1.96 * ibs0_sd, color="gray", linestyle="--", linewidth=0.7,
            label=f"-1.96 SD = {ibs0_mean_signed - 1.96 * ibs0_sd:.4f}")
ax4.set_xlabel("Mean IBS0 [(WGS + UKB) / 2]")
ax4.set_ylabel("Difference (WGS - UKB)")
ax4.set_title("Bland-Altman Plot: IBS0 Agreement")
ax4.legend(fontsize=8)
plt.tight_layout()
plt.savefig("ibs0_bland_altman.png", dpi=200)
print("  ibs0_bland_altman.png saved.")

print("\nDone — all kinship QC outputs written.")
