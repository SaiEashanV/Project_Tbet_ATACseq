#!/bin/bash
set -euo pipefail

TF_INPUT="${1:-}"
FLANK="${2:-100}"

if [[ -z "$TF_INPUT" ]]; then
  echo "USAGE: bash $0 <TF_NAME> [FLANK]"
  echo "  Example: bash $0 TBX21 100"
  exit 1
fi

BINBASE="/pl/active/getahun_lab/Project_tbet_ATACseq/TOBIAS/BINDetect"
CORRDIR="/pl/active/getahun_lab/Project_tbet_ATACseq/TOBIAS/ATACorrect_Output"
OUTDIR="/pl/active/getahun_lab/Project_tbet_ATACseq/TOBIAS/FootprintPlots"

mkdir -p "$OUTDIR"

BW_L_WT="${CORRDIR}/L_WT.merged/L_WT.merged_corrected.bw"
BW_L_KO="${CORRDIR}/L_KO.merged/L_KO.merged_corrected.bw"
BW_LG_WT="${CORRDIR}/LG_WT.merged/LG_WT.merged_corrected.bw"
BW_LG_KO="${CORRDIR}/LG_KO.merged/LG_KO.merged_corrected.bw"
BW_LA21_WT="${CORRDIR}/LA21_WT.merged/LA21_WT.merged_corrected.bw"
BW_LA21_KO="${CORRDIR}/LA21_KO.merged/LA21_KO.merged_corrected.bw"
BW_LA21G_WT="${CORRDIR}/LA21G_WT.merged/LA21G_WT.merged_corrected.bw"
BW_LA21G_KO="${CORRDIR}/LA21G_KO.merged/LA21G_KO.merged_corrected.bw"

for f in "$BW_L_WT" "$BW_L_KO" "$BW_LG_WT" "$BW_LG_KO" \
         "$BW_LA21_WT" "$BW_LA21_KO" "$BW_LA21G_WT" "$BW_LA21G_KO"; do
  if [[ ! -s "$f" ]]; then
    echo "[ERROR] Missing/empty bigwig: $f"
    exit 1
  fi
done

TF=$(echo "$TF_INPUT" | xargs)

find_bed () {
  local condition="$1"
  local MOTIFDIR
  MOTIFDIR=$(find "${BINBASE}/${condition}_KO_vs_WT" -maxdepth 1 -type d -iname "*${TF}*" | head -n 1)
  if [[ -z "$MOTIFDIR" ]]; then
    echo "[ERROR] Could not find motif folder for TF=${TF} in ${condition}_KO_vs_WT" >&2
    exit 1
  fi
  local BED
  BED=$(find "${MOTIFDIR}/beds" -type f -iname "*_WT_bound.bed" | head -n 1)
  if [[ -z "$BED" || ! -s "$BED" ]]; then
    echo "[ERROR] Missing/empty BED for TF=${TF} in ${condition}: ${BED:-NONE}" >&2
    exit 1
  fi
  echo "$BED"
}

BED_L=$(find_bed "L")
BED_LG=$(find_bed "LG")
BED_LA21=$(find_bed "LA21")
BED_LA21G=$(find_bed "LA21G")

echo "[INFO] TF: ${TF}"
echo "[INFO] BED L:     ${BED_L}"
echo "[INFO] BED LG:    ${BED_LG}"
echo "[INFO] BED LA21:  ${BED_LA21}"
echo "[INFO] BED LA21G: ${BED_LA21G}"

OUTPNG="${OUTDIR}/${TF}_4conditions_KO_vs_WT_flank${FLANK}.png"
TMPPY="${OUTDIR}/plot_${TF}_tmp.py"

cat > "$TMPPY" << 'PYEOF'
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pyBigWig

tf      = sys.argv[1]
beds    = sys.argv[2].split(",")
bw_wt   = sys.argv[3].split(",")
bw_ko   = sys.argv[4].split(",")
outpng  = sys.argv[5]
flank   = int(sys.argv[6])

CONDITIONS = ["L", "LG", "LA21", "LA21G"]
COLORS = {"WT": "#2166AC", "KO": "#D7191C"}
OUTLIER_PERCENTILE = 99
SMOOTH_SIGMA = 2

def gaussian_smooth(y, sigma=SMOOTH_SIGMA):
    radius = max(1, int(4 * sigma + 0.5))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-(x**2) / (2 * sigma**2))
    kernel /= kernel.sum()
    return np.convolve(y, kernel, mode="same")

def load_sites(bedfile):
    sites = []
    with open(bedfile) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end   = int(parts[2])
            strand = parts[5] if len(parts) > 5 else "+"
            center = (start + end) // 2
            motif_width = end - start
            sites.append((chrom, center, strand, start, end, motif_width))
    return sites

def average_profile(bwfile, sites, flank, label=""):
    bw = pyBigWig.open(bwfile)
    chrom_sizes = bw.chroms()
    profiles = []
    n_boundary_skip = 0
    n_runtime_skip  = 0
    for chrom, center, strand, _, _, _ in sites:
        ext_start = center - flank
        ext_end   = center + flank
        if ext_start < 0 or chrom not in chrom_sizes or ext_end > chrom_sizes[chrom]:
            n_boundary_skip += 1
            continue
        try:
            vals = bw.values(chrom, ext_start, ext_end, numpy=True)
        except RuntimeError:
            n_runtime_skip += 1
            continue
        if vals is None or len(vals) != 2 * flank:
            continue
        vals = np.array(vals, dtype=float)
        if strand == "-":
            vals = vals[::-1]
        profiles.append(vals)
    bw.close()
    print(f"[INFO] {label}: {n_boundary_skip} boundary skip, {n_runtime_skip} runtime skip", file=sys.stderr)
    if len(profiles) == 0:
        raise ValueError(f"No valid signal from {bwfile}")
    arr = np.vstack(profiles)
    row_maxes = np.nanmax(arr, axis=1)
    upper = np.percentile(row_maxes, OUTLIER_PERCENTILE)
    keep = row_maxes <= upper
    arr = arr[keep]
    print(f"[INFO] {label}: {len(arr)} retained, {np.sum(~keep)} outliers removed", file=sys.stderr)
    mean_p = np.nanmean(arr, axis=0)
    n_valid = np.sum(~np.isnan(arr), axis=0)
    sem_p  = np.nanstd(arr, axis=0, ddof=1) / np.sqrt(n_valid)
    x = np.arange(-flank, flank)
    return x, mean_p, sem_p

def calc_fpd(y, flank, motif_width):
    half_motif = max(1, motif_width // 2)
    center_idx = list(range(flank - half_motif, flank + half_motif))
    flank_idx  = list(range(0, flank - half_motif)) + list(range(flank + half_motif, len(y)))
    return np.mean(y[flank_idx]) - np.mean(y[center_idx])

fig, axes = plt.subplots(1, 4, figsize=(14, 3.5), sharey=True)
fig.suptitle(tf, fontsize=14, fontweight="bold", y=1.02)

for idx, (cond, bed, bw_wt_f, bw_ko_f) in enumerate(zip(CONDITIONS, beds, bw_wt, bw_ko)):
    ax = axes[idx]
    sites = load_sites(bed)
    motif_width = sites[0][5] if sites else 20
    motif_half  = motif_width // 2

    for geno, bwf in [("WT", bw_wt_f), ("KO", bw_ko_f)]:
        x, y, sem = average_profile(bwf, sites, flank, label=f"{cond}/{geno}")
        y   = gaussian_smooth(y)
        sem = gaussian_smooth(sem)
        fpd = calc_fpd(y, flank, motif_width)
        ax.plot(x, y, linewidth=1.5, color=COLORS[geno],
                label=f"{geno} (FPD={fpd:.3f})", zorder=3)
        ax.fill_between(x, y - sem, y + sem, color=COLORS[geno], alpha=0.15, zorder=2)

    ax.axvspan(-motif_half, motif_half, color="grey", alpha=0.08, zorder=1)
    ax.axvline(-motif_half, color="grey", linewidth=0.8, linestyle="--", zorder=1)
    ax.axvline( motif_half, color="grey", linewidth=0.8, linestyle="--", zorder=1)

    ax.set_title(cond, fontsize=12, fontweight="bold")
    ax.tick_params(labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.set_xlim(-flank, flank)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(50))
    ax.set_xlabel("bp from centre", fontsize=9)
    ax.legend(fontsize=8, frameon=False, loc="upper right")

    if idx == 0:
        ax.set_ylabel("mean corrected cutsite signal", fontsize=10)
    else:
        ax.spines["left"].set_visible(False)
        ax.tick_params(left=False)

plt.tight_layout()
plt.savefig(outpng, dpi=300, bbox_inches="tight", pad_inches=0.15)
print(f"[DONE] Saved: {outpng}")
PYEOF

python "$TMPPY" \
  "$TF" \
  "${BED_L},${BED_LG},${BED_LA21},${BED_LA21G}" \
  "${BW_L_WT},${BW_LG_WT},${BW_LA21_WT},${BW_LA21G_WT}" \
  "${BW_L_KO},${BW_LG_KO},${BW_LA21_KO},${BW_LA21G_KO}" \
  "$OUTPNG" \
  "$FLANK"

rm -f "$TMPPY"
echo "[DONE] $OUTPNG"