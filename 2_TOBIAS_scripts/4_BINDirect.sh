#!/bin/bash
#SBATCH --account=amc-general
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --job-name=TOBIAS_BINDetect_tbet
#SBATCH --output=TOBIAS_BINDetect_%j.log
#SBATCH --error=TOBIAS_BINDetect_%j.err
#SBATCH --mail-user=saieashan.vankamamidi@cuanschutz.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load miniforge
conda activate atac_seq

# ====== INPUTS ======
GENOME="/pl/active/getahun_lab/atacseq_refdata/genomes/mm10_primary/mm10.fa"
PEAKS="/pl/active/getahun_lab/Project_tbet_ATACseq/consensus_peaks/ALL8_union_peaks.bed"
FPDIR="/pl/active/getahun_lab/Project_tbet_ATACseq/TOBIAS/FootprintScores_Output"
MOTIFS="/pl/active/getahun_lab/atacseq_refdata/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme/JASPAR2022_CORE_vertebrates.meme"
OUTBASE="/pl/active/getahun_lab/Project_tbet_ATACseq/TOBIAS/BINDetect"

run_bindetect () {
  local outdir="$1"; shift
  local cond1="$1"; shift
  local cond2="$1"; shift
  local sig1="$1"; shift
  local sig2="$1"; shift
  mkdir -p "$outdir"
  TOBIAS BINDetect \
    --motifs $MOTIFS \
    --signals "$sig1" "$sig2" \
    --genome "$GENOME" \
    --peaks "$PEAKS" \
    --outdir "$outdir" \
    --cond_names "$cond1" "$cond2" \
    --cores "${SLURM_CPUS_PER_TASK}"
}

# ====== 1) L: KO vs WT ======
run_bindetect \
  "${OUTBASE}/L_KO_vs_WT" \
  "KO" "WT" \
  "${FPDIR}/L_KO.merged_footprints.bw" \
  "${FPDIR}/L_WT.merged_footprints.bw"

# ====== 2) LG: KO vs WT ======
run_bindetect \
  "${OUTBASE}/LG_KO_vs_WT" \
  "KO" "WT" \
  "${FPDIR}/LG_KO.merged_footprints.bw" \
  "${FPDIR}/LG_WT.merged_footprints.bw"

# ====== 3) LA21: KO vs WT ======
run_bindetect \
  "${OUTBASE}/LA21_KO_vs_WT" \
  "KO" "WT" \
  "${FPDIR}/LA21_KO.merged_footprints.bw" \
  "${FPDIR}/LA21_WT.merged_footprints.bw"

# ====== 4) LA21G: KO vs WT ======
run_bindetect \
  "${OUTBASE}/LA21G_KO_vs_WT" \
  "KO" "WT" \
  "${FPDIR}/LA21G_KO.merged_footprints.bw" \
  "${FPDIR}/LA21G_WT.merged_footprints.bw"

echo "DONE: BINDetect completed for all four conditions."