#!/bin/bash
#SBATCH --account=amc-general
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --mem=48G
#SBATCH --job-name=TOBIAS_Footprints_tbet
#SBATCH --output=TOBIAS_Footprint_scores_%j.log
#SBATCH --error=TOBIAS_Footprint_scores_%j.err
#SBATCH --mail-user=saieashan.vankamamidi@cuanschutz.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load miniforge
conda activate atac_seq

# Inputs
PEAKS="/pl/active/getahun_lab/Project_tbet_ATACseq/consensus_peaks/ALL8_union_peaks.bed"
AC_BASE="/pl/active/getahun_lab/Project_tbet_ATACseq/TOBIAS/ATACorrect_Output"

# Output folder
OUTDIR="/pl/active/getahun_lab/Project_tbet_ATACseq/TOBIAS/FootprintScores_Output"
mkdir -p $OUTDIR

SAMPLES=(
  L_WT.merged
  L_KO.merged
  LG_WT.merged
  LG_KO.merged
  LA21_WT.merged
  LA21_KO.merged
  LA21G_WT.merged
  LA21G_KO.merged
)

for S in "${SAMPLES[@]}"; do
  CORR_BW="${AC_BASE}/${S}/${S}_corrected.bw"
  OUT_BW="${OUTDIR}/${S}_footprints.bw"
  [[ -s "$CORR_BW" ]] || { echo "Missing corrected bigWig: $CORR_BW" >&2; exit 1; }
  echo "Running FootprintScores for $S"
  TOBIAS FootprintScores \
    --signal "$CORR_BW" \
    --regions "$PEAKS" \
    --output "$OUT_BW" \
    --cores "${SLURM_CPUS_PER_TASK}"
done

echo "DONE: FootprintScores finished for all samples."