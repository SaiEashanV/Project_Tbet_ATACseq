#!/bin/bash
#SBATCH --account=amc-general
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --job-name=TOBIAS_ATACorrect_tbet
#SBATCH --output=TOBIAS_ATACorrect_%j.log
#SBATCH --error=TOBIAS_ATACorrect_%j.err
#SBATCH --mail-user=saieashan.vankamamidi@cuanschutz.edu
#SBATCH --mail-type=END,FAIL

# ATACorrect models and removes Tn5 sequence insertion bias from ATAC-seq reads
set -euo pipefail
echo "Starting TOBIAS ATACorrect"
date

module load miniforge
conda activate atac_seq

# ===== PATHS =====
BAMDIR="/pl/active/getahun_lab/Project_tbet_ATACseq/merged_bams"
GENOME="/pl/active/getahun_lab/atacseq_refdata/genomes/mm10_primary/mm10.fa"
BLACKLIST="/pl/active/getahun_lab/atacseq_refdata/genomes/mm10_primary/mm10-blacklist.v2.bed"
PEAKS="/pl/active/getahun_lab/Project_tbet_ATACseq/consensus_peaks/ALL8_union_peaks.bed"
OUTBASE="/pl/active/getahun_lab/Project_tbet_ATACseq/TOBIAS/ATACorrect_Output"
mkdir -p $OUTBASE

# ===== SAMPLES =====
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

# ===== RUN =====
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Processing $SAMPLE"
  TOBIAS ATACorrect \
    --bam ${BAMDIR}/${SAMPLE}.bam \
    --genome $GENOME \
    --peaks $PEAKS \
    --blacklist $BLACKLIST \
    --outdir ${OUTBASE}/${SAMPLE} \
    --cores 8
done

echo "Finished ATACorrect"
date