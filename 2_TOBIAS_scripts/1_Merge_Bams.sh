#!/bin/bash
#SBATCH --account=amc-general
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --job-name=merge_bams_tbet
#SBATCH --output=merge_bams_%j.log
#SBATCH --error=merge_bams_%j.err
#SBATCH --mail-user=saieashan.vankamamidi@cuanschutz.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "Starting BAM merge"
date

module load miniforge
conda activate atac_seq

BAMDIR="/pl/active/getahun_lab/Project_tbet_ATACseq/bam_files"
OUTDIR="/pl/active/getahun_lab/Project_tbet_ATACseq/merged_bams"
mkdir -p $OUTDIR

merge_and_index() {
  GROUP=$1
  shift
  BAMS=("$@")

  echo "Merging $GROUP"
  samtools merge -f -@ 8 \
    ${OUTDIR}/${GROUP}.merged.bam \
    "${BAMS[@]}"

  echo "Indexing $GROUP"
  samtools index -@ 8 -b \
    ${OUTDIR}/${GROUP}.merged.bam \
    ${OUTDIR}/${GROUP}.merged.bam.bai

  echo "Done: $GROUP"
}

merge_and_index L_WT \
  ${BAMDIR}/LWT1.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LWT2.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LWT3.sorted.rmdup.qc.removebl.bam

merge_and_index L_KO \
  ${BAMDIR}/LKO1.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LKO2.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LKO3.sorted.rmdup.qc.removebl.bam

merge_and_index LG_WT \
  ${BAMDIR}/LGWT1.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LGWT2.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LGWT3.sorted.rmdup.qc.removebl.bam

merge_and_index LG_KO \
  ${BAMDIR}/LGKO1.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LGKO2.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LGKO3.sorted.rmdup.qc.removebl.bam

merge_and_index LA21_WT \
  ${BAMDIR}/LA21WT1.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LA21WT2.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LA21WT3.sorted.rmdup.qc.removebl.bam

merge_and_index LA21_KO \
  ${BAMDIR}/LA21KO1.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LA21KO2.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LA21KO3.sorted.rmdup.qc.removebl.bam

merge_and_index LA21G_WT \
  ${BAMDIR}/LA21GWT1.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LA21GWT2.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LA21GWT3.sorted.rmdup.qc.removebl.bam

merge_and_index LA21G_KO \
  ${BAMDIR}/LA21GKO1.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LA21GKO2.sorted.rmdup.qc.removebl.bam \
  ${BAMDIR}/LA21GKO3.sorted.rmdup.qc.removebl.bam

echo "All done"
date