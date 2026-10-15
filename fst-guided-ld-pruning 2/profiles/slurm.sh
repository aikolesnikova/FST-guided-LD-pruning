#!/usr/bin/env bash
#SBATCH --job-name=fst_ld_prune
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
set -euo pipefail
# Adapt paths, modules and SLURM settings to your HPC environment.
python scripts/fst_ld_pruning.py \
  --bfile "${BFILE:?Set BFILE}" \
  --scores "${SCORES:?Set SCORES}" \
  --outdir "${OUTDIR:?Set OUTDIR}" \
  --threads "${SLURM_CPUS_PER_TASK:-1}"
