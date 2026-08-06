#!/bin/bash
# phipmake2 SLURM submission. Keeps the same --export interface as the
# runphipmake.sh it replaces.
#
#   sbatch --export=wd="/path/to/screen/",plan="Counts-FoldChange-Enrichment-Hits" \
#          runphipmake2.sh
#
# Variables:
#   wd        screen directory containing the input matrices (required)
#   plan      dash-separated stages. Default Counts-FoldChange-Enrichment-Hits;
#             append -Polyclonal to include polyclonal scoring
#   metadata  peptide library metadata directory
#   threads   data.table threads. Default: --cpus-per-task
#
# Memory is bounded by one sub-library at a time, so these limits do not need
# to scale with screen size.

#SBATCH --job-name=phipmake2
#SBATCH --time=4:00:00
#SBATCH --partition=parallel
#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=%J-phipmake2.out
#SBATCH --mail-type=END

set -euo pipefail

date
echo "wd=${wd:?wd is required}"
echo "plan=${plan:=Counts-FoldChange-Enrichment-Hits}"
echo "threads=${threads:=${SLURM_CPUS_PER_TASK:-4}}"

umask 000
ml r/4.3.0 || module load r/4.3.0

ARGS=(--wd "$wd" --stages "$plan" --threads "$threads")
if [ -n "${metadata:-}" ]; then ARGS+=(--metadata "$metadata"); fi

CLI=$(Rscript -e 'cat(system.file("cli/phipmake2.R", package = "phipmake2"))')
if [ -z "$CLI" ]; then
  echo "phipmake2 is not installed in this R library path." >&2
  exit 1
fi
Rscript "$CLI" "${ARGS[@]}"

mkdir -p "${wd}/logs"
if [ -n "${SLURM_JOB_ID:-}" ] && [ -f "${SLURM_JOB_ID}-phipmake2.out" ]; then
  cp "${SLURM_JOB_ID}-phipmake2.out" "${wd}/logs/${SLURM_JOB_ID}-phipmake2.txt"
fi

chmod -R 777 "$wd" || true
echo "phipmake2 done"
date
