#!/bin/bash
#SBATCH --job-name=AAR-012-cupin1
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --output=/scratch/p318738/SSA/AAR-COMP-012/AAR-COMP-012-Step3/logs/cupin1_%j.log
#SBATCH --error=/scratch/p318738/SSA/AAR-COMP-012/AAR-COMP-012-Step3/logs/cupin1_%j.err

# AAR-COMP-012 Step 3 — Cupin domain structure predictions, batch 1 (9 sequences)
# Spb40, Tri28, Aza12, NngM, DobE, PyrN, Afn9, ThzN, Apy9

set -euo pipefail

YAML_DIR="/scratch/p318738/SSA/AAR-COMP-012/AAR-COMP-012-Step3/yaml/cupin_batch1"
OUT_DIR="/scratch/p318738/SSA/AAR-COMP-012/AAR-COMP-012-Step3/cupin"
CACHE_DIR="/scratch/p318738/boltz2/cache"

module purge
module load CUDA/12.1.1

source /cvmfs/hpc.rug.nl/versions/2023.01/rocky8/x86_64/amd/zen3/software/Miniconda3/22.11.1-1/etc/profile.d/conda.sh
conda activate boltz2

export BOLTZ_CACHE="$CACHE_DIR"
mkdir -p "$OUT_DIR" "$CACHE_DIR"

echo "Host:      $(hostname)"
echo "Job ID:    $SLURM_JOB_ID"
echo "YAML dir:  $YAML_DIR"
echo "Output:    $OUT_DIR"
echo "Sequences: $(ls $YAML_DIR/*.yaml | wc -l)"
nvidia-smi || true
python -c "import torch; print('CUDA available:', torch.cuda.is_available())"

echo ""
echo "Starting boltz predict — cupin batch 1"
echo "======================================="

srun boltz predict "$YAML_DIR" \
    --cache "$CACHE_DIR" \
    --out_dir "$OUT_DIR" \
    --use_msa_server \
    --recycling_steps 3 \
    --diffusion_samples 3 \
    --use_potentials

echo ""
echo "cupin batch 1 complete."
