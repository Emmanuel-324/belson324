DATA_DIR="/home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_03"
SCRIPT="/home/emmanuel324/projects/belson324/aspect_ratio.py"

cat > "$DATA_DIR/submit.sh" <<'SH'
#!/bin/bash
#SBATCH -A amcorrosion
#SBATCH -p normal_q
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J es03_aspect
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --array=1-567%40

set -euo pipefail
ROOT="/home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_03"
FILELIST="$ROOT/files.txt"
OUTDIR="$ROOT/out"
SCRIPT="/home/emmanuel324/projects/belson324/aspect_ratio.py"
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VTK_SMP_MAX_THREADS=${SLURM_CPUS_PER_TASK}
mkdir -p "$OUTDIR"
python "$SCRIPT" --file "$ROOT/$FILE" --outdir "$OUTDIR"
SH

chmod +x "$DATA_DIR/submit.sh"
