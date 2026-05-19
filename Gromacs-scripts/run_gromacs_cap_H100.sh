#!/bin/bash -l
#SBATCH --job-name=capH100
#SBATCH --nodes=1
#SBATCH --gres=gpu:h100:1
#SBATCH --reservation=powercapped-ihpc161h-h100
#SBATCH --time=24:00:00
##### SBATCH --output=FreqA100-%j.out
#SBATCH --export=ALL

# === Load modules ===
module load gromacs/2025.2-mkl-cuda12 #for helma 

# === GROMACS GPU environment settings ===
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_USE_GPU_BUFFER_OPS=1
export GMX_DISABLE_GPU_TIMING=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export CUDA_VISIBLE_DEVICES=0
export NSTEPS=200000
export OPTIONS="-maxh 0.2 -ntomp 16 -bonded gpu -update gpu -pme gpu -nb gpu -ntmpi 1 -pin on -pinstride 1"
export INPUTFILE_LOCATION="/home/hpc/ihpc/ihpc161h/GROMACS-BAthesis/inputs/"

for ((J=200; J<=700; J+=20)); do 
POWER_LIMIT=$J
if (( J > 700)); then
  echo "edge case for over 700"
  POWER_LIMIT=700
fi

# === Setup output directory ===
DATE=$(date +%Y%m%d)
OUTDIR=/home/hpc/ihpc/ihpc161h/GROMACS-BAthesis/runs/$1/${POWER_LIMIT}
mkdir -p $OUTDIR
cd $OUTDIR

# === Define benchmark list ===
BENCHMARKS=(2md_start0 FL_md1_berendsen PI_large_test eag1 rnanvt stmv_pme_nvt)

# === Set GPU clocks ===
echo "Setting GPU Powercap to ${POWER_LIMIT}"
sudo /usr/bin/nvidia-smi --power-limit=${POWER_LIMIT}

# === Run each benchmark ===
for i in "${BENCHMARKS[@]}"; do
    echo ">>> Running benchmark: $i"

    if [ ! -f "${INPUTFILE_LOCATION}$i.tpr" ]; then
        echo "ERROR: ${INPUTFILE_LOCATION}$i.tpr not found! Skipping $i."
        continue
    fi
    
    # Power logging
    nvidia-smi -i $CUDA_VISIBLE_DEVICES --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > ${i}_${POWER_LIMIT}_${HOSTNAME}_powerlog.csv & POWER_PID=$!

    # Main GROMACS run
    gmx mdrun -nsteps $NSTEPS $OPTIONS -s ${INPUTFILE_LOCATION}$i.tpr -deffnm ${i}_${POWER_LIMIT}_${HOSTNAME}_perflog  -g ${i}_${POWER_LIMIT}_${HOSTNAME}_perflog.log -gpu_id $CUDA_VISIBLE_DEVICES > ${i}_${POWER_LIMIT}_${HOSTNAME}-${SLURM_JOB_ID}.out 2>&1 > ${i}_${POWER_LIMIT}_${HOSTNAME}_perflog.stdout

    echo "End: $(date) / $(date +%s)"

    kill $POWER_PID
    sleep 5
done
done

