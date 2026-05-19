#!/bin/bash -l
#SBATCH --job-name=freqA100
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --reservation=powercapped-ihpc161h-a100
#SBATCH --gres=gpu:a100:1
#SBATCH --partition=a100
#SBATCH -C el8
#SBATCH --output=FreqA100-%j.out

# === Load modules ===
module load gromacs/2024.4-gcc11.2.0-mkl-cuda # for alex

# === GROMACS GPU environment settings ===
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_USE_GPU_BUFFER_OPS=1
export GMX_DISABLE_GPU_TIMING=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export CUDA_VISIBLE_DEVICES=0
export NSTEPS=200000
export OPTIONS="-maxh 0.2 -ntomp 16 -bonded gpu -update gpu -pme gpu -nb gpu -ntmpi 1 -pin on -pinstride 1"
export INPUTFILE_LOCATION="/home/hpc/ihpc/ihpc161h/GROMACS-BAthesis/inputs/"

for ((J=210; J<=1454; J+=45)); do
if (( J > 1410)); then 
	echo "edge case 1410"
	J=1410
fi

# === Parse GPU clock settings ===
GPU_MEM_CLOCK="1251"
GPU_GRAPHICS_CLOCK="$J"
FREQ_TAG="1215-${GPU_GRAPHICS_CLOCK}"

# === Setup output directory ===
DATE=$(date +%Y%m%d)
OUTDIR=/home/hpc/ihpc/ihpc161h/GROMACS-BAthesis/runs/$1/${FREQ_TAG}
mkdir -p $OUTDIR
cd $OUTDIR

# === Define benchmark list ===
BENCHMARKS=(2md_start0 FL_md1_berendsen PI_large_test eag1 rnanvt stmv_pme_nvt)

# === Set GPU clocks ===
echo "Setting GPU clocks to GRAPHICS=$GPU_GRAPHICS_CLOCK"
sudo /usr/bin/nvidia-smi --lock-gpu-clocks="$GPU_GRAPHICS_CLOCK"

# === Signal cleanup for power logging ===
trap "kill 0" SIGINT SIGTERM EXIT

# === Run each benchmark ===
for i in "${BENCHMARKS[@]}"; do
    echo ">>> Running benchmark: $i"

    if [ ! -f "${INPUTFILE_LOCATION}$i.tpr" ]; then
        echo "ERROR: ${INPUTFILE_LOCATION}$i.tpr not found! Skipping $i."
        continue
    fi
    
    # Power logging
    nvidia-smi -i $CUDA_VISIBLE_DEVICES --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > ${i}_${FREQ_TAG}_${HOSTNAME}_powerlog.csv & POWER_PID=$!

    # Main GROMACS run
    gmx mdrun -nsteps $NSTEPS $OPTIONS -s ${INPUTFILE_LOCATION}$i.tpr -deffnm ${i}_${FREQ_TAG}_${HOSTNAME}_perflog  -g ${i}_${FREQ_TAG}_${HOSTNAME}_perflog.log -gpu_id $CUDA_VISIBLE_DEVICES > ${i}_${FREQ_TAG}_${HOSTNAME}-${SLURM_JOB_ID}.out 2>&1 > ${i}_${FREQ_TAG}_${HOSTNAME}_perflog.stdout

    echo "End: $(date) / $(date +%s)"

    kill $POWER_PID
    sleep 5
done
done

