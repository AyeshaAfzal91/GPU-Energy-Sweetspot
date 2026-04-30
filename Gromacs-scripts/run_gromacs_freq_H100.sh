#!/bin/bash -l
#SBATCH --job-name=gmx-multi-bench
#SBATCH --nodes=1
##### SBATCH -w a0805
#SBATCH --gres=gpu:h100:1
##### SBATCH --reservation=powerbench-A100
#SBATCH --time=24:00:00
#SBATCH --output=gmx-bench-%j.out
#SBATCH --export=ALL

# === Load modules ===
# module load gromacs/2024.2-gcc-mkl-cuda11.7 # for testcluster
module load gromacs/2024.4-gcc11.2.0-mkl-cuda # for alex

# === GROMACS GPU environment settings ===
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_USE_GPU_BUFFER_OPS=1
export GMX_DISABLE_GPU_TIMING=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export CUDA_VISIBLE_DEVICES=0
export NSTEPS=200000
export OPTIONS="-maxh 0.2 -ntomp 16 -bonded gpu -update gpu -pme gpu -nb gpu -ntmpi 1 -pin on -pinstride 1"
export INPUTFILE_LOCATION="/lustre/ihpc/ihpc040h/GROMACS/"

# === Validate CLOCK_REF ===
if [ -z "$CLOCK_REF" ]; then
    echo "ERROR: CLOCK_REF not set!"
    exit 1
fi

# === Parse GPU clock settings ===
GPU_MEM_CLOCK="${CLOCK_REF%%,*}"
GPU_GRAPHICS_CLOCK="${CLOCK_REF#*,}"
FREQ_TAG="${GPU_MEM_CLOCK}-${GPU_GRAPHICS_CLOCK}"

cat /proc/cpuinfo

# === Setup output directory ===
DATE=$(date +%Y%m%d)
HOSTNAME=H100-$(hostname -s)
OUTDIR=Host${HOSTNAME}-Freq${FREQ_TAG}-ID${SLURM_JOB_ID}
mkdir -p $OUTDIR
cd $OUTDIR

# === Define benchmark list ===
BENCHMARKS=(2md_start0 FL_md1_berendsen PI_large_test eag1 rnanvt stmv_pme_nvt)

# === Set GPU clocks ===
echo "Setting GPU clocks to MEM=$GPU_MEM_CLOCK, GRAPHICS=$GPU_GRAPHICS_CLOCK"
sudo /usr/bin/nvidia-smi --applications-clocks=$CLOCK_REF

# === gives fan speeds and power usage as reported by the nodes management card  ===
sudo /usr/sbin/ipmi-sensors --sensor-types=Fan,Power_Supply,Current

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
    gmx mdrun -nsteps $NSTEPS $OPTIONS -s ${INPUTFILE_LOCATION}$i.tpr -deffnm ${i}_${FREQ_TAG}_${HOSTNAME}_perflog  -g ${i}_${FREQ_TAG}_${HOSTNAME}_perflog.log > ${i}_${FREQ_TAG}_${HOSTNAME}_perflog.stdout

    echo "End: $(date) / $(date +%s)"

    kill $POWER_PID
    sleep 5
done

# === Reset GPU clocks ===
# sudo /usr/bin/nvidia-smi --reset-applications-clocks
