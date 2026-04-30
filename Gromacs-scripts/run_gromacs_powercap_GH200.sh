#!/bin/bash -l
#SBATCH --job-name=powercap-gracehop1-GH200
#SBATCH --nodes=1
#SBATCH -w gracehop1                
#SBATCH --reservation=powerbench-GPU
#SBATCH --time=24:00:00
##### SBATCH --output=gmx-bench-%j.out
#SBATCH --export=ALL

# === Load modules ===
module load cuda nvhpc gromacs #/2024.3-gcc-mkl-cuda11.7 (testcluster) (no module for gracehop1) # module load gromacs/2024.2-gcc-mkl-cuda11.7 # for testcluster gromacs/2024.4-gcc11.2.0-mkl-cuda (alex)

# === List GPUs, manually expose the GH200 GPU and confirm which GPU is visible  ===
nvidia-smi --query-gpu=index,name,uuid,pci.bus_id --format=csv,noheader
nvidia-smi -L
export CUDA_VISIBLE_DEVICES=0
echo "Using GPU: $CUDA_VISIBLE_DEVICES"
nvidia-smi

# === GROMACS GPU environment settings ===
export GMX_GPU_PME_DECOMPOSITION=1
export GMX_USE_GPU_BUFFER_OPS=1
export GMX_DISABLE_GPU_TIMING=1
export GMX_ENABLE_DIRECT_GPU_COMM=1
export NSTEPS=200000
export OPTIONS="-maxh 0.2 -ntomp 16 -bonded gpu -update gpu -pme gpu -nb gpu -ntmpi 1 -pin on -pinstride 1"
export INPUTFILE_LOCATION="/home/vault/ihpc/ihpc040h/GROMACS/"

# === Validate CLOCK_REF ===
#if [ -z "$CLOCK_REF" ]; then
#    echo "ERROR: CLOCK_REF not set!"
#    exit 1
#fi

# === Validate POWER_LIMIT ===
if [ -z "$POWER_LIMIT" ]; then
    echo "ERROR: POWER_LIMIT not set!"
    exit 1
fi

# === Parse GPU clock settings ===
#GPU_MEM_CLOCK="${CLOCK_REF%%,*}"
#GPU_GRAPHICS_CLOCK="${CLOCK_REF#*,}"
#FREQ_TAG="${GPU_MEM_CLOCK}-${GPU_GRAPHICS_CLOCK}"

cat /proc/cpuinfo

# === Setup output directory ===
DATE=$(date +%Y%m%d)
HOSTNAME=${NODE}-${CUDA_VISIBLE_DEVICES}-$(hostname -s)
OUTDIR=${HOSTNAME}-PowerCAP${POWER_LIMIT}W-ID${SLURM_JOB_ID}
mkdir -p $OUTDIR
cd $OUTDIR

# === Define benchmark list ===
BENCHMARKS=(2md_start0 FL_md1_berendsen PI_large_test eag1 rnanvt stmv_pme_nvt system7_KP)

# === Set GPU clocks ===
# echo "Setting GPU clocks to MEM=$GPU_MEM_CLOCK, GRAPHICS=$GPU_GRAPHICS_CLOCK"
# sudo /usr/bin/nvidia-smi --applications-clocks=$CLOCK_REF

# === Set Power Limit ===
sudo /usr/bin/nvidia-smi --power-limit=${POWER_LIMIT}

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
    nvidia-smi -i $CUDA_VISIBLE_DEVICES --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > ${i}_${POWER_LIMIT}_${HOSTNAME}_powerlog.csv & POWER_PID=$!

    # Main GROMACS run
    gmx mdrun -nsteps $NSTEPS $OPTIONS -s ${INPUTFILE_LOCATION}$i.tpr -deffnm ${i}_${POWER_LIMIT}_${HOSTNAME}_perflog  -g ${i}_${POWER_LIMIT}_${HOSTNAME}_perflog.log -gpu_id $CUDA_VISIBLE_DEVICES${i}_${POWER_LIMIT}_${HOSTNAME}-${SLURM_JOB_ID}.out 2>&1
# for gracehop: /home/atuin/unrz/unrz007h/new_spack/opt/spack/linux-ubuntu22.04-neoverse_v2/gcc-12.3.0/gromacs-2024.1-gnqvkzl5p7vjlwhttzsrr63sqlu264yt/bin/gmx mdrun -nsteps $NSTEPS $OPTIONS -s ${INPUTFILE_LOCATION}$i.tpr -deffnm ${i}_${FREQ_TAG}_${HOSTNAME}_perflog  -g ${i}_${FREQ_TAG}_${HOSTNAME}_perflog.log > ${i}_${FREQ_TAG}_${HOSTNAME}_perflog.stdout
# gmx mdrun -nsteps $NSTEPS $OPTIONS -s ${INPUTFILE_LOCATION}$i.tpr -deffnm ${i}_${FREQ_TAG}_${HOSTNAME}_perflog2  -g ${i}_${FREQ_TAG}_${HOSTNAME}_perflog2.log

    echo "End: $(date) / $(date +%s)"

    kill $POWER_PID
    sleep 5
done

# === Reset GPU clocks ===
# sudo /usr/bin/nvidia-smi --reset-applications-clocks
