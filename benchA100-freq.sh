#!/bin/bash -l
#SBATCH --job-name=benchA100-freq
#SBATCH --nodes=1
#SBATCH --gres=gpu:a100:1
#SBATCH --time=24:00:00

if [ $# -ne 1 ]; then
	echo "Usage: $0 comment-string"
	exit 1
fi

echo "Saving bench log also to log file ..."
echo "Please be patient ..."
exec > >( tee amber24-bench-`date +%Y%m%d`-`hostname -s`-$1.$$.log) 2>&1

# save script for reference :)
echo "### This is $0"
cat $0

echo "Start: `date`  / `date +%s`"

echo "########################################"

# testcluster, tinygpu
#module load amber-gpu/24p02-at24p04-gnu-cuda12.4

# alex
module load amber/24p02-at24p05-gnu-cuda12.4  

#cat /proc/cpuinfo

RESULT=$1
mkdir -p "$RESULT"

export CUDA_VISIBLE_DEVICES=0
nvidia-smi -i $CUDA_VISIBLE_DEVICES

for ((J=210; J<=1410; J+=600)); do

# === PARSE GPU CLOCK SETTINGS ===
#A100 only supports 1215
GPU_MEM_CLOCK="1215"
GPU_GRAPHICS_CLOCK="$J"
FREQ_TAG="1215-$J"

# === SET GPU CLOCKS ===  
echo "SETTING GPU CLOCKS TO MEM=$GPU_MEM_CLOCK, GRAPHICS=$GPU_GRAPHICS_CLOCK"
# example clock ref: "1215,1410"
sudo /usr/bin/nvidia-smi --applications-clocks=$CLOCK_REF

OUTDIR=$FREQ_TAG
mkdir -p "$1/$OUTDIR"

for i in Cellulose JAC STMV FactorIX; do
	for n in NPT NVE; do

#nvidia-smi -q -i $CUDA_VISIBLE_DEVICES -d POWER,UTILIZATION --loop > "$RESULT/$OUTDIR/$i-$n-`date +%Y%m%d`-`hostname -s`-$1.$$".powerlog & PID=$!

nvidia-smi -i $CUDA_VISIBLE_DEVICES --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > $RESULT/$OUTDIR/${i}_${n}_${FREQ_TAG}_powerlog.csv & POWER_PID=$!

echo "!!!! SMT IST AKTIV! ALSO NUR 8 PHYSIKALISCHE KERNE !!!!"

pmemd.cuda -O -i PME/${i}_production_${n}_4fs/mdin.GPU -o "$RESULT/$OUTDIR/${i}_${n}_$1.$$".mdout -p PME/Topologies/$i.prmtop -c PME/Coordinates/$i.inpcrd

echo "End: `date`  / `date +%s`"

kill $POWER_PID
	done
done

for i in TRPCage myoglobin nucleosome; do 

#nvidia-smi -q -i $CUDA_VISIBLE_DEVICES -d POWER,UTILIZATION --loop > "$RESULT/$OUTDIR/$i-`date +%Y%m%d`-`hostname -s`-$1.$$".powerlog &PID=$!

nvidia-smi -i $CUDA_VISIBLE_DEVICES --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > $RESULT/$OUTDIR/${i}_${FREQ_TAG}_powerlog.csv & POWER_PID=$!

echo "!!!! SMT ist aktiv! Also nur 8 physikalische Kerne !!!!"

pmemd.cuda -O -i GB/$i/mdin.GPU -o "$RESULT/$OUTDIR/${i}_$1.$$".mdout -p GB/$i/prmtop -c GB/$i/inpcrd

echo "End: `date`  / `date +%s`"

kill $POWER_PID

done
done 
