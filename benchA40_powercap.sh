#!/bin/bash -l
#SBATCH --job-name=benchA40-powercap
#SBATCH --nodes=1
#SBATCH --gres=gpu:a40:1
#SBATCH --time=24:00:00

if [ $# -ne 1 ]; then
	echo "Usage: $0 comment-string"
	exit 1
fi

# example clock ref: "1215,1410"

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

for ((J=100; J<=300; J+=100)); do  #in the paper the range was set from 100 to 300 for the A40 with a TDP OF 300
POWER_LIMIT=$J

# === SET GPU CLOCKS ===  
echo "SETTING GPU POWER LIMIT TO ${POWER_LIMIT}W"
sudo /usr/bin/nvidia-smi --power-limit=${POWER_LIMIT}

OUTDIR=$POWER_LIMIT
mkdir -p "$1/$OUTDIR"

for i in Cellulose JAC STMV FactorIX; do
	for n in NPT NVE; do

#nvidia-smi -q -i $CUDA_VISIBLE_DEVICES -d POWER,UTILIZATION --loop > "$RESULT/$OUTDIR/$i-$n-`date +%Y%m%d`-`hostname -s`-$1.$$".powerlog & PID=$!

nvidia-smi -i $CUDA_VISIBLE_DEVICES --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > $RESULT/$OUTDIR/${i}_${n}_${POWER_LIMIT}_powerlog.csv & POWER_PID=$!

echo "!!!! SMT IST AKTIV! ALSO NUR 8 PHYSIKALISCHE KERNE !!!!"

pmemd.cuda -O -i PME/${i}_production_${n}_4fs/mdin.GPU -o "$RESULT/$OUTDIR/${i}_${n}_$1.$$".mdout -p PME/Topologies/$i.prmtop -c PME/Coordinates/$i.inpcrd

echo "End: `date`  / `date +%s`"

kill $POWER_PID
	done
done

for i in TRPCage myoglobin nucleosome; do 

#nvidia-smi -q -i $CUDA_VISIBLE_DEVICES -d POWER,UTILIZATION --loop > "$RESULT/$OUTDIR/$i-`date +%Y%m%d`-`hostname -s`-$1.$$".powerlog &PID=$!

nvidia-smi -i $CUDA_VISIBLE_DEVICES --loop-ms=100 --query-gpu=timestamp,power.draw,utilization.gpu --format=csv > $RESULT/$OUTDIR/${i}_${POWER_LIMIT}_powerlog.csv & POWER_PID=$!

echo "!!!! SMT ist aktiv! Also nur 8 physikalische Kerne !!!!"

pmemd.cuda -O -i GB/$i/mdin.GPU -o "$RESULT/$OUTDIR/${i}_$1.$$".mdout -p GB/$i/prmtop -c GB/$i/inpcrd

echo "End: `date`  / `date +%s`"

kill $POWER_PID

done
done 
