#!/bin/bash -l
### running on 2 nodes with 72 cores per node = 2*72 = 144 mpi tasks
#SBATCH --job-name=20p12-at21p11-ompi-intel
##SBATCH --output=gmx.o%j
##SBATCH --error=gmx.e%j
##SBATCH --partition=work
#SBATCH -N 3
#SBATCH --tasks-per-node=72
#SBATCH --partition=multinode
#SBATCH --time=24:00:00
#SBATCH --export=NONE
##SBATCH --mail-type=end
##SBATCH --mail-user=tobias.kloeffel@fau.de
#SBATCH --get-user-env
echo "Starting: $(date)"
unset SLURM_EXPORT_ENV

module load amber/20p12-at21p11-ompi-intel

export I_MPI_PIN_RESPECT_CPUSET=0
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
export I_MPI_PIN=1

numactl -H
num=$((SLURM_NNODES*72))
space=1
runtime=2.0

cd ${SLURM_SUBMIT_DIR}

work=./${SLURM_JOBID}

unset KMP_AFFINITY

cd ${SLURM_SUBMIT_DIR}
run=n004
mkdir -p $run
cd $run
root=$PWD

export OMP_PROC_BIND=close
export OMP_PLACES=cores

input="factorix jac"

binary=pmemd.MPI
ldd $binary

for inputnum in $input
do
    export OMP_NUM_THREADS=1
    tmpi=$((num/1/space))
    tnode=$((tmpi/SLURM_NNODES))
    cp ${SLURM_SUBMIT_DIR}/${inputnum}* .

	mkdir -p $work
	rm -rf $work/*
	cd $work
	
	mpi="srun"
	opts=" -O -i ../${inputnum}.mdin -c ../${inputnum}.inpcrd -p ../${inputnum}.prmtop -o ${inputnum}.output "

	$mpi $binary $opts &
	pid=$!
        (sleep ${runtime}h; kill -s 9 $pid) &
        pid1=$!
	wait $pid
	kill -s 9 $pid1
	
	mkdir -p $root/${SLURM_NNODES}/pmemd/${inputnum}
	cp -a ./* $root/${SLURM_NNODES}/pmemd/${inputnum}
	
	cd $root
done

input="jac dhfr"

binary=sander.MPI
ldd $binary

for inputnum in $input
do
    export OMP_NUM_THREADS=1
    tmpi=$((num/1/space))
    tnode=$((tmpi/SLURM_NNODES))
    cp ${SLURM_SUBMIT_DIR}/${inputnum}* .

	mkdir -p $work
	rm -rf $work/*
	cd $work
	
	mpi="srun"
	opts=" -O -i ../${inputnum}-sander.mdin -c ../${inputnum}.inpcrd -p ../${inputnum}.prmtop -o ${inputnum}.output "

	$mpi $binary $opts &
	pid=$!
        (sleep ${runtime}h; kill -s 9 $pid) &
        pid1=$!
	wait $pid
	kill -s 9 $pid1
	
        grep -q Performance ${inputnum}.output
        if [ $? == 0 ]
        then
   	   mkdir -p $root/${SLURM_NNODES}/sander/${inputnum}
   	   cp -a ./* $root/${SLURM_NNODES}/sander/${inputnum}
        fi
	
	cd $root
done

input="trx"

binary=sander.MPI
ldd $binary

for inputnum in $input
do
    export OMP_NUM_THREADS=1
    tmpi=$((num/1/space))
    tnode=$((tmpi/SLURM_NNODES))
    cp ${SLURM_SUBMIT_DIR}/${inputnum}* .

	mkdir -p $work
	rm -rf $work/*
	cd $work
	
	mpi="srun"
	opts=" -O -i ../${inputnum}.mdin -c ../${inputnum}.x -idip ../${inputnum}.dip -p ../${inputnum}.prmtop -o ${inputnum}.output "

	$mpi $binary $opts &
	pid=$!
        (sleep ${runtime}h; kill -s 9 $pid) &
        pid1=$!
	wait $pid
	kill -s 9 $pid1
	
        grep -q Performance ${inputnum}.output
        if [ $? == 0 ]
        then
   	   mkdir -p $root/${SLURM_NNODES}/sander/${inputnum}
   	   cp -a ./* $root/${SLURM_NNODES}/sander/${inputnum}
        fi
	
	cd $root
done

rm -rf $work
