This is Satellite Tobacco Mosaic Virus (STMV) in a TIP3P Water box - 1,067,095
atoms.

Run on multiple CPUs with e.g.:

mpirun -np 4 $AMBERHOME/exe/pmemd -O -i mdin.CPU -o mdout -p prmtop -c inpcrd

Run on a GPU with e.g.:

$AMBERHOME/exe/pmemd.cuda -O -i mdin.GPU -o mdout -p prmtop -c inpcrd

This is a more realistic version of the STMV benchmark run with what would be
typical NPT production parameters.  It uses shake with a 4fs timestep.  It has
a 9 Angstrom cutoff and writes output and trajectory files every 1000 steps
(4ps).  Additionally, it writes the restart file only at the end of this short
run.  For performance reasons, it is set to produce NETCDF binary trajectory
files.
