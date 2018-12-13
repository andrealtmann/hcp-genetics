#! /bin/bash

#SGE_parms
#$ -cwd
#$ -S /bin/bash
#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=23:59:59

#pgs="education"
pgs=$1
cvmeth="famaware"
#cvmeth="standard"

R --no-save -q < cluster_call_benchmark.R --pgs=$pgs --nparc=300 --cv=$cvmeth


