#! /bin/bash

#SGE_parms
#$ -cwd
#$ -S /bin/bash
#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=9:59:59

#cvmeth="famaware"
#cvmeth="standard"
cvmeth=$1
nparc=$2

#feature to be worked on based on sge task id in array job (1-68)
fsfeat=$SGE_TASK_ID

R --no-save -q < cluster_call_FSbenchmark.R --ICVreg=true --sex=both --CEU=0.0 --nperm=50 --fsfn=$fsfeat --nparc=$nparc --cv=$cvmeth
