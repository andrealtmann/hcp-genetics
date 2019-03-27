#! /bin/bash

#SGE_parms
#$ -cwd
#$ -S /bin/bash
#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=23:59:59

#pgs="education"
pgs=$1
nparc=$2

R --no-save -q < cluster_call.R --pgs=$pgs --nparc=$nparc --nperm=1000


