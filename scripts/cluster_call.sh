#! /bin/bash

#SGE_parms
#$ -cwd
#$ -S /bin/bash
#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=23:59:59

#pgs="education"
pgs=$1

R --no-save -q < cluster_call.R --pgs=$pgs --nparc=300 --nperm=250


