#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=4:00:00
#PBS -j oe
#PBS -V
#PBS -e wx_Ncpu.err
#PBS -o wx_Ncpu.out
#PBS -N wx_Ncpu
cd $PBS_O_WORKDIR

JOBID=(`echo $PBS_JOBID | tr '.' ' ' `)
LOGFILE="$JOBID.log"


module load warpx/warpx-21.11-intel2020

mpirun -np 1  warpx  slab_WeiHou.wx3 >> $LOGFILE 2>&1

