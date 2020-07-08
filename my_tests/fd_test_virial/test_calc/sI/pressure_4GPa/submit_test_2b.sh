#!/bin/bash
#$ -pe smp 8 
#$ -l h_rt=24:00:00
#$ -q orinoco
#$ -S /bin/bash
#$ -N fd_virial
#$ -j yes
#$ -cwd

echo "-- New Job --"
export  OMP_NUM_THREADS=${NSLOTS}
source /home/es732/miniconda3/etc/profile.d/conda.sh

conda activate mbx-quippy
python virial_test_2b.py

