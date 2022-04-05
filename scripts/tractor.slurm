#!/bin/bash
#SBATCH -J tractor_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=8G
#SBATCH --time=11:59:00
#SBATCH --export=ALL
#SBATCH -o tractor_test.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"

module purge
. /home/jiaxuanl/Research/Packages/kuaizi/diezi/setup_env.sh
export OMP_NUM_THREADS=1

python3 tractor_fitting.py

now=$(date +"%T")
echo "end time ... $now"
