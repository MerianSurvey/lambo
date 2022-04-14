#!/bin/bash
#SBATCH -J tractor_b1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4G
#SBATCH --time=11:59:00
#SBATCH --export=ALL
#SBATCH -o tractor_b1.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"

module purge
. /home/jiaxuanl/Research/Merian/setup_env.sh
export OMP_NUM_THREADS=1

python tractor_fitting.py --njobs=42 --low=0 --high=500
# python tractor_fitting.py --njobs=1 --low=96 --high=97

now=$(date +"%T")
echo "end time ... $now"
