#!/bin/bash
#SBATCH -J gaap23
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=5:30:00
#SBATCH --export=ALL

#SBATCH -o ./log/slurm_gaap23.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"

module purge
. /home/jiaxuanl/Research/Merian/merian_tractor/scripts/setup_env_w40.sh

nice -n 0 python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[8] --patch_jobs=4 --hsc_type="S20A"

now=$(date +"%T")
echo "end time ... $now"
