#!/bin/bash
#SBATCH -J NDE_gaapCosmos20_25_20_25
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24G
#SBATCH --time=5:30:00
#SBATCH --export=ALL

#SBATCH -o ./log/gaap_gaapCosmos20_25_20_25.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"
echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

module purge
. /home/jiaxuanl/Research/Merian/merian_tractor/scripts/setup_env_gen3.sh

python run_gaap.py --seed_low=20 --seed_high=25 --n_jobs=5


now=$(date +"%T")
echo "end time ... $now"
