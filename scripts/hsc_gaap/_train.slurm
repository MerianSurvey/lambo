#!/bin/bash
#SBATCH -J gaapCosmos_9813_[0]
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=50G
#SBATCH --time=12:30:00
#SBATCH --export=ALL
#SBATCH -o ./log/gaapCosmos_9813_[0].o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"
echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

module purge
. /home/jiaxuanl/Research/Merian/merian_tractor/scripts/setup_env_w40.sh

python run_gaap.py --tract=9813 --patch_cols='[3]' --patch_rows='[0]' --patch_jobs=None --filter_jobs=5 --hsc_type='S20A'


now=$(date +"%T")
echo "end time ... $now"
