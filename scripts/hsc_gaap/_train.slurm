#!/bin/bash
#SBATCH -J NDE_gaapCosmos_[4, 8]
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --mem=20G
#SBATCH --time=12:30:00
#SBATCH --export=ALL
#SBATCH -o ./log/gaap_gaapCosmos_[4, 8].o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"
echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

module purge
. /home/jiaxuanl/Research/Merian/merian_tractor/scripts/setup_env_w40.sh

python run_gaap.py --patch_cols='[7]' --patch_rows='[4, 8]' --patch_jobs=None --filter_jobs=2 --hsc_type='S20A'


now=$(date +"%T")
echo "end time ... $now"
