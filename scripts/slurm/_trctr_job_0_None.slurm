#!/bin/bash
#SBATCH -J trctr_job_0_None
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --mem-per-cpu=2G
#SBATCH --time=03:59:00
#SBATCH --export=ALL
#SBATCH -o trctr_job_0_None.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"

module purge
source /home/jiaxuanl/Research/Merian/merian_tractor/scripts/setup_env.sh
export OMP_NUM_THREADS=1

python tractor_fitting_native.py /scratch/gpfs/jiaxuanl/Data/Merian/Cutout/stars/stars-2022-06-26.fits --suffix stars_native --njobs 54 --point_source True \
--low 0 --high None --ind_list None --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
--CUTOUT_SUBDIR ./Cutout/stars/ --CATALOG_SUBDIR ./Catalogs/stars/ \



now=$(date +"%T")
echo "end time ... $now"
