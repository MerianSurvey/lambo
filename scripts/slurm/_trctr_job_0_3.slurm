#!/bin/bash
#SBATCH -J trctr_job_0_3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=11:59:00
#SBATCH --export=ALL
#SBATCH -o trctr_job_0_3.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"

module purge
. /home/jiaxuanl/Research/Merian/setup_env.sh
export OMP_NUM_THREADS=1

python tractor_fitting.py /scratch/gpfs/jiaxuanl/Data/Merian/Cutout/stars/stars-2022-04-25.fits --suffix stars --njobs 3 --point_source True \
--low 0 --high 3 --ind_list None --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
--CUTOUT_SUBDIR ./Cutout/stars/ --CATALOG_SUBDIR ./Catalogs/stars/ \



now=$(date +"%T")
echo "end time ... $now"
