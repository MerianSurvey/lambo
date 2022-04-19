#!/bin/bash
#SBATCH -J trctr_job2_500_1000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --mem-per-cpu=2G
#SBATCH --time=11:59:00
#SBATCH --export=ALL
#SBATCH -o trctr_job2_500_1000.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"

module purge
. /home/jiaxuanl/Research/Merian/setup_env.sh
export OMP_NUM_THREADS=1

python tractor_fitting.py /scratch/gpfs/jiaxuanl/Data/Merian/Catalogs/COSMOS_cutouts_tractor_gaap_lowz_calexp.fits --suffix lowz_calexp --njobs 54 \
--low 500 --high 1000 --ind_list None --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
--CUTOUT_SUBDIR ./Cutout/_new/ --CATALOG_SUBDIR ./Catalogs/ \



now=$(date +"%T")
echo "end time ... $now"
