#!/bin/bash
#SBATCH -J trctr_job_0_10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=11:59:00
#SBATCH --export=ALL
#SBATCH -o trctr_job_0_10.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"

module purge
. /home/jiaxuanl/Research/Merian/setup_env.sh
export OMP_NUM_THREADS=1

python tractor_fitting.py /scratch/gpfs/jiaxuanl/Data/Merian/Catalogs/magellan/magellan_spec_obj_cat.fits --suffix magellan --njobs 3 \
--low 0 --high 10 --ind_list None --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
--CUTOUT_SUBDIR ./Cutout/magellan_spec/ --CATALOG_SUBDIR ./Catalogs/magellan/ \



now=$(date +"%T")
echo "end time ... $now"
