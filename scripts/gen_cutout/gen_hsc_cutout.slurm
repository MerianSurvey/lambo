#!/bin/bash
#SBATCH -J hsc_cutout
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=11:59:00
#SBATCH --export=ALL
#SBATCH -o hsc_cutout.o
#SBATCH --mail-type=all
#SBATCH --mail-user=jiaxuanl@princeton.edu

now=$(date +"%T")
echo "start time ... $now"

. /home/jiaxuanl/Research/Packages/Merian/setup_env.sh

export OMP_NUM_THREADS=1
export DATADIR="/scratch/gpfs/jiaxuanl/Data/" # Directory of all data
export LSBGDIR="/scratch/gpfs/jiaxuanl/Data/Merian/"

python3 ./s18a_batch_cutout.py \
    $DATADIR\
    $LSBGDIR"/Cutout/stars/stars-2022-04-25_test.fits" \
    --ra_name ra --dec_name dec \
    --low 0 --high None --size half_size --unit arcsec \
    --bands grizy --name name --prefix s18a_wide \
    --njobs 20 --overwrite False


now=$(date +"%T")
echo "end time ... $now"
