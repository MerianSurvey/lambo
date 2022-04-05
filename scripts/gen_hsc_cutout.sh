#!/bin/sh

. /home/jiaxuanl/Research/Packages/kuaizi/diezi/setup_env.sh

export DATADIR="/scratch/gpfs/jiaxuanl/Data/" # Directory of all data
export LSBGDIR="/scratch/gpfs/jiaxuanl/Data/Merian/"

# everything should be downloaded to /scratch/gpfs/jiaxuanl/Data/HSC/LSBG/Cutout
# cat['radius'] is always the cutout size we should use!!!!

python3 ./s18a_batch_cutout.py \
    $DATADIR\
    $LSBGDIR"Catalogs/COSMOS_cutouts_tractor_gaap_test.fits" \
    --bands grizy --ra_name ra --dec_name dec \
    --name "prefix" --output $LSBGDIR"/Cutout/NSA/z001_002" \
    --catalog_dir $LSBGDIR"/Catalogs" --catalog_suffix "" \
    --size 0.3 --prefix "bug" \
    --njobs 1 --psf