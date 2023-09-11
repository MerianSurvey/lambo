import os, sys
import glob
import numpy as np
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from astropy.table import Table, vstack, hstack, Column
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.time import Time

import time
import pandas as pd

#from gaap import consolidateMerianCats
#from find_patches_to_reduce import findReducedPatches

import fire

cat_dir = '/scratch/gpfs/sd8758/merian/catalog/S20A/'

def get_summary_cat():
	
	merian_id = []
	merian_ra = []
	merian_dec = []
	hsc_id = []
	hsc_ra = []
	hsc_dec = []
	tracts = []

	print(len(glob.glob(cat_dir + '*/*use*')))
	cat_array = glob.glob(cat_dir + '*/*use*')

	for f in cat_array:
		tract_num = f.split('_')[2]
		print('Adding tract #' + str(tract_num) + ' to the summary catalog')
		# read a single catalog
		dat = Table.read(f)
		df = dat.to_pandas()
		print('This tract has ' + str(len(df)) + ' sources')
		print('') 
		tracts = tracts + [int(tract_num)] * len(df)
		mer_id_single = df['objectId_Merian'].values.tolist()
		merian_id = merian_id + mer_id_single 
		mer_ra_single = df['coord_ra_Merian'].values.tolist()
		merian_ra = merian_ra + mer_ra_single
		mer_dec_single = df['coord_dec_Merian'].values.tolist()
		merian_dec = merian_dec + mer_dec_single
		hsc_id_single = df['object_id_HSCS20A'].values.tolist()
		hsc_id = hsc_id + hsc_id_single
		hsc_ra_single = df['ra_HSCS20A'].values.tolist()
		hsc_ra = hsc_ra + hsc_ra_single
		hsc_dec_single = df['dec_HSCS20A'].values.tolist()
		hsc_dec = hsc_dec + hsc_dec_single
		
	# generate a summary table
	t = Table()
	t['objectId_Merian'] = merian_id
	t['coord_ra_Merian'] = merian_ra
	t['coord_dec_Merian'] = merian_dec
	t['object_id_HSCS20A'] = hsc_id
	t['ra_HSCS20A'] = hsc_ra
	t['dec_HSCS20A'] = hsc_dec
	t['tract'] = tracts

	t.write(cat_dir + 'meriandr1_summary_catalog.fits', overwrite=True)

	return None

if __name__ == '__main__':
    start_time = time.time()

    fire.Fire(get_summary_cat)

    print("--- %s seconds ---" % (time.time() - start_time))

