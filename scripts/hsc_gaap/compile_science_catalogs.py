import os, sys
import glob
import numpy as np
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from astropy.table import Table, vstack, hstack, Column
from astropy.coordinates import SkyCoord
import astropy.units as u
import lsst.daf.butler as dafButler
import numpy.ma as ma
#from gaap import consolidateMerianCats
#from find_patches_to_reduce import findReducedPatches
#from hsc_gaap.gaap import findReducedPatches, consolidateMerianCats
import time
import re
from rtree import index
import pandas as pd
import time

#from hsc_gaap.gaap import consolidateMerianCats
#from hsc_gaap.find_patches_to_reduce import findReducedPatches
from gaap import consolidateMerianCats
from find_patches_to_reduce import findReducedPatches

from unagi import hsc

import fire

keep_merian_file = 'keep_table_columns_merian.txt'
keep_gaap_file = 'keep_table_columns_gaap.txt'
repo = '/scratch/gpfs/am2907/Merian/gaap/'
repo_out = '/scratch/gpfs/sd8758/merian/catalog/'
s20a = hsc.Hsc(dr='dr3', rerun='s20a_wide')

def load_param_file(filename='cat.param'):
	dat = pd.read_csv(filename, delimiter='\t', header=0)
	
	return dat

def merge_merian_catalogs(tracts=[], repo=repo, alltracts=False, rewrite=False):
    keepColumns_merian = list(np.genfromtxt(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/hsc_gaap/', keep_merian_file), dtype=None, encoding="ascii"))
    keepColumns_gaap = list(np.genfromtxt(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/hsc_gaap/', keep_gaap_file), dtype=None, encoding="ascii"))

    keepColumns_merian_540 = [c for c in keepColumns_merian if "540" in c]
    keepColumns_merian_708 = [c for c in keepColumns_merian if "708" in c]
    keepColumns_merian_all = [c for c in keepColumns_merian if ("708" not in c) & ("540" not in c)]

    if alltracts:
        tracts = glob.glob(os.path.join(repo, "log/*"))
        tracts = [int(t.split("/")[-1]) for t in tracts]
        tracts.sort()

    for tract in tracts:
        catDir = os.path.join(repo, f"S20A/gaapTable/{tract}/")
        outCatFile = os.path.join(catDir, f'objectTable_{tract}_S20A.fits')

        if os.path.isfile(outCatFile):
            if not rewrite:
                print(f"CATALOG FOR TRACT {tract} ALREADY WRITTEN - SKIPPING")
                continue
            else:
                print(f"CATALOG FOR TRACT {tract} ALREADY WRITTEN - REWRITING")
                
        objectTables = np.array(glob.glob(catDir + "*/objectTable*"))
        
        patches_708 = findReducedPatches(tract, band="N708")
        patches_540 = findReducedPatches(tract, band="N540")
        patches = np.unique(np.concatenate([patches_708, patches_540]))

        ot_patches = np.array([ot.split("_")[2] for ot in objectTables])
        ot_patches = np.array([int(otp[0]) + int(otp[2])*9 for otp in ot_patches])
        objectTables, ot_patches = objectTables[ot_patches.argsort()], ot_patches[ot_patches.argsort()]

        if not np.all(patches == ot_patches):
            print(f"Object tables do not match expected patches for tract {tract}")

        if len(ot_patches)==0:
            print(f"No catalogs for tract {tract} - skipping")
            continue

        print(f'COMPILING CATALOG FOR TRACT {tract} WITH {len(ot_patches)} PATCHES')

        tableList = []
        for patchno, patchpath in zip(ot_patches, objectTables):
            try:
                if (patchno in patches_708) & (patchno in patches_540):
                    merian = consolidateMerianCats([patchno], tract)[keepColumns_merian]
                elif patchno in patches_708:
                    merian = consolidateMerianCats([patchno], tract)[keepColumns_merian_all + keepColumns_merian_708]
                elif patchno in patches_540:
                    merian = consolidateMerianCats([patchno], tract)[keepColumns_merian_all + keepColumns_merian_540]
            except:
                print(f"NO MERIAN OBJECT TABLES FOR TRACT {tract} PATCH {patchno} - SKIPPING")
                continue
            try:
                gaapTable = Table.read(patchpath)[keepColumns_gaap]
            except:
                print(f"GAAP COLUMN EXCEPTION FOR {tract}")
                gaapTable = Table.read(patchpath)
                keepColumns_exception = list(set(keepColumns_gaap).intersection(set(gaapTable.colnames)))
                gaapTable = gaapTable[keepColumns_exception]

            assert np.all(merian['objectId'] == gaapTable['objectId']) , "GAaP table does not match Merian table!"
            gaapTable.remove_columns(["objectId", "coord_ra", "coord_dec", "ebv"])
            
            stack_table = hstack([merian, gaapTable])
            tableList.append(stack_table)

        if len (tableList) == 0:
            continue

        fullTable = vstack(tableList)
        print(f"COMPILED TABLE OF {len(fullTable)} ROWS and {len(fullTable.colnames)} COLUMNS")
        fullTable.write(outCatFile, overwrite=True)
        print(f"WROTE TABLE TO {outCatFile}")

def select_unique_objs(tracts, repo=repo, alltracts=False):

	for tract in tracts:
		catDir = os.path.join(repo, f"S20A/gaapTable/{tract}/")
		outCatFile = os.path.join(catDir, f'objectTable_{tract}_S20A.fits')

		tractTable = Table.read(outCatFile)	
		unique_flag = (tractTable['detect_isPrimary']==True) #Select unique objects using detect_isPrimary
		tablePrimary = tractTable[unique_flag]
		print(f"COMPILED TABLE OF UNIQUE SCIENCE OBJECTS WITH {len(tablePrimary)} ROWS and {len(tablePrimary.colnames)} COLUMNS")
			
		outCatDir = os.path.join(repo_out, f"S20A/{tract}/")
		isExist = os.path.exists(outCatDir)
		if not isExist:
			os.makedirs(outCatDir)

		outSciCatFile = os.path.join(outCatDir, f'meriandr1_unique_{tract}_S20A.fits') 
		tablePrimary.write(outSciCatFile, overwrite=True)


def parse_circle(line, tract_center):
    tract_cen_ra, tract_cen_dec = tract_center
    
    # Extract the coordinates and radius from a circle line
    match = re.search(r'circle\(([-+]?[0-9]*\.?[0-9]+), ([-+]?[0-9]*\.?[0-9]+), ([-+]?[0-9]*\.?[0-9]+)d\)', line)
    if match:
        ra = float(match.group(1))
        dec = float(match.group(2))
        radius = float(match.group(3))
        
        distance = ((ra - tract_cen_ra) ** 2 + (dec - tract_cen_dec) ** 2) ** 0.5
        if distance <= 5:
            return ra, dec, radius
    return None


def check_object_in_circles(objects, circles):
    objects_within_circles = []
    p = index.Property()
    p.dimension = 2
    idx = index.Index(properties=p)
    
    # Build R-tree index for circles
    for i, circle in enumerate(circles):
        circle_ra, circle_dec, circle_radius = circle
        left = circle_ra - circle_radius
        bottom = circle_dec - circle_radius
        right = circle_ra + circle_radius
        top = circle_dec + circle_radius
        idx.insert(i, (left, bottom, right, top))
        
    # Check objects within circles
    for j, obj in enumerate(objects):
        obj_ra, obj_dec = obj['ra'], obj['dec']
        result = list(idx.intersection((obj_ra, obj_dec, obj_ra, obj_dec)))

            circle_ra, circle_dec, circle_radius = circle
            distance = ((obj_ra - circle_ra) ** 2 + (obj_dec - circle_dec) ** 2) ** 0.5
            if distance <= circle_radius:
                objects_within_circles.append({'obj': obj, 'index': j})
                break
    return objects_within_circles


def apply_bright_star_mask(tracts, repo=repo_out, alltracts=False):

	maskDir = '/projects/MERIAN/starmask_s20a/'
	maskFile = os.path.join(maskDir, f'updated_S20A_mask.reg')
	for tract in tracts:
		catDir = os.path.join(repo, f"S20A/{tract}/")
		catFile = os.path.join(catDir, f'meriandr1_unique_{tract}_S20A.fits')


		# Read catalog
		tractTable = Table.read(catFile)
		objects = []
		for ra, dec in zip(tractTable['coord_ra'], tractTable['coord_dec']):
    			objects.append({'ra': ra, 'dec': dec})
		# Get tract ~center
		objects_ra = list(map(lambda x: x['ra'], objects))
		objects_dec = list(map(lambda x: x['dec'], objects))
		tract_center = (np.median(objects_ra), np.median(objects_dec))
		
		# Read the circles from the file, choose only circles at the vicinity of the tract
		circles = []
		with open(maskFile, 'r') as file:
    			for line in file:
        			if line.startswith('circle'):
            				circle = parse_circle(line, tract_center)
            				if circle:
                				circles.append(circle)
	

		# Check objects within circles
		start_time = time.time()
		objects_within_circles = check_object_in_circles(objects, circles)
		elapsed_time = time.time() - start_time
		print(f"Elapsed time: {elapsed_time} seconds")

		# Add a new column (IsMask) to table with 0 for not masked by a bright star and 1 for being mask 
		indices = list(map(lambda x: x['index'], objects_within_circles))
		objects_list_len = len(objects)
		IsMask = np.zeros(objects_list_len)
		IsMask[indices] = 1
		tractTable['IsMask'] = IsMask

		# Write to a new fits file
		outCatDir = os.path.join(repo_out, f"S20A/{tract}/")
		outCatFile = os.path.join(outCatDir, f'meriandr1_mask__{tract}_S20A.fits')
		tractTable.write(outCatFile, overwrite=True)
		

def download_s20a(tracts, save=False, outdir='/projects/MERIAN/repo/'):

	sql_test = open(
        os.path.join(os.getenv("LAMBO_HOME"), 'lambo/scripts/hsc_gaap/HSC_S20A_tract.SQL'), 'r').read()

	for tract in tracts:

		outCatDir = os.path.join(repo_out, f"S20A/{tract}/")		

		sql_test = sql_test.replace('8524', str(tract))	
		print(
		f'# SQL QUERY S20A FROM HSC DATABASE (s20a_wide.summary) FOR TRACT = {tract}')

		result_test = s20a.sql_query(sql_test, from_file=False, preview=False, verbose=True)
		
		outCatFile = os.path.join(outCatDir, f'hsc_gaap_{tract}_S20A.fits')
		result_test.write(outCatFile, overwrite=True)


def merge_merian_hscS20A(tracts):
	"""
	Crossmatch Merian gaap catalog with S20A gaap catalog.
	Matched objects are selected by being separated by less than matchdist arcseconds.
	"""
	
	# read the matching radius from the cat.param file
	df_param = load_param_file(filename='cat.param')
	matchdist = df_param[df_param['parameter']=='matchdist']['value'].values[0] 

	for tract in tracts:
		catDir = os.path.join(repo_out, f"S20A/{tract}/")
		hscGaap_file = f'hsc_gaap_{tract}_S20A.fits'
		merianGaap_file = f'meriandr1_mask__{tract}_S20A.fits'
	
		if not os.path.isfile(os.path.join(catDir, hscGaap_file)):
			raise ValueError(f"No HSC S20A gaap catalog at {os.path.join(catDir, hscGaap_file)}")

		if not os.path.isfile(os.path.join(catDir, merianGaap_file)):
			raise ValueError(f"No HSC Merian gaap catalog at {os.path.join(catDir, merianGaap_file)}")

		hscCat = Table.read(os.path.join(catDir, hscGaap_file))
		# add a HSCS20A prefix to hsc catalog column names
		hsc_cols_old = hscCat.colnames
		hsc_cols_new = [x+'_HSCS20A' for x in hsc_cols_old]
		hscCat.rename_columns(hsc_cols_old, hsc_cols_new)


		merianCat = Table.read(os.path.join(catDir, merianGaap_file))
		# add a Merian prefix to merian catalog column names
		merian_cols_old = merianCat.colnames
		merian_cols_new = [x+'_Merian' for x in merian_cols_old]
		merianCat.rename_columns(merian_cols_old, merian_cols_new)


		_hsc = SkyCoord(hscCat['ra_HSCS20A'], hscCat['dec_HSCS20A'], unit='deg')
		_merian = SkyCoord(merianCat['coord_ra_Merian'], merianCat['coord_dec_Merian'], unit='deg')
		idx, d2d, _ = _merian.match_to_catalog_sky(_hsc)

		match = (d2d < matchdist * u.arcsec)

		combined_table = merianCat
	
		for i in range(len(hscCat.colnames)):
			current_col = hscCat.colnames[i]
			c = Column(length=len(merianCat), name=hscCat.colnames[i], dtype=hscCat[current_col].dtype)
			combined_table.add_column(c)
		


		hsc_match = np.zeros(len(combined_table))
		for i in range(len(match)):
			if match[i] == True:
				hsc_match[i] = 1
				hsc_idx = idx[i]
				for j in range(len(hscCat.colnames)):
					col = hscCat.colnames[j]
					hsc_col_value = hscCat[col][hsc_idx]
					if ma.is_masked(hsc_col_value) == False:
						combined_table[col][i] = hsc_col_value
		combined_table['hsc_match'] = hsc_match
	
		print(type(combined_table['object_id_HSCS20A'][0]))
		print(type(combined_table['z_apertureflux_10_flux_HSCS20A'][0]))
	
		# Write to a new fits file
		outCatDir = os.path.join(repo_out, f"S20A/{tract}/")
		outCatFile = os.path.join(outCatDir, f'meriandr1_hscmerged_{tract}_S20A.fits')
		combined_table.write(outCatFile, overwrite=True)
		print(f"MERGED HSC S20A + MERIAN TABLE OF SCIENCE OBJECTS WITH {len(combined_table)} ROWS and {len(combined_table.colnames)} COLUMNS")


if __name__ == '__main__':
    start_time = time.time()
#    fire.Fire(merge_merian_catalogs)

    fire.Fire(select_unique_objs)
    fire.Fire(apply_bright_star_mask)
    fire.Fire(download_s20a)
    fire.Fire(merge_merian_hscS20A)

    print("--- %s seconds ---" % (time.time() - start_time))
