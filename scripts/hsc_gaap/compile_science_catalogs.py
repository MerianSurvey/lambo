import os, sys
import glob
import numpy as np
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from astropy.table import Table, vstack, hstack
import lsst.daf.butler as dafButler
<<<<<<< HEAD
from gaap import findReducedPatches, consolidateMerianCats
#from hsc_gaap.gaap import findReducedPatches, consolidateMerianCats
import time
import re
from rtree import index

=======
from hsc_gaap.gaap import consolidateMerianCats
from hsc_gaap.find_patches_to_reduce import findReducedPatches
>>>>>>> 2cf4b8335fab198c80db45db8982e3b958684cc5

import fire

keep_merian_file = 'keep_table_columns_merian.txt'
keep_gaap_file = 'keep_table_columns_gaap.txt'
repo = '/scratch/gpfs/am2907/Merian/gaap/'
repo_out = '/scratch/gpfs/sd8758/merian/catalog/'

def merge_merian_catalogs(tracts, repo=repo, alltracts=False, rewrite=False):
    keepColumns_merian = list(np.genfromtxt(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/hsc_gaap/', keep_merian_file), dtype=None, encoding="ascii"))
    keepColumns_gaap = list(np.genfromtxt(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/hsc_gaap/', keep_gaap_file), dtype=None, encoding="ascii"))

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
                print(f"CATALOG FOR TRACT {tract} ALREADY WRITTEN - REWRITTING")
                
        objectTables = np.array(glob.glob(catDir + "*/objectTable*"))
        
        patches = findReducedPatches(tract)
        patches.sort()

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
                merian = consolidateMerianCats([patchno], tract)[keepColumns_merian]
            except:
                try:
                    merian = consolidateMerianCats([patchno], tract)[keepColumns_merian[:15]]
                except:
                    print(f"NO MERIAN OBJECT TABLES FOR TRACT {tract} PATCH {patchno} - SKIPPING")
                    continue
            try:
                gaapTable = Table.read(patchpath)[keepColumns_gaap]
            except:
                gaapTable = Table.read(patchpath)
                keepColumns_exception = list(set(keepColumns_gaap).intersection(set(gaapTable.colnames)))
                gaapTable = gaapTable[keepColumns_exception]

            assert np.all(merian['objectId'] == gaapTable['objectId']) , "GAaP table does not match Merian table!"
            gaapTable.remove_columns(["objectId", "coord_ra", "coord_dec", "ebv"])
            
            stack_table = hstack([merian, gaapTable])
            tableList.append()

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

		outSciCatFile = os.path.join(outCatDir, f'SciObjectTable_{tract}_S20A.fits') 
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
        for i in result:
            circle = circles[i]
            circle_ra, circle_dec, circle_radius = circle
            distance = ((obj_ra - circle_ra) ** 2 + (obj_dec - circle_dec) ** 2) ** 0.5
            if distance <= circle_radius:
                objects_within_circles.append({'obj': obj, 'index': j})
                break
    return objects_within_circles


def apply_bright_star_mask(tracts, repo=repo_out, alltracts=False):

	maskDir = '/projects/MERIAN/starmask_s19a/'
	maskFile = os.path.join(maskDir, f'BrightObjectMask-99999-9,9-HSC-I.reg') # I-band, it this what we want?

	for tract in tracts:
		catDir = os.path.join(repo, f"S20A/{tract}/")
		catFile = os.path.join(catDir, f'SciObjectTable_{tract}_S20A.fits')


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
		outCatFile = os.path.join(outCatDir, f'SciObjectTableIsMask_{tract}_S20A.fits')
		tractTable.write(outCatFile, overwrite=True)
		


if __name__ == '__main__':
    fire.Fire(merge_merian_catalogs)
    fire.Fire(select_unique_objs)
    fire.Fire(apply_bright_star_mask)
	
	

