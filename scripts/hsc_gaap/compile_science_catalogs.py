import os, sys
import glob
import numpy as np
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from astropy.table import Table, vstack, hstack
import lsst.daf.butler as dafButler
from hsc_gaap.gaap import consolidateMerianCats
from hsc_gaap.find_patches_to_reduce import findReducedPatches

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

if __name__ == '__main__':
    fire.Fire(merge_merian_catalogs)
    fire.Fire(select_unique_objs)
	
	

