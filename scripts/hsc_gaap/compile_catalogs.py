import os, sys
import glob
import numpy as np
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from astropy.table import Table, vstack, hstack
import lsst.daf.butler as dafButler
from hsc_gaap.gaap import findReducedPatches, consolidateMerianCats

import fire

def compileCatalogs(tracts, repo= "/scratch/gpfs/am2907/Merian/gaap/", alltracts=False):
    keepColumns_merian = list(np.load(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/hsc_gaap/', "keep_table_columns_merian.npy")))
    keepColumns_gaap = list(np.load(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/hsc_gaap/', "keep_table_columns_gaap.npy")))

    for tract in tracts:
        patches = findReducedPatches(tract)
        patches.sort()
        catDir = os.path.join(repo, f"S20A/gaapTable/{tract}/")
        objectTables = np.array(glob.glob(catDir + "*/objectTable*"))

        ot_patches = np.array([ot.split("_")[2] for ot in objectTables])
        ot_patches = np.array([int(otp[0]) + int(otp[2])*9 for otp in ot_patches])
        objectTables, ot_patches = objectTables[ot_patches.argsort()], ot_patches[ot_patches.argsort()]

        assert np.all(patches == ot_patches), "Object tables do not match expected patches"

        print(f'COMPILING CATALOG FOR TRACT {tract} WITH {len(patches)} PATCHES')

        tableList = []
        for patchno, patchpath in zip(ot_patches, objectTables):
            merian = consolidateMerianCats([patchno], tract)[keepColumns_merian]
            gaapTable = Table.read(patchpath)[keepColumns_gaap]   

            assert np.all(merian['objectId'] == gaapTable['objectId']) , "GAaP table does not match Merian table!"
            gaapTable.remove_columns(["objectId", "coord_ra", "coord_dec", "ebv"])
            tableList.append(hstack([merian, gaapTable]))

        fullTable = vstack(tableList)
        print(f"COMPILED TABLE OF {len(fullTable)} ROWS and {len(fullTable.colnames)} COLUMNS")
        outCatFile = os.path.join(catDir, f'objectTable_{tract}_S20A.fits')
        fullTable.write(outCatFile, overwrite=True)
        print(f"WROTE TABLE TO {outCatFile}")

if __name__ == '__main__':
    fire.Fire(compileCatalogs)
