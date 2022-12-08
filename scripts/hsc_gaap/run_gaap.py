import numpy as np
from astropy.table import QTable, Table, hstack, vstack
import lsst.pex.config
import lsst.afw.table
import lsst.meas.algorithms
import lsst.pex.exceptions
import lsst.meas.extensions.gaap
from functools import partial
from itertools import product
import fire
import sys
import os
import gc
sys.path.append('/home/jiaxuanl/Research/Merian/merian_tractor/scripts/')
import multiprocess as mp
mp.freeze_support()

old_patches = [name for name in os.listdir(
    "/projects/MERIAN/repo/S20A/deepCoadd_calexp/9813/")]
new_patches = [int(name[0]) + int(name[2]) * 9 for name in old_patches]
merian_patches = [int(name) for name in os.listdir(
    "/projects/MERIAN/repo/DECam/runs/merian/dr1_wide/20220921T193246Z/deepCoadd_forced_src/9813")]
common_patches = np.intersect1d(new_patches, merian_patches)

from hsc_gaap.gaap import GaapTask


def runGaap(patch, bands='gri', hsc_type='w_2022_40'):
    try:
        patch
        assert patch in common_patches, "Patch not in common patches"
        print('### Processing patch =', patch, ', bands =', bands)
        if hsc_type == 'S20A':
            gaap = GaapTask(9813, patch, bands, hsc_type='S20A',
                            repo='/projects/MERIAN/repo/',
                            collections='S20A/deepCoadd_calexp')
            gaap._checkHSCfile()
        elif hsc_type == 'w_2022_40':
            gaap = GaapTask(9813, patch, bands, hsc_type='w_2022_40',
                            repo='/projects/HSC/repo/main',
                            collections='HSC/runs/RC2/w_2022_40/DM-36151')
        elif hsc_type == 'w_2022_04':
            gaap = GaapTask(9813, patch, bands, hsc_type='w_2022_04',
                            repo='/projects/MERIAN/repo/',
                            collections='HSC/runs/RC2/w_2022_04/DM-33402')

        gaap.load_merian_reference(band='N708',
                                   repo='/projects/MERIAN/repo/',
                                   collections='DECam/runs/merian/dr1_wide'
                                   )
        gaap.setDefaultMeasureConfig()
        gaap.runAll()
        gaap.writeObjectTable()
        gaap.transformObjectCatalog(
            functorFile='/home/jiaxuanl/Research/Merian/merian_tractor/scripts/hsc_gaap/Object.yaml')
        gaap.saveObjectTable()

        del gaap
        gc.collect()
        print('\n')
    except Exception as e:
        print(e)


def runGaapMultiJobs(patch_low, patch_high, bands='gri', njobs=4, hsc_type='w_2022_40'):
    patches = np.arange(patch_low, patch_high + 1, 1)
    iterables = product(patches, list(bands))
    pool = mp.Pool(njobs)
    pool.map(partial(runGaap, hsc_type=hsc_type), iterables)
    pool.close()
    pool.join()


def runGaapRowColumn(patch_cols, patch_rows, bands='grizy', njobs=4, hsc_type='S20A'):
    """
    Parameters
    ----------
    patch_cols : list. 
        The column numbers of the patches to be processed, following the old patch pattern.
        patch_new = int(patch_old[0]) + int(patch_old[2]) * 9
    patch_rows : list.
        The row numbers of the patches to be processed, following the old patch pattern.

    """
    import itertools
    patches_old = list(itertools.product(patch_cols, patch_rows))
    patches = [item[0] + item[1] * 9 for item in patches_old]
    pool = mp.Pool(njobs)
    pool.map(partial(runGaap, bands=bands, hsc_type=hsc_type), patches)
    pool.close()
    pool.join()


if __name__ == '__main__':
    fire.Fire(runGaapRowColumn)


########## S20A ##########
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[0,1] --njobs=6 --hsc_type="S20A" # done in gaap1
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[2,3] --njobs=6 --hsc_type="S20A" # gaap2
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[4,5] --njobs=6 --hsc_type="S20A" # gaap3
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[6,7] --njobs=6 --hsc_type="S20A" # gaap4
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[8] --njobs=2 --hsc_type="S20A" # not done yet


########## w_2022_40 ##########
# python run_gaap.py --patch_cols=[3,4,5,6,7] --patch_rows=[2] --njobs=2 --hsc_type="w_2022_40" --bands='ri' # gaap_w40
