import numpy as np
from astropy.table import QTable, Table, hstack, vstack
import lsst.pex.config
import lsst.afw.table
import lsst.meas.algorithms
import lsst.pex.exceptions
import lsst.meas.extensions.gaap
from functools import partial
from itertools import product
import traceback
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

from hsc_gaap.gaap import GaapTask, NaiveLogger
import gc

def runGaap(patch, tract=9813, bands='gri', hsc_type='w_2022_40', logger=None, filter_pool=None):
    try:
        if patch not in common_patches:
            if logger is not None:
                logger.error(f'Patch {patch} not in common patches')
                raise ValueError(f'Patch {patch} not in common patches')

        if logger is not None:
            logger.info(f'### Processing patch = {patch}, bands = {bands}')
        print('### Processing patch =', patch, ', bands =', bands)

        if hsc_type == 'S20A':
            gaap = GaapTask(tract, patch, bands, hsc_type='S20A',
                            repo='/projects/MERIAN/repo/',
                            collections='S20A/deepCoadd_calexp',
                            logger=logger)
            gaap._checkHSCfile()
        elif hsc_type == 'w_2022_40':
            gaap = GaapTask(tract, patch, bands, hsc_type='w_2022_40',
                            repo='/projects/HSC/repo/main',
                            collections='HSC/runs/RC2/w_2022_40/DM-36151',
                            logger=logger)
        elif hsc_type == 'w_2022_04':
            gaap = GaapTask(tract, patch, bands, hsc_type='w_2022_04',
                            repo='/projects/MERIAN/repo/',
                            collections='HSC/runs/RC2/w_2022_04/DM-33402',
                            logger=logger)

        gaap.load_merian_reference(band='N708',
                                   repo='/projects/MERIAN/repo/',
                                   collections='DECam/runs/merian/dr1_wide'
                                   )
        gaap.setDefaultMeasureConfig()
        gaap.runAll(filter_pool)
        gaap.writeObjectTable()
        gaap.transformObjectCatalog(
            functorFile='/home/jiaxuanl/Research/Merian/merian_tractor/scripts/hsc_gaap/Object.yaml')
        gaap.saveObjectTable()

        if logger is not None:
            logger.info(
                f'    - Patch {patch}: Succeeded for tract = {tract}, patch = {patch}, bands = {bands}')
        del gaap
        gc.collect()
        print('\n')
    except Exception as e:
        print('ERROR:', e)
        print(traceback.format_exc())


def runGaapRowColumn(patch_cols, patch_rows, bands='grizy', patch_jobs=5, filter_jobs=None, hsc_type='S20A'):
    """
    This function uses ``multiprocessing`` to run ``runGaap`` in parallel.
    The "parallel" is in the sense that the patches are processed in parallel, not the bands.
    For a given patch, all bands are processed sequentially, which takes ~2.5 hours. 

    Parameters
    ----------
    patch_cols : list. 
        The column numbers of the patches to be processed, following the old patch pattern.
        ``patch_new = int(patch_old[0]) + int(patch_old[2]) * 9``
    patch_rows : list.
        The row numbers of the patches to be processed, following the old patch pattern.

    """
    import itertools
    if os.path.isdir('./log') is False:
        os.mkdir('./log')
    logger = NaiveLogger(f'./log/gaap_{patch_rows}.log')
    patches_old = list(itertools.product(patch_cols, patch_rows))
    patches = [item[0] + item[1] * 9 for item in patches_old]

    # if patch_jobs is not None:
    #     pool = mp.Pool(patch_jobs)
    # else:
    #     pool = None

    # if pool is not None:
    #     pool.map(partial(runGaap, bands=bands,
    #                      hsc_type=hsc_type, logger=logger, filter_pool=filter_pool), patches)
    #     pool.close()
    #     pool.join()
    # else:
    for patch in patches:
        if filter_jobs is not None:
            print('Using filter pool: ', filter_jobs, 'jobs')
            filter_pool = mp.Pool(filter_jobs)
        else:
            filter_pool = None
        runGaap(patch, bands=bands,
                hsc_type=hsc_type, logger=logger, filter_pool=filter_pool)
        gc.collect()
    

if __name__ == '__main__':
    fire.Fire(runGaapRowColumn)


########## S20A ##########
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[0,1] --patch_jobs=6 --hsc_type="S20A" # gaap1
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[2,3] --patch_jobs=6 --hsc_type="S20A" # gaap2
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[4,5] --patch_jobs=6 --hsc_type="S20A" # gaap3
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[6,7] --patch_jobs=6 --hsc_type="S20A" # gaap4
# python run_gaap.py --patch_cols=[0,1,2,3,4,5,6,7,8] --patch_rows=[8] --patch_jobs=6 --hsc_type="S20A" # gaap5

########## w_2022_40 ##########
# python run_gaap.py --patch_cols=[3,4,5,6,7] --patch_rows=[2] --njobs=2 --hsc_type="w_2022_40" --bands='ri' # gaap_w40

########## test ##########
# python run_gaap.py --patch_cols=[2] --patch_rows=[4] --patch_jobs=None --filter_jobs=5 --hsc_type="S20A"
