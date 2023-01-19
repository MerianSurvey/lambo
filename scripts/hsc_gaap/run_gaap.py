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
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
import multiprocess as mp
mp.freeze_support()
from astropy.coordinates import SkyCoord
import astropy.units as u

from hsc_gaap.gaap import GaapTask, NaiveLogger, findReducedPatches
import lsst.daf.butler as dafButler


def runGaap(patch, tract=9813, bands='gri', hsc_type='w_2022_40', logger=None, filter_jobs=None, repo='/projects/MERIAN/repo/'):
    old_patches = [name for name in os.listdir(
        f"{repo}/S20A/deepCoadd_calexp/{tract}/")]
    new_patches = [int(name[0]) + int(name[2]) * 9 for name in old_patches]

    print("Starting butler")
    merian_butler = dafButler.Butler("/projects/MERIAN/repo")
    merian_patches = findReducedPatches(tract)
    common_patches = np.intersect1d(new_patches, merian_patches)
    logger.info(
        f'In tract = {tract}, there are {len(common_patches)} common patches')
        
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
                            repo=repo,
                            collections='S20A/deepCoadd_calexp',
                            logger=logger, merian_butler=merian_butler)
        elif hsc_type == 'w_2022_40':
            gaap = GaapTask(tract, patch, bands, hsc_type='w_2022_40',
                            repo=repo,
                            collections='HSC/runs/RC2/w_2022_40/DM-36151',
                            logger=logger)
        elif hsc_type == 'w_2022_04':
            gaap = GaapTask(tract, patch, bands, hsc_type='w_2022_04',
                            repo=repo,
                            collections='HSC/runs/RC2/w_2022_04/DM-33402',
                            logger=logger)

        gaap.load_merian_reference(band='N708',
                                   repo='/projects/MERIAN/repo/',
                                   collections='DECam/runs/merian/dr1_wide'
                                   )
        del merian_butler
        gaap.setDefaultMeasureConfig()

        if filter_jobs is not None:
            filter_pool = mp.Pool(min(len(gaap.bands), filter_jobs))
        else:
            filter_pool=None

        gaap.runAll(filter_pool)
        print(f'BANDS WRITING TO CATALOG: {gaap.bands}')
        gaap.writeObjectTable()
        gaap.transformObjectCatalog(
            functorFile=os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/hsc_gaap/Object.yaml'))
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


def runGaapRowColumn(tract, patch_cols, patch_rows, bands='grizy', patch_jobs=5, filter_jobs=None, hsc_type='S20A', repo='/projects/MERIAN/repo/'):
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
    if os.path.isdir(f'./log/{tract}') is False:
        os.makedirs(f'./log/{tract}')
    if os.path.isfile(f'./log/{tract}//{patch_cols[0] + patch_rows[0] * 9}.log'):
        os.system(f'rm ./log/{tract}//{patch_cols[0] + patch_rows[0] * 9}.log')
    logger = NaiveLogger(f'./log/{tract}/{patch_cols[0] + patch_rows[0] * 9}.log')
    patches_old = list(product(patch_cols, patch_rows))
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

        runGaap(patch, tract, bands=bands,
                hsc_type=hsc_type, logger=logger, filter_jobs=filter_jobs, repo=repo)

        matchBlendedness(tract, patch_cols, patch_rows, repo=repo)

        gc.collect()


def matchBlendedness(tract, patch_cols, patch_rows, repo='/projects/MERIAN/repo/', matchdist=0.5, matchmag=10):
    """
    Crossmatch gaap catalog with S20A blendedness catalog.
    Matches objects that are separated by less than matchdist arcseconds and 
    have 3-pixel-aperture g-band flux with at most matchmag nJy disagreement.

     Parameters
    ----------
    matchdist : float
        Maximum matching distance in arcseconds
    matchmag : float.
        Maximum 3-pixel-aperture flux difference for matching.

    """
    patches = list(product(patch_cols, patch_rows))

    for patch in patches:


        gaapCat_dir = os.path.join(
            repo, "S20A/gaapTable/", f"{tract}", f"{patch[0]},{patch[1]}")
        gaapCat_file = f"objectTable_{tract}_{patch[0]},{patch[1]}_S20A.fits"
        blendCat_file = f"S20A_blendedness_{patch[0]},{patch[1]}.fits"
        
        if not os.path.isfile(os.path.join(gaapCat_dir, gaapCat_file)):
            raise ValueError(f"No gaap catalog at {os.path.join(gaapCat_dir, gaapCat_file)}")

        if not os.path.isfile(os.path.join(gaapCat_dir, gaapCat_file)):
             raise ValueError(f"No blendedness catalog at {os.path.join(gaapCat_dir, blendCat_file)}")

        gaapCat = Table.read(os.path.join(gaapCat_dir, gaapCat_file))
        blendCat = Table.read(os.path.join(
            gaapCat_dir, blendCat_file))

        old_cols = [f'm_{filt}_blendedness_abs' for filt in 'grizy'] + \
            [f'm_{filt}_blendedness_flag' for filt in 'grizy']
        old_cols += [f'{filt}_apertureflux_10_flux' for filt in 'grizy'] + \
            [f'{filt}_apertureflux_10_flag' for filt in 'grizy']
        new_cols = [f'{filt}_blendedness' for filt in 'grizy'] + \
            [f'{filt}_blendedness_flag' for filt in 'grizy']
        new_cols += [f'{filt}_apertureflux10_S20A' for filt in 'grizy'] + \
            [f'{filt}_apertureflux10_flag_S20A' for filt in 'grizy']
        blendCat.rename_columns(old_cols, new_cols)

        _gaap = SkyCoord(gaapCat['coord_ra'], gaapCat['coord_dec'], unit='deg')
        _blend = SkyCoord(blendCat['ra'], blendCat['dec'], unit='deg')
        ind, dist, _ = _gaap.match_to_catalog_sky(_blend)

        match_flux_col = np.array(gaapCat.colnames)[["ap03Flux" in col for col in gaapCat.colnames]][0]
    
        match = (dist < .5 * u.arcsec) & (
            abs(gaapCat[match_flux_col] - blendCat["g_apertureflux10_S20A"][ind]) < 10)

        for col in new_cols:
            gaapCat[col] = np.ones(len(gaapCat)) * np.nan
            gaapCat[col][match] = blendCat[ind[match]][col]

        gaapCat.write(os.path.join(gaapCat_dir, gaapCat_file), overwrite=True)
        print(f'Matched GAaP table with S20A blendedness based on {match_flux_col[0]} band.')



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
# python run_gaap.py --patch_cols=[2,3,4] --patch_rows=[4] --patch_jobs=None --filter_jobs=5 --hsc_type="S20A" --tract=9813
