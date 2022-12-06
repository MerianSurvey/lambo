import numpy as np
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

from hsc_gaap.gaap import GaapTask, joinMerianCatPatches


def runGaap(patch_filt, hsc_type='w_2022_40'):
    patch, filt = patch_filt
    assert patch in common_patches, "Patch not in common patches"
    print('### Processing patch', patch, 'band', filt)
    if hsc_type == 'S20A':
        gaap = GaapTask(9813, patch, filt, hsc_type='S20A',
                        repo='/projects/MERIAN/repo/',
                        collections='S20A/deepCoadd_calexp')
        gaap._checkHSCfile()
    elif hsc_type == 'w_2022_40':
        gaap = GaapTask(9813, patch, filt, hsc_type='w_2022_40',
                        repo='/projects/HSC/repo/main',
                        collections='HSC/runs/RC2/w_2022_40/DM-36151')
    elif hsc_type == 'w_2022_04':
        gaap = GaapTask(9813, patch, filt, hsc_type='w_2022_04',
                        repo='/projects/MERIAN/repo/',
                        collections='HSC/runs/RC2/w_2022_04/DM-33402')

    gaap.load_merian_reference(band='N708',
                               repo='/projects/MERIAN/repo/',
                               collections='DECam/runs/merian/dr1_wide',
                               )
    gaap.setDefaultMeasureConfig()
    gaap.run()
    cat2 = gaap.writeObjectTable(save=True)
    cat_ref = joinMerianCatPatches([23])
    cat_ref = cat_ref[np.in1d(cat_ref['objectId'], cat2['id'])]
    # cat2.remove_columns(['id', 'coord_ra', 'coord_dec'])
    # cat = hstack([cat_ref, cat2])
    cat_ref.write(
        f'/projects/MERIAN/repo/S20A/gaapTable/9813/5,2/MerianTable_{gaap.band.upper()}_{gaap.tract}_{gaap.patch_old}.fits')

    _ = gaap.writeObjectTable()
    del gaap
    gc.collect()
    print('\n')


def runGaapMultiJobs(seed_low, seed_high, bands='ri', njobs=4, hsc_type='w_2022_40'):
    patches = common_patches[seed_low:seed_high]
    iterables = product(patches, list(bands))
    pool = mp.Pool(njobs)
    pool.map(partial(runGaap, hsc_type=hsc_type), iterables)
    pool.close()
    pool.join()


if __name__ == '__main__':
    fire.Fire(runGaapMultiJobs)

# python run_gaap.py --seed_low=20 --seed_high=23 --njobs=9 --hsc_type="S20A"
# python run_gaap.py --seed_low=20 --seed_high=23 --njobs=6 --hsc_type="w_2022_04"
