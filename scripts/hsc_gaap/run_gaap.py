import numpy as np
import lsst.pex.config

import lsst.afw.table
import lsst.meas.algorithms
import lsst.pex.exceptions
import lsst.meas.extensions.gaap

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


def runGaap(patch):
    assert patch in common_patches, "Patch not in common patches"
    for filt in ['i']:  # list('igrzy'):
        print('### Processing patch', patch, 'band', filt)
        gaap = GaapTask(9813, patch, filt,
                        repo='/projects/MERIAN/repo/',
                        collections='S20A/deepCoadd_calexp')
        gaap._checkHSCfile()
        gaap.load_merian_reference(band='N540',
                                   repo='/projects/MERIAN/repo/',
                                   collections='DECam/runs/merian/dr1_wide',
                                   #    range=(0, 200)
                                   )
        gaap.setDefaultMeasureConfig()
        gaap.run()
        _ = gaap.writeObjectTable()
        del gaap
        gc.collect()
        print('\n')


def runGaapMultiJobs(seed_low, seed_high, njobs=4):
    patches = common_patches[seed_low:seed_high]
    pool = mp.Pool(njobs)
    pool.map(runGaap, patches)
    pool.close()
    pool.join()


if __name__ == '__main__':
    fire.Fire(runGaapMultiJobs)

# python run_gaap.py --seed_low=20 --seed_high=23 --n_jobs=3
