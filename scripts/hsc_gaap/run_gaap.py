import numpy as np
import time
import matplotlib.pyplot as plt
import lsst.meas.base
import lsst.pex.config

import lsst.afw.table
import lsst.meas.algorithms
import lsst.pex.exceptions
import lsst.meas.extensions.gaap

from astropy.io import fits
import sys
import os
import gc
sys.path.append('/home/jiaxuanl/Research/Merian/merian_tractor/scripts/')

old_patches = [name for name in os.listdir(
    "/projects/MERIAN/repo/S20A/deepCoadd_calexp/9813/")]
new_patches = [int(name[0]) + int(name[2]) * 9 for name in old_patches]
merian_patches = [int(name) for name in os.listdir(
    "/projects/MERIAN/repo/DECam/runs/merian/dr1_wide/20220921T193246Z/deepCoadd_forced_src/9813")]
common_patches = np.intersect1d(new_patches, merian_patches)

from hsc_gaap.gaap import GaapTask

for patch in common_patches[:1]:
    for filt in list('grizy'):
        print('### Processing patch', patch, 'band', filt)
        gaap = GaapTask(9813, patch, filt,
                        repo='/projects/MERIAN/repo/',
                        collections='S20A/deepCoadd_calexp')
        gaap._checkHSCfile()
        gaap.load_merian_reference(band='N708',
                                   repo='/projects/MERIAN/repo/',
                                   collections='DECam/runs/merian/dr1_wide',
                                   range=None  # (0, 9020)
                                   )
        gaap.setDefaultMeasureConfig()
        gaap.run()
        outcat = gaap.writeObjectTable()
        del gaap
        gc.collect()
