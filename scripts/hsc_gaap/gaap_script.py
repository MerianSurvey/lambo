import numpy as np
import time
import matplotlib.pyplot as plt
import lsst.meas.base
import lsst.pex.config
import lsst.afw.display as afwDisplay
import lsst.afw.geom as afwGeom
import lsst.afw.table
import lsst.meas.algorithms
import lsst.meas.deblender
import lsst.pex.exceptions
import lsst.meas.extensions.gaap

# from kuaizi.display import display_single
# from astropy.io import fits
# from astropy.table import Table, hstack


# Load Merian catalog as reference
# filt = 'i'
filt = 'N708'
tract = 9813
patch = 23
patch_old = f'{patch % 9},{patch // 9}'

# repo = '/projects/HSC/repo/main'
# collection = 'HSC/runs/RC2/w_2022_40/DM-36151'
# instr = 'HSC'
repo = '/projects/MERIAN/repo/'
collection = 'DECam/runs/merian/dr1_wide'
instr = 'DECam'

import lsst.daf.butler as dafButler
butler = dafButler.Butler(repo)
dataId = dict(tract=tract, patch=patch, band=filt)

refCat = butler.get(
    'deepCoadd_forced_src',
    collections=collection,
    dataId=dataId,
    instrument=instr,
    skymap='hsc_rings_v1',
)[:10]  # [0:8358:8357]

refExposure = butler.get(
    'deepCoadd_calexp',
    collections=collection,
    dataId=dataId,
    instrument=instr,
    skymap='hsc_rings_v1',
)

refCatInBand = butler.get(
    'deepCoadd_ref',
    collections=collection,
    dataId=dataId,
    instrument=instr,
    skymap='hsc_rings_v1',
)[:10]  # [0:8358:8357]

expID = dict(tract=tract, patch=patch, band=filt, skymap='hsc_rings_v1')


hsc_type = 'S20A'

if hsc_type == 'w40':
    exposure = lsst.afw.image.ExposureF(
        f'/scratch/arunj/deepCoadd_calexp_9813_{patch}_i_hsc_rings_v1_HSC_runs_RC2_w_2022_40_DM-36356_20221006T193226Z.fits'
    )
elif hsc_type == 'w04':
    exposure = lsst.afw.image.ExposureF(
        f"/projects/MERIAN/repo/HSC/runs/RC2/w_2022_04/DM-33402/20220128T212035Z/deepCoadd_calexp/9813/{patch}/i/deepCoadd_calexp_9813_{patch}_i_hsc_rings_v1_HSC_runs_RC2_w_2022_04_DM-33402_20220128T212035Z.fits"
    )
elif hsc_type == 'S20A':
    exposure = lsst.afw.image.ExposureF(
        f"/projects/MERIAN/repo/S20A/deepCoadd_calexp/9813/{patch_old}/calexp-HSC-I-{tract}-{patch_old}.fits"
    )


measureConfig = lsst.meas.base.ForcedPhotCoaddConfig()
measureConfig.footprintDatasetName = 'ScarletModelData'  # 'DeblendedFlux'

# measureConfig.measurement.plugins.names.add("ext_gaap_GaapFlux")
measureConfig.measurement.plugins.names.add("base_SdssShape")
measureConfig.measurement.plugins.names.add("base_GaussianFlux")
measureConfig.measurement.plugins.names.add("base_PsfFlux")
measureConfig.measurement.plugins.names.add("base_SdssCentroid")
# measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doMeasure = True  # Set it to False for timing comparison
# measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doPsfPhotometry = True
# measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doOptimalPhotometry = True
# measureConfig.measurement.plugins["ext_gaap_GaapFlux"].sigmas = [0.5, 0.75, 1.0, 1.5, 2.0]

photTask = lsst.meas.base.ForcedPhotCoaddTask(
    refSchema=refCat.schema, config=measureConfig)
photTask.log.setLevel('TRACE')


measCat, exposureID = photTask.generateMeasCat(exposureDataId=butler.registry.expandDataId(expID),
                                               exposure=exposure,
                                               refCat=refCat,
                                               refCatInBand=refCatInBand,
                                               refWcs=refExposure.wcs,
                                               idPackerName='tract_patch',
                                               footprintData=refCatInBand)


print("# Starting the measureTask at ", time.ctime())
t1 = time.time()
photTask.run(measCat, exposure, refCat=refCat,
             refWcs=refExposure.wcs, exposureId=exposureID)
t2 = time.time()
print("# Finished measureTask in %.2f seconds." % (t2 - t1))
cat2 = measCat.copy(deep=True).asAstropy()
