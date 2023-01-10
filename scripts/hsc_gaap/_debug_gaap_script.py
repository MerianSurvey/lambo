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


fig, axes = plt.subplots(1, 1, figsize=(14, 6), sharey=True)

plt.sca(axes)  # [0])
aper = '1p0'
flag = cat['spec_z']
sct = plt.scatter(cat['z_cosmos'],
                  - 2.5 *
                  np.log10(cat[f'i_gaap{aper}Flux'] /
                           cat[f'N708_gaap{aper}Flux']),
                  c=-2.5 *
                  np.log10(cat[f'r_gaap{aper}Flux'] /
                           cat[f'N540_gaap{aper}Flux']),
                  cmap='Spectral', vmin=-1.0, vmax=0.3, edgecolors='none', s=30)

sct = plt.scatter(cat['z_cosmos'],
                  - 2.5 *
                  np.log10(cat[f'i_gaap{aper}Flux'] /
                           cat[f'N708_gaap{aper}Flux']),
                  c=-2.5 *
                  np.log10(cat[f'r_gaap{aper}Flux'] /
                           cat[f'N540_gaap{aper}Flux']),
                  cmap='Spectral', vmin=-1.0, vmax=0.3, edgecolors='none', s=30)

plt.axvline(0.058, ls='--', alpha=0.6)
plt.axvline(0.10, ls='--', alpha=0.6)
plt.ylim(-1, 1.0)
plt.xlim(0, 1.5)
plt.title(r'Manual GAAP $\texttt{1p0}$')
plt.xlabel('$z_\mathrm{COSMOS}$', fontsize=18)
plt.ylabel('$i$ - N708', fontsize=18)

# plt.sca(axes[1])
# sct = plt.scatter(obj_cat['z_cosmos'],
#             meas_cat['mag'][:, 2] - meas_cat['mag'][:, 5],
#             c=meas_cat['mag'][:, 1] - meas_cat['mag'][:, 6],
#             cmap='Spectral', vmin=-0.8, vmax=0.4)
# plt.axvline(0.058, ls='--', alpha=0.6)
# plt.axvline(0.10, ls='--', alpha=0.6)
# plt.ylim(-0.4, 1.0)
# plt.title('Tractor')
# plt.xlabel('$z_\mathrm{COSMOS}$', fontsize=18)

fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.88, 0.13, 0.02, 0.74])
cbar = fig.colorbar(sct, cax=cbar_ax)
cbar.set_label('$r$ - N540', fontsize=19)
# for t in cbar.ax.get_label():
#     t.set_fontsize(18)

# # fig.colorbar(im, cax=cax, orientation='vertical')
# fig.colorbar(sct, cax=cax, label='$r$ - N540', )
plt.subplots_adjust(wspace=0.05)
