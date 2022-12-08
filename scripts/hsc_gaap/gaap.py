import numpy as np
import time
import os
import lsst.meas.base
import lsst.pex.config
import lsst.afw.display as afwDisplay
import lsst.afw.geom as afwGeom
import lsst.afw.table
import lsst.meas.algorithms
import lsst.meas.deblender
import lsst.pex.exceptions
import lsst.meas.extensions.gaap

import lsst.meas.modelfit.psf.psfContinued
import lsst.meas.modelfit.optimizer.optimizer
import lsst.meas.modelfit
import lsst.meas.modelfit.cmodel.cmodelContinued
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegion
import lsst.meas.modelfit.cmodel.cmodel

import lsst.meas.base.forcedPhotCoadd
import lsst.meas.base.catalogCalculation
import lsst.meas.base.scaledApertureFlux
import lsst.meas.base.sdssShape
import lsst.ip.diffim.psfMatch
import lsst.meas.base.psfFlux
import lsst.meas.algorithms.subtractBackground
import lsst.meas.modelfit.psf.psfContinued
import lsst.meas.base.applyApCorr
import lsst.meas.base.noiseReplacer
import lsst.meas.base.blendedness
import lsst.meas.extensions.convolved.convolved
import lsst.meas.base.baseMeasurement
import lsst.meas.extensions.shapeHSM.hsmShapeControl
import lsst.meas.modelfit.optimizer.optimizer
import lsst.meas.modelfit
import lsst.afw.math._warper
import lsst.meas.base.forcedMeasurement
import lsst.meas.modelfit.cmodel.cmodelContinued
import lsst.ip.diffim._dipoleAlgorithms
import lsst.meas.base.footprintArea
import lsst.meas.modelfit.priors.priors
import lsst.meas.base.classification
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegion
import lsst.meas.base.localBackground
import lsst.meas.base.wrappers
import lsst.meas.extensions.photometryKron.photometryKron
import lsst.meas.base.pixelFlags
import lsst.pipe.base.config
import lsst.meas.extensions.gaap._gaap
import lsst.meas.extensions.gaap._gaussianizePsf
import lsst.meas.extensions.shapeHSM.hsmMomentsControl
import lsst.meas.base.plugins
import lsst.meas.base.apertureFlux
import lsst.meas.base.naiveCentroid
import lsst.meas.base.sdssCentroid
import lsst.meas.base.peakLikelihoodFlux
import lsst.meas.modelfit.cmodel.cmodel
import lsst.meas.base.gaussianFlux

import lsst.daf.butler as dafButler
import astropy.units as u
from astropy.table import QTable, Table, hstack, vstack


class GaapTask(object):
    def __init__(self, tract, patch, band, hsc_type='S20A',
                 repo='/projects/MERIAN/repo/', collections='S20A/deepCoadd_calexp',
                 is_merian=False, log_level='ERROR'):
        self.tract = tract
        self.patch = patch
        self.band = band
        self.repo = repo
        self.collections = collections
        self.log_level = log_level
        self.patch_old = f'{self.patch % 9},{self.patch // 9}'
        if not is_merian:
            self.hsc_type = hsc_type
            self._get_exposure()

    def _checkHSCfile(self):
        self.filename = os.path.join(self.repo, self.collections, str(self.tract), str(self.patch_old),
                                     f'calexp-HSC-{self.band.upper()}-{self.tract}-{self.patch_old}.fits')
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f'File {self.filename} not found')

    def load_merian_reference(self, band='N540', repo='/projects/MERIAN/repo/',
                              collections='DECam/runs/merian/dr1_wide', range=None):
        """
        Here the range is to indicate the range of PARENTS! 
        The children will be included automatically.

        """
        self.merian_refBand = band
        self.merian = GaapTask(self.tract, self.patch, None,
                               band, repo, collections,
                               is_merian=True)
        self.merian.butler = dafButler.Butler(repo)
        self.merian.dataId = dict(tract=self.tract, patch=self.patch,
                                  band=band, skymap='hsc_rings_v1')
        self.refCat = self.merian.butler.get(
            'deepCoadd_ref',
            collections=self.merian.collections,
            dataId=self.merian.dataId,
            instrument='DECam',
        )
        self.refCatInBand = self.merian.butler.get(
            'deepCoadd_meas',
            collections='DECam/runs/merian/dr1_wide',
            dataId=self.merian.dataId,
            instrument='DECam',
        )
        if range is not None:
            parents = np.arange(range[0], range[1], 1)
            temp = self.refCat[parents[0]:parents[-1] + 1]
            for parent in parents:
                childrens = np.where(
                    self.refCat['parent'] == self.refCat['id'][parent])[0]
                temp.extend(self.refCat[childrens[0]:childrens[-1] + 1])
            self.refCat = temp.copy()

            temp = self.refCatInBand[parents[0]:parents[-1] + 1]
            for parent in parents:
                childrens = np.where(
                    self.refCatInBand['parent'] == self.refCatInBand['id'][parent])[0]
                temp.extend(self.refCatInBand[childrens[0]:childrens[-1] + 1])
            self.refCatInBand = temp.copy()

        # We use the scarlet model data to get the footprint
        self.footprintCatInBand = self.merian.butler.get(
            'deepCoadd_scarletModelData',
            collections='DECam/runs/merian/dr1_wide',
            dataId=self.merian.dataId,
            instrument='DECam',
        )

        self.refExposure = self.merian.butler.get(
            'deepCoadd_calexp',
            collections=self.merian.collections,
            dataId=self.merian.dataId,
            instrument='DECam',
        )

        # There can be very small differences in the WCS, so we need to make sure they are the same
        # If they are the same, then we overwrite the WCS in the exposure with the reference WCS
        if checkWcsEqual(self.exposure.getWcs(), self.refExposure.getWcs()):
            self.exposure.setWcs(self.refExposure.getWcs())

        print('Loaded Merian reference catalog and image')

    def _get_exposure(self):
        if self.hsc_type == 'S20A':
            self._checkHSCfile()
            self.exposure = lsst.afw.image.ExposureF(self.filename)
        elif self.hsc_type == 'w_2022_04':
            assert self.hsc_type in self.collections, 'w04 data not in repo'
            butler = dafButler.Butler(self.repo)
            dataId = dict(tract=self.tract, patch=self.patch, band=self.band)
            self.exposure = butler.get(
                'deepCoadd_calexp',
                collections=self.collections,
                dataId=dataId,
                instrument='HSC',
                skymap='hsc_rings_v1',
            )
        elif self.hsc_type == 'w_2022_40':
            assert self.hsc_type in self.collections, 'w40 data not in repo'
            butler = dafButler.Butler(self.repo)
            dataId = dict(tract=self.tract, patch=self.patch, band=self.band)
            self.exposure = butler.get(
                'deepCoadd_calexp',
                collections=self.collections,
                dataId=dataId,
                instrument='HSC',
                skymap='hsc_rings_v1',
            )
        print(f'Loaded HSC {self.hsc_type} deepCoadd_calexp image')

    def setDefaultMeasureConfig(self):
        measureConfig = lsst.meas.base.ForcedPhotCoaddConfig()
        measureConfig.footprintDatasetName = 'ScarletModelData'
        # measureConfig.footprintDatasetName = 'DeblendedFlux'
        # measureConfig.measurement.doReplaceWithNoise = False

        measureConfig.measurement.plugins.names.add("base_PsfFlux")
        measureConfig.measurement.plugins.names.add(
            "base_CircularApertureFlux")
        measureConfig.measurement.plugins.names.add("base_SdssShape")
        measureConfig.measurement.plugins.names.add("base_SdssCentroid")
        measureConfig.measurement.plugins.names.add("base_Blendedness")
        measureConfig.measurement.plugins.names.add(
            'modelfit_DoubleShapeletPsfApprox')
        measureConfig.measurement.plugins.names.add('modelfit_CModel')

        measureConfig.measurement.plugins.names.add("ext_gaap_GaapFlux")
        measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doMeasure = True
        measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doPsfPhotometry = True
        measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doOptimalPhotometry = True
        measureConfig.measurement.plugins["ext_gaap_GaapFlux"].sigmas = [
            0.5, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.5, 2.0, 2.5, 3.5]
        # measureConfig.measurement.plugins["ext_gaap_GaapFlux"].scalingFactors = [
        #     1.5]
        # [1.15, 1.25, 1.5]
        self.measureConfig = measureConfig

    def setLogger(self, task):
        import logging
        task.log.setLevel(self.log_level)
        if not task.log.hasHandlers():
            logger = logging.getLogger('simple_example')
            logger.setLevel(logging.DEBUG)
            ch = logging.StreamHandler()
            ch.setLevel(logging.DEBUG)
            task.log.addHandler(ch)

    def run(self):
        measureTask = lsst.meas.base.ForcedPhotCoaddTask(
            refSchema=self.refCat.schema, config=self.measureConfig)
        self.setLogger(measureTask)
        measCat, exposureID = measureTask.generateMeasCat(exposureDataId=self.merian.butler.registry.expandDataId(self.merian.dataId),
                                                          exposure=self.exposure,
                                                          refCat=self.refCat,
                                                          refCatInBand=self.refCatInBand,
                                                          refWcs=self.refExposure.wcs,
                                                          idPackerName='tract_patch',
                                                          footprintData=self.footprintCatInBand)
        self.measCat = measCat
        # return
        print("# Starting the GAaP measureTask at ", time.ctime())
        t1 = time.time()
        measureTask.run(measCat,
                        self.exposure,
                        refCat=self.refCat,
                        refWcs=self.refExposure.wcs,
                        exposureId=exposureID)
        t2 = time.time()
        print("# Finished the GAaP measureTask in %.2f seconds." % (t2 - t1))
        self.measCat = measCat

    def writeObjectTable(self, save=True):
        outCat = self.measCat.copy(deep=True).asAstropy()
        outCat['coord_ra'] = outCat['coord_ra'].to(u.deg)
        outCat['coord_dec'] = outCat['coord_dec'].to(u.deg)
        rawCat = outCat.copy()  # First write a table with raw columns

        old_gaap_cols = [
            item for item in self.measCat.schema.getNames() if 'gaap' in item and 'apCorr' not in item]
        old_gaap_cols += ['base_PsfFlux_instFlux',
                          'base_PsfFlux_instFluxErr', 'base_PsfFlux_flag']

        outCat = outCat[['id', 'coord_ra', 'coord_dec'] + old_gaap_cols]

        # PhotCalib
        for col in old_gaap_cols:
            if 'instFlux' in col:
                outCat[col] = outCat[col].value * \
                    self.exposure.getPhotoCalib().instFluxToNanojansky(1) * u.nanomaggy

        new_gaap_cols = []
        for col in old_gaap_cols:
            name = col.replace('base', f'{self.band}')
            name = name.replace('ext_gaap_GaapFlux', f'{self.band}_gaap')
            name = name.replace('_instFlux', 'Flux').replace('PsfFlux', 'Psf')
            if 'Flux' in name:
                aper = name.split(
                    "x_")[-1].replace('FluxErr', '').replace('Flux', '')
                name = name.replace('_1_15x', '')
                name = name.replace('_' + aper, aper.replace('_', 'p'))

            if 'flag' in name:
                aper = name.split(
                    "x_")[-1].replace('_flag_bigPsf', '').replace('_flag', '')
                name = name.replace('_1_15x', '')
                if not 'gauss' in name:
                    name = name.replace(
                        '_' + aper, aper.replace('_', 'p') + 'Flux')
            new_gaap_cols.append(name)

        outCat.rename_columns(old_gaap_cols, new_gaap_cols)
        self.outCatDir = os.path.join('/projects/MERIAN/repo/', 'S20A', 'gaapTable',
                                      str(self.tract), str(
                                          self.patch_old))
        self.outCatFileName = os.path.join(self.outCatDir,
                                           f'gaapTable_{self.band.upper()}_{self.tract}_{self.patch_old}_{self.hsc_type}.fits')
        self.outCat = QTable(outCat)

        if save:
            if not os.path.isdir(self.outCatDir):
                os.makedirs(self.outCatDir)
            self.outCat.write(self.outCatFileName, overwrite=True)
            print('Wrote GAaP table to', self.outCatFileName)
            rawCat.write(os.path.join(self.outCatDir,
                                      f'_raw_gaapTable_{self.band.upper()}_{self.tract}_{self.patch_old}_{self.hsc_type}.fits'))
        return self.outCat


def joinCatBands(patch=23, filters='griz', tract=9813, hsc_type='w_2022_40'):
    """
    Merge catalogs in multiple bands
    """
    cats = []
    patch_old = f'{patch % 9},{patch // 9}'
    for i, filt in enumerate(list(filters)):
        temp = Table.read(
            f'/projects/MERIAN/repo/S20A/gaapTable/{tract}/{patch_old}/gaapTable_{filt.upper()}_{tract}_{patch_old}_{hsc_type}.fits')
        if i > 0:
            temp.remove_columns(['id', 'coord_ra', 'coord_dec'])
        cats.append(temp)
    return hstack(cats)


def joinCatPatches(patches, filters='griz', tract=9813, hsc_type='w_2022_40'):
    """
    Concatenate catalogs in multiple patches
    """
    cats = []
    for i in patches:
        cats.append(joinCatBands(patch=i, filters=filters,
                    tract=tract, hsc_type=hsc_type))
    return vstack(cats)


def joinMerianCatPatches(patches, tract=9813):
    import lsst.daf.butler as dafButler
    butler = dafButler.Butler('/projects/MERIAN/repo/')

    cats = []
    for patch in patches:
        patch_old = f'{patch % 9},{patch // 9}'
        dataId = dict(tract=tract, patch=patch)
        refCat = butler.get(
            'objectTable',
            collections='DECam/runs/merian/dr1_wide',
            dataId=dataId,
            instrument='DECam',
            skymap='hsc_rings_v1',
        )
        refCat = Table.from_pandas(refCat, index=True)
        cats.append(refCat)
    return vstack(cats)


def checkWcsEqual(wcs1, wcs2):
    # Check origin pixel
    d = wcs1.getPixelOrigin() - wcs2.getPixelOrigin()
    flag = (np.abs(d.x) < 1e-5) & (np.abs(d.y) < 1e-5)
    flag &= np.all(np.abs(wcs1.getCdMatrix() - wcs2.getCdMatrix()) < 1e-8)
    return flag


def transformObjectCatalog():
    pass
