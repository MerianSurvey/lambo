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
import lsst.daf.butler as dafButler
import astropy.units as u

class GaapTask(object):
    def __init__(self, tract, patch, band, repo='/projects/MERIAN/repo/', collections='S20A/deepCoadd_calexp', is_merian=False):
        self.tract = tract
        self.patch = patch
        self.band = band
        self.repo = repo
        self.collections = collections
        self.patch_old = f'{self.patch % 9},{self.patch // 9}'
        if not is_merian:
            self._get_exposure()

    def _checkHSCfile(self):
        self.filename = os.path.join(self.repo, self.collections, str(self.tract), str(self.patch_old),
                                     f'calexp-HSC-{self.band.upper()}-{self.tract}-{self.patch_old}.fits')
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f'File {self.filename} not found')

    def load_merian_reference(self, band='N708', repo='/projects/MERIAN/repo/',
                              collections='DECam/runs/merian/dr1_wide', range=None):
        self.merian = GaapTask(self.tract, self.patch,
                               band, repo, collections, is_merian=True)
        self.merian.butler = dafButler.Butler(repo)
        dataId = dict(tract=self.tract, patch=self.patch, band=band)
        self.refCat = self.merian.butler.get(
            'deepCoadd_forced_src',
            collections=self.merian.collections,
            dataId=dataId,
            instrument='DECam',
            skymap='hsc_rings_v1',
        )
        if range is not None:
            self.refCat = self.refCat[range[0]:range[1]]
        self.refExposure = self.merian.butler.get(
            'deepCoadd_calexp',
            collections=self.merian.collections,
            dataId=dataId,
            instrument='DECam',
            skymap='hsc_rings_v1',
        )
        print('Loaded Merian reference catalog and image')

    def _get_exposure(self):
        self._checkHSCfile()
        self.exposure = lsst.afw.image.ExposureF(self.filename)
        print('Loaded HSC deepCoadd_calexp image')

    def setDefaultMeasureConfig(self):
        measureConfig = lsst.meas.base.ForcedMeasurementConfig()
        measureConfig.plugins.names.add("ext_gaap_GaapFlux")
        measureConfig.plugins.names.add("base_SdssShape")
        measureConfig.plugins.names.add("base_SdssCentroid")
        measureConfig.plugins.names.add("ext_gaap_GaapFlux")
        measureConfig.plugins["ext_gaap_GaapFlux"].doMeasure = True
        measureConfig.plugins["ext_gaap_GaapFlux"].doPsfPhotometry = True
        measureConfig.plugins["ext_gaap_GaapFlux"].doOptimalPhotometry = True
        measureConfig.plugins["ext_gaap_GaapFlux"].sigmas = [
            0.3, 0.4, 0.5, 0.6, 0.75, 1.0, 1.5, 2.0]
        self.measureConfig = measureConfig

    def run(self):
        measureTask = lsst.meas.base.ForcedMeasurementTask(
            refSchema=self.refCat.schema, config=self.measureConfig)
        measCat = measureTask.generateMeasCat(
            self.refExposure, self.refCat, self.refExposure.wcs, self.refCat.getIdFactory())

        print("# Starting the GAaP measureTask at ", time.ctime())
        t1 = time.time()
        measureTask.run(measCat, self.exposure,
                        refCat=self.refCat, refWcs=self.refExposure.wcs)
        t2 = time.time()
        print("# Finished the GAaP measureTask in %.2f seconds." % (t2 - t1))
        self.measCat = measCat

    def writeObjectTable(self):
        outCat = self.measCat.copy(deep=True).asAstropy()
        outCat['coord_ra'] = outCat['coord_ra'].to(u.deg)
        outCat['coord_dec'] = outCat['coord_dec'].to(u.deg)
        # self.measCat.writeFits(filename)
