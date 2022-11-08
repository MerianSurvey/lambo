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
from astropy.table import QTable


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
        self.merian.dataId = dict(tract=self.tract, patch=self.patch,
                                  band=band, skymap='hsc_rings_v1')
        self.refCat = self.merian.butler.get(
            'deepCoadd_forced_src',
            collections=self.merian.collections,
            dataId=self.merian.dataId,
            instrument='DECam',
        )
        self.refCatInBand = self.merian.butler.get(
            'deepCoadd_ref',
            collections='DECam/runs/merian/dr1_wide',
            dataId=self.merian.dataId,
            instrument='DECam',
        )
        if range is not None:
            self.refCat = self.refCat[range[0]:range[1]]
            self.refCatInBand = self.refCatInBand[range[0]:range[1]]

        self.refExposure = self.merian.butler.get(
            'deepCoadd_calexp',
            collections=self.merian.collections,
            dataId=self.merian.dataId,
            instrument='DECam',
        )
        print('Loaded Merian reference catalog and image')

    def _get_exposure(self):
        self._checkHSCfile()
        self.exposure = lsst.afw.image.ExposureF(self.filename)
        print('Loaded HSC deepCoadd_calexp image')

    def setDefaultMeasureConfig(self):
        measureConfig = lsst.meas.base.ForcedPhotCoaddConfig()
        measureConfig.footprintDatasetName = 'DeblendedFlux'
        measureConfig.measurement.plugins.names.add("ext_gaap_GaapFlux")
        measureConfig.measurement.plugins.names.add("base_SdssShape")
        measureConfig.measurement.plugins.names.add("base_SdssCentroid")
        measureConfig.measurement.plugins.names.add("ext_gaap_GaapFlux")
        measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doMeasure = True
        measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doPsfPhotometry = True
        measureConfig.measurement.plugins["ext_gaap_GaapFlux"].doOptimalPhotometry = True
        measureConfig.measurement.plugins["ext_gaap_GaapFlux"].sigmas = [
            0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.5, 2.0]
        # measureConfig.measurement.plugins["ext_gaap_GaapFlux"].scalingFactors = [
        #     1.5]
        # [1.15, 1.25, 1.5]
        self.measureConfig = measureConfig

    def run(self):
        measureTask = lsst.meas.base.ForcedPhotCoaddTask(
            refSchema=self.refCat.schema, config=self.measureConfig)
        measCat, exposureID = measureTask.generateMeasCat(exposureDataId=self.merian.butler.registry.expandDataId(self.merian.dataId),
                                                          exposure=self.exposure,
                                                          refCat=self.refCat,
                                                          refCatInBand=self.refCatInBand,
                                                          refWcs=self.refExposure.wcs,
                                                          idPackerName='tract_patch',
                                                          footprintData=self.refCatInBand)
        self.measCat = measCat
        # return
        print("# Starting the GAaP measureTask at ", time.ctime())
        t1 = time.time()
        measureTask.run(measCat, self.exposure, refCat=self.refCat,
                        refWcs=self.refExposure.wcs, exposureId=exposureID)
        t2 = time.time()
        print("# Finished the GAaP measureTask in %.2f seconds." % (t2 - t1))
        self.measCat = measCat

    def writeObjectTable(self):
        outCat = self.measCat.copy(deep=True).asAstropy()
        outCat['coord_ra'] = outCat['coord_ra'].to(u.deg)
        outCat['coord_dec'] = outCat['coord_dec'].to(u.deg)

        old_gaap_cols = [
            item for item in self.measCat.schema.getNames() if 'gaap' in item and 'apCorr' not in item]
        outCat = outCat[['id', 'coord_ra', 'coord_dec'] + old_gaap_cols]

        # PhotCalib
        for col in old_gaap_cols:
            if 'instFlux' in col:
                outCat[col] = outCat[col].value * \
                    self.exposure.getPhotoCalib().instFluxToNanojansky(1) * u.nanomaggy
                # outCat[col] = outCat[col].value * \
                #     self.exposure.photoCalib.instFluxToMagnitude(
                #         1) * u.nanomaggy

        new_gaap_cols = []
        for col in old_gaap_cols:
            name = col.replace('ext_gaap_GaapFlux', f'{self.band}_gaap')
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
        self.outCatDir = os.path.join(self.repo, self.collections.split('/')[0], 'gaapTable',
                                      str(self.tract), str(
            self.patch_old))
        self.outCatFileName = os.path.join(self.outCatDir,
                                           f'gaapTable_{self.band.upper()}_{self.tract}_{self.patch_old}.fits')
        self.outCat = QTable(outCat)

        if not os.path.isdir(self.outCatDir):
            os.makedirs(self.outCatDir)
        self.outCat.write(self.outCatFileName, overwrite=True)
        print('Wrote GAaP table to', self.outCatFileName)
        return self.outCat
