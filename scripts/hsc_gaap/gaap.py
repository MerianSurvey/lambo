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
from lsst.geom import Box2I, Extent2I, Point2I, Point2D

# Write catalog
from lsst.pipe.tasks.postprocess import WriteObjectTableTask, TransformObjectCatalogTask, TransformObjectCatalogConfig
from lsst.pipe.tasks.parquetTable import MultilevelParquetTable
from lsst.afw.table import SourceCatalog

import lsst.daf.butler as dafButler
import astropy.units as u
from astropy.table import QTable, Table, vstack
from multiprocessing import Process, Manager


class GaapTask(object):
    def __init__(self, tract, patch, bands, hsc_type='S20A',
                 repo='/projects/MERIAN/repo/', collections='S20A/deepCoadd_calexp',
                 is_merian=False, logger=None, log_level='INFO', merian_butler=None):
        """
        Run GAaP on one patch of one tract of HSC data.

        Parameters
        ----------
        tract : int.
            Tract number. E.g., 9813 (COSMOS field).
        patch : int.
            Patch number. E.g., 23. There are 81 patches in a tract.
        bands : list of str.
            List of bands to run GAaP on. E.g., ``['g', 'r', 'i', 'z', 'y']`` or ``'grizy'``.
        hsc_type : str, optional.
            The version of HSC data. E.g., 'S20A' or 'w_2022_40'. The default is 'S20A'.
        repo : str, optional.
            The path to the data repo. The default is '/projects/MERIAN/repo/'.
            If you use ``w_2022_40`` data, you need to specify the repo to be ``'/projects/HSC/repo/main'``.
        collections : str, optional.
            The collection of the data. The default is ``'S20A/deepCoadd_calexp'``.
            If you use ``w_2022_40`` data, you need to specify the collection to be ``'HSC/runs/RC2/w_2022_40/DM-36151'``.
        is_merian : bool, optional.
            Whether the data is from MERIAN. The default is False.
        logger : logging.Logger, optional.
            External logger other than the LSSTpipe logger. The default is None.
        log_level : str, optional.
            The log level. The default is ``'ERROR'``.
        """
        self.tract = tract
        self.patch = patch
        self.bands = bands
        self.repo = repo
        self.collections = collections
        self.log_level = log_level
        self.patch_old = f'{self.patch % 9},{self.patch // 9}'
        self.logger = logger
        self.merian_butler = merian_butler

        if self.logger is not None:
            logger.setLevel(self.log_level)

        if not is_merian:
            self.hsc_type = hsc_type
            self.exposures = {}
            self._get_exposure()
            manager = Manager()
            self.forcedSrcCats = manager.dict()

    def _get_exposure(self):
        """
        Load HSC exposures.
        """
        if self.hsc_type == 'S20A':
            self._checkHSCfile()
            for band in self.bands:
                self.exposures[band] = lsst.afw.image.ExposureF(
                    self.filenames[band])
        elif self.hsc_type == 'w_2022_04':
            assert self.hsc_type in self.collections, 'w04 data not in repo'
            butler = dafButler.Butler(self.repo)
            for band in self.bands:
                dataId = dict(tract=self.tract, patch=self.patch, band=band)
                self.exposures[band] = butler.get(
                    'deepCoadd_calexp',
                    collections=self.collections,
                    dataId=dataId,
                    instrument='HSC',
                    skymap='hsc_rings_v1',
                )
        elif self.hsc_type == 'w_2022_40':
            assert self.hsc_type in self.collections, 'w40 data not in repo'
            butler = dafButler.Butler(self.repo)
            for band in self.bands:
                dataId = dict(tract=self.tract, patch=self.patch, band=band)
                self.exposures[band] = butler.get(
                    'deepCoadd_calexp',
                    collections=self.collections,
                    dataId=dataId,
                    instrument='HSC',
                    skymap='hsc_rings_v1',
                )
        if self.logger is not None:
            self.logger.debug(
                'Loaded HSC {} deepCoadd_calexp images in {} bands'.format(self.hsc_type, self.bands))

    def _checkHSCfile(self):
        """
        Check if HSC S20A data for the given tract and patch exist.
        """
        self.filenames = {}
        foundbands = ""
        for band in self.bands:
            self.filenames[band] = os.path.join(self.repo, self.collections, str(self.tract), str(self.patch_old),
                                                f'calexp-HSC-{band.upper()}-{self.tract}-{self.patch_old}.fits')
            if os.path.isfile(self.filenames[band]):
                foundbands += band
        notfoundbands = self.bands.replace(foundbands, "")

        if foundbands == "":       
            raise FileNotFoundError(
                f'File {self.filenames[notfoundbands[0]]} not found for all filters')
        elif foundbands != self.bands:
            print(f'File {self.filenames[notfoundbands[0]]} not found for filters {notfoundbands}')
            print(f'Running GAaP only on filters {foundbands}')
            self.bands = foundbands

    def load_merian_reference(self, band='N708', repo='/projects/MERIAN/repo/',
                              collections='DECam/runs/merian/dr1_wide', range=None):
        """
        Load Merian data as reference for forced photometry. 
        We use the scarlet footprint from Merian.

        Parameters
        ----------
        band: str, optional.
            The band of the reference catalog. The default is 'N708'.
        repo: str, optional.
            The path to the Merian repo. The default is '/projects/MERIAN/repo/'.
        collections: str, optional.
            The collection of the Merian data. The default is 'DECam/runs/merian/dr1_wide'.
            You may change it if you want to use other Merian data release.
        range: list of int, optional.
            This is used to specify the range of PARENTS in the reference catalog.
            In the debug mode, we only use a small range of PARENTS to save time.
            The children of these parents will be included automatically.
            In production mode, please set ``range=None``.
        """
        self.merian_refBand = band
        self.merian = GaapTask(self.tract, self.patch, None,
                               [band], repo, collections,
                               is_merian=True)

        if self.merian_butler is None:
            print("New butler")
            self.merian_butler = dafButler.Butler(repo)

        self.merian.dataId = dict(tract=self.tract, patch=self.patch,
                                  band=band, skymap='hsc_rings_v1')
        # self.ref_deepCoadd_obj = self.merian_butler.get(
        #     'deepCoadd_obj',
        #     collections=self.merian.collections,
        #     dataId=self.merian.dataId,
        #     instrument='DECam',
        # )
        self.refCat = self.merian_butler.get(
            'deepCoadd_ref',
            collections=self.merian.collections,
            dataId=self.merian.dataId,
            instrument='DECam',
        )
        self.refCatInBand = self.merian_butler.get(
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
        self.footprintCatInBand = self.merian_butler.get(
            'deepCoadd_scarletModelData',
            collections='DECam/runs/merian/dr1_wide',
            dataId=self.merian.dataId,
            instrument='DECam',
        )

        self.refExposure = self.merian_butler.get(
            'deepCoadd_calexp',
            collections=self.merian.collections,
            dataId=self.merian.dataId,
            instrument='DECam',
        )

        self.exposureDataId = self.merian_butler.registry.expandDataId(
            self.merian.dataId)

        print("Deleting merian_butler")
        del self.merian_butler

        # There can be very small differences in the WCS, so we need to make sure they are the same
        # If they are the same, then we overwrite the WCS in the exposure with the reference WCS
        for band in self.bands:
            if checkWcsEqual(self.exposures[band].getWcs(), self.refExposure.getWcs()) and self.hsc_type == 'S20A':
                self.exposures[band].setWcs(self.refExposure.getWcs())

        if self.logger is not None:
            self.logger.debug(
                f'Loaded Merian reference catalog and image in {self.merian_refBand} band')

    def setDefaultMeasureConfig(self):
        """
        Set default measurement config for ``lsst.meas.base.ForcedPhotCoaddTask``.
        For a complete list of the measurement plugins, please see
        ``/projects/MERIAN/repo/DECam/runs/merian/dr1_wide/20220921T193246Z/forcedPhotCoadd_config``.
        """
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
        self.measureConfig = measureConfig

    def setLogger(self, task):
        """
        Set logger for ``lsst.meas.base.ForcedPhotCoaddTask``.
        Typically used for debugging.
        """
        import logging
        task.log.setLevel(self.log_level)
        if not task.log.hasHandlers():
            logger = logging.getLogger('simple_example')
            logger.setLevel(logging.DEBUG)
            ch = logging.StreamHandler()
            ch.setLevel(logging.DEBUG)
            task.log.addHandler(ch)

    def runAll(self, pool=None):
        """
        Run ``gaap`` photometry on all bands.

        Parameters
        ----------
        pool : multiprocessing.Pool, optional.
            The multiprocessing pool to parallelize the measurements on different bands.
        """
        if pool is not None:
            if self.logger is not None:
                self.logger.debug(
                    'Running ForcedPhotCoaddTask on all bands using multiprocessing')
            successbands = pool.map(self.run, self.bands)
            pool.close()
            pool.join()
            self.bands = "".join(successbands)
            
        else:
            for band in self.bands:
                self.run(band)

    def run(self, band):
        """
        Run ``gaap`` photometry on a single band.
        """
        
        measureTask = lsst.meas.base.ForcedPhotCoaddTask(
            refSchema=self.refCat.schema, config=self.measureConfig)
        self.setLogger(measureTask)
        if self.logger is not None:
            logger = self.logger
        else:
            logger = measureTask.log

        logger.info(
            f'    - Patch {self.patch}: Generating ForcedPhotCoadd catalog for {band} band')
        
        # Get rid of parents that will crash the photometry because
        # Scarlet doesn't deal properly with weird patches in images
        catalog = measureTask.measurement.generateMeasCat(self.exposures[band], self.refCat, self.refExposure.wcs)
        parents = catalog[catalog["parent"] == 0]
        for parent in parents:
            parentId = parent.getId()
            try:
                self.footprintCatInBand.blends[parentId]   
                try:
                    position = Point2D(self.footprintCatInBand.blends[parentId].psfCenter)
                    self.exposures[band].getPsf().computeKernelImage(position).array[None, :, :]
                except:
                    del self.footprintCatInBand.blends[parentId]   

            except:
                continue


        measCat, exposureID = measureTask.generateMeasCat(exposureDataId=self.exposureDataId,
                                                        exposure=self.exposures[band],
                                                        refCat=self.refCat,
                                                        refCatInBand=self.refCatInBand,
                                                        refWcs=self.refExposure.wcs,
                                                        idPackerName='tract_patch',
                                                        footprintData=self.footprintCatInBand)
        logger.info(
            f'    - Patch {self.patch}: Running ForcedPhotCoaddTask on {band} band')
        print(
            f"# Starting the GAaP measureTask for patch={self.patch}, band={band} at", time.ctime())
        t1 = time.time()
        measureTask.run(measCat,
                        self.exposures[band],
                        refCat=self.refCat,
                        refWcs=self.refExposure.wcs,
                        exposureId=exposureID)
        t2 = time.time()
        logger.info(
        f'    - Patch {self.patch}: band {band}, Finished ForcedPhotCoaddTask in {(t2 - t1):.2f} seconds')
        print(
            f"# Finished the GAaP measureTask for band {band} in %.2f seconds." % (t2 - t1))
        self.forcedSrcCats[band] = measCat
        return(band)
        
    def writeObjectTable(self):
        """
        This step turns the ``measCat`` (which is effectively ``deepCoadd_forced_src``) into a ``deepCoadd_obj`` table, 
        where measurements of all filters are merged into one table. It returns a ``pandas.DataFrame``.

        References
        ----------
        https://github.com/lsst/pipe_tasks/blob/eee7ff78d7d7b7bcbdcc59fce9cbef2d184e5a8c/python/lsst/pipe/tasks/postprocess.py
        """
        writeTask = WriteObjectTableTask()
        catalogs = {}
        for band in self.bands:
            catalogs[band] = {'forced_src': self.forcedSrcCats[band],
                              'ref': SourceCatalog(),  # empty
                              'meas': SourceCatalog()  # empty
                              }
        self.deepCoadd_obj = writeTask.run(catalogs, self.tract, self.patch)
        print(self.logger)
        if self.logger is not None:
            self.logger.info(
                f'    - Patch {self.patch}: Finished WriteObjectTableTask')

    def transformObjectCatalog(self, functorFile=None):
        """
        Produce a flattened Object Table, in the same format of ``lsstpipe`` output ``objectTable``.
        The mapping between ``deepCoadd_obj`` and ``objectTable`` is specified in a functor file.
        The default functor file is ``'$PIPE_TASKS_DIR/schemas/Object.yaml'``. However, it doesn't
        contain customized ``gaap`` apertures. So I modified the original functor file and put it in
        ``scripts/hsc_gaap/Object.yaml``.

        Parameters
        ----------
        functorFile : str.
            If not specified, use the default functor file.

        References
        ----------
        https://github.com/lsst/pipe_tasks/blob/eee7ff78d7d7b7bcbdcc59fce9cbef2d184e5a8c/python/lsst/pipe/tasks/postprocess.py
        """
        import warnings
        warnings.filterwarnings("ignore")
        from contextlib import suppress

        parq = MultilevelParquetTable(dataFrame=self.deepCoadd_obj)
        transConfig = TransformObjectCatalogConfig()
        if functorFile is not None:
            transConfig.functorFile = functorFile
        with suppress(NotImplementedError):
            transTask = TransformObjectCatalogTask(config=transConfig)
            transTask.funcs.log.setLevel('FATAL')
            self.objectTable = transTask.run(parq)

        if self.logger is not None:
            self.logger.info(
                f'    - Patch {self.patch}: Finished TransformObjectCatalogTask')

    def saveObjectTable(self):
        """
        Save the ``objectTable`` to a parquet file.

        Parameters
        ----------
        fileName : str
            If not specified, use the default file name.
        """
        self.outCat = Table.from_pandas(self.objectTable, index='objectId')
        self.outCatDir = os.path.join(self.repo, 'S20A', 'gaapTable',
                                      str(self.tract), str(
                                          self.patch_old))
        self.outCatFileName = os.path.join(self.outCatDir,
                                           f'objectTable_{self.tract}_{self.patch_old}_{self.hsc_type}.fits')
        self.outCat = QTable(self.outCat)
        if not os.path.isdir(self.outCatDir):
            os.makedirs(self.outCatDir)
        self.outCat.write(self.outCatFileName, overwrite=True)
        if self.logger is not None:
            self.logger.info(
                f'    - Patch {self.patch}: Wrote GAaP table to {self.outCatFileName}')
        print('Wrote GAaP table to', self.outCatFileName)


def consolidateObjectTables(patches, tract=9813, hsc_type='S20A', repo='/projects/MERIAN/repo/'):
    """
    Concatenate catalogs in multiple patches.

    Parameters
    ----------
    patches : list of int.
        List of patches, in modern format. If given an empty list, all Merian-reduced patches
        in the tract are used. 
    tract : int.
        Tract number.
    hsc_type : str.
        The version of HSC data. E.g., ``S20A`` or ``w_2020_40``.
    repo : str, optional.
        The path to the data repo. The default is '/projects/MERIAN/repo/'.    """

    if len(patches) == 0:

        output_collection = "DECam/runs/merian/dr1_wide"
        skymap = "hsc_rings_v1"
        patches = findReducedPatches(tract, band="N708", output_collection=output_collection, skymap=skymap)

    cats = []
    for patch in patches:
        patch_old = f'{patch % 9},{patch // 9}'
        inCatDir = os.path.join(repo, 'S20A', 'gaapTable',
                                str(tract), str(patch_old))
        inCatFileName = os.path.join(inCatDir,
                                     f'objectTable_{tract}_{patch_old}_{hsc_type}.fits')
        cats.append(Table.read(inCatFileName))
    return vstack(cats)


def consolidateMerianCats(patches, tract=9813, collections='DECam/runs/merian/dr1_wide'):
    """
    Concatenate Merian catalogs in multiple patches.

    Parameters
    ----------
    patches : list of int.
        List of patches, in modern format.
    tract : int.
        Tract number.
    """
    import lsst.daf.butler as dafButler
    butler = dafButler.Butler('/projects/MERIAN/repo/')

    cats = []
    for patch in patches:
        patch_old = f'{patch % 9},{patch // 9}'
        dataId = dict(tract=tract, patch=patch)
        refCat = butler.get(
            'objectTable',
            collections=collections,
            dataId=dataId,
            instrument='DECam',
            skymap='hsc_rings_v1',
        )
        refCat = Table.from_pandas(refCat, index=True)
        cats.append(refCat)
    return vstack(cats)


def checkWcsEqual(wcs1, wcs2):
    """
    Check whether two WCS are equal (within some tolerance).

    Parameters
    ----------
    wcs1, wcs2 : `lsst.afw.geom.SkyWcs`.
    """
    # Check origin pixel
    d = wcs1.getPixelOrigin() - wcs2.getPixelOrigin()
    flag = (np.abs(d.x) < 1e-6) & (np.abs(d.y) < 1e-6)
    # Check transfom matrix
    flag &= np.all(np.abs(wcs1.getCdMatrix() - wcs2.getCdMatrix()) < 1e-8)
    return flag


def NaiveLogger(filename):
    import logging
    logging.basicConfig(filename=filename,
                        filemode='a',
                        format='%(asctime)s, %(name)s: %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.INFO)
    logger = logging.getLogger('GaapMaster')

    return logger

def findReducedPatches(tract, band="N708", output_collection='DECam/runs/merian/dr1_wide', skymap='hsc_rings_v1'):
    """Find all patches that have the necessary Merian data products for gaap photometry"""

    butler = dafButler.Butler('/projects/MERIAN/repo/')

    dataId = dict(tract=tract, band=band, skymap=skymap)
    
    deepCoadd_ref_patches = set([item.dataId["patch"] for item in 
                                butler.registry.queryDatasets('deepCoadd_ref', 
                                dataId=dataId, collections=output_collection,
                                skymap=skymap)])

    deepCoadd_meas_patches = set([item.dataId["patch"] for item in 
                                butler.registry.queryDatasets('deepCoadd_meas', 
                                dataId=dataId, collections=output_collection,
                                skymap=skymap)])

    deepCoadd_scarletModelData_patches = set([item.dataId["patch"] for item in 
                                butler.registry.queryDatasets('deepCoadd_scarletModelData', 
                                dataId=dataId, collections=output_collection,
                                skymap=skymap)])

    deepCoadd_calexp_patches = set([item.dataId["patch"] for item in 
                                butler.registry.queryDatasets('deepCoadd_calexp', 
                                dataId=dataId, collections=output_collection,
                                skymap=skymap)])

    del butler

    patches = deepCoadd_ref_patches.intersection(deepCoadd_meas_patches, 
                                                deepCoadd_scarletModelData_patches, 
                                                deepCoadd_calexp_patches)
    return(np.array(list(patches)))