import lsst.daf.butler as dafButler
import numpy as np
import glob
import os, sys


def findReducedPatches(tract, band="N708", output_collection='DECam/runs/merian/dr1_wide', skymap='hsc_rings_v1'):
        """Returns all merian patches for a given tract that have the necessary data products to run gaap"""

        butler = dafButler.Butler('/projects/MERIAN/repo/')
        dataId = dict(tract=tract, band=band, skymap=skymap)
        
        deepCoadd_ref_patches = set([item.dataId["patch"] for item in butler.registry.queryDatasets('deepCoadd_ref', 
                dataId=dataId, collections=output_collection,
                skymap=skymap)])

        deepCoadd_meas_patches = set([item.dataId["patch"] for item in butler.registry.queryDatasets('deepCoadd_meas', 
                dataId=dataId, collections=output_collection,
                skymap=skymap)])

        deepCoadd_scarletModelData_patches = set([item.dataId["patch"] for item in butler.registry.queryDatasets('deepCoadd_scarletModelData', 
                dataId=dataId, collections=output_collection,
                skymap=skymap)])

        deepCoadd_calexp_patches = set([item.dataId["patch"] for item in butler.registry.queryDatasets('deepCoadd_calexp', 
                dataId=dataId, collections=output_collection,
                skymap=skymap)])

        del butler

        patches = deepCoadd_ref_patches.intersection(deepCoadd_meas_patches, deepCoadd_scarletModelData_patches, deepCoadd_calexp_patches)
        return(np.array(list(patches)))


def findGaapReducedPatches(tract, repo = '/scratch/gpfs/am2907/Merian/gaap/'):
    """Returns all merian patches in a given tract that have already 
       been gaap reduced (have gaap catalogs in appropriate directory)"""

    dir = os.path.join(repo, f"S20A/gaapTable/{tract}/")
    cats = glob.glob(os.path.join(dir, "**/*objectTable*"), recursive=True)
    cat_tract = glob.glob(os.path.join(dir, "objectTable*"), recursive=False)

    if len(cats) > 0:
        if len(cat_tract) > 0 :
            cats.remove(cat_tract[0])

        patches = np.sum(np.array([i.split(f"objectTable_{tract}_")[1].split("_S20A")[0].split(",")
                                    for i in cats], dtype=int)*[1,9],1)
        patches.sort()

    else:
        patches = []
    return (patches)


def hasCompiledCatalog(tract, repo='/scratch/gpfs/am2907/Merian/gaap'):
    dir = os.path.join(repo, f"S20A/gaapTable/{tract}/")
    cat_tract = glob.glob(os.path.join(dir, "objectTable*"), recursive=False)

    if len(cat_tract) > 0 :
        return (cat_tract[0])
    else:
         return(False)


def saveMerianReducedPatchList(tracts, outfile, band="N708"):
     f = open(outfile, "w")
     for tract in tracts:
        patches = ",".join(findReducedPatches(tract, band=band).astype(str))
        if len(patches) > 0:
            f.write(f'{tract},{patches}\n')
     f.close()
     print(f"Saved file to {outfile}.")


def saveGaapReducedPatchList(tracts, outfile):
     f = open(outfile, "w")
     for tract in tracts:
        patches = ",".join(findGaapReducedPatches(tract).astype(str))
        if len(patches) > 0:
            f.write(f'{tract},{patches}\n')
     f.close()
     print(f"Saved file to {outfile}.")


def saveGaapNotReducedPatchList(tracts, outfile, patches = [], band="N708", notionformat=False):
     f = open(outfile, "w")
     if notionformat:
         f.write("Tract,Tract #,# patches,Band\n")

     for tract in tracts:
        if len(patches) ==0 :
            patches_mer  = findReducedPatches(tract, band=band)
        else:
            patches_mer = patches[np.where(tracts==tract)[0][0]]
        patches_gaap = findGaapReducedPatches(tract)
        patches_i = np.array(list(set(patches_mer) - set(patches_gaap)))
        npatches = len(patches_i)
        patches_i = ",".join(patches_i.astype(str))
        if npatches > 0:
            if not notionformat:
                f.write(f'{tract},{patches_i}\n')
            else:
                f.write(f'{tract},{tract},{npatches},{band}\n')
     f.close()
     print(f"Saved file to {outfile}.")
