import os, sys
import fire
import numpy as np
import lsst.daf.butler as dafButler
from unagi import hsc
s20a = hsc.Hsc(dr='dr3', rerun='s20a_wide')
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from hsc_gaap.gaap import findReducedPatches

def download_hsc_coadd(tract=9813, patch='8,6', filt='i', outdir='/projects/MERIAN/repo/'):
    outdir = os.path.join(outdir, 'S20A/deepCoadd_calexp', str(tract), str(patch))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    filename = f'calexp-HSC-{filt.upper()}-{tract}-{patch}.fits'
    filename = os.path.join(outdir, filename)
    if not os.path.isfile(filename):
        s20a.download_patch(
            tract, patch, filt=f'HSC-{filt.upper()}', output_file=filename, verbose=False)
        print('# Downloading ' + filename + ' finished! ')


def download_blendedness(tract=9813, patch='8,6', save=False, outdir='/projects/MERIAN/repo/'):
    outdir = os.path.join(outdir, 'S20A/gaapTable')

    patch_x = int(patch[0])
    patch_y = int(patch[2])
    skymap_id = tract * 10000 + patch_x * 100 + patch_y

    sql_test = open(
        os.path.join(os.getenv("LAMBO_HOME"), 'lambo/scripts/hsc_gaap/blendedness.SQL'), 'r').read()
    sql_test = sql_test.replace('10000000', str(skymap_id))
    print(
        f'# SQL query blendedness from HSC database for tract = {tract}, patch = {patch}')

    result_test = s20a.sql_query(
        sql_test, from_file=False, preview=False, verbose=True)

    if save:
        if not os.path.isdir(f'{outdir}/{tract}/{patch}'):
            os.makedirs(
                f'{outdir}/{tract}/{patch}')
        result_test.write(
            f'{outdir}/{tract}/{patch}/S20A_blendedness_{patch}.fits',
            overwrite=True)
    return result_test


def runDownload(tract=9813, outdir='/projects/MERIAN/repo/', only_merian=True, alltracts=False):
    # download all patches for all tracts with merian reduced data
    if alltracts:
        output_collection = "DECam/runs/merian/dr1_wide"
        data_type = "deepCoadd_calexp"
        skymap = "hsc_rings_v1"
        butler = dafButler.Butler('/projects/MERIAN/repo/', collections=output_collection, skymap=skymap)

        patches = np.array([[data_id['tract'], data_id["patch"]] for data_id in butler.registry.queryDataIds (['tract','patch'], datasets=data_type, 
                                                        collections=output_collection, skymap=skymap)])
        tracts = np.unique(patches[:,0])
        del butler 

        for t in tracts:
            runDownload(tract=t, outdir=outdir, only_merian=only_merian, alltracts=False)
        return()

    if type(tract) is list:
        for t in tract:
            runDownload(tract=t, outdir=outdir, only_merian=only_merian, alltracts=False)
        return()

    if only_merian:
        patches = findReducedPatches(tract)
        patches = [[p % 9, int(p/9)] for p in patches]

    else:
        from tqdm import tqdm
        import itertools
        patches = tqdm(itertools.product(range(0, 9), range(0, 9)))

    for i, j in patches:
        for filt in list('grizy'):
            try:
                download_hsc_coadd(tract=tract, patch=f'{i},{j}', filt=filt,
                                   outdir=outdir)
            except Exception as e:
                print(e)
                print('Failed for tract = ', tract, ' patch = ',
                      f'{i},{j}', ' filter = ', filt)
        if not os.path.isfile(os.path.join(outdir, f'S20A/gaapTable/{tract}/{i},{j}/S20A_blendedness_{i},{j}.fits')):
            download_blendedness(tract=tract, patch=f'{i},{j}', save=True, outdir=outdir)

    print('Finished for tract = ', tract)



if __name__ == '__main__':
    fire.Fire(runDownload)

# python download_S20A.py --tract=9812
