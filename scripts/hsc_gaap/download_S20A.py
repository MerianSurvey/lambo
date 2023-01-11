import os
import fire

from unagi import hsc
s20a = hsc.Hsc(dr='dr3', rerun='s20a_wide')


def download_hsc_coadd(tract=9813, patch='8,6', filt='i', outdir='/projects/MERIAN/repo/S20A/deepCoadd_calexp'):
    outdir = os.path.join(outdir, str(tract), str(patch))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    filename = f'calexp-HSC-{filt.upper()}-{tract}-{patch}.fits'
    filename = os.path.join(outdir, filename)
    if not os.path.isfile(filename):
        s20a.download_patch(
            tract, patch, filt=f'HSC-{filt.upper()}', output_file=filename, verbose=False)
        print('# Downloading ' + filename + ' finished! ')


def download_blendedness(tract=9813, patch='8,6', save=False):
    patch_x = int(patch[0])
    patch_y = int(patch[2])
    skymap_id = tract * 10000 + patch_x * 100 + patch_y

    sql_test = open(
        '/home/jiaxuanl/Research/Merian/merian_tractor/scripts/hsc_gaap/blendedness.SQL', 'r').read()
    sql_test = sql_test.replace('10000000', str(skymap_id))
    print(
        f'# SQL query blendedness from HSC database for tract = {tract}, patch = {patch}')

    result_test = s20a.sql_query(
        sql_test, from_file=False, preview=False, verbose=True)

    if save:
        if not os.path.isdir(f'/projects/MERIAN/repo/S20A/gaapTable/{tract}/{patch}'):
            os.makedirs(
                f'/projects/MERIAN/repo/S20A/gaapTable/{tract}/{patch}')
        result_test.write(
            f'/projects/MERIAN/repo/S20A/gaapTable/{tract}/{patch}/S20A_blendedness_{patch}.fits',
            overwrite=True)
    return result_test


def runDownload(tract=9813):
    from tqdm import tqdm
    import itertools
    for i, j in tqdm(itertools.product(range(0, 8), range(0, 8))):
        for filt in list('grizy'):
            try:
                download_hsc_coadd(tract=tract, patch=f'{i},{j}', filt=filt,
                                   outdir='/projects/MERIAN/repo/S20A/deepCoadd_calexp/')
            except Exception as e:
                print(e)
                print('Failed for tract = ', tract, ' patch = ',
                      f'{i},{j}', ' filter = ', filt)
        download_blendedness(tract=tract, patch=f'{i},{j}', save=True)

    print('Finished for tract = ', tract)


if __name__ == '__main__':
    fire.Fire(runDownload)

# python download_S20A.py --tract=9812