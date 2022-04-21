"""
This file runs tractor fitting for objects in a give catalog.
"""
import fire
import os
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import pickle
import multiprocess as mp
mp.freeze_support()
from functools import partial

import kuaizi
from kuaizi.detection import Data
from kuaizi.utils import padding_PSF
from kuaizi.tractor.utils import initialize_meas_cat, _write_to_row
from kuaizi.tractor.fit import tractor_hsc_sep_blob_by_blob

import sys
sys.path.append('/home/jiaxuanl/software/astrometry.net-0.89')
sys.path.append('/home/jiaxuanl/Research/Packages/tractor/')


def fitting_obj(index, obj_cat, meas_cat,
                channels=['g', 'r', 'i', 'z', 'N708', 'N540'], ref_filt='i'):
    """
    Wrapper for fitting one object.
    """
    row = meas_cat[index]
    obj = obj_cat[index]

    ID_name = 'ID' if 'ID' in obj_cat.colnames else 'id'
    PREFIX_name = 'PREFIX' if 'PREFIX' in obj_cat.colnames else 'prefix'

    row['ID'] = obj[ID_name]

    forced_channels = [filt for filt in channels if filt != ref_filt]

    print(f'### Tractor modeling for obj {index}, ID = {obj[ID_name]}')

    try:
        # Load data
        cutout = []
        for filt in channels:
            if 'N' in filt:
                cutout.append(fits.open(os.path.join(
                    obj['dir'], f"{obj[PREFIX_name]}_{filt}_deepCoadd_calexp.fits")))
            else:
                cutout.append(fits.open(os.path.join(
                    obj['dir'], f's18a_wide_{obj[ID_name]}_{filt}.fits')))

        psf_list = []
        for filt in channels:
            if 'N' in filt:
                psf_list.append(fits.open(os.path.join(
                    obj['dir'], f"{obj[PREFIX_name]}_{filt}_psf.fits")))
            else:
                psf_list.append(fits.open(os.path.join(
                    obj['dir'], f's18a_wide_{obj[ID_name]}_{filt}_psf.fits')))

        # Reconstruct data
        images = np.array([hdu[1].data for hdu in cutout if len(hdu) > 1])
        # note: all bands share the same WCS here
        w = wcs.WCS(cutout[0][1].header)
        weights = 1 / np.array([hdu[3].data for hdu in cutout])
        psf_pad = padding_PSF(psf_list)  # Padding PSF cutouts from HSC
        data = Data(images=images, weights=weights, wcs=w,
                    psfs=psf_pad, channels=channels)
    except Exception as e:
        print(f'IMCOMPLETE FILE FOR {index} = obj {obj[ID_name]}', e)

    # Start fitting
    try:
        model_dict = {}
        # fitting in the i-band first: then pass the i-band parameters of target galaxy to other bands
        model_dict[ref_filt], _obj_cat_i = tractor_hsc_sep_blob_by_blob(
            obj, ref_filt, data.channels, data,
            freeze_dict={'pos': False, 'shape': False, 'shape.re': False, 'shape.ab': False, 'shape.phi': False,
                         'sersicindex': False},  # don't fix shape/sersic
            verbose=False)

        for filt in forced_channels:
            pos = True
            fix_all = True
            ref_catalog = model_dict[ref_filt].catalog.copy()
            model_dict[filt], _ = tractor_hsc_sep_blob_by_blob(
                obj, filt, data.channels, data,
                fix_all=fix_all, tractor_cat=ref_catalog,
                obj_cat=_obj_cat_i,
                freeze_dict={'pos': pos, 'shape': True, 'shape.re': True, 'shape.ab': True, 'shape.phi': True,
                             'sersicindex': True},  # don't fix shape/sersic
                verbose=False)
        # with open(os.path.join(obj['dir'], f'cosmos_{obj[ID_name]}_tractor.pkl'), 'wb') as f:
        #     pickle.dump(model_dict, f)
    except Exception as e:
        print(f'FITTING ERROR FOR {index} = obj {obj[ID_name]}', e)
    try:
        row = _write_to_row(row, model_dict, channels=channels)
    except Exception as e:
        print(f'MEAUREMENT ERROR FOR {index} = obj {obj[ID_name]}', e)
    return row


# cat_dir = f'./Catalogs/COSMOS_cutouts_tractor_gaap_{cat_prefix}.fits'
def multiprocess_fitting(cat_dir, suffix='midz', njobs=16, low=0, high=250, ind_list=None,
                         DATADIR='/scratch/gpfs/jiaxuanl/Data/Merian',
                         CUTOUT_SUBDIR='./Cutout/_new/',
                         CATALOG_SUBDIR='./Catalogs/'):
    print('SET ENVIRONMENT')
    os.chdir(DATADIR)
    print("CURRENT WORKING DIRECTORY:", os.getcwd())
    print('Number of processor to use:', njobs)

    # Load catalog
    print('Load catalog')

    obj_cat = Table.read(cat_dir)
    obj_cat['dir'] = [file.replace('/projects/MERIAN/poststamps/g09_broadcut_new',
                                   os.path.join(DATADIR, CUTOUT_SUBDIR),
                                   ) for file in obj_cat['dir']]

    if low is None:
        low = 0
    if high is None:
        high = len(obj_cat)

    # Setup fitting
    channels = list('grizy') + ['N708', 'N540']
    ref_filt = 'i'

    # Initialize output catalog
    meas_cat = initialize_meas_cat(obj_cat, channels=channels)

    # Multiprocessing
    pwel = mp.Pool(njobs)
    iterable = range(low, high) if ind_list is None else ind_list
    meas_cat_tractor = pwel.map(
        partial(fitting_obj, obj_cat=obj_cat, meas_cat=meas_cat,
                channels=channels, ref_filt=ref_filt), iterable)

    if ind_list is not None:
        output_filename = f'{DATADIR}/{CATALOG_SUBDIR}/tractor_{suffix}_output_ind_list.pkl'
    else:
        output_filename = f'{DATADIR}/{CATALOG_SUBDIR}/tractor_{suffix}_output_{low}_{high}.pkl'

    with open(output_filename, 'wb') as f:
        pickle.dump(meas_cat_tractor, f)

    meas_cat_tractor = vstack(meas_cat_tractor)
    output_filename = output_filename.replace('.pkl', '.fits')
    print('#############')
    print('Writing to file:', output_filename)
    meas_cat_tractor.write(output_filename, overwrite=True)
    print('#############')


# cat_dir = './Catalogs/magellan/magellan_spec_obj_cat.fits'
def multiprocess_fitting_magellan(cat_dir, suffix='', njobs=16, low=0, high=250, ind_list=None,
                                  DATADIR='/scratch/gpfs/jiaxuanl/Data/Merian',
                                  CUTOUT_SUBDIR='./Cutout/magellan_spec/',
                                  CATALOG_SUBDIR='./Catalogs/magellan/'):
    print('SET ENVIRONMENT')
    os.chdir(DATADIR)
    print("CURRENT WORKING DIRECTORY:", os.getcwd())
    print('Number of processor to use:', njobs)

    # Load catalog
    print('Load catalog')

    obj_cat = Table.read(cat_dir)
    # obj_cat['dir'] = [file.replace('/projects/MERIAN/poststamps/g09_broadcut_new',
    #                                os.path.join(DATADIR, CUTOUT_SUBDIR),
    #                                ) for file in obj_cat['dir']]

    if low is None:
        low = 0
    if high is None:
        high = len(obj_cat)

    # Setup fitting
    channels = list('grizy') + ['N708', 'N540']
    ref_filt = 'i'

    # Initialize output catalog
    meas_cat = initialize_meas_cat(obj_cat, channels=channels)

    # Multiprocessing
    pwel = mp.Pool(njobs)
    iterable = range(low, high) if ind_list is None else ind_list
    meas_cat_tractor = pwel.map(
        partial(fitting_obj, obj_cat=obj_cat, meas_cat=meas_cat,
                channels=channels, ref_filt=ref_filt), iterable)

    if ind_list is not None:
        output_filename = f'{DATADIR}/{CATALOG_SUBDIR}/tractor_{suffix}_output_ind_list.pkl'
    else:
        output_filename = f'{DATADIR}/{CATALOG_SUBDIR}/tractor_{suffix}_output_{low}_{high}.pkl'

    with open(output_filename, 'wb') as f:
        pickle.dump(meas_cat_tractor, f)

    meas_cat_tractor = vstack(meas_cat_tractor)
    output_filename = output_filename.replace('.pkl', '.fits')
    print('#############')
    print('Writing to file:', output_filename)
    meas_cat_tractor.write(output_filename, overwrite=True)
    print('#############')


if __name__ == '__main__':
    fire.Fire(multiprocess_fitting_magellan)

# python tractor_fitting.py /scratch/gpfs/jiaxuanl/Data/Merian/Catalogs/magellan/magellan_spec_obj_cat.fits --suffix "magellan" --njobs 3 \
# --low 0 --high 10 --ind_list None --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
# --CUTOUT_SUBDIR './Cutout/magellan_spec/' --CATALOG_SUBDIR './Catalogs/magellan/' \
