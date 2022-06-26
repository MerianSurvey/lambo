"""
This file runs tractor fitting for objects in a give catalog.
"""
import fire
import os
import numpy as np
import copy
from astropy import wcs
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import pickle
import multiprocess as mp
mp.freeze_support()
from functools import partial

from kuaizi.detection import Data
# from kuaizi.utils import padding_PSF
from kuaizi.tractor.utils import initialize_meas_cat, _write_to_row
from kuaizi.tractor.fit import tractor_hsc_sep_blob_by_blob

import sys
sys.path.append('/home/jiaxuanl/software/astrometry.net-0.89')
sys.path.append('/home/jiaxuanl/Research/Packages/tractor/')
sys.path.append('/home/jiaxuanl/Research/Packages/carpenter/src/')

import matplotlib.pyplot as plt

import kuaizi
kuaizi.set_matplotlib(style='JL', usetex=True, dpi=70)
from carpenter.display import display_merian_cutout_rgb


def padding_PSF(psf_hdu):
    '''
    If the sizes of HSC PSF in all bands are not the same, this function pads the smaller PSFs.

    The output size must be odd.


    Parameters:
        psf_list: a list returned by `unagi.task.hsc_psf` function

    Returns:
        psf_pad: a list including padded PSFs. They now share the same size.
    '''
    # Padding PSF cutouts from HSC
    max_len = max(psf_hdu[0].data.shape)
    if max_len % 2 == 0:
        max_len += 1

    y_len, x_len = psf_hdu[0].data.shape
    dy = ((max_len - y_len) // 2, (max_len - y_len) // 2)
    dx = ((max_len - x_len) // 2, (max_len - x_len) // 2)

    if (max_len - y_len) == 1:
        dy = (1, 0)
    if (max_len - x_len) == 1:
        dx = (1, 0)

    temp = np.pad(psf_hdu[0].data.astype('float'),
                  (dy, dx), 'constant', constant_values=0)

    if temp.shape == (max_len, max_len):
        return temp
    else:
        raise ValueError('Wrong size!')


def fitting_obj(index, obj_cat, meas_cat, point_source=False,
                channels=['g', 'r', 'i', 'z', 'N708', 'N540'], ref_filt='r'):
    """
    Wrapper for fitting one object.
    """
    row = meas_cat[index]
    obj = obj_cat[index]

    ID_name = 'ID' if 'ID' in obj_cat.colnames else 'id'
    PREFIX_name = 'PREFIX' if 'PREFIX' in obj_cat.colnames else 'prefix'

    row['ID'] = obj[ID_name]

    # [filt for filt in channels if filt != ref_filt]
    forced_channels = channels

    if isinstance(point_source, str):
        point_source = (point_source.lower() == 'true')

    print(f'### Tractor modeling for obj {index}, ID = {obj[ID_name]}')
    if point_source:
        print('    Point source mode')

    try:
        # Load data
        cutout = []
        for filt in channels:
            # if 'N' in filt:
            #     cutout.append(fits.open(os.path.join(
            #         obj['dir'], f"{obj[PREFIX_name]}_{filt}_deepCoadd_calexp.fits")))
            # else:
            #     cutout.append(fits.open(os.path.join(
            #         obj['dir'], f's18a_wide_{obj[ID_name]}_{filt}.fits')))
            cutout.append(fits.open(os.path.join(
                obj['dir'], f"{obj['prefix']}_{filt}_deepCoadd_calexp.fits")))

        psf_list = []
        for filt in channels:
            # if 'N' in filt:
            #     psf_list.append(fits.open(os.path.join(
            #         obj['dir'], f"{obj[PREFIX_name]}_{filt}_psf.fits")))
            # else:
            #     psf_list.append(fits.open(os.path.join(
            #         obj['dir'], f's18a_wide_{obj[ID_name]}_{filt}_psf.fits')))
            psf_list.append(fits.open(os.path.join(
                obj['dir'], f"{obj['prefix']}_{filt}_psf.fits")))

        # Reconstruct data
        images = np.array([hdu[1].data for hdu in cutout if len(hdu) > 1])
        # note: all bands share the same WCS here
        w = wcs.WCS(cutout[0][1].header)
        weights = 1 / np.array([hdu[3].data for hdu in cutout])
        psf_pad = [padding_PSF(_psf) for _psf in psf_list]
        # padding_PSF(psf_list)  # Padding PSF cutouts from HSC
        data = Data(images=images, weights=weights, wcs=w,
                    psfs=psf_pad, channels=channels)
    except Exception as e:
        print(f'IMCOMPLETE FILE FOR {index} = obj {obj[ID_name]}', e)

    # Start fitting
    try:
        model_dict = {}
        # fitting in the i-band first: then pass the i-band parameters of target galaxy to other bands
        model_dict[ref_filt], _obj_cat_i = tractor_hsc_sep_blob_by_blob(
            obj, ref_filt, data.channels, data, point_source=point_source,
            freeze_dict={'pos': False, 'shape': False, 'shape.re': False, 'shape.ab': False, 'shape.phi': False,
                         'sersicindex': False},  # don't fix shape/sersic
            verbose=False)

        temp = model_dict[ref_filt].catalog.copy()
        temp = temp[model_dict[ref_filt].target_ind:model_dict[ref_filt].target_ind + 1]

        for filt in forced_channels:
            pos = True
            fix_all = True
            # ref_catalog = model_dict[ref_filt].catalog.copy(
            # )[model_dict[ref_filt].target_ind:model_dict[ref_filt].target_ind + 1]
            model_dict[filt], _ = tractor_hsc_sep_blob_by_blob(
                obj, filt, data.channels, data,
                fix_all=fix_all, tractor_cat=copy.deepcopy(temp),
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

    try:
        # Visualize our model
        channels = list('griz') + ['N708', 'N540']
        from kuaizi.tractor.utils import HiddenPrints
        with HiddenPrints():
            # model_img = np.asarray([
            #     model_dict[key].getModelImage(
            #         0, srcs=model_dict[key].catalog[model_dict[key].target_ind:model_dict[key].target_ind + 1]
            #     ) for key in channels])
            model_img = np.asarray([
                model_dict[key].getModelImage(
                    0, srcs=model_dict[key].catalog
                ) for key in channels])
        fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(11, 5))
        Q = 1.5
        stretch = 0.7
        _img = np.vstack([data.images[0:4], data.images[-2:]])
        _, img_rgb = display_merian_cutout_rgb(_img, filters=list('griz') + ['N708', 'N540'],
                                               ax=ax1, color_norm=None, Q=Q, stretch=stretch,
                                               channel_map=None, N708_strength=3)

        _, _ = display_merian_cutout_rgb(model_img, filters=list('griz') + ['N708', 'N540'],
                                         ax=ax2, Q=Q, stretch=stretch,
                                         color_norm=None,
                                         channel_map=None)

        _, _ = display_merian_cutout_rgb(_img - model_img, filters=list('griz') + ['N708', 'N540'],
                                         ax=ax3, Q=Q, stretch=stretch,
                                         color_norm=None,
                                         channel_map=None)
        plt.subplots_adjust(wspace=0.03)
        plt.tight_layout()
        plt.savefig(
            f"/tigress/jiaxuanl/public_html/Merian/stars/{obj['prefix']}_tractor.png")
        plt.close()
    except Exception as e:
        print(f'VISUALIZATION ERROR FOR {index} = obj {obj[ID_name]}', e)

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
def multiprocess_fitting_magellan(cat_dir, suffix='', point_source=True, njobs=16, low=0, high=250, ind_list=None,
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
    obj_cat['id'] = obj_cat['name']
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
                channels=channels, ref_filt=ref_filt, point_source=point_source), iterable)

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

# python tractor_fitting.py /scratch/gpfs/jiaxuanl/Data/Merian/Cutout/stars/stars-2022-05-29.fits \
# --point_source "True" --suffix "stars" --njobs 1 \
# --low 420 --high 421 --ind_list None --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
# --CUTOUT_SUBDIR './Cutout/stars/' --CATALOG_SUBDIR './Catalogs/stars/' \
