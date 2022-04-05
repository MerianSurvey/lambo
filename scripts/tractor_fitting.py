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
from kuaizi.tractor.utils import tractor_hsc_sep_blob_by_blob, initialize_meas_cat, _write_to_row
kuaizi.set_env(project='Merian', name='',
               data_dir='/scratch/gpfs/jiaxuanl/Data')


import sys
sys.path.append('/home/jiaxuanl/software/astrometry.net-0.89')
sys.path.append('/home/jiaxuanl/Research/Packages/tractor/')

# Load catalog
obj_cat = Table.read('./Catalogs/COSMOS_cutouts_tractor_gaap_test.fits')
obj_cat['dir'] = [file.replace('/projects/MERIAN/poststamps/g09_broadcut',
                               '/scratch/gpfs/jiaxuanl/Data/Merian/Cutout'
                               ) for file in obj_cat['dir']]
# Initialize output catalog
meas_cat = initialize_meas_cat(obj_cat)

# Setup fitting
channels = 'grizyN'
ref_filt = 'r'
forced_channels = [filt for filt in channels if filt != ref_filt]


def fitting_obj(index, obj_cat, meas_cat,):
    """
    Wrapper for fitting one object.
    """
    row = meas_cat[index]
    obj = obj_cat[index]
    row['ID'] = obj['ID']

    print(f'### Tractor modeling for obj {index}, ID = {obj["ID"]}')
    try:
        cutout = [fits.open(os.path.join(
            obj['dir'], f's18a_wide_{obj["ID"]}_{filt}.fits')) for filt in 'grizy']
        cutout += [fits.open(os.path.join(obj['dir'],
                                          f"{obj['PREFIX']}_N708_deepCoadd.fits"))]

        psf_list = [fits.open(os.path.join(
            obj['dir'], f's18a_wide_{obj["ID"]}_{filt}_psf.fits')) for filt in 'grizy']
        psf_list += [fits.open(os.path.join(obj['dir'],
                                            f"{obj['PREFIX']}_N708_psf.fits"))]
        # Reconstruct data
        from kuaizi.detection import Data
        from kuaizi.utils import padding_PSF

        images = np.array([hdu[1].data for hdu in cutout if len(hdu) > 1])
        # note: all bands share the same WCS here
        w = wcs.WCS(cutout[0][1].header)
        weights = 1 / np.array([hdu[3].data for hdu in cutout])
        psf_pad = padding_PSF(psf_list)  # Padding PSF cutouts from HSC
        data = Data(images=images, weights=weights, wcs=w,
                    psfs=psf_pad, channels=channels)
    except Exception as e:
        print(f'IMCOMPLETE FILE FOR {index} = obj {obj["ID"]}', e)

    # Start fitting
    try:
        model_dict = {}
        # fitting in the i-band first: then pass the i-band parameters of target galaxy to other bands
        model_dict[ref_filt] = tractor_hsc_sep_blob_by_blob(
            obj, ref_filt, data.channels, data,
            freeze_dict={'pos': False, 'shape': False,
                         'sersicindex': False},  # don't fix shape/sersic
            verbose=False)

        for filt in forced_channels:
            pos = True  # False if filt == 'N' else True
            model_dict[filt] = tractor_hsc_sep_blob_by_blob(
                obj, filt, data.channels, data,
                ref_source=model_dict[ref_filt].catalog[model_dict[ref_filt].target_ind],
                freeze_dict={'pos': pos, 'shape': True,
                             'sersicindex': True},  # don't fix shape/sersic
                verbose=False,)
        with open(os.path.join(obj['dir'], f'cosmos_{obj["ID"]}_tractor.pkl'), 'wb') as f:
            pickle.dump(model_dict, f)

        # Write to catalog
        row = _write_to_row(row, model_dict)
        return row

    except Exception as e:
        print(f'FITTING ERROR FOR {index} = obj {obj["ID"]}', e)
        return


# Main
ncpu = 12
low = 0
high = 200
pwel = mp.Pool(ncpu)
meas_cat_tractor = pwel.map(
    partial(fitting_obj, obj_cat=obj_cat, meas_cat=meas_cat), range(low, high))
meas_cat_tractor = vstack(meas_cat_tractor)

meas_cat_tractor.write(
    f'./Catalogs/tractor_test_output_{low}_{high}.fits', overwrite=True)
