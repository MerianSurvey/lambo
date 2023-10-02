import os, sys
import glob
import numpy as np
from astropy.table import Table
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))

from tqdm import tqdm
from hsc_gaap.gaap import findReducedPatches
from hsc_gaap.check_gaap_run import checkRun


if __name__ == '__main__':

    blendcols = ['g_blendedness',
            'r_blendedness',
            'i_blendedness',
            'z_blendedness',
            'y_blendedness',
            'g_blendedness_flag',
            'r_blendedness_flag',
            'i_blendedness_flag',
            'z_blendedness_flag',
            'y_blendedness_flag',
            'g_apertureflux10_S20A',
            'r_apertureflux10_S20A',
            'i_apertureflux10_S20A',
            'z_apertureflux10_S20A',
            'y_apertureflux10_S20A',
            'g_apertureflux10_flag_S20A',
            'r_apertureflux10_flag_S20A',
            'i_apertureflux10_flag_S20A',
            'z_apertureflux10_flag_S20A',
            'y_apertureflux10_flag_S20A']
    
    repo = "/scratch/gpfs/am2907/Merian/gaap"

    cat_paths = glob.glob(os.path.join(repo, "S20A/gaapTable/", "**/objectTable*,*"), recursive=True)


    for path in tqdm(cat_paths):
        try:
            cat = Table.read(path)
        except:
            print(path)
        try:
            cat.remove_columns(blendcols)
            cat.write(path, overwrite=True)
        except:
            continue
