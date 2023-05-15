import os, sys
import glob
import numpy as np
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from astropy.table import Table, vstack, hstack
import lsst.daf.butler as dafButler
from hsc_gaap.gaap import findReducedPatches, consolidateMerianCats
from hsc_gaap.compile_catalogs import merge_merian_catalogs
import fire


#keep_merian_file = 'keep_table_columns_merian.txt'
#keep_gaap_file = 'keep_table_columns_gaap.txt'
#repo = '/scratch/gpfs/am2907/Merian/gaap/'


if __name__ == '__main__':
    fire.Fire(merge_merian_catalogs)
		
