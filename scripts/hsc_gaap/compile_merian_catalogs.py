import os, sys
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from hsc_gaap.compile_science_catalogs import merge_merian_catalogs
import fire


#keep_merian_file = 'keep_table_columns_merian.txt'
#keep_gaap_file = 'keep_table_columns_gaap.txt'
#repo = '/scratch/gpfs/am2907/Merian/gaap/'


if __name__ == '__main__':
    fire.Fire(merge_merian_catalogs)
		
