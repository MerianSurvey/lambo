import sys, os
import fire
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from hsc_gaap.gaap import consolidateObjectTables

def saveTractTable(tract, repo="/scratch/gpfs/am2907/Merian/gaap/",  hsc_type = "S20A",
                   copy_to_projects=False, delete=False):

    """Stack and save all of the gaap catalogs for a given tract"""
    stacked_tables = consolidateObjectTables(patches=[], tract=tract, repo=repo)
    print("Stacked")
    filename = os.path.join(repo, f'{hsc_type}/gaapTable/{tract}/objectTable_{tract}_{hsc_type}.fits')
    if not delete:
         stacked_tables.write(filename, overWrite=True)
         print ("Wrote table to {filename}")
    if copy_to_projects:
        filename = os.path.join("/projects/MERIAN/repo/", f'{hsc_type}/gaapTable/{tract}/objectTable_{tract}_{hsc_type}.fits')
        stacked_tables.write(filename, overWrite=True)
        print ("Wrote table to {filename}")
        


def stackAndClean(tract, stack=True, copy_to_projects=False, clean=False, 
                    repo="/scratch/gpfs/am2907/Merian/gaap/", hsc_type="S20A"):

    """Stack and save all of the gaap catalogs for a given tract and erase all the
    used and moved images and catalogs"""
    if stack:
        saveTractTable(tract, repo, hsc_type, copy_to_projects, delete=clean)

    if clean:
        os.system(f'rm -r {repo}/{hsc_type}/deepCoadd_calexp/{tract}')
        os.system(f'rm -r {repo}/{hsc_type}/gaapTable/{tract}')

if __name__ == '__main__':
    fire.Fire(stackAndClean)
