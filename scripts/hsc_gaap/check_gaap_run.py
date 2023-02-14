import os, sys
import glob
import numpy as np
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
import fire

from hsc_gaap.gaap import findReducedPatches

def checkRun(tract, repo = '/scratch/gpfs/am2907/Merian/gaap/', output=True):
    if output:
        print(f'TRACT: {tract}')
    problem_patches=[]

    patches = findReducedPatches(tract)
    dir = os.path.join(repo, f"log/{tract}/")
    logs = np.array(glob.glob(dir + "*.o"))
    logs_patches = np.array([log.split("/")[-1].split(".")[0] for log in logs]).astype(int)
    logs, logs_patches = logs[logs_patches.argsort()], logs_patches[logs_patches.argsort()]

    if len(logs) != len(patches):
        print("MISSING PATCHES")

    error1 = "Wrote GAaP table to "
    error2 = "Running GAaP only on filters"
    error3 = "Matched GAaP table with S20A blendedness"
    error4 = "not found for all filters"
    for i, log in enumerate(logs):
        logfile = open(log, "r").read()
        
        if error4 in logfile:
            if output:
                print(f'PROBLEM IN PATCH {logs_patches[i]}: Missing all filter images ')
            problem_patches.append(logs_patches[i])
            continue

        if error1 not in logfile:
            if output:
                print (f'PROBLEM IN PATCH {logs_patches[i]}: Catalog not saved')
            problem_patches.append(logs_patches[i])

        if (error2 in logfile):
            if output:
                print (f'PROBLEM IN PATCH {logs_patches[i]}: Missing some filter images')
            problem_patches.append(logs_patches[i])
        
        elif (logfile.count("Finished the GAaP measureTask for band") < 5):
            if output:
                print (f'PROBLEM IN PATCH {logs_patches[i]}: Failed for {5 - logfile.count("Finished the GAaP measureTask for band")} bands')
            problem_patches.append(logs_patches[i])

        if error3 not in logfile:
            if output:
                print (f'PROBLEM IN PATCH {logs_patches[i]}: Blend matching failed')
            problem_patches.append(logs_patches[i])

    if (len(problem_patches)==0) & (len(logs) == len(patches)):
        print ("NO PROBLEMS")
    print()

    return np.unique(problem_patches)

if __name__ == '__main__':
    fire.Fire(checkRun)