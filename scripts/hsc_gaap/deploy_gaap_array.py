'''
python script to deploy slurm jobs for running gaap
'''
import os
import sys
import fire
import numpy as np
import lsst.daf.butler as dafButler
sys.path.append(os.path.join(os.getenv('LAMBO_HOME'), 'lambo/scripts/'))
from hsc_gaap.gaap import findReducedPatches
from hsc_gaap.check_gaap_run import checkRun

def deploy_training_job(tract, filter_jobs=5,
                        python_file='lambo/scripts/hsc_gaap/run_gaap.py',
                        name='gaapCosmos', multijobs=False, email="am2907@princeton.edu", outname = None, 
                        repo='/projects/MERIAN/repo/', submit=False, fixpatches=False):
        
    ''' Create slurm script to process all patches already reduced in Merian for a given tract.
    '''

    # Get a list of reduced patches for a given tract  
    if not fixpatches:
        patches = findReducedPatches(tract)
    else:
        patches = checkRun(tract, repo=repo, output=False)

    time = "1:30:00"
    # name = name
    python_file = os.path.join(os.getenv("LAMBO_HOME"), python_file)
    cntnt = '\n'.join([
        "#!/bin/bash",
        f"#SBATCH -J {tract}_gaapCosmos",
        "#SBATCH --nodes=1",
        f"#SBATCH --ntasks-per-node={filter_jobs+1}",  # {njobs}
        f"#SBATCH --mem={int(filter_jobs * 10)}G",
        "#SBATCH --time=%s" % time,
        f"#SBATCH --array={','.join(patches.astype(str))}",
        "#SBATCH --export=ALL",
        f"#SBATCH -o ./log/{tract}/%a.o",
        "#SBATCH --mail-type=all",
        f"#SBATCH --mail-user={email}",
        "",
        'now=$(date +"%T")',
        'echo "start time ... $now"',
        'echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."',
        'echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"',
        'echo "Executing on the machine:" $(hostname)',
        "",
        "module purge",
        f". {os.getenv('LAMBO_HOME')}/lambo/scripts/setup_env_w40.sh",
        f"export LAMBO_HOME='{os.getenv('LAMBO_HOME')}'",
        "",
        f"python {python_file} --tract={tract} --patch_cols=\"[`expr $SLURM_ARRAY_TASK_ID % 9`]\" --patch_rows=\"[`expr $SLURM_ARRAY_TASK_ID / 9`]\" --filter_jobs={filter_jobs} --hsc_type='S20A' --repo={repo}",
        "",
        "",
        'now=$(date +"%T")',
        'echo "end time ... $now"',
        ""])

    # create the slurm script execute it and remove it
    if outname is None:
        outname = str(tract)
        if fixpatches:
            outname += "_fixing"
    if not os.path.isdir('./slurmscripts'):
        os.mkdir("./slurmscripts")
    if not os.path.isdir(f'./log/{tract}'):
        os.makedirs(f'./log/{tract}')
    f = open(f'slurmscripts/{outname}.slurm', 'w')
    f.write(cntnt)
    f.close()
    if submit:
        os.system(f'sbatch ./slurmscripts/{outname}.slurm')
    return None
    
if __name__ == '__main__':
    fire.Fire(deploy_training_job)


# Examples
# python3 deploy_gaap_array.py --tract=8520  --repo="/scratch/gpfs/am2907/Merian/gaap" --submit=True