'''
python script to deploy slurm jobs for running gaap
'''
import os
import sys
import fire
import numpy as np
import lsst.daf.butler as dafButler

def deploy_training_job(tract, filter_jobs=5,
                        python_file='lambo/scripts/hsc_gaap/run_gaap.py',
                        name='gaapCosmos', multijobs=False, email="am2907@princeton.edu", outname = None, 
                        repo='/projects/MERIAN/repo/', submit=False):
        
    ''' Create slurm script to process all patches already reduced in Merian for a given tract.
    '''

    # Get a list of reduced patches for a given tract

    output_collection = "DECam/runs/merian/dr1_wide"
    data_type = "deepCoadd_calexp"
    skymap = "hsc_rings_v1"

    butler = dafButler.Butler('/projects/MERIAN/repo/', collections=output_collection, skymap=skymap)
    
    patches = np.array([[data_id['tract'], data_id["patch"]] for data_id in butler.registry.queryDataIds (['tract','patch'], datasets=data_type, 
                                                    collections=output_collection, skymap=skymap)])
    patches = patches[patches[:, 0].argsort()]
    del butler 

    patches = np.unique(patches[patches[:,0] == tract][:,1])

    time = "1:30:00"
    # name = name
    python_file = os.path.join(os.getenv("LAMBO_HOME"), python_file)
    cntnt = '\n'.join([
        "#!/bin/bash",
        f"#SBATCH -J {name}_array_{tract}",
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
        outname = tract
    if not os.path.isdir('./slurmscripts'):
        os.mkdir("./slurmscripts")
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