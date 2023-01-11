'''
python script to deploy slurm jobs for running gaap
'''
import os
import sys
import fire


def deploy_training_job(tract, patch_min, patch_max, patch_jobs=None, filter_jobs=5,
                        python_file='lambo/scripts/hsc_gaap/run_gaap.py',
                        name='gaapCosmos', multijobs=False, email="am2907@princeton.edu", outname = "_train", 
                        repo='/projects/MERIAN/repo/'):
    ''' create slurm script and then submit 
    '''
    time = "1:30:00"
    # name = name
    python_file = os.path.join(os.getenv("LAMBO_HOME"), python_file)
    cntnt = '\n'.join([
        "#!/bin/bash",
        f"#SBATCH -J {name}_array_{tract}_{patch_min}_{patch_max}",
        "#SBATCH --nodes=1",
        f"#SBATCH --ntasks-per-node={filter_jobs+1}",  # {njobs}
        f"#SBATCH --mem={int(filter_jobs * 10)}G",
        "#SBATCH --time=%s" % time,
        f"#SBATCH --array={patch_min}-{patch_max}%5",
        "#SBATCH --export=ALL",
        f"#SBATCH -o ./log/{name}_array_{tract}_%a.o",
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
        f"python {python_file} --tract={tract} --patch_cols=\"[`expr $SLURM_ARRAY_TASK_ID % 9`]\" --patch_rows=\"[`expr $SLURM_ARRAY_TASK_ID / 9`]\" --patch_jobs={patch_jobs} --filter_jobs={filter_jobs} --hsc_type='S20A' --repo={repo}",
        "",
        "",
        'now=$(date +"%T")',
        'echo "end time ... $now"',
        ""])

    # create the slurm script execute it and remove it
    f = open(f'{outname}.slurm', 'w')
    f.write(cntnt)
    f.close()
    # os.system(f'sbatch {outname}.slurm')
    # os.system('rm _train.slurm')
    return None

if __name__ == '__main__':
    fire.Fire(deploy_training_job)

# Examples
# python deploy_gaap.py --tract=9813 --patch_cols="[0,1,2,3,4,5,6,7,8]" --patch_rows="[0]" --patch_jobs=None --filter_jobs=5
