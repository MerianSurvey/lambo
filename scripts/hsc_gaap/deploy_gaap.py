'''
python script to deploy slurm jobs for running gaap
'''
import os
import sys
import fire


def deploy_training_job(patch_cols, patch_rows, patch_jobs=None, filter_jobs=5,
                        python_file='run_gaap.py',
                        name='gaapCosmos', multijobs=False):
    ''' create slurm script and then submit 
    '''
    time = "12:30:00"
    # name = name

    cntnt = '\n'.join([
        "#!/bin/bash",
        f"#SBATCH -J NDE_%s_{patch_rows}" % (name),
        "#SBATCH --nodes=1",
        f"#SBATCH --ntasks-per-node={filter_jobs+1}",  # {njobs}
        f"#SBATCH --mem={int(filter_jobs * 10)}G",
        "#SBATCH --time=%s" % time,
        "#SBATCH --export=ALL",
        # f"#SBATCH --array={seed_low}-{seed_high}" if multijobs else "",
        f"#SBATCH -o ./log/gaap_{name}_{patch_rows}.o",
        "#SBATCH --mail-type=all",
        "#SBATCH --mail-user=jiaxuanl@princeton.edu",
        "",
        'now=$(date +"%T")',
        'echo "start time ... $now"',
        'echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."',
        'echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"',
        'echo "Executing on the machine:" $(hostname)',
        "",
        "module purge",
        ". /home/jiaxuanl/Research/Merian/merian_tractor/scripts/setup_env_w40.sh",
        "",
        f"python {python_file} --patch_cols='{patch_cols}' --patch_rows='{patch_rows}' --patch_jobs={patch_jobs} --filter_jobs={filter_jobs} --hsc_type='S20A'",
        "",
        "",
        'now=$(date +"%T")',
        'echo "end time ... $now"',
        ""])

    # create the slurm script execute it and remove it
    f = open('_train.slurm', 'w')
    f.write(cntnt)
    f.close()
    os.system('sbatch _train.slurm')
    # os.system('rm _train.slurm')
    return None


if __name__ == '__main__':
    fire.Fire(deploy_training_job)

# 22.09.14
# python deploy_gaap.py --patch_cols="[0,1,2,3,4,5,6,7,8]" --patch_rows="[7]" --patch_jobs=None --filter_jobs=5
# python deploy_gaap.py --patch_cols="[0,1,2,3,4,5,6,7,8]" --patch_rows="[8]" --patch_jobs=None --filter_jobs=5
# python deploy_gaap.py --patch_cols="[3,4,5,6,7,8]" --patch_rows="[3]" --patch_jobs=None --filter_jobs=5
# python deploy_gaap.py --patch_cols="[3,4,5,6,7,8]" --patch_rows="[4]" --patch_jobs=None --filter_jobs=5
# python deploy_gaap.py --patch_cols="[7]" --patch_rows="[4,8]" --patch_jobs=None --filter_jobs=2
