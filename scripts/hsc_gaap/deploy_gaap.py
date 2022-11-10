'''
python script to deploy slurm jobs for running gaap
'''
import os
import sys
import fire


def deploy_training_job(seed_low, seed_high, njobs=5, python_file='run_gaap.py',
                        name='gaapCosmos', multijobs=False):
    ''' create slurm script and then submit 
    '''
    time = "5:30:00"
    name = name + str(seed_low) + '_' + str(seed_high)

    cntnt = '\n'.join([
        "#!/bin/bash",
        f"#SBATCH -J NDE_%s_{seed_low}_{seed_high}" % (name),
        "#SBATCH --nodes=1",
        f"#SBATCH --ntasks-per-node=1",  # {njobs}
        "#SBATCH --mem=24G",
        "#SBATCH --time=%s" % time,
        "#SBATCH --export=ALL",
        f"#SBATCH --array={seed_low}-{seed_high}" if multijobs else "",
        f"#SBATCH -o ./log/gaap_{name}_{seed_low}_{seed_high}.o",
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
        ". /home/jiaxuanl/Research/Merian/merian_tractor/scripts/setup_env_gen3.sh",
        "",
        f"python {python_file} --seed_low={seed_low} --seed_high={seed_high} --n_jobs={njobs}",
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
    #os.system('rm _train.slurm')
    return None


if __name__ == '__main__':
    fire.Fire(deploy_training_job)

# 22.09.14

# python deploy_gaap.py --seed_low=20 --seed_high=25 --njobs=5
