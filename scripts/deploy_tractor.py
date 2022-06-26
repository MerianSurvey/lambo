import os
import sys
import fire

SUFFIX = 'stars'  # 'lowz_calexp'
CAT_DIR = "/scratch/gpfs/jiaxuanl/Data/Merian/Cutout/stars/stars-2022-05-29.fits"


def deploy_modeling_job(low=0, high=100, ind_list=None, name='trctr_job', ncpu=12, njobs=44,
                        point_source=False,
                        DATADIR='/scratch/gpfs/jiaxuanl/Data/Merian',
                        CUTOUT_SUBDIR='./Cutout/stars/',
                        CATALOG_SUBDIR='./Catalogs/'):
    ''' 
    create slurm script and then submit 
    '''
    time = "11:59:00"
    run_tractor_content = '\n'.join([
        f"python tractor_fitting.py {CAT_DIR} --suffix {SUFFIX} --njobs {njobs} --point_source {point_source} \\",
        f"--low {low} --high {high} --ind_list {ind_list} --DATADIR {DATADIR} \\",
        f"--CUTOUT_SUBDIR {CUTOUT_SUBDIR} --CATALOG_SUBDIR {CATALOG_SUBDIR} \\"
    ])

    cntnt = '\n'.join([
        "#!/bin/bash",
        f"#SBATCH -J {name}_{low}_{high}" if ind_list is None else f"#SBATCH -J {name}_ind_list",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks-per-node=%i" % ncpu,
        "#SBATCH --mem-per-cpu=2G",
        "#SBATCH --time=%s" % time,
        "#SBATCH --export=ALL",
        f"#SBATCH -o {name}_{low}_{high}.o" if ind_list is None else f"#SBATCH -o {name}_ind_list.o",
        "#SBATCH --mail-type=all",
        "#SBATCH --mail-user=jiaxuanl@princeton.edu",
        "",
        'now=$(date +"%T")',
        'echo "start time ... $now"',
        "",
        "module purge",
        "source /home/jiaxuanl/Research/Merian/merian_tractor/scripts/setup_env.sh",
        "export OMP_NUM_THREADS=1",
        "",
        run_tractor_content,
        "",
        "",
        "",
        'now=$(date +"%T")',
        'echo "end time ... $now"',
        ""])

    # create the slurm script execute it and remove it
    if not os.path.isdir('./slurm'):
        os.mkdir('./slurm')
    f = open(f'./slurm/_{name}_{low}_{high}.slurm', 'w') if ind_list is None else open(
        f'./slurm/_{name}_ind_list.slurm', 'w')
    f.write(cntnt)
    f.close()
    os.system(f'sbatch ./slurm/_{name}_{low}_{high}.slurm') if ind_list is None else os.system(
        f'sbatch ./slurm/_{name}_ind_list.slurm')
    return None


if __name__ == '__main__':
    fire.Fire(deploy_modeling_job)


############ EXAMPLE #############
# python deploy_tractor.py --low 0 --high 500 --ncpu 15 --njobs 54 --name trctr_job1
# python deploy_tractor.py --low 500 --high 1000 --ncpu 15 --njobs 54 --name trctr_job2
# python deploy_tractor.py --low 1000 --high 1500 --ncpu 15 --njobs 54 --name trctr_job3
# python deploy_tractor.py --low 1500 --high None --ncpu 15 --njobs 54 --name trctr_job4

# python deploy_tractor.py --low 500 --high 510 --ncpu 1 --njobs 3 --name trctr_job

# python deploy_tractor.py --low 0 --high 400 --ncpu 15 --njobs 54 --name trctr_job --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
# --CUTOUT_SUBDIR './Cutout/magellan_spec/' --CATALOG_SUBDIR './Catalogs/magellan/'

# python deploy_tractor.py --low 400 --high None --ncpu 15 --njobs 54 --name trctr_job --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
# --CUTOUT_SUBDIR './Cutout/magellan_spec/' --CATALOG_SUBDIR './Catalogs/magellan/'

# python deploy_tractor.py --low 0 --high None --ncpu 15 --njobs 54 --name trctr_job --DATADIR /scratch/gpfs/jiaxuanl/Data/Merian \
# --CUTOUT_SUBDIR './Cutout/stars/' --CATALOG_SUBDIR './Catalogs/stars/' --point_source "True"
