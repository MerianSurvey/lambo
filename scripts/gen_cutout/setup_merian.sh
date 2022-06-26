# setup environment for generating Merian cutouts
# ssh tiger2-sumire

module purge
module load rh/devtoolset/6

LSST_CONDA_ENV_NAME=lsst-scipipe-3.0.0
source "/projects/HSC/LSST/stack/loadLSST.sh"

# source /projects/HSC/LSST/stack/loadLSST.bash
setup lsst_distrib
setup lsst_apps
setup obs_subaru