Scripts for integrating BFM in FABM


To use these with in seamless notebooks 
clone bfmforfabm in folder in extern

git clone --recurse-submodules git@github.com:inogs/bfmforfabm.git ogs
git checkout dev_1D_diverse


go to fabm submodule and checkout the neccton branch
git checkout neccton
go to gotm submodule and check that you are using branch v6.0
idifferently from before no modifications of gotm are needed

+++++++++++++++++++++++++++++++
+++++++++++++++++++++++++++++++


use the following my_install

# This script is intended to be source'd, not executed

set -e

REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#export FFLAGS='-fcheck=all -O0 -Wall -Wextra -g -cpp -DDEBUG -DQUICK'

# Add additional FABM "institutes" (bfm, ecosmo, etc.) and their base directory on the line below.
FABM_ARGS="-DFABM_INSTITUTES=gotm;msi;ersem;pisces;ogs -DFABM_ERSEM_BASE=${REPO_DIR}/extern/ersem -DFABM_PISCES_BASE=${REPO_DIR}/extern/pisces -DFABM_OGS_BASE=${REPO_DIR}/extern/ogs -DCMAKE_Fortran_COMPILER
=gfortran  -DCMAKE_BUILD_TYPE=release"

# Build pyfabm
#WORK_DIR=`mktemp -d`
#cd ${WORK_DIR}
#cmake ${REPO_DIR}/extern/fabm/src/drivers/python $FABM_ARGS
#make -j4 install
#cd -
#rm -rf ${WORK_DIR}

# Build gotm-fabm
WORK_DIR=`mktemp -d`
#WORK_DIR=$CINECA_SCRATCH/TEST_NECCTON
#mkdir $WORK_DIR
cd ${WORK_DIR}
###cmake ${REPO_DIR}/extern/gotm -DFABM_BASE=${REPO_DIR}/extern/fabm $FABM_ARGS

### For building the executable if everything (gotm,fabm,ogs) is fine in the repo dir (dev_1D): 
cmake ${REPO_DIR}/extern/gotm -DFABM_BASE=${REPO_DIR}/extern/fabm $FABM_ARGS 

### For building the executable that works (better) for dev_1D_diverse:
###cmake /g100_work/OGS_devC/V10C/RUNS_SETUP/SEAMLESS/eat/extern/gotm -DFABM_BASE=/g100/home/userexternal/ealvarez/seamless-notebooks_div/extern/fabm $FABM_ARGS -DADJ_BIOPTIMOD_LIBRARIES="/g100_scratch/user
external/ealvarez/seamless-notebooks/extern/Forward_Adjoint/lib"


make -j4 VERBOSE=1
cp -v gotm /g100_work/OGS23_PRACE_IT/plazzari/seamless-notebooks-neccton-test/bin/
cd -
rm -rf ${WORK_DIR}

## Build eat
#WORK_DIR=`mktemp -d`
#cd ${WORK_DIR}
#cmake ${REPO_DIR}/extern/eat -DFABM_BASE=${REPO_DIR}/extern/fabm $FABM_ARGS
#make -j4 install
#cd -
#rm -rf ${WORK_DIR}

#cd ../..
