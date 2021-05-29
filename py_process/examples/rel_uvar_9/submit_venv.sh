#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hashtag574@gmail.com
#SBATCH --account=rrg-straub-ab
#SBATCH --job-name="py-ruvar9-cnt"
#SBATCH --nodes=1
#SBATCH --output=results
#SBATCH --ntasks-per-node=1
#SBATCH --mem=96GB   # increase as needed
#SBATCH --time=4:00:00            # time (DD-HH:MM:SS)

# Set Path
ENVDIR=/home/tpeng/ENV/bin/
SOURCEDIR=/home/tpeng/MITgcm_post/py_process/


GROUPNAME=/
CASENAME=rel_uvar_9
START_IND=480000
END_IND=550000
LOC_STR=cnt

# Load python
module load python/3.8
# Prepare virtualenv
source ${ENVDIR}/activate

# run python and scripts
python ${SOURCEDIR}/read_ext_plt_rty_spec.py  $GROUPNAME $CASENAME $START_IND $END_IND $LOC_STR 

