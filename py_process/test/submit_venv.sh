#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hashtag574@gmail.com
#SBATCH --account=rrg-straub-ab
#SBATCH --job-name="py-test"
#SBATCH --nodes=1
#SBATCH --output=results
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0      # increase as needed
#SBATCH --time=0:15:00            # time (DD-HH:MM:SS)

# Set Path
SOURCEDIR=./

# Load python
module load python/3.8
# Prepare virtualenv
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip3 install --no-index --upgrade pip
# install packages using requirements script
pip3 install --no-index -r ~/backup/py_requirements.txt

# install independent packages
#pip install --no-index  ~/backup/py_packs/mpl_scatter_density-0.7-py3-none-any.whl
pip install   ~/backup/py_packs/MITgcmutils-0.1.2-py3-none-any.whl 

# run python and scripts
python $SOURCEDIR/py_test.py

