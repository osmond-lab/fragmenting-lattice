#!/bin/sh
#
# description: get everything set up after logging on
#

# load python
module load NiaEnv/2022a python/3.11.5

# make directory for virtual environments if havent already
#mkdir ~/.virtualenvs

# create virtual environment, if havent already
myenv="frag"
#virtualenv --system-site-packages ~/.virtualenvs/$myenv

# activate virtual environment
source ~/.virtualenvs/$myenv/bin/activate 

# set up for jupyter
#pip install ipykernel
#python -m ipykernel install --name $myenv --user
#venv2jup

# install snakemake in the environment if havent already
#pip install snakemake==8.10.8

# snakemake write to write-able location
export XDG_CACHE_HOME=$SCRATCH

# libraries for Snakefile
#pip install snakemake-executor-plugin-slurm
# note: i had to comment out four lines of this plugin to prevent mem requests, as niagara does not allow them
