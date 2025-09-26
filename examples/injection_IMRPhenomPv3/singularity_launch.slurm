#!/bin/sh

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --job-name=inject_ex
#SBATCH --partition=GravityTheory
#SBATCH --output=%x-%j.output
#SBATCH --error=%x-%j.error
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=<EMAIL>

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load singularity
cd $HOME/scratch                        # change to scratch
mkdir -p $SLURM_JOB_NAME/$SLURM_JOB_ID  # create job unique directory
cd $SLURM_JOB_NAME/$SLURM_JOB_ID        # switch to it

mkdir bin
mkdir build
mkdir data
mkdir src

# Copy in the data

gwatparent=$HOME/scratch/gwat
dir=$gwatparent/gw_analysis_tools/examples/injection_IMRPhenomPv3
dockerimage=$gwatparent/gwat_pv3.tar.gz

cp $dir/src/injection.cpp ./src/
cp $dir/makefile ./
cp $dir/runscript_for_slurm.sh ./ 

# Some people may need to put the cache in a specific place, as home directories don't always have enough space
export SINGULARITY_CACHEDIR=$HOME/scratch/.singularity/

# Launch the command. Note the bind option to mount the current directory (scratch project folder) to the container. The "pwd" option tells the container to run the following command from that directory in the container.
singularity exec --no-home  --bind $(pwd):/opt/project --pwd /opt/project docker://$dockerimage /bin/bash ./runscript_for_slurm.sh
