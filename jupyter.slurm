#!/bin/bash -login
#SBATCH -p med2                # partition, or queue, to assign to
#SBATCH -J Jupyter            # name for job
#SBATCH -N 1                  # one "node", or computer
#SBATCH -n 1                  # one task for this node
#SBATCH -c 8                 # cores per task
#SBATCH -t 02-00:00:00        # ask for 2 days
#SBATCH --mem=20000           # memory (20,000 mb = 20gb)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mabuelanin@gmail.com

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate cattle_sv2

# fail on weird errors
set -o nounset
set -o errexit
set -x

#Run Jupyter
jupyter lab --no-browser --ip "*" --notebook-dir /home/mshokrof/workshop_12Jan_2023/SV_calling_LR/

# print out various information about the job
env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
