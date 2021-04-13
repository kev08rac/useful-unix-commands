#!/bin/bash
#
#SBATCH --job-name=new-job
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=username@mail.edu
#SBATCH -n 1
#SBATCH -t 0-00:60 # Runtime in D-HH:MM
#SBATCH --output=new-job-%A_%a.out
#SBATCH --array=1-10

srun -N 1 -c 8 -t 2:00:00 --pty bash # Request an interactive node with 8 cores for 2 hours

# This is a very basic template script for bash on a remote server using parallelization

trap '(read -p "[$BASH_SOURCE:$LINE0] $BASH_COMMAND?")' DEBUG # FOR DEBUGGING: Turn off set -u for this to work
set -euo pipefail
module load gcc/9.2.0

path = /path/to/directory
job_num = ${SLURM_ARRAY_TASK_ID}

echo $path
echo $job_num
