#!/bin/bash
#SBATCH -J cpujob
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=3
#SBATCH --output=test.out
#SBATCH --error=test.out
#SBATCH --array=0

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
cd /home/yutong/software/
root=/home/yutong/data/MGH_ICH/results/
singularity exec --bind $root:/opt/data /bin/bash /opt/preprocess.sh \
    $SLURM_ARRAY_TASK_COUNT $SLURM_ARRAY_TASK_ID $SLURMS_CPU_PER_TASKS
