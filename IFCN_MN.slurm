#!/bin/bash
#
#SBATCH --mem=10G
#SBATCH --job-name=IFCN_MN
#SBATCH --array=1-250
#SBATCH --output=/IFCN_MN/%a.out
#SBATCH --error=/IFCN_MN/%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH -c 1

uname -n 1>&1

# Adjust path as needed
#cd /TMS_Thresholding

Batch_id=$(($SLURM_ARRAY_TASK_ID))

# Call MATLAB function, modify path for MATLAB as needed
/usr/bin/time -v /admin/apps/rhel7/matlabR2020a/bin/matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "batch_IFCN_MN(${Batch_id}, 100, 200); quit;"
