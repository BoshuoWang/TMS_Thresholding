#!/bin/bash
#
#SBATCH --mem=20G
#SBATCH --job-name=RM1
#SBATCH --array=1-250
#SBATCH --output=/RM1/%a.out
#SBATCH --error=/RM1/%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH -c 1

uname -n 1>&1

# Adjust path as needed
#cd /TMS_Thresholding

Batch_id=$(($SLURM_ARRAY_TASK_ID))

# Call MATLAB function, modify path for MATLAB as needed
/usr/bin/time -v /admin/apps/rhel8/matlabR2021a/bin/matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "batch_RM(${Batch_id}, 100, 200, 1); quit;"	
