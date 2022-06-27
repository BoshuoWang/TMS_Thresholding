#!/bin/bash
#
#SBATCH -p wmglab
#SBATCH --mem=10G
#SBATCH --job-name=MLE1
#SBATCH --array=1-250
#SBATCH --output=/work/wmglab/bw152/TMS_Thresholding_v36/out_err/MLE1/%a.out
#SBATCH --error=/work/wmglab/bw152/TMS_Thresholding_v36/out_err/MLE1/%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH -c 1

uname -n 1>&1

# Adjust path as needed
#cd /TMS_Thresholding

Batch_id=$(($SLURM_ARRAY_TASK_ID))

# Call MATLAB function, modify path for MATLAB as needed
/usr/bin/time -v /admin/apps/rhel7/matlabR2020a/bin/matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "batch_MLE(${Batch_id}, 100, 200, 1); quit;"
