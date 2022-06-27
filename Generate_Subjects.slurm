#!/bin/bash
#
#SBATCH --mem=10G
#SBATCH --job-name=TMS_IO_Generate_Subjs
#SBATCH --output=subjs.out
#SBATCH --error=subjs.err
#SBATCH --mail-type=END,FAIL
#SBATCH -c 51

uname -n 1>&1

# Adjust path as needed
#cd /TMS_Thresholding

# Call MATLAB function, modify path for MATLAB as needed
/usr/bin/time -v /admin/apps/rhel8/matlabR2021a/bin/matlab -nodisplay -nodesktop -nosplash -singleCompThread -r "GenerateSubjects; quit;"
