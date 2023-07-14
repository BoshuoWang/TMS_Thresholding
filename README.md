# TMS_Thresholding

Repository of MATLAB code related to manuscript: B.Wang, A.V.Peterchev, and S.M.Goetz, Analysis and Comparison of Methods for Determining Transcranial Magnetic Stimulation Thresholds, 2022.
bioRxiv preprint DOI https://doi.org/10.1101/2022.06.26.495134

Functions:
batch_Bayesian.m & compile_Bayesian.m: submit jobs running Bayesian methods on cluster and compile results from different batches.
batch_IFCN_MN.m & compile_IFCN_MN.m: submit jobs running IFCN relative frequency and Mills-Nithi methods on cluster and compile results from different batches. 
batch_MLE.m, compile_MLE_ver.m, & compile_Bayesian.m: submit jobs running MLE/MAP/Bayesian methods on cluster and compile results from different batches.
batch_RM.m, compile_RM_all.m, compile_RM_nonparam.m, & compile_RM_param.m: submit jobs running stochastic root-finding methods (Robbins-Monro) on cluster and compile results from different batches.

Folder "Statistical-MEP-Model" is updated version of code from repository: https://github.com/sgoetzduke/Statistical-MEP-Model
	For details, see readme within folder.

Folder "Functions" contains TMS thresholding methods and supporting functions and data.
	BayesianAdaptiveEstimation.m: Bayesian methods.
	fsolve_diffser.m: used for finding point on noiseless response curve.
	generateLUT.m: generate lookup table for BayesianAdaptiveEstimation.
	IFCNmethod.m: relative frequency methods.
	loglikelyhood.m: calculate log-likelihood.
	logMAP.m: calculate log-a-posterior-likelihood.
	LUT.mat: lookup table for Bayesian methods
	MillsNithiMethod.m: Mills-Nithi two-threshold method.
	MLEstimation.m: maximum likelihood method and variations (including MAP).
	p_ab.m: fitted prior distribution of parameters for Gaussian model used in MAP and Bayesian methods.
	pMAP.m: a-posterior likelihood
	pMLE.m: likelihood.
	StochasticApproximation.m stochastic root-finding methods (all variants).
	
Folder "DataStorage" contains .mat file (parameters of 25000 virtual subjects) and stores output data of TMS thresholding methods.
	Subjects_25000.mat: subject parameters generated for v36 using GenerateSubjects.m. Threshold defined as amplitude generating 50 uV MEP using noiseless response model. Starting amplitude defined as amplitude generating 500 uV MEP using noiseless response model.
	Subjects_25000_reprocessed.mat: same subject parameters, with threshold and starting amplitude reprocessed for v37 using reGenerateSubjects.m, compile_reSubjects.m,and reprocess_startamp.m. Threshold defined as amplitude generating suprathreshold MEP (greater or equal than 50 uV) with 50% probablity. Starting amplitude defined as amplitude generating suprathreshold MEP with at least 95% probablity (90% probablity version available). Contains additional information.


Example slurm files are provided to run TMS thresholding methods on a cluster. Modify paths and parameters as needed.

For running individual methods, load Subjects_25000_reprocessed.mat, select one element of the Subjects variable, and setup params variable (see e.g., batch_MLE.m for details). Run function for specific method with arguments (params, version). Details on versions can be found within the function for a specific method. 
-------------------------------------------------------------------------------------------------------------------------
The copyrights of this software are owned by Duke University. As such, two licenses for this software are offered:

1. An open-source license under the GPLv2 license for non-commercial use.
2. A custom license with Duke University, for commercial use without the GPLv2 license restrictions. 
 
As a recipient of this software, you may choose which license to receive the code under. Outside contributions to the Duke-owned code base cannot be accepted unless the contributor transfers the copyright to those changes over to Duke University.

To enter a custom license agreement without the GPLv2 license restrictions, please contact the Digital Innovations department at the Duke Office for Translation & Commercialization (OTC) (https://olv.duke.edu/software/) at otcquestions@duke.edu with reference to “OTC File No. 8073” in your email.

Please note that this software is distributed AS IS, WITHOUT ANY WARRANTY; and without the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.