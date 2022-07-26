This repository contains codes used for statistical analysis utilizing Statistics and Machine Learning Toolbox
<br>

* CYJ_linearModel_brain.m: <br>
performing voxel-wise GLM/LMM in volume data, parpool is used to parallelize jobs.

* plot_lm_adjusted_version2.m:<br>
adjusted response plot for linear regression model fitted by fitlm.

* plot_lmeLevel1_adjustedCV.m:<br>
adjusted response plot for level1 variables of linear mixed model fitted by fitlme. See Figure 2 in https://doi.org/10.1101/2022.07.07.499164

* CYJ_lvOneOut_valida_corr.m:<br>
performing validation of correlation using leave-one-out cross validation.

* CYJ_CCA_adjustedCV_NORM.m:<br>
performing CCA between two sets of variables with covariates being controled and p value detemined via permutation tests.
