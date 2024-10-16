# glioma-Tcell-MDSC-treatment

<a href="https://github.com/stepien-lab/glioma-Tcell-MDSC-treatment/"><img src="https://img.shields.io/badge/GitHub-stepien--lab%2Fglioma--Tcell--MDSC--treatment-blue" /></a> <a href="https://doi.org/10.1016/j.jtbi.2024.111951"><img src="https://img.shields.io/badge/doi-10.1016%2Fj.jtbi.2024.111951-orange"></a> <a href="https://doi.org/10.1101/2024.04.29.591725"><img src="https://img.shields.io/badge/bioRxiv-2024.04.29.591725-orange" /></a> <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" /></a>

The code contained in the glioma-Tcell-MDSC-treatment project was developed to numerically simulate the interactions of glioma cells, T cells, and myeloid-derived suppressor cells (MDSCs) in a model of brain cancer with treatment. It is described in:
>[Hannah G. Anderson](https://github.com/HannahGrace314), [Gregory P. Takacs](https://www.linkedin.com/in/gregoryptakacs/), [Jeffrey K. Harrison](https://pharmacology.med.ufl.edu/profile/harrison-jeffrey/), [Libin Rong](https://people.clas.ufl.edu/libinrong/), and [Tracy L. Stepien](https://github.com/tstepien/), Optimal control of combination immunotherapy for a virtual murine cohort in a glioblastoma-immune dynamics model, _Journal of Theoretical Biology_, 595, 111951 (2024) doi: [10.1016/j.jtbi.2024.111951](https://doi.org/10.1016/j.jtbi.2024.111951)

## Description of Folders
+ [Optimal Control](Optimal_Control): MATLAB code to optimize treatment regimens for the combination immunotherapy in glioblastoma. Contains code to create Figures 2-6. 
+ [Parameter Identifiability](Parameter_Identifiability): MATLAB code to perform practical identifiability analysis on the treatment-free version of the model. Also contains Mathematica file with information regarding structural identifiability analysis. Contains code to create Figure D7.

## Programs
[Optimal Control](Optimal_Control) Folder:
+ [OCHiPerGator.m](Optimal_Control/OCHiPerGator.m): This script samples n parameter sets to represent the tumor microenvironments of n gliomas, and then calls GBMoptimalcontrol.m to determine optimized, personalized treatment regimens for the virtual cohort. Since this code takes some time, we used the UF supercomputer, HiPerGator, for the treatment personalization of 10,000 virtual subjects.
+ [OCHiPerGatorSlurmscript.sh](Optimal_Control/OCHiPerGatorSlurmscript.sh): Example slurm script for submitting OCHiPerGator.m to a supercomputer, in our case, HiPerGator. To make this slurm script your own, be sure to modify mail-user, account, and qos. Likely, you'll also need to increase the time depending on the sample size.
+ [GBMoptimalcontrol.m](Optimal_Control/GBMoptimalcontrol.m): Function for finding the optimal treatment regimen using GPOPS-II. Called by: OCHiPerGator.m. Uses: GBMFuncoptimal.m, gliomaImmunotherapyContinuous.m, gliomaImmunotherapyEndpoint.m.
+ [GBMFuncoptimal.m](Optimal_Control/GBMFuncoptimal.m): Function containing the model of treatment-free glioblastoma-immune dynamics to be used in conjunction with the optimal control code. The purpose of this code is to determine system dynamics pre- and post-treatment.
+ [gliomaImmunotherapyContinuous.m](Optimal_Control/gliomaImmunotherapyContinuous.m): Glioblastoma-immune dynamics system (with treatment) and objective functional for optimal control using GPOPS-II. 
+ [gliomaImmunotherapyEndpoint.m](Optimal_Control/gliomaImmunotherapyEndpoint.m): Endpoint function needed for optimal control using GPOPS-II.
+ [OCoptimizedexample.m](Optimal_Control/OCoptimizedexample.m): This script outputs Figure 2 representing the optimal regimen for an average mouse with glioblastoma. Uses: GBMFuncoptimal.m (and output from OCHiPerGator.m + GBMoptimalcontrol.m).
+ [virtualcohortFigures.m](Optimal_Control/virtualcohortFigures.m): This script analyzes the virtual murine cohort from the OCHiPerGator.m output and generates Figures 3-6.

[Parameter Identifiability](Parameter_Identifiability) Folder:
+ [GBM_identifability_main.m](Parameter_Identifiability/GBM_identifiability_main.m): Main code for practical identifiability, edited from Marisa Eisenberg's original code (marisae@umich.edu). This code loads treatment-free glioblastoma data from Anderson et al. (2023), finds the parameter set of best fit using GBMCost.m along with the gradient descent method, generates the Fisher information matrix using GBMMiniFisher.m, and determines the profile likelihoods using GBMProfLike.m. The code also uses GBMFuncidentifiable.m, which contains the treatment-free glioblastoma (GBM) model. Use this script to generate Figure D7.
+ [GBMFuncidentifiable.m](Parameter_Identifiability/GBMFuncidentifiable.m): Function containing GBM model to be used for practical identifiability analysis (GBM_identifiability_main.m). The code allows certain parameters to be fixed in GBM_identifiability_main.m while others can be varied in order to be fitted to data.
+ [GBMCost.m](Parameter_Identifiability/GBMCost.m): Cost function for the GBM model, edited from Marisa Eisenberg's original code. We used the weighted ordinary least squares (OLS) error, but this script contains multiple other error functions that can be used: Poisson ML, root mean square error (RMSE), mean square error (MSE), and relative error.
+ [GBMMiniFisher.m](Parameter_Identifiability/GBMMiniFisher.m): Function to calculate the simplified Fisher information matrix. Original code from Marisa Eisenberg. 
+ [GBMProfLike.m](Parameter_Identifiability/GBMProfLike.m): Profile Likelihood Generator. Original code from Marisa Eisenberg.
+ [Differential_Algebra_Approach_for_Global_Structural_Identifiability_for_GBM_model.nb](Parameter_Identifiability/Differential_Algebra_Approach_for_Global_Structural_Identifiability_for_GBM_model.nb): Wolfram Mathematica document calculating the input-output relation / polynomial for the GBM model. The input-output relation is used to determine the global structural identifiability of the model using the differential algebra approach.

## Lead Developer
The lead developer of this code is [Hannah G. Anderson](https://github.com/HannahGrace314).

## Licensing
Copyright 2022-2024 [Tracy Stepien's Research Lab Group](https://github.com/stepien-lab/). This is free software made available under the MIT License. For details see the [LICENSE](LICENSE) file.
