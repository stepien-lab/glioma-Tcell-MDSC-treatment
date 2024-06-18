# glioma-Tcell-MDSC-treatment

<a href="https://github.com/stepien-lab/glioma-Tcell-MDSC-treatment/"><img src="https://img.shields.io/badge/GitHub-stepien--lab%2Fglioma--Tcell--MDSC--treatment-blue" /></a> <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" /></a>

The code contained in the glioma-Tcell-MDSC-treatment project was developed to numerically simulate the interactions of glioma cells, T cells, and myeloid-derived suppressor cells (MDSCs) in a model of brain cancer with treatment. It is described in:
>[Hannah G. Anderson](https://github.com/HannahGrace314), Gregory P. Takacs, Jeffrey K. Harrison, Libin Rong, and [Tracy L. Stepien](https://github.com/tstepien/), "Optimal control of combination immunotherapy in a glioblastoma-immune dynamics model"

## Description of Folders
+ [Optimal Control](Optimal_Control): Code to create Figures... in the paper
+ [Parameter Identifiability](Parameter_Identifiability): Code to create Figures... in the paper
+ [foldername](foldername): Code to create Figures... in the paper

## Programs
Optimal Control Folder:
+ [OCHiPerGator.m](Optimal_Control/OCHiPerGator.m): This script samples n parameter sets to represent the tumor microenvironments of n gliomas, and then calls GBMoptimalcontrol.m to determine optimized, personalized treatment regimens for the virtual cohort. Since this code takes some time, we used the UF supercomputer, HiPerGator, for the treatment personalization of 10,000 virtual subjects.
+ [OCHiPerGatorSlurmscript.sh](Optimal_Control/OCHiPerGatorSlurmscript.sh): Example slurm script for submitting OCHiPerGator.m to a supercomputer, in our case, the UF supercomputer HiPerGator. To make this slurmscript your own, be sure to modify mail-user, account, and qos. Likely, you'll also need to increase the time.
+ [GBMoptimalcontrol.m](Optimal_Control/GBMoptimalcontrol.m): Function for finding the optimal treatment regimen using GPOPS-II. Called by: OCHiPerGator.m. Uses: gliomaImmunotherapyContinuous.m, gliomaImmunotherapyEndpoint.m.
+ [GBMFuncoptimal.m](Optimal_Control/GBMFuncoptimal.m): Function containing the model of treatment-free glioblastoma-immune dynamics to be used in conjunction with the optimal control code. The purpose of this code is to determine system dynamics pre- and post-treatment.
+ [gliomaImmunotherapyContinuous.m](Optimal_Control/gliomaImmunotherapyContinuous.m): 
+ [gliomaImmunotherapyEndpoint.m](Optimal_Control/gliomaImmunotherapyEndpoint.m): 
+ [OCoptimizedexample.m](Optimal_Control/OCoptimizedexample.m): 
+ [virtualcohortFigures.m](Optimal_Control/virtualcohortFigures.m): 

Parameter Identifiability Folder:
+ [GBM_identifability_main.m](Parameter_Identifiability/GBM_identifiability_main.m): 
+ [GBMFuncidentifiable.m](Parameter_Identifiability/GBMFuncidentifiable.m):
+ [GBMCost.m](Parameter_Identifiability/GBMCost.m): 
+ [GBMMiniFisher.m](Parameter_Identifiability/GBMMiniFisher.m): 
+ [GBMProfLike.m](Parameter_Identifiability/GBMProfLike.m): 
+ [Differential_Algebra_Approach_for_Global_Structural_Identifiability_for_GBM_model.nb](Parameter_Identifiability/Differential_Algebra_Approach_for_Global_Structural_Identifiability_for_GBM_model.nb): 

+ [filename.m](filename.m): run this program to solve the model equations

## Lead Developer
The lead developer of this code is [Hannah G. Anderson](https://github.com/HannahGrace314).

## Licensing
Copyright 2022 [Tracy Stepien's Research Lab Group](https://github.com/stepien-lab/). This is free software made available under the MIT License. For details see the [LICENSE](LICENSE) file.
