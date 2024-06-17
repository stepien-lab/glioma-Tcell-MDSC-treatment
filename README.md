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
+ [OCHiPerGator.m](Optimal_Control/OCHiPerGator.m): 
+ [OCHiPerGatorSlurmscript.sh](Optimal_Control/OCHiPerGatorSlurmscript.sh): 
+ [GBMoptimalcontrol.m](Optimal_Control/GBMoptimalcontrol.m): 
+ [GBMFuncoptimal.m](Optimal_Control/GBMFuncoptimal.m): 
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
