function dxdt = GBMFuncoptimal(t,x,p)

% Anderson et al. (2024) model of treatment-free glioblastoma-immune dynamics 
% to be used in conjunction with the optimal control code.

% The purpose of this code is to determine system dynamics pre- and
% post-treatment.

%--------------------------------------------------------------------------

% INPUTS: 

C = x(1);
T = x(2);
M = x(3);

% x is 3x1 vector with:
% C = number of tumor cells
% T = number of T cells
% M = number of MDSCs 

%--------------------------------------------------------------------------

%Parameter values for tumor growth:

lambdaC = p(1);
Cmax = p(2); 
eta = p(3); 

%   lambdaC = tumor cell growth rate
%   Cmax = carrying capacity of tumor cells
%   eta = kill rate of tumor cells by T cells

%--------------------------------------------------------------------------

%Parameter values for T cells activation/suppression:

aT = p(4);
sT = p(5);
rho = p(6);
eC = p(7);
r = p(8); 
dT = p(9);

%   aT = activation rate of T cells
%   sT = stimulation rate of T cells
%   rho = inhibition rate of T cells by PD-L1-PD-1
%   eC = expression of PD-L1 in tumor cells vs. T cells
%   r = inhibition of T cells by MDSCs
%   dT = death rate of T cells

%--------------------------------------------------------------------------

% Parameter values for MDSCs
sM = p(10);  
dM = p(11);

% sM = stimulation by CCL2/CCL7
% dM = death rate of MDSCs

%--------------------------------------------------------------------------

% ODE FUNCTIONS: 
dC = lambdaC*C*(1-(C/Cmax))-eta*T*C; 
dT = ((aT+(sT*T*C))/(1+(rho*T*(T+eC*C)))) -r*T*M -dT*T;
dM = (sM*C)-dM*M;
    
% COMBINE RESULTS:
dxdt = [dC; dT; dM];
