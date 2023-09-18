function dxdt = GBMFuncoptimal(t,x,p)
% function dxdt = GBMFuncoptimal(t,x,p)

% This is our model of treatment free tumor-immune dynamics 
% without the MDSC expansion term, to fit the model used for optimal
% control.
% INPUTS: 
C = x(1);
T = x(2);
M = x(3);
dxdt = zeros(3,1);

% x is 3x1 vector with:
% C = number of tumor cells
% T = number of T cells
% M = number of MDSCs 


%Parameter values for tumor growth:
lambdaC = p(1);
Cmax = p(2); 
eta = p(3); 

%   lambdaC = tumor cell growth rate
%   Cmax = carrying capacity of tumor cells
%   eta = kill rate of tumor cells by T cells


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


% Parameter values for MDSCs
sM = p(10);  
dM = p(11);

% sM = stimulation by CCL2/CCL7
% dM = death rate of MDSCs


%  output:
%   dxdt(1) = derivative of the tumor cell equation
%   dxdt(2) = derivative of the T cell equation
%   dxdt(3) = derivative of the MDSC equation


dxdt(1) = lambdaC*C*(1-(C/Cmax))-eta*T*C; 
%dxdt(2) = (aT+(sT*T*C))-dT*T;
dxdt(2) = ((aT+(sT*T*C))/(1+(rho*T*(T+eC*C))))...
  -r*T*M-dT*T;
dxdt(3)= (sM*C)-dM*M;
