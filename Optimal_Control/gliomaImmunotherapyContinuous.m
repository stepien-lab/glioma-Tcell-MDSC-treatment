function phaseout = gliomaImmunotherapyContinuous(input)
% glioblastoma-immune dynamics system and objective functional for optimal control using GPOPS-II 
% used in the code GBMoptimalcontrol.m

% Note: if you decide to use a one-term integral constraint (Intconstraint) instead of a three-term 
% constraint in  OCHiPerGator.m, you'll have to modify the phaseout.integrand at the bottom of this script.


%Parameter values for tumor cells:
lambdaC = input.auxdata.lambdaC;
Cmax = input.auxdata.Cmax; 
eta = input.auxdata.eta; 
%   lambdaC = tumor cell growth rate
%   Ck = carrying capacity of tumor cells
%   eta = kill rate of tumor cells by T cells


%Parameter values for T cells:
aT = input.auxdata.aT;
sT = input.auxdata.sT;
rho = input.auxdata.rho; 
eC = input.auxdata.eC;
r = input.auxdata.r; 
dT = input.auxdata.dT;
%   aT = T cell activation rate
%   sT = rate of tumor cell-mediated proliferation of T cells
%   rho = inhibition rate of T cells by PD-L1-PD-1
%   eC = expression of PD-L1 in tumor cells vs. T cells
%   r = inhibition of T cells by MDSCs
%   dT = T cell death rate


% Parameter values for MDSCs
sM = input.auxdata.sM; 
dM = input.auxdata.dM;
% sM = MDSC recruitment rate by tumor production of chemokines
% dM = MDSC death rate


gammaA = input.auxdata.gammaA;
gammaR = input.auxdata.gammaR;
% gammaA = efficacy of therapeutic dose of anti-PD-1
% gammaR = efficacy of therapeutic dose of CCR2 antagonist


weight1 = input.auxdata.weight1;
weight2 = input.auxdata.weight2;
weight3 = input.auxdata.weight3;


t = input.phase.time;
x = input.phase.state; % cell populations
u = input.phase.control; % immunotherapies


C = x(:,1); %tumor
T = x(:,2); %T cell
M = x(:,3); %MDSCs
A = u(:,1); %anti-PD-1 
R = u(:,2); %CCR2 antagonist
Cdot = lambdaC*C.*(1-(C/Cmax))-eta*T.*C; 
Tdot = ((aT+(sT*T.*C))./(1+(rho*(1-gammaA.*A).*T.*(T+eC*C))))...
    -(r*T.*M)-(dT*T);
Mdot = sM*C.*(1 - gammaR.* R) -dM*M;


phaseout.dynamics = [Cdot, Tdot, Mdot];
phaseout.integrand = [weight1.*C + weight2.*A.^2 + weight3.*R.^2, A, R];%integrand of the objective functional 
