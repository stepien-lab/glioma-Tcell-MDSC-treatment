function phaseout = gliomaImmunotherapyContinuous(input)

%disp('starting cts')

%Parameter values for tumor growth:
lambdaC = input.auxdata.lambdaC;
Cmax = input.auxdata.Cmax; 
eta = input.auxdata.eta; 

%   lambdaC = tumor cell growth rate
%   Ck = carrying capacity of tumor cells
%   eta = kill rate of tumor cells by T cells


%Parameter values for T cells activation/suppression:
aT = input.auxdata.aT;
sT = input.auxdata.sT;
rho = input.auxdata.rho; %may need to change value of rho since now including gamma/A in T eqn?
eC = input.auxdata.eC;
r = input.auxdata.r; 
dT = input.auxdata.dT;
%sC = input.auxdata.sC;

%   aT = activation rate of T cells
%   sT = stimulation rate of T cells
%   rho = inhibition rate of T cells by PD-L1-PD-1
%   eC = expression of PD-L1 in tumor cells vs. T cells
%   r = inhibition of T cells by MDSCs
%   dT = death rate of T cells
%   sC = steepness coefficient of the tumor recruitment of T cells curve


% Parameter values for MDSCs
sM = input.auxdata.sM; 
%alpha = input.auxdata.alpha;
%q = input.auxdata.q;
%Mhat = input.auxdata.Mhat;
dM = input.auxdata.dM;


% sM = stimulation by CCL2/CCL7
% alpha = MDSC expansion coefficient
% q = steepness coefficient of the MDSC production curve
% Mhat = rate of change of MDSC outside brain in absence of tumor
% dM = death rate of MDSCs

%mu = input.auxdata.mu; % values range from 1.566 -175.743 cm^3/g
gammaA = input.auxdata.gammaA; % values range from 1.178e7 -3.135e7 cm^3/g
gammaR = input.auxdata.gammaR;

% mu = blocking rate of PD-1 by anti-PD-1
% gamma = source of anti-PD-1 or CCR2 antagonist

weight1 = input.auxdata.weight1;
weight2 = input.auxdata.weight2;

t = input.phase.time;
x = input.phase.state;
u = input.phase.control; %anti-PD-1 and CCR2 Antagonist (how to input these??)
% ODE for anti-PD-1 and CCR2 Antagonist not in yet

C = x(:,1);
T = x(:,2);
M = x(:,3);
A = u(:,1);
R = u(:,2);
Cdot = lambdaC*C.*(1-(C/Cmax))-eta*T.*C; %Logistic
Tdot = ((aT+(sT*T.*C))./(1+(rho*(1-gammaA.*A).*T.*(T+eC*C))))...
    -(r*T.*M)-(dT*T);
%Tdot = ((aT+(sT*T.*C))./(1+(rho.*T.*(T+eC*C))))...
%    -(r*T.*M)-(dT*T);
%Mdot = (sM*C)-dM*M;
%Mdot = ((sM*C)./(1+ gammaR.*R))-dM*M;
Mdot = sM*C.*(1 - gammaR.* R) -dM*M;

phaseout.dynamics = [Cdot, Tdot, Mdot];
phaseout.integrand = [weight1.*C + weight2.*A.^2 + R.^2]; %, A, R]; %, A, R]; %, A, R]; %OBJECTIVE FUNCTIONAL GOES HERE(just innards, no integral sign) %%%%
% check to see if they have way to include whole objective functional
% choose different values for weights
% have weights depend on Cmax to balance terms

%disp('end cts')