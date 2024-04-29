% OCHiPerGator = Optimal Control using HiPerGator

% This script samples n parameter sets to represent the tumor
% microenvironments of n gliomas, and then calls GBMoptimalcontrol.m to
% determine optimized, personalized treatment regimens for a combination
% immunotherapy. Since this code takes some time, we used the UF
% supercomputer, HiPerGator, for the treatment personalization of 10,000 
% virtual subjects.

% Save file?
    savefile = 0; % 0 is no, 1 is yes
    filename = 'OCtest.mat';

n = 5; % number of parameter sets to vary
maxdose = [0.8,0.9]; %anti-PD-1 and CCR2 antagonist max percent reduction

multiple = 1;   % if multiple = 0, there is only 1 intconstraint
                % if multiple = 1, there are multiple (namely 3)

if multiple == 0
    Intconstraint = 50;
    % if using a single integral constraint, need to change 
    % phaseout.integrand to [weight1.*C + weight2.*A.^2 + weight3.*R.^2]
    % in gliomaImmunotherapyContinuous.m
end

if multiple == 1
    Intconstraint = [200, 15, 30];
end

% choose parameters to vary in sampling:
chooseparam = [1;     % 1  - lambdaC
               0;     % 2  - Cmax
               1;     % 3  - η 
               0;     % 4  - a_T 
               0;     % 5  - s_T 
               1;     % 6  - ρ        
               0;     % 7  - ε_C 
               1;     % 8  - r 
               0;     % 9  - d_T
               0;     % 10 - s_M 
               1];    % 11 - d_M
           
param = [0.431;      % 1  - lambdaC
         3.04e6;     % 2  - Cmax
         2.57e-8;    % 3  - η 
         3.26e6;     % 4  - a_T 
         8.56e6;     % 5  - s_T 
         0.107;      % 6  - ρ        
         16.2;       % 7  - ε_C 
         6.92e-6;    % 8  - r 
         0.0221;     % 9  - d_T
         0.0466;     % 10 - s_M  
         0.251];     % 11 - d_M
                 
 mpr = [0            0.5;        % 1  - lambdaC
        1e6          5e7;        % 2  - Cmax
        0            1e-6;       % 3  - η 
        5e1          5e6;        % 4  - a_T 
        1e2          1e7;        % 5  - s_T 
        0            0.5;        % 6  - ρ     
        1            100;        % 7  - ε_C 
        0            1e-4;       % 8  - r 
        0            0.75;       % 9  - d_T
        0            0.1;        % 10 - s_M 
        0            0.5];       % 11 - d_M  


icp = find(chooseparam == 1); %index of chosen parameters 


%% Sample parameter sets 
Psample = (param*ones(1,n))';

% OPTION 1: SAMPLE PARAMETERS USING LATIN HYPERCUBE

rng default
LHSMatrix = lhsdesign(n,length(icp));

for j = 1:length(icp)
     Psample(:,icp(j)) = (mpr(icp(j),2) - mpr(icp(j),1)).*LHSMatrix(:,j) + mpr(icp(j),1);
end


% OPTION 2: SAMPLE PRACTICALLY IDENTIFIABLE PARAMETERS ACCORDING TO 
% PROBABILITY DISTRIBUTIONS FROM ANDERSON ET AL (2023)

% Psample(:,1) = random('Gamma',5.40,0.050,n,1);      %lambdaC
% Psample(:,3) = random('Weibull',1.84e-7,0.999,n,1); %eta
% Psample(:,6) = random('Exponential',0.223,n,1);     %rho
% Psample(:,8) = random('Gamma',0.801,3.56e-5,n,1);   %r
% Psample(:,11) = random('Normal',0.258,0.143,n,1);   %dM


%% Setup for Optimal Control
Cgrow = 35000;
Tgrow = 100;
Mgrow = 0; 

initialcondition = [Cgrow;Tgrow;Mgrow];

%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%

CMax = 2e8; % humane endpoint for mice (our maximum cell count was 41,395,783)
% Humane endpont:   https://policies.unc.edu/TDClient/2833/Portal/KB/ArticleDet?ID=132214
CMin = 0;
TMax = CMax; 
TMin = 0; 
MMax = CMax; 
MMin = 0;
AMax = maxdose(1);
AMin = 0;
RMax = maxdose(2);
RMin = 0;
t0Max = 7;   %change to 0 if you want treatment to start immediately
t0Min = 7;   %changed to 7 since Flores-Toro (2020) starts treatment at day 7
tfMax = 200; %adjust to allow longer treatment period
tfMin = 10;  %adjust to force longer treatment period

bounds.phase.initialtime.lower = t0Min;
bounds.phase.initialtime.upper = t0Max;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.state.lower = [CMin, TMin, MMin];
bounds.phase.state.upper = [CMax, TMax, MMax]; 
bounds.phase.finalstate.lower = [CMin, TMin, MMin];
bounds.phase.control.lower = [AMin, RMin];
bounds.phase.control.upper = [AMax, RMax]; 


if multiple == 0
    bounds.phase.integral.lower = (0); 
    bounds.phase.integral.upper = (Intconstraint); 
end

if multiple == 1
    bounds.phase.integral.lower = [0,0,0]; 
    bounds.phase.integral.upper = Intconstraint; 
end

bounds.phase.duration.lower = 10; %according to Flores-Toro et al. (2020) Fig.2, 100% of mice survive to this point in treatment
bounds.phase.duration.upper = 43; % Adam et al. (2021) experiment on toxicity of anti-PD-1 -- treated mice for 6 weeks (42 days)

guess.phase.time    = [0; 10]; 
guess.phase.control = [[AMax/4; AMax/4],[RMax/4;RMax/4]]; 
guess.phase.integral = (Intconstraint/2); 


%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre'; % mesh refinement method
mesh.tolerance       = 1e-3; % desired accuracy tolerance of the mesh
mesh.maxiterations   = 10; %default is 10
mesh.phase.colpoints = 4*ones(1,10); % row vector of integers in interal (1,10)
mesh.phase.fraction  = 0.1*ones(1,10);


%% Evaluate Optimal Control Problem

output = cell(n,3);
preoutput = cell(n,1);


parfor i = 1:n
    p = Psample(i,:);
    output{i,1} = p;
    preoutput{i} = GBMoptimalcontrol(p,initialcondition,bounds,guess,mesh);
end


parfor i = 1:n
    output{i,2} = preoutput{i};
end

output{1,3} = {maxdose, multiple, Intconstraint,initialcondition};


%% Save File

if savefile == 1
    save(filename, "output")
end
