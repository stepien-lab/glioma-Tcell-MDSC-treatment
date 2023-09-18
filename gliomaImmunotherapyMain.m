%-------------- Glioma Combination Immunotherapy Problem -----------------%

clear
clc

% Parameters:
%-------------------------------------------------------------------%
%-------------------- Data Required by Problem ---------------------%
%-------------------------------------------------------------------%
disp('startingrun')
auxdata.lambdaC = 0.431; %0.431; %bifurcation LB: 0.3448, UB: 0.5172 (shift 20% up/down)
auxdata.Cmax = 3.04e6; %3.04e6; % LB: 2.432e6, UB: 3.648e6
auxdata.eta = (2.57e-8); %2.57e-8; %bifurcation LB: 2.06e-8, UB: 3.08e-8
auxdata.aT = 3.26e6; %bifurcation ? not sensitive though, but science paper
auxdata.sT = 8.56e6; %bifurcation
auxdata.rho = 0.107; %0.107; % LB: 0.0856, UB: 0.1284
auxdata.eC = 16.2; %16.2; % LB: 12.96, UB: 19.44
auxdata.r = (6.92e-6); %6.92e-6; % LB: 5.536e-6, UB: 8.304e-6
auxdata.dT = 0.0221;
auxdata.sM = 0.0466; %0.0466; % LB: 0.0373, UB: 0.0559
auxdata.dM = 0.251; %0.251; % LB: 0.201, UB: 0.301
auxdata.gammaA = 1; %2e7;
% I input the parameter set of minimal error according ABC output for average data (except for gamma)
auxdata.gammaR = 1; %% WHAT IS THIS??

auxdata.weight1 = 2/auxdata.Cmax;
auxdata.weight2 = 2;


% HIGHEST INDIVIDUAL DOSE FOR ANTI-PD-1: 1000 ug (ie 1 mg) 
% (double highest dose in Joe's paper 2020)
% HIGHEST OVERALL DOSE FOR ANTI-PD-1: 4800 ug (ie 4.8 mg)
%(double highest overall dose in Adam et al. 2021)

% HIGH INDIVIDUAL DOSE FOR CCR2 ANTAGONIST (needs more work):
% FROM JOE'S PAPER: roughly 2 mg (so highest: 4 mg)
% HIGH OVERALL DOSE FOR CCR2 ANTAGONIST (needs more work):
% FROM JOE'S PAPER: 1.98 * twice daily * 21 days = roughly 84mg (so
% highest: 168 mg)


%both treatments in terms of day^-1 cm^3/g (where the density^-1 is for the
%drug not the cells)
%a = 4.5e-8; %1e-2; % was 5e-8 if this were = 1/gammaA (should be < 1/gammaA)
%is this the bound on an individual treatment?? yes
% I'm pretty sure, that based on the boundary conditions, a is the bound
% that the treatment can be in the body at any time t
%r = 1e8;

multiple = 0;   % if multiple = 0, there is only 1 intconstraint
                % if multiple = 1, there are multiple (namely 3)


a = 0.7;
r = 0.9; %perhaps change back to 1 or 0.9

if multiple == 0
    Intconstraint = 50;
end

if multiple == 1
    Intconstraint = [35, 15, 15];
end

%Intconstraint = 1e15; % 1e-4; % is this the bound on the integral constraint? yes
% is the integral constraint a bound on the integral of the objective
% functional???
%%
%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%

% need to consider boundary conditions %%%%
CMax = 2e8; % humane endpoint for mice (our maximum cell count was 41,395,783)
% Humane endpont:   https://policies.unc.edu/TDClient/2833/Portal/KB/ArticleDet?ID=132214
CMin = 0;
TMax = CMax; 
TMin = 0; %?
MMax = CMax; 
MMin = 0;
yMax = (Intconstraint); % y is the integral constraint on the treatment
%yMin = [0,0,0];
if multiple == 0
    yMin = 0;
end

if multiple == 1
    yMin = [0,0,0];
end
AMax = a; %probably need to change when add in CCR2 antagonist
AMin = 0;
RMax = r;
RMin = 0;
t0Max = 7;% 7; %change to 0 if you want treatment to start immediately
t0Min = 7; % 7; %changed to 7 since Flores-Toro starts treatment at day 7
tfMax = 200; %adjust to allow longer treatment period
tfMin = 10; %adjust to force longer treatment period

% solve for initial conditions (not collected in best way--would be better
% if we could just refer to gliomaImmunotherapyContinuous.m instead of
% GBMFunc.m)

Cgrow = 35000;
Tgrow = 100;
Mgrow = 0; 

p = [auxdata.lambdaC, auxdata.Cmax,auxdata.eta, auxdata.aT, auxdata.sT, auxdata.rho, auxdata.eC, auxdata.r, auxdata.dT, auxdata.sM,auxdata.dM];

tspan = [0:1:t0Min];

addpath('/Users/hannahanderson/Downloads')
[T,X] = ode45(@(t,x) GBMFuncoptimal(t,x,p),tspan,[Cgrow;Tgrow;Mgrow]);

C0 = X(length(tspan),1);
T0 = X(length(tspan),2);
M0 = X(length(tspan),3);


%C0 = 408039; %either grab day 7 time point, or run simulation and grab day 7
%maybe run my simulation and grab a later point? or use data to start at a later time point so T0, M0 not 0?
%T0 = 81221;
%M0 = 56353; % for these initial conditions, I picked a random tumor at day 24
y0 = 0; %[0,0,0];

bounds.phase.initialtime.lower = t0Min;
bounds.phase.initialtime.upper = t0Max;
bounds.phase.finaltime.lower = tfMin;
bounds.phase.finaltime.upper = tfMax;
bounds.phase.initialstate.lower = [C0, T0, M0]; % could leave wiggle room for initial state?
bounds.phase.initialstate.upper = [C0, T0, M0];
bounds.phase.state.lower = [CMin, TMin, MMin];
bounds.phase.state.upper = [CMax, TMax, MMax]; % max should never be above tumor size that would not be survivable to a mouse
bounds.phase.finalstate.lower = [CMin, TMin, MMin];
bounds.phase.finalstate.upper = [CMax/100, TMax, MMax];
bounds.phase.control.lower = [AMin, RMin]; % change this when figure out treatment
bounds.phase.control.upper = [AMax, RMax]; %
% MAYBE GET RID OF INTEGRAL CONSTRAINT BOUNDS? DO I REALLY CARE ABOUT IT?
%bounds.phase.integral.lower = [0,0,0]; %probably will stay the same
%bounds.phase.integral.upper = Intconstraint; %

if multiple == 0
    bounds.phase.integral.lower = (0); %probably will stay the same
    bounds.phase.integral.upper = (Intconstraint); %
end

if multiple == 1
    bounds.phase.integral.lower = [0,0,0]; %probably will stay the same
    bounds.phase.integral.upper = Intconstraint; %
end

bounds.phase.duration.lower = 10; %according to Flores-Toro Fig.2, 100% of mice survive to this point in treatment
bounds.phase.duration.upper = 50; %also tried 43, also test 60; % treatment in Flores-Toro lasts 21 days (starts 7 days after implantation) 
% also reasonable duration UB according to Adam et al. toxicity anti-PD-1
% -- treated mice for 6 weeks (42 days)
% bounds.parameter.lower = 
% bounds.parameter.upper =
guess.phase.time    = [0; 10]; %?
guess.phase.state   = [[C0; CMax],[T0; TMax],[M0; MMax]];
guess.phase.control = [[AMax/4; AMax/4],[RMax/4;RMax/4]]; %
guess.phase.integral = (Intconstraint/2); %
%guess.parameter = vertcat(auxdata); %% DO I EVEN NEED THIS?




%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre'; % mesh refinement method
mesh.tolerance       = 1e-3; % desired accuracy tolerance of the mesh
mesh.maxiterations   = 50; %default is 10
% mesh.colpointsmin
% mesh.colpointsmax % max allowable # of collocation points in a mesh
% interval
% mesh.splitmult % only used with method "hp-DarbyRao"
% mesh.curveratio % only used with method "hp-DarbyRao"
% mesh.R % only used wtih method "hp-LiuRao"
% mesh.sigma % only used with method "hp-LiuRao-Legendre"
mesh.phase.colpoints = 4*ones(1,10); % row vector of integers in interal (1,10)
mesh.phase.fraction  = 0.1*ones(1,10);

% collocation method is a way of approximating the solution of the
% derivative by creations n +1 conditions (for n collocation points) in
% order to solve a polynomial of degree n, which is an approximation of the
% solution (example: trapezoid method is a collocation method using 2
% collocation points). These collocation methods are in fact implicit
% Runge-Kutta methods.
% An orthogonal collocation method uses orthogonal poly basis to approx, such as
% the Legendre polynomials

%-------------------------------------------------------------------%
%--------------------------- Problem Setup -------------------------%
%-------------------------------------------------------------------%
setup.name                        = 'Glioma-Immunotherapy-Problem';
setup.functions.continuous        = @gliomaImmunotherapyContinuous;
setup.functions.endpoint          = @gliomaImmunotherapyEndpoint;
setup.displaylevel                = 2;
setup.auxdata                     = auxdata;
setup.bounds                      = bounds;
setup.guess                       = guess;
setup.mesh                        = mesh;
setup.nlp.solver                  = 'ipopt';
setup.derivatives.supplier        = 'sparseCD'; %type of derivative approximation
setup.derivatives.derivativelevel = 'second'; %could set to first order
setup.scales.method               = 'automatic-hybridUpdate';
setup.method                      = 'RPM-Differentiation';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
% HOW WILL WE USE gpops2 in HiPerGator????

%% Final outputs of interest

solution = output.result.solution;
time = solution.phase(1).time;
state = solution.phase(1).state;
control = solution.phase(1).control;

Cfinal = solution.phase(1).state(end,1);
Csizemax = max(solution.phase(1).state(:,1));
tf = solution.phase.time(end);
cumulativedosageantiPD1 = trapz(solution.phase(1).time, solution.phase(1).control(:,1));
cumulativedosageCCR2antagonist = trapz(solution.phase(1).time, solution.phase(1).control(:,2));


% CHECK IF THESE WORK
Cthreshold = 20000;
timesaboveCthreshold = (solution.phase(1).time).*(solution.phase(1).state(:,1) > Cthreshold);

antiPD1threshold = 0.5;
timesaboveantiPD1threshold = (solution.phase(1).time).*(solution.phase(1).control(:,1) > antiPD1threshold);

CCR2threshold = 0.5;
timesaboveCCR2threshold = (solution.phase(1).time).*(solution.phase(1).control(:,2) > CCR2threshold);

