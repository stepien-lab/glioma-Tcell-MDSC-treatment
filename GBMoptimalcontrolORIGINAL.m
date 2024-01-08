% OPTIMAL CONTROL FUNCTION

function GBMoutput = GBMoptimalcontrolORIGINAL(p, Cthreshold, antiPD1threshold, CCR2threshold, maxdose, multiple, Intconstraint)

% Parameters:
%-------------------------------------------------------------------%
%-------------------- Data Required by Problem ---------------------%
%-------------------------------------------------------------------%
disp('startingrun')
auxdata.lambdaC = p(1);
auxdata.Cmax = p(2);
auxdata.eta = p(3);
auxdata.aT = p(4);
auxdata.sT = p(5);
auxdata.rho = p(6);
auxdata.eC = p(7);
auxdata.r = p(8);
auxdata.dT = p(9);
auxdata.sM = p(10);
auxdata.dM = p(11);
auxdata.gammaA = 1; %2e7;
% I input the parameter set of minimal error according ABC output for average data (except for gamma)
auxdata.gammaR = 1; %% WHAT IS THIS??

auxdata.weight1 = 2/auxdata.Cmax;
auxdata.weight2 = 2;

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
if multiple == 0
    yMin = 0;
end

if multiple == 1
    yMin = [0,0,0];
end
AMax = maxdose(1);
AMin = 0;
RMax = maxdose(2);
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

%-------------------------------------------------------------------%

solution = output.result.solution;

Cfinal = solution.phase(1).state(end,1);
Csizemax = max(solution.phase(1).state(:,1));
tf = solution.phase.time(end);
tmax = max(solution.phase.time(:) .*(solution.phase(1).state(:,1) == Csizemax));
cumulativedosageantiPD1 = trapz(solution.phase(1).time, solution.phase(1).control(:,1));
cumulativedosageCCR2antagonist = trapz(solution.phase(1).time, solution.phase(1).control(:,2));


% Note: times do not include the first week (or specified time) of no
% treatment
timesaboveCthreshold = (solution.phase(1).time).*(solution.phase(1).state(:,1) > Cthreshold);

timesaboveantiPD1threshold = (solution.phase(1).time).*(solution.phase(1).control(:,1) > antiPD1threshold);

timesaboveCCR2threshold = (solution.phase(1).time).*(solution.phase(1).control(:,2) > CCR2threshold);



preoutput = {[Cfinal;tf], [Csizemax; tmax], [cumulativedosageantiPD1; cumulativedosageCCR2antagonist], timesaboveCthreshold, [timesaboveantiPD1threshold, timesaboveCCR2threshold]};

% final tumor size and final time DONE
% max tumor size and time of this DONE
% total cumulative dosage (did trapezoid integral) DONE
% total time tumor is above a certain level DONE 
% total time treatment is above a certain dosage DONE 


GBMoutput = preoutput;
end
