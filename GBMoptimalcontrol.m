function GBMoutput = GBMoptimalcontrol(p,initialcondition,bounds,guess,mesh)
% function to find the optimal treatment regimen using GPOPS-II
% Called by: OCHiPerGator.m
% Uses: gliomaImmunotherapyContinuous.m, gliomaImmunotherapyEndpoint.m


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
auxdata.gammaA = 1; % gammaA and gammaR determined from Tu et al. 2020 data
auxdata.gammaR = 1; 

auxdata.weight1 = 2/auxdata.Cmax;
auxdata.weight2 = 2;
auxdata.weight3 = 1;

%% Get initial value at treatment start time

tspan = [0:1:bounds.phase.initialtime.lower];

[T,X] = ode45(@(t,x) GBMFuncoptimal(t,x,p),tspan,initialcondition);

C0 = X(end,1);
T0 = X(end,2);
M0 = X(end,3);

bounds.phase.finalstate.upper = [C0, bounds.phase.state.upper(2), bounds.phase.state.upper(3)]; 
bounds.phase.initialstate.lower = [C0, T0, M0]; 
bounds.phase.initialstate.upper = [C0, T0, M0];
guess.phase.state   = [[C0; bounds.phase.state.upper(1)],[T0; bounds.phase.state.upper(2)],[M0; bounds.phase.state.upper(3)]];

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

%-------------------------------------------------------------------%

solution = output.result.solution;

state = solution.phase(1).state;
time = solution.phase(1).time;
control = solution.phase(1).control;


GBMoutput = {state, control, time, T, X};

end
