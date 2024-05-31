function output = GliomaImmunotherapyEndpoint(input)
% For optimal control using GPOPS-II 
% used in the code GBMoptimalcontrol.m

Cf = input.phase.finalstate(1);
output.objective = Cf; 
