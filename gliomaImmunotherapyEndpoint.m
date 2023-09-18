function output = GliomaImmunotherapyEndpoint(input)

Cf = input.phase.finalstate(1);
output.objective = Cf; 
%idk if I should change more, just changed pf to Cf
%I know it is the objective functional that we're seeking to minimize