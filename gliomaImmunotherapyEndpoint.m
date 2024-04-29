function output = GliomaImmunotherapyEndpoint(input)

Cf = input.phase.finalstate(1);
output.objective = Cf; 
