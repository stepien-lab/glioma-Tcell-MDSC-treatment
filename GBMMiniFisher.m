% Simplified FIM for GBM
% Original Source: Marisa Eisenberg (marisae@umich.edu)

function FIM = GBMMiniFisher(tspan,fixed,params,indexchooseparam,initialcondition,variablenum)

delta = 0.001; %percent change we'll use, default 0.1%
X = [];


for j=1:length(params)
    params1 = params;
    params2 = params;
    params1(j) = (1+delta)*params(j);
    params2(j) = (1-delta)*params(j);
    
    [t, x1] = ode45(@GBMFuncidentifiable,tspan,initialcondition,[],fixed,params1,indexchooseparam);
    [t, x2] = ode45(@GBMFuncidentifiable,tspan,initialcondition,[],fixed,params2,indexchooseparam);
    
    X = [X, (x1(:,variablenum) - x2(:,variablenum))./(2*delta*params(j))]; % approximation of the derivative of the output wrt parameter
    %this fills in the jth column of the design matrix with the sensitivities to parameter j 
    %at each time point.
end

%In case the initial conditions are fixed
X = X(2:end,:);

%FIM (simplified w/o weighting term)
FIM = X'*X;


