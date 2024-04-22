function dxdt = GBMFuncidentifiable(t,x,fixed,params,indexchooseparam) 
    % This is our model of treatment free glioblastoma-immune dynamics to
    % be used for practical identifiability analysis (GBM_identifiability_main.m). 
    % The code allows certain parameters to be fixed while others are 
    % varied in order to be fitted to data.
    
    
    % EXTRACT VARIABLES:
    C = x(1); % Cancer cells
    T = x(2); % T-Cells
    M = x(3); % MDSCs
    
    
    % EXTRACT PARAMETERS: 
    mp = fixed;
    mp(indexchooseparam(:)) = params(:);  
    
    lambdaC = mp(1); 
    Cmax = mp(2);
    eta = mp(3);
    aT = mp(4);
    sT = mp(5);
    rho = mp(6);
    eC = mp(7);
    r = mp(8);
    dT = mp(9);
    sM = mp(10);
    dM = mp(11);
    
    % ODE FUNCTIONS: 
    dC = lambdaC*C*(1-(C/Cmax))-eta*T*C; 
    dT = ((aT+(sT*T*C))/(1+(rho*T*(T+eC*C)))) -r*T*M -dT*T;
    dM = (sM*C)-dM*M;
    
    % COMBINE RESULTS:
    dxdt = [dC; dT; dM];
    
end