% Code for HiPerGator run

% Save file?
    savefile = 0; % 0 is no, 1 is yes
    filename = 'OCtest.mat';

n = 2; % number of parameter sets to vary
Cthreshold = 20000;
antiPD1threshold = 0.5;
CCR2threshold = 0.5;
maxdose = [0.7,0.9]; %anti-PD-1 max dose, CCR2 antagonist max dose

multiple = 0;   % if multiple = 0, there is only 1 intconstraint
                % if multiple = 1, there are multiple (namely 3)

if multiple == 0
    Intconstraint = 50;
end

if multiple == 1
    Intconstraint = [35, 15, 15];
end

chooseparam = [1;     % 1  - lambdaC
               0;     % 2  - Cmax
               0;     % 3  - η 
               0;     % 4  - a_T 
               0;     % 5  - s_T 
               0;     % 6  - ρ        
               0;     % 7  - ε_C 
               0;     % 8  - r 
               0;     % 9  - d_T
               0;     % 10 - s_M 
               0];    % 11 - d_M
           
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
    
%% Sample parameter sets
   
rng(0, 'twister'); % to uniformly generate random numbers 

Psample = zeros(n,11);

for j = 1:length(param)
    if chooseparam(j) == 1
        Psample(:,j) = (mpr(j,2) - mpr(j,1)).*rand(n,1) + mpr(j,1);
    else
        Psample(:,j) = param(j);
    end
end

%% Evaluate Optimal Control Problem

output = cell(n,6);


for i = 1:n
    p = Psample(i,:);
    output{i,1} = p;
    %output{i,2:7} = GBMoptimalcontrol(p, Cthreshold, antiPD1threshold, CCR2threshold, maxdose, multiple, Intconstraint);
    preoutput = GBMoptimalcontrolORIGINAL(p, Cthreshold, antiPD1threshold, CCR2threshold, maxdose, multiple, Intconstraint);
    for j = 2:6
        output{i,j} = preoutput{j-1};
    end
end

%% Save File

if savefile == 1
    save(filename)
end
