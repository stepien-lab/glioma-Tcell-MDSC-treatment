% 
clc; close all; clear;

%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
% TIME PARAMETERS:
t0 = 0;       % Initial Time
tb = 20000;   % Final Time for Bifurcation


% INITIAL CONDITIONS:
u0 = [35000;   % (C) Tumor Cells 
      0;  % (T) T-Cells  
      0]; % (M) MDSCs
  
  
% BIFURCATION OPTIONS:
index = 3;  % Index Number for Bifurcation Parameter % TRY 6 AGAIN
N = 100;    % Resolution
tmt = 0.95; % Percentage of Solution to Throw Out
  

% VARIABLE AND PARAMETER NAMES:
vn = ["C"; "T"; "M"]; 
pn = ["\lambda_C";
      "C_{max}";
      "\eta";
      "a_T";
      "s_T";
      "\rho";
      "\epsilon_C";
      "r";
      "d_T";
      "s_M";
      "\alpha";
      "q";
      "d_M"];
  

% MODEL PARAMETERS:
mp = [0.431;         %  1 - λ_C   - Maximum Growth Rate of Tumor 
      3.04e6;        %  2 - C_max - Carrying Capacity of Tumor Cells 
      2.57e-8;       %  3 - η     - kill rate of tumor cells by T cells                      
      3.26e6;        %  4 - a_T   - activation rate by IL-12
      8.56e6;        %  5 - s_T   - stimulation rate by IL-2
      0.107;         %  6 - ρ     - inhibition of T cells by PD-1-PD-L1 complex
      16.2;          %  7 - ε_C   - expression of PD-L1 in tumor cells vs T cells
      6.92e-6;       %  8 - r     - suppression rate of T cells by MDSCs
      0.0221;        %  9 - d_T   - death rate of T cells
      0.0372;        % 10 - s_M   - stimulation by chemokines CCL2 and/or CCL7
      3.60e8;        % 11 - α     - MDSC expansion coefficient
      3.51e10;       % 12 - q     - steepness coefficient of the MDSCs production curve
      0.251];        % 13 - d_M   - death rate of MDSCs
  
  
% MODEL PARAMETER RANGES: (Min/Max)
 mpr = [0            0.5;        % 1  - lambdaC
        1e6          5e7;        % 2  - Cmax
        0            1e-6;       % 3  - η 
        5e1          5e6;        % 4  - a_T 
        1e2          1e7;        % 5  - s_T 
        1e-16        0.5;        % 6  - ρ     
        1            100;        % 7  - ε_C 
        0            1e-4;       % 8  - r 
        0            0.75;       % 9  - d_T
        0            0.1;        % 10 - s_M 
        1e7          5e8;        % 11 - α
        1e9          1e11;       % 12 - q
        0            0.5];       % 13 - d_M


% INITIALIZE SIZE PARAMETERS:
ne = length(u0); % Number of Equations
np = length(mp); % Number of Parameters


%% BIFURCATION DIAGRAM

vs = zeros(1,N);
     
% CREATE BIFURACTION PARAMETER VECTOR:
vs(:) = linspace(mpr(index,1),mpr(index,2),N);     
    
    
% SOLVE ODE WITH EACH PARAMETER VALUE:
tFm = tmt*tb;
U = zeros(N,ne,2);
for i=1:N    
        
    % PRINT STATUS:
    clc; fprintf('Iteration: %3.0f / %3.0f\n',i,N);
        
    % SOLVE ODE WITH PARAMETER VALUE i:
    mp(index) = vs(i);
    [t,u] = SolveODE(@f,u0,mp,t0,tb);
        
    I = find(t>=tFm);
    U(i,:,:) = [min(u(I,1)) max(u(I,1))
                min(u(I,2)) max(u(I,2))
                min(u(I,3)) max(u(I,3))];
end
      

% PLOT RESULTS:
figure;
    for i = 1:ne       
        subplot(3,1,i);
        plot(vs(:),U(:,i,1),vs(:),U(:,i,2),'--','linewidth',2);
        xlim(mpr(index,:)); 
        legend('Min','Max','location','northeast');
        if i == ne
            xlabel(pn(index),'fontsize',20);
        end
        ylabel(vn(i),'fontsize',20);
    end
    
 
%% ODE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dudt = f(~,u,mp)  


    % EXTRACT VARIABLES:
    C = u(1); % Cancer Cells
    T = u(2); % T-Cells
    M = u(3); % MDSCs

    
    % EXTRACT PARAMETERS:
    rc  = mp(1);  % λ_C
    kc  = mp(2);  % C_max
    dc  = mp(3);  % η
    at  = mp(4);  % a_T
    st  = mp(5);  % s_T
    rho  = mp(6); % ρ
    ec  = mp(7); % ε_C
    r   = mp(8);  % r
    dt  = mp(9);  % d_T
    sm  = mp(10); % s_M
    a   = mp(11); % α
    q   = mp(12); % q
    dm  = mp(13); % d_M   
    
    
    % PRELIMINARY FUNCTIONS:
    Q = T*(T+ec*C);
    F1 = 1/(1+rho*Q);
    F2 = C;
    
    
    % ODE FUNCTIONS:
    dC = C*(rc*(1-(C/kc)) - dc*T);           % (C) Cancer Cells
    dT = (at+st*T*C)*F1 - r*T*M - dt*T;      % (T) T Cells
    dM = sm*F2 + a*C/(q+C) - dm*M;           % (M) MDSCs
    
    
    % COMBINE RESULTS:
    dudt = [dC; dT; dM];
  
    
end


%% FUNCTION TO SOLVE ODE
function [tout,uout] = SolveODE(F,u0,mp,t0,tf)

    [tout,uout] = ode45(F,[t0 tf],u0,[],mp);
    
    
end