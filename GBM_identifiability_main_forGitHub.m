% GBM Code for Practical Identifiability 
% Original Source: Marisa Eisenberg (marisae@umich.edu) 

% This code uses: GBMFuncidentifiable.m, GBMCost.m, GBMMiniFisher.m, GBMProfLife.m

clear
                
%% Load Data

addpath('/Users/hannahanderson/MATLAB-Drive')

data = NaN(6,2,3);

% TUMOR
data(:,:,1) =  table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',1,'Range','E3:F8'));
% two columns of data containing the average number of tumor cells (col 2)
% on days 7,13, 20, 24, 27, 34 (in that order) (col 1)

% T CELL
data(1:4,:,2) =  table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',2,'Range','E3:F6'));
% two columns of data containing the average number of T cells (column 2)
% on days 7, 24, 27, 34 (in that order)

% MDSC
data(:,:,3) =  table2array(readtable('GBMtreatmentfreedata.xlsx','Sheet',3,'Range','E3:F8'));
% two columns of data containing the average number of MDSCs (column 2)
% on days 7, 13, 20, 24, 27, 34 (in that order) (column 1)

times = [7; 13; 20; 24; 27; 34]; % same as times for tumor and MDSC (T cell has two less data points)


%% Setup
initialcondition = [35000;100;0];

% choose parameters to vary:
chooseparam = [1;     % 1  - lambdaC
               1;     % 2  - Cmax
               1;     % 3  - η 
               1;     % 4  - a_T 
               1;     % 5  - s_T 
               1;     % 6  - ρ        
               1;     % 7  - ε_C 
               1;     % 8  - r 
               1;     % 9  - d_T 
               1;     % 10 - s_M 
               1];    % 11 - d_M
           
indexchooseparam = find(chooseparam == 1);
           
% Starting parameter values:
fixed = [0.431;      % 1  - lambdaC
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

params = nonzeros(chooseparam.*fixed);  

paramnames = {'\lambda_C','C_{max}','\eta','a_T','s_T', '\rho','\epsilon_C','r','d_T','s_M','d_M'};

np = length(indexchooseparam); % number of parameters


%% Simulate and Plot the Model

[t,x] = ode45(@GBMFuncidentifiable,times,initialcondition,[],fixed,params,indexchooseparam);

figure(1)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(t,x(:,1),'b','LineWidth',2);
    plot(data(:,1,1),data(:,2,1),'ko','LineWidth',2);
    ylabel('Tumor Population');  
    xlabel('Time (days)');
    
figure(2)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(t,x(:,2),'b','LineWidth',2);
    plot(data(1:4,1,2),data(1:4,2,2),'ko','LineWidth',2);
    ylabel('T cell Population');  
    xlabel('Time (days)');
    
figure(3)
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(t,x(:,3),'b','LineWidth',2);
    plot(data(:,1,3),data(:,2,3),'ko','LineWidth',2);
    ylabel('MDSC Population');  
    xlabel('Time (days)');
    
    
%% Parameter Estimation

% Estimate the model parameters 
[paramests, fval, exitflag] = fminsearch(@(p) GBMCost(data(:,1,1:3),fixed,p,indexchooseparam,data(:,2,1:3),initialcondition),params,optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000));

% Note: if exitflag = 1, the function converged to a solution
%       if exitflag = 0, the number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.MaxFunEvals.
%       if exitflag = -1, the algorithm was terminated by the output function.

paramests = abs(paramests);
% Re-simulate the model with the final parameter estimates
[test,xest] = ode45(@GBMFuncidentifiable,times,initialcondition,[],fixed,paramests,indexchooseparam);
  

%% Plot results

% Model Fit
figure(1) 
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
     hold on
    plot(test(:,1),xest(:,1),'k','LineWidth',2);
    plot(data(:,1,1),data(:,2,1),'ko','LineWidth',2);
    legend('Initial Simulation with Starting Parameters','Data','Model with Parameter Estimates','Location','se'); 
    ylabel('Tumor Population');  
    xlabel('Time (days)');
    
figure(2) 
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
     hold on
    plot(test(:,1),xest(:,2),'k','LineWidth',2);
    plot(data(1:4,1,2),data(1:4,2,2),'ko','LineWidth',2);
    legend('Initial Simulation with Starting Parameters','Data','Model with Parameter Estimates','Location','se'); 
    ylabel('T cell Population');  
    xlabel('Time (days)');
    
figure(3)  
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
     hold on
    plot(test(:,1),xest(:,3),'k','LineWidth',2);
    plot(data(:,1,3),data(:,2,3),'ko','LineWidth',2);
    legend('Initial Simulation with Starting Parameters','Data','Model with Parameter Estimates','Location','se'); 
    ylabel('MDSC Population');  
    xlabel('Time (days)');

% Residuals of the fit (shows how off the model is to the data)    
figure(4)  
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(test(:,1),(data(:,2,1))-xest(:,1),'o');
    ylabel('Residuals (tumor)');  
    xlabel('Time (days)'); 
    
figure(5)  
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(test(:,1),(data(:,2,2))-xest(:,2),'o');
    ylabel('Residuals (T cells)');  
    xlabel('Time (days)'); 
    
figure(6)  
    set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial')
    hold on
    plot(test(:,1),(data(:,2,3)-xest(:,3)),'o');
    ylabel('Residuals (MDSCs)');  
    xlabel('Time (days)'); 

    
%% Fisher Information Matrix
FIM = NaN(np,np,3);
completeFIM = zeros(np,np);

% Numerically approximate the FIM for each variable individually
for variablenum = 1:3
    FIM(:,:,variablenum) = GBMMiniFisher(data(:,1,1),fixed,paramests,indexchooseparam,initialcondition, variablenum);
end

tumorrref = rref(FIM(:,:,1));
tcellrref = rref(FIM(:,:,2)); 
MDSCrref = rref(FIM(:,:,3));

%Calculate & print rank of FIM
ranktumor = rank(tumorrref)
ranktcell = rank(tcellrref)
rankMDSC = rank(MDSCrref)


%% Generate Profile Likelihoods

% Wrapper function for parameter estimation
costfun = @(p) GBMCost(data(:,1,1:3),fixed,p,indexchooseparam,data(:,2,1:3),initialcondition);

threshold = chi2inv(0.95,np)/2 + fval;

numpoints = 10; % increasing numpoints should result in a smoother profile
profrange = 0.75; %percent of range for profile to run across

profiles = NaN(numpoints*2,3+np,np);

parfor i=1:np
    [i]
    %Generate a profile for parameter i, using paramests as the starting
    %value and the fitter to do the parameter estimation:
    profiles(:,:,i) = GBMProfLike(paramests,i,costfun,profrange,numpoints);
        % each profile has columns: profiled parameter value, resulting
        % cost-function (e.g. RSS) value, any flags from the optimizer, and
        % then columns for each of the other parameter estimates.
end


%% Plot Profiles

numofrows = ceil(length(indexchooseparam)/3);

figure
profilefig = tiledlayout(numofrows,3,'TileSpacing','compact');
for i = 1:length(indexchooseparam)
    nexttile(i);
    set(gca,'LineWidth',1,'FontSize',18,'FontName','Arial')
    hold on
    plot(profiles(:,1,i),profiles(:,2,i),'k','LineWidth',2) %identifiability curve
    %plot(paramests(i),fval,'r*','LineWidth',2) % technically min point found by fminsearch (ie best parameter)
    plot(profiles(10,1,i),profiles(10,2,i),'r*','LineWidth',2) % better at capturing min point found by fminsearch
    plot([profiles(1,1,i) profiles(end,1,i)],[threshold threshold],'r--') % red dotted threshhold to compare identifiability curve to
    xlabel(paramnames{indexchooseparam(i)})
end
profilefig.YLabel.String = 'Error in best fit';
profilefig.YLabel.FontSize = 20;
set(gcf,'Position',[0 1000 1000 1000],'PaperPositionMode','auto');   
title(profilefig,'Profile Likelihoods','FontSize',18)

