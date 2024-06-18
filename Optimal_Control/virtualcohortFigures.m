% This script analyzes the virtual murine cohort from the OCHiPerGator.m 
% output and generates figures 3-6 in Anderson et al. "Optimal control 
% of combination immunotherapy in a glioblastoma-immune dynamics model"
% (code does not generate the pi charts, which were created using
% PowerPoint)

% Uses output from OCHiPerGator.m

% To get Figure 3 and 4, set figmode = 0
% To get Figure 5, set figmode = 1
% - For 5a, set yearsofsurvival = 2
% - For 5b, set yearsofsurvival = 5

%% Figure mode

figmode = 0; % 0 - during treatment (mortal, quality of life issues, healthy)
             % 1 - after treatment (dfs, pfs, recurrence, failure)

             
%% Thresholds / pre-specified values

Carryingcap = 3.04e6; % tumor carrying capacity
initialtumor = 35000; % initial tumor size

% Cancer thresholds
Cmorbid =  Carryingcap*0.5; % above this threshold, subject experience quality of life issues due to tumor size
Cmortal = Carryingcap*0.9; % above this threshold, subject dies
Cmorbidtime = 22+7; % (time above Cthreshold) (approximately two mouse "year"--bit more) -- need to add 7 days because first 7 days don't have treatment
Cthreshold = initialtumor; % if subject is above Cthreshold for a time period of Cmorbidtime, subject is categorized as having quality of life issues

% Anti-PD-1 thresholds 
antiPD1morbiddose = 15*0.8; % 15 days at top efficacy (gives wiggle room, since could be not as dangerous for longer but at lower dose)
antiPD1morbidtime = 12; % Adam et al. 2021 (time above antiPD1threshold)
antiPD1threshold = 0.75; % checking to see how many days approx at max dosage (0.8)

% CCR2 antagonist thresholds 
CCR2morbiddose = 25*0.9; % 25 days at top efficacy (gives wiggle room)
CCR2morbidtime = 25; % Flores-Toro et al. 2020 (time above CCR2threshold)
CCR2threshold = 0.85; % checking to see how many days approx at max dosage (0.9)

% Values for survival plots (DFS, PFS, recurrence, failure)
yearsofsurvival = 5; % years
tpforstd = 10; % number of time points you want to check standard deviation on for PFS
tumorcapfordfs = 1; % anything under tumorcapfordfs is considered a disease-free mouse

%% Choose which parameters were tested in OCHiPerGator.m

chooseparam = [1;     % 1  - lambdaC
               0;     % 2  - Cmax
               1;     % 3  - η 
               0;     % 4  - a_T 
               0;     % 5  - s_T 
               1;     % 6  - ρ        
               0;     % 7  - ε_C 
               1;     % 8  - r 
               0;     % 9  - d_T 
               0;     % 10 - s_M 
               1];    % 11 - d_M
     
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

axislabels = {'0',            '0.5';            % 1  - lambdaC
              '10^6',         '5 \times 10^7';  % 2  - Cmax
              '0',            '10^{-6}';        % 3  - η 
              '50',           '5 \times 10^6';  % 4  - a_T 
              '100',          '10^7';           % 5  - s_T
              '0',            '0.5';            % 6  - ρ        
              '1',            '100';            % 7  - ε_C 
              '0',            '10^{-4}';        % 8  - r 
              '0',            '0.75';           % 9  - d_T
              '0',            '0.1';            % 10 - s_M 
              '0',            '0.5'};           % 11 - d_M


np = length(chooseparam);
icp = find(chooseparam == 1);
paramnames = {'\lambda_C','C_{max}','\eta','a_T','s_T', '\rho','\epsilon_C','r','d_T','s_M','d_M'};

%% Load output
% This section loads a variable called "output."
% output{:,1} contains the sampled parameter sets (11x1).
% output{:,2} contains the output from GBMoptimalcontrol: {state, control, time, T, X}.
% output{1,3} contains the initial conditions: {maxdose, multiple, Intconstraint, initialcondition}.

load('virtualcohortOutput.mat')
% ^ here is where you load the matrix you saved from running OCHiPerGator.m
% This matrix represents the optimized regimens and cell dynamics of the
% virtual murine cohort.

% Recall: In the Anderson et al. paper, Figures 3-5 used a uniform sampling 
% of the virtual cohort, while Figure 6 sampled the cohort according to
% GBM-specific parameter distributions. (This means that, while you can 
% get all figures with one matrix, you should run through this code once
% with the uniform matrix for Figures 3-5 and a second time with the
% GBM-specific matrix to get Figure 6).


%% Extracting numerical data and determining thresholds
n = length(output(:,1)); % number of parameter sets
Psample = zeros(n,11);

% max tumor size and time of this 
Csizemax = zeros(n,1);
tmax = zeros(n,1);

% final tumor size and final time 
Cfinal = zeros(n,1);
tf = zeros(n,1);

% total cumulative dosage 
cumulativedosageantiPD1 = zeros(n,1);
cumulativedosageCCR2antagonist = zeros(n,1);

% total time tumor is above Cthreshold in size
totaltimeaboveCthreshold = zeros(n,1);

% total time treatment is above drug thresholds
totaltimeaboveantiPD1threshold = zeros(n,1);
totaltimeaboveCCR2threshold = zeros(n,1);

% final condition of the states at treatment endpoint
statesinitial = zeros(n,3);
statesmid = zeros(n,3);
statesfinal = zeros(n,3);

% Tumor size on days throughout treatment
tumoronday = ones(n,50); % tumor size on days 1 to 50

for i = 1:n
    [i]
    GBMoutput = output{i,2};
    Psample(i,:)= output{i,1};
    state = GBMoutput{1,1};
    control = GBMoutput{1,2};
    time = GBMoutput{1,3};
    T = GBMoutput{1,4};
    X = GBMoutput{1,5};
    
    
    % Maximum / Final tumor size and times
    Csizemaxpre = max(X(:,1));
    tmaxpre = max(T.*(X(:,1) == Csizemaxpre));
    Csizemaxpost = max(state(:,1));
    tmaxpost = max(time.*(state(:,1) == Csizemaxpost));
    Csizemax(i) = max(Csizemaxpre,Csizemaxpost);
    if Csizemax(i) == Csizemaxpre
        tmax(i) = tmaxpre;
    else
        tmax(i) = tmaxpost;
    end

    Cfinal(i) = state(end,1);
    tf(i) = time(end);
    cumulativedosageantiPD1(i) = trapz(time, control(:,1)); % trapezoidal integral
    cumulativedosageCCR2antagonist(i) = trapz(time, control(:,2));

    
    % Times above thresholds
    statestumor = [X(1:end-1,1); state(:,1)];
    times = [T(1:end-1,1);time];
    
    timesaboveCthreshold = (times+1).*(statestumor > Cthreshold); %need to add one in case initial value at day 0 is above threshold (function will read that day as not above the threshold if it stays as day "0")--note, if threshold is larger than initial value, the initial value will be set back to 0 anyway
    
    % Note: times below do not include the first week (or specified time) of no
    % treatment
    timesaboveantiPD1threshold = time.*(control(:,1) > antiPD1threshold);
    timesaboveCCR2threshold = time.*(control(:,2) > CCR2threshold);
    
    totaltimeaboveCthreshold(i) = totaltime(timesaboveCthreshold);
    totaltimeaboveantiPD1threshold(i) = totaltime(timesaboveantiPD1threshold);
    totaltimeaboveCCR2threshold(i) = totaltime(timesaboveCCR2threshold);
    
    
    % States at implantation, treatment start, and treatment end
    statesinitial(i,:) = X(1,:); % states at implantation (day 0)
    statesmid(i,:) = X(end,:); % states at the start of treatment
    statesfinal(i,:) = state(end,1:3); % states at the end of treatment
    
    
    % Tumor sizes throughout treatment
    for j = 1:floor(tf(i))
        indexofday = find(times>=j);
        tumoronday(i,j) = statestumor(indexofday(1));
    end

end

%% All tumor sizes at all time points from implantation to natural death

endptcheck = 800; % since the mouse is implanted with a tumor at day 0 when they are at least several weeks old, day 800 after implantation is essentially when a mouse would naturally die
tumoralldays = zeros(n,endptcheck);


for i=1:n
    [i]
    tumoralldays(i,1:floor(tf(i))) = tumoronday(i,1:floor(tf(i)));
    [Tafter,Xafter] = ode45(@GBMFuncoptimal,[floor(tf(i)):1:endptcheck],statesfinal(i,:),[],Psample(i,:));
    tumoralldays(i,floor(tf(i)):endptcheck) = Xafter(:,1);
end

% This section takes a long time, so I'd suggest doing it once, saving
% the output, and then commenting out this code and loading the saved 
% matrix for future runs
% ex: load('micealltumordays.mat')

%% Determine if therapy was a failure, PFS, or DFS

mouselifespan = 836; % days, C57BL/6 mice, 794 days for females, 878 days for males (Kunstyr 1975)
humanlifespan = 76.1; % years, according to 2022 CDC report on US lifespans
survivaldays = ceil(mouselifespan*(yearsofsurvival/humanlifespan)); 

% Failed therapy (failed therapy, ie death, if tumor surpasses Cmortal) 
failedtherapy = zeros(n,1); % if failedtherapy(i) = 1, patient i has failed therapy (bad)
                            % if failedtherapy(i) = 0, patient i has NOT failed therapy (good)

% Progression-free survival (PFS)
pfs = zeros(n,1); % if pfs(i) = 0, patient i does NOT have progression-free survival (bad)
                  % if pfs(i) = 1, patient i has progression-free survival (good)
standarddeviation = zeros(n,1);
meanend = zeros(n,1);

% Disease-free survival (DFS)
dfs = zeros(n,1); % if dfs(i) = 0, patient i does NOT have disease-free survival (bad)
                  % if dfs(i) = 1, patient i has disease-free survival (good)
                     
for i = 1:n
    [i]
      
    [Tafter,Xafter] = ode45(@GBMFuncoptimal,[tf(i):1:tf(i)+survivaldays],statesfinal(i,:),[],Psample(i,:));
    
    if max(tumoralldays(i,1:floor(tf(i)+survivaldays))) >= Cmortal
        failedtherapy(i) = 1;
    end
 
    % PFS
    gatherpfs = 0;
    if Xafter(1,1) >= Xafter(end,1) % ie decreasing after treatment
        pfs(i) = 1;
    else
        meanend(i) = mean(Xafter(end -tpforstd +1:end,1)); % last 10 timespoints
        for j =1:tpforstd
            if abs(Xafter(end -j +1,1) - meanend(i)) <= 100 % if each is within 100 cells of the mean value
                gatherpfs = gatherpfs+1;
            end
        end
        if gatherpfs == tpforstd
            pfs(i) = 1;
        end
    end
    
    % DFS 
    if Xafter(end,1) < tumorcapfordfs 
        dfs(i) = 1;
    end
end

% all disease-free patients are automatically progression-free
for i=1:n
    if dfs(i) ==1 && pfs(i)==0 % all disease-free are progression-free
        pfs(i) = 1;
    end
    if failedtherapy(i) ==1 && pfs(i)==1 % all who die are not progression-free
        pfs(i) = 0;
    end
end


% Plot dfs,pfs,failure
plotdfs = dfs;
plotfail = failedtherapy;
plotallpfs = pfs;
plotpfs = pfs; % separate vector from "plotallpfs" to show patients who are progression-free only (ie not pfs and dfs)
plotrecur = zeros(n,1); % cancer recurrence

for i=1:n
    if plotdfs(i) ==1 || plotfail(i) == 1 % ie the patient is first and foremost dfs or failure, pfs takes second place in terms of plotting
        plotpfs(i) = 0;
    end
end

checkrecur = plotdfs + plotpfs + plotfail;
for i =1:n
    if checkrecur(i) == 0 % ie not dfs, pfs, or fail (death)
        plotrecur(i) = 1;
    end
end

%% Analysis of morbidity versus mortality
% only thing that can cause mortality is tumor size being above a certain
% threshold--everything else morbidity

Cmortality = zeros(n,1);
Cmorbidity = zeros(n,1);
Csizemorbidity = zeros(n,1);
Ctimemorbidity = zeros(n,1);
antiPD1morbidity = zeros(n,1);
antiPD1timemorbidity = zeros(n,1);
antiPD1dosemorbidity = zeros(n,1);
CCR2morbidity = zeros(n,1);
CCR2timemorbidity = zeros(n,1);
CCR2dosemorbidity = zeros(n,1);
mortality = zeros(n,1);
morbidity = zeros(n,1);
drugmorbidity = zeros(n,1);


for i = 1:n
    
% TUMOR

% mortality
if Csizemax(i) > Cmortal
    Cmortality(i) = 1;
end

% morbidity
if Csizemax(i) > Cmorbid
    Cmorbidity(i) = 1;
    Csizemorbidity(i) = 1;
end

if totaltimeaboveCthreshold(i) > Cmorbidtime
    Cmorbidity(i) = 1; 
    Ctimemorbidity(i) = 1;
end

% ANTI PD 1 morbidity
    
if totaltimeaboveantiPD1threshold(i) > antiPD1morbidtime
    antiPD1morbidity(i) = 1; 
    antiPD1timemorbidity(i) = 1;
end

if cumulativedosageantiPD1(i) > antiPD1morbiddose
    antiPD1morbidity(i) = 1;
    antiPD1dosemorbidity(i) = 1;
end


% CCR2 ANTAGONIST morbidity
     
if totaltimeaboveCCR2threshold(i) > CCR2morbidtime
    CCR2morbidity(i) = 1; 
    CCR2timemorbidity(i) = 1;
end

if cumulativedosageCCR2antagonist(i) > CCR2morbiddose
    CCR2morbidity(i) = 1;
    CCR2dosemorbidity(i) = 1;
end

end

premortality = Cmortality;
premorbidity = Cmorbidity + antiPD1morbidity + CCR2morbidity;
predrugmorbidity = antiPD1morbidity + CCR2morbidity;

mortality(find(premortality~=0))=1;
morbidity(find(premorbidity~=0))=1;
drugmorbidity(find(predrugmorbidity~=0))=1;


% Compile mortality and morbidity vectors for plotting
plotmortal = mortality;
plotmorbid = morbidity;
plotantiPD1morbidity = antiPD1morbidity;
plotCCR2morbidity = CCR2morbidity;
plotCmorbidity = Cmorbidity;

plothealthy = zeros(n,1);
plotmorbidboth = zeros(n,1);
plotmorbidConly = zeros(n,1);
plotmorbiddrugonly = zeros(n,1);

for i=1:n
    if mortality(i) == 1
        plotmorbid(i) = 0; %ie mortality takes precendence over morbidity for the plots
    end
    if plotmortal(i) == 0 && plotmorbid(i) == 0 % if neither mortal nor morbid, then treatment has no issues
        plothealthy(i) = 1;
    end
    if Cmorbidity(i) == 1 && drugmorbidity(i) == 1
        plotmorbidboth(i) = 1;
    end
    if Cmorbidity(i) == 1 && drugmorbidity(i) == 0
        plotmorbidConly(i) = 1;
    end
    if Cmorbidity(i) == 0 && drugmorbidity(i) == 1
        plotmorbiddrugonly(i) = 1;
    end
end

if length(nonzeros(plotmortal + plotmorbid + plothealthy))<n || length(nonzeros(plotmortal + plotmorbidboth + plotmorbidConly + plotmorbiddrugonly + plothealthy))<n
    keyboard % you have not captured all your data points, n
end


%% Figure set up for figure 4 (as well as various other outcomes tested but not included in paper)
scatcolor = ["","","","",""];

if figmode == 0
    plottype = zeros(n,5);
    legt = ["Mortality concerns","QoL impacted by tumor and drug","QoL impacted by tumor size only","QoL impacted by drug only","No concerns"];
    
    % mortality 
    mortk = 5;
    plottype(:,mortk) = plotmortal;
    scatcolor(mortk) = 'k';

    % morbidity for both tumor and drug 
    bmorbk = 4;
    plottype(:,bmorbk) = plotmorbidboth;
    scatcolor(bmorbk) = "#A2142F";

    % morbidity for tumor only 
    tmorbk = 3; 
    plottype(:,tmorbk) = plotmorbidConly;
    scatcolor(tmorbk) = "#D95319";

    % morbidity for drug only 
    dmorbk = 2; 
    plottype(:,dmorbk) = plotmorbiddrugonly;
    scatcolor(dmorbk) = "#EDB120";
 
    % "healthy" 
    healthk = 1;
    plottype(:,healthk) = plothealthy;
    scatcolor(healthk) = "#4DBEEE";
    
end


if figmode == 1
    plottype = zeros(n,4);
    legt = ["Disease-free","Progression-free","Recurrence","Death"];
    
    % disease-free survival 
    dfsk = 1; 
    plottype(:,dfsk) = plotdfs;
    scatcolor(dfsk)="#77AC30";
    
    % progression-free survival 
    pfsk = 2; 
    plottype(:,pfsk) = plotpfs;
    scatcolor(pfsk) = "#4DBEEE";

    
    % recurrence 
    reck=3; 
    plottype(:,reck) = plotrecur;
    scatcolor(reck) = 'r';    
    
    % death 
    dk = 4;
    plottype(:,dk) = plotfail;
    scatcolor(dk) = 'k';
end


ptl = length(plottype(1,:));

Psampleplot = cell(ptl,1);
Cfinalplot = cell(ptl,1);
tfplot = cell(ptl,1);
Csizemaxplot = cell(ptl,1);
tmaxplot = cell(ptl,1);
cumulativedosageantiPD1plot = cell(ptl,1);
cumulativedosageCCR2antagonistplot = cell(ptl,1);
totaltimeaboveCthresholdplot = cell(ptl,1);
totaltimeaboveantiPD1thresholdplot = cell(ptl,1);
totaltimeaboveCCR2thresholdplot = cell(ptl,1);

for i = 1:length(plottype(1,:))
    prePsample= zeros(n,np);
    for j = 1:n
        if plottype(j,i) == 1
            prePsample(j,:) = Psample(j,:);
        end
    end
    Psampleplot{i} = reshape(nonzeros(prePsample),length(nonzeros(prePsample))/np,np);
    Cfinalplot{i} = pp(Cfinal,plottype(:,i));
    tfplot{i} = pp(tf,plottype(:,i));
    Csizemaxplot{i} = pp(Csizemax,plottype(:,i));
    tmaxplot{i} = pp(tmax,plottype(:,i));
    cumulativedosageantiPD1plot{i} = pp(cumulativedosageantiPD1,plottype(:,i));
    cumulativedosageCCR2antagonistplot{i} = pp(cumulativedosageCCR2antagonist,plottype(:,i));
    totaltimeaboveCthresholdplot{i} = pp(totaltimeaboveCthreshold,plottype(:,i));
    totaltimeaboveantiPD1thresholdplot{i} = pp(totaltimeaboveantiPD1threshold,plottype(:,i));
    totaltimeaboveCCR2thresholdplot{i} = pp(totaltimeaboveCCR2threshold,plottype(:,i));
end

%% 2D-figures of 1 parameter compared to specific outcomes (testing results--none of these figures included in paper)

fsl = 12; % font size for figure legends

numofrows = 3;
numofcols = 3;
ptitle = {'Tumor growth rate, \lambda_C';
    'Tumor carrying capacity, C_{max}';
    'T cell kill rate by tumor cells, \eta';
    'T cell activation rate, a_T';
    'T cell stimulation rate, s_T';
    'T cell inhibition rate by PD-L1-PD-1, \rho';
    'Tumor upregulation of PD-L1, \epsilon_C';
    'T cell inhibition rate by MDSCs, r';
    'T cell death rate, d_T';
    'MDSC recruitment rate, s_M';
    'MDSC death rate, d_M'};

for i = 1:length(icp)
figure(i)
pfig(i) = tiledlayout(numofrows,numofcols,'TileSpacing','compact');

% Final tumor size
nexttile(1);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on
for j = 1:length(Psampleplot)
    scatter(Psampleplot{j}(:,icp(i)),Cfinalplot{j},'filled',"MarkerFaceColor",scatcolor(j))
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('tumor cells (C)')
title('Final tumor size')

% Endpoint of treatment (time)
nexttile(2);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on
for j = 1:length(Psampleplot)
    scatter(Psampleplot{j}(:,icp(i)),tfplot{j},'filled',"MarkerFaceColor",scatcolor(j))
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('time (days)')
title('Time of treatment conclusion') 

% Maximum tumor size   
ax1 = nexttile(3);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on 
for j = 1:length(Psampleplot)
    m1(j) = scatter(Psampleplot{j}(:,icp(i)),Csizemaxplot{j},'filled',"MarkerFaceColor",scatcolor(j));
end
if figmode == 0
    leg1 = legend(ax1,[m1(mortk),m1(bmorbk),m1(tmorbk),m1(dmorbk),m1(healthk)],[legt,'Example patient'],'FontSize',fsl);
    leg1.Layout.Tile = 'south';
else
    leg1 = legend(ax1,[m1(dfsk),m1(pfsk),m1(reck),m1(dk)],[legt,'Example patient'],'FontSize',fsl);
    leg1.Layout.Tile = 'south';
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('tumor cells (C)')
title('Maximum tumor size')


% Time of maximum tumor size
nexttile(4);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on
for j = 1:length(Psampleplot)
    scatter(Psampleplot{j}(:,icp(i)),tmaxplot{j},'filled',"MarkerFaceColor",scatcolor(j))
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('time (days)')
title('Time of maximum tumor size')

% Total cumulative dosage of anti-PD-1
nexttile(5);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on
for j = 1:length(Psampleplot)
    scatter(Psampleplot{j}(:,icp(i)),cumulativedosageantiPD1plot{j},'filled',"MarkerFaceColor",scatcolor(j))
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('cumulative % reduction') 
title('Anti-PD1-1 cumulative % reduction')

% Total cumulative dosage of CCR2 antagonist
nexttile(6);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on
for j = 1:length(Psampleplot)
    scatter(Psampleplot{j}(:,icp(i)),cumulativedosageCCR2antagonistplot{j},'filled',"MarkerFaceColor",scatcolor(j))
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('cumulative % reduction')
title('CCR2 antagonist cumulative % reduction')

% Time above tumor threshold 
nexttile(7);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on
for j = 1:length(Psampleplot)
    scatter(Psampleplot{j}(:,icp(i)),totaltimeaboveCthresholdplot{j},'filled',"MarkerFaceColor",scatcolor(j))
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('time (days)')
title('Time above tumor threshold')

% Time above anti-PD-1 threshold 
nexttile(8);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on
for j = 1:length(Psampleplot)
    scatter(Psampleplot{j}(:,icp(i)),totaltimeaboveantiPD1thresholdplot{j},'filled',"MarkerFaceColor",scatcolor(j))
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('time (days)')
title('Time above anti-PD-1 threshold')

% Time above CCR2 antagonist threshold 
nexttile(9);
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
hold on
for j = 1:length(Psampleplot)
    scatter(Psampleplot{j}(:,icp(i)),totaltimeaboveCCR2thresholdplot{j},'filled',"MarkerFaceColor",scatcolor(j))
end
hold off
xlim(mpr(icp(i),:))
xlabel(paramnames{icp(i)})
ylabel('time (days)')
title('Time above CCR2 inhibitor threshold')


set(gcf,'Position',[0 1000 950 1000],'PaperPositionMode','auto');   
title(pfig(i),ptitle{icp(i)},'FontSize',20)
end

%% 2D-figures of specific outcomes (includes Figure 4 regarding cumulative % reduction of both drugs)
fsl = 16; % font size for legend

% Final tumor size and time
ax2 = figure(20);
m21 = scatter(Cfinalplot{1},tfplot{1},'filled',"MarkerFaceColor",scatcolor(1));
hold on
m22 = scatter(Cfinalplot{2},tfplot{2},'filled',"MarkerFaceColor",scatcolor(2));
m23 = scatter(Cfinalplot{3},tfplot{3},'filled',"MarkerFaceColor",scatcolor(3));
m24 = scatter(Cfinalplot{4},tfplot{4},'filled',"MarkerFaceColor",scatcolor(4));
if figmode == 0
    m25 = scatter(Cfinalplot{5},tfplot{5},'filled',"MarkerFaceColor",scatcolor(5));
    leg2 = legend([m25,m24,m23,m22,m21],legt,'FontSize',fsl);
else
    leg2 = legend([m21,m22,m23,m24],legt,'FontSize',fsl);
end
hold off
set(leg2,'Location','southoutside') 
xlabel('tumor cells (C)')
ylabel('time (days)')
title('Final tumor size and time')
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
set(gcf,'Position',[0 500 550 500],'PaperPositionMode','auto');  

% Maximum tumor size and time
ax3 = figure(21);
m31 = scatter(Csizemaxplot{1},tmaxplot{1},'filled',"MarkerFaceColor",scatcolor(1));
hold on
m32 = scatter(Csizemaxplot{2},tmaxplot{2},'filled',"MarkerFaceColor",scatcolor(2));
m33 = scatter(Csizemaxplot{3},tmaxplot{3},'filled',"MarkerFaceColor",scatcolor(3));
m34 = scatter(Csizemaxplot{4},tmaxplot{4},'filled',"MarkerFaceColor",scatcolor(4));
if figmode == 0
    m35 = scatter(Csizemaxplot{5},tmaxplot{5},'filled',"MarkerFaceColor",scatcolor(5));
    leg3 = legend([m35,m34,m33,m32,m31],legt,'FontSize',fsl);
else
    leg3 = legend([m31,m32,m33,m34],legt,'FontSize',fsl);
end
hold off
set(leg3,'Location','southoutside');
xlabel('maximum tumor size (cells)')
ylabel('time of maximum tumor size (days)')
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
set(gcf,'Position',[0 500 550 500],'PaperPositionMode','auto');  

% Total cumulative dosage of drugs (FIGURE 4)
ax4 = figure(22);
m41 = scatter(cumulativedosageantiPD1plot{1},cumulativedosageCCR2antagonistplot{1},'filled',"MarkerFaceColor",scatcolor(1));
hold on
m42 = scatter(cumulativedosageantiPD1plot{2},cumulativedosageCCR2antagonistplot{2},'filled',"MarkerFaceColor",scatcolor(2));
m43 = scatter(cumulativedosageantiPD1plot{3},cumulativedosageCCR2antagonistplot{3},'filled',"MarkerFaceColor",scatcolor(3));
m44 = scatter(cumulativedosageantiPD1plot{4},cumulativedosageCCR2antagonistplot{4},'filled',"MarkerFaceColor",scatcolor(4));
if figmode == 0
    m45 = scatter(cumulativedosageantiPD1plot{5},cumulativedosageCCR2antagonistplot{5},'filled',"MarkerFaceColor",scatcolor(5));
    leg4 = legend([m45,m44,m43,m42,m41],legt,'FontSize',fsl);
else
    leg4 = legend([m41,m42,m43,m44],legt,'FontSize',fsl);
end
hold off
set(leg4,'Location','northwest') 
xlabel('anti-PD-1 cumulative % reduction')
ylabel('CCR2 antagonist cumulative % reduction') 
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
set(gcf,'Position',[0 500 550 500],'PaperPositionMode','auto');  

% Time above drug thresholds
ax5 = figure(23);
m51 = scatter(totaltimeaboveantiPD1thresholdplot{1},totaltimeaboveCCR2thresholdplot{1},'filled',"MarkerFaceColor",scatcolor(1));
hold on
m52 = scatter(totaltimeaboveantiPD1thresholdplot{2},totaltimeaboveCCR2thresholdplot{2},'filled',"MarkerFaceColor",scatcolor(2));
m53 = scatter(totaltimeaboveantiPD1thresholdplot{3},totaltimeaboveCCR2thresholdplot{3},'filled',"MarkerFaceColor",scatcolor(3));
m54 = scatter(totaltimeaboveantiPD1thresholdplot{4},totaltimeaboveCCR2thresholdplot{4},'filled',"MarkerFaceColor",scatcolor(4));
if figmode == 0
    m55 = scatter(totaltimeaboveantiPD1thresholdplot{5},totaltimeaboveCCR2thresholdplot{5},'filled',"MarkerFaceColor",scatcolor(5));
    leg5 = legend([m55,m54,m53,m52,m51],legt,'FontSize',fsl);
else
    leg5 = legend([m51,m52,m53,m54],legt,'FontSize',fsl);
end
hold off
set(leg5,'Location','southoutside') 
xlabel('time (days) - anti-PD-1') % can move these to be in correct position with arrow button after plotting
ylabel('time (days) - CCR2 antagonist')
title('Time above drug thresholds')
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
set(gcf,'Position',[0 500 550 500],'PaperPositionMode','auto');  

% Time above initial tumor size versus final time
ax6 = figure(24);
m61 = scatter(totaltimeaboveCthresholdplot{1},tfplot{1},'filled',"MarkerFaceColor",scatcolor(1));
hold on
m62 = scatter(totaltimeaboveCthresholdplot{2},tfplot{2},'filled',"MarkerFaceColor",scatcolor(2));
m63 = scatter(totaltimeaboveCthresholdplot{3},tfplot{3},'filled',"MarkerFaceColor",scatcolor(3));
m64 = scatter(totaltimeaboveCthresholdplot{4},tfplot{4},'filled',"MarkerFaceColor",scatcolor(4));
if figmode == 0
    m65 = scatter(totaltimeaboveCthresholdplot{5},tfplot{5},'filled',"MarkerFaceColor",scatcolor(5));
    leg2 = legend([m65,m64,m63,m62,m61],legt,'FontSize',fsl);
else
    leg2 = legend([m61,m62,m63,m64],legt,'FontSize',fsl);
end
hold off
set(leg2,'Location','southoutside') 
xlabel('time above initial tumor size (days)')
ylabel('final time (days)')
title('Final tumor size and time')
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
set(gcf,'Position',[0 500 550 500],'PaperPositionMode','auto');  

% Maximum tumor size and time
ax3 = figure(25);
m31 = scatter(Csizemaxplot{1},tmaxplot{1},'filled',"MarkerFaceColor",scatcolor(1));
hold on
m32 = scatter(Csizemaxplot{2},tmaxplot{2},'filled',"MarkerFaceColor",scatcolor(2));
m33 = scatter(Csizemaxplot{3},tmaxplot{3},'filled',"MarkerFaceColor",scatcolor(3));
m34 = scatter(Csizemaxplot{4},tmaxplot{4},'filled',"MarkerFaceColor",scatcolor(4));
if figmode == 0
    m35 = scatter(Csizemaxplot{5},tmaxplot{5},'filled',"MarkerFaceColor",scatcolor(5));
    leg3 = legend([m35,m34,m33,m32,m31],legt,'FontSize',fsl);
else
    leg3 = legend([m31,m32,m33,m34],legt,'FontSize',fsl);
end
hold off
set(gca,'LineWidth',1,'FontSize',fsl+2,'FontName','Arial')
set(leg3,'Location','southoutside','FontSize',14);
xlabel('tumor cells (C)')
ylabel('time (days)')
title('Maximum tumor size and time')


%% Figure set up for Figure 3 and 5

if figmode == 0
    plottype = zeros(n,5);
    legt = ["Mortality","QoL impacted by tumor","QoL impacted by anti-PD-1","QoL impacted by CCR2 antagonist","No concerns"];
    
    % mortality 
    plottype(:,5) = plotmortal;
    scatcolor(5) = 'k';

    % morbidity for tumor
    plottype(:,4) = plotCmorbidity;
    scatcolor(4) = "#A2142F";

    % morbidity for anti-PD-1
    plottype(:,3) = plotantiPD1morbidity;
    scatcolor(3) = "#EDB120";

    % morbidity for CCR2
    plottype(:,2) = plotCCR2morbidity;
    scatcolor(2) = "#77AC30";
 
    % "healthy" 
    plottype(:,1) = plothealthy;
    scatcolor(1) = "#4DBEEE";
    
end


if figmode == 1
    plottype = zeros(n,4);
    legt = ["Disease-free","Progression-free","Recurrence","Death"];

    % disease-free survival 
    plottype(:,2) = plotdfs;
    scatcolor(2)="#77AC30";
    
    % progression-free survival 
    plottype(:,1) = plotallpfs;
    scatcolor(1) = "#4DBEEE";
    
    % recurrence 
    plottype(:,3) = plotrecur;
    scatcolor(3) = 'r';    
    
    % death 
    plottype(:,4) = plotfail;
    scatcolor(4) = 'k';
end


ptl = length(plottype(1,:));

Psampleplot = cell(ptl,1);
Cfinalplot = cell(ptl,1);
tfplot = cell(ptl,1);
Csizemaxplot = cell(ptl,1);
tmaxplot = cell(ptl,1);
cumulativedosageantiPD1plot = cell(ptl,1);
cumulativedosageCCR2antagonistplot = cell(ptl,1);
totaltimeaboveCthresholdplot = cell(ptl,1);
totaltimeaboveantiPD1thresholdplot = cell(ptl,1);
totaltimeaboveCCR2thresholdplot = cell(ptl,1);

for i = 1:length(plottype(1,:))
    prePsample= zeros(n,np);
    for j = 1:n
        if plottype(j,i) == 1
            prePsample(j,:) = Psample(j,:);
        end
    end
    Psampleplot{i} = reshape(nonzeros(prePsample),length(nonzeros(prePsample))/np,np);
    Cfinalplot{i} = pp(Cfinal,plottype(:,i));
    tfplot{i} = pp(tf,plottype(:,i));
    Csizemaxplot{i} = pp(Csizemax,plottype(:,i));
    tmaxplot{i} = pp(tmax,plottype(:,i));
    cumulativedosageantiPD1plot{i} = pp(cumulativedosageantiPD1,plottype(:,i));
    cumulativedosageCCR2antagonistplot{i} = pp(cumulativedosageCCR2antagonist,plottype(:,i));
    totaltimeaboveCthresholdplot{i} = pp(totaltimeaboveCthreshold,plottype(:,i));
    totaltimeaboveantiPD1thresholdplot{i} = pp(totaltimeaboveantiPD1threshold,plottype(:,i));
    totaltimeaboveCCR2thresholdplot{i} = pp(totaltimeaboveCCR2threshold,plottype(:,i));
end

%% Scatter plots and distributions for figure 3 and 5
fsl = 18; % font size
sdim = length(icp);

figure
scatterdist = tiledlayout(sdim,sdim,'TileSpacing','compact');
for i = 1:sdim
    for j = 1:sdim
        if icp(i) > icp(j)
            ax2 = nexttile(sdim*(i-1) + j); %position in block
            
            m1 = scatter(Psampleplot{1}(:,icp(j)),Psampleplot{1}(:,icp(i)),'filled',"MarkerFaceColor",scatcolor(1));
            hold on 
            m2 = scatter(Psampleplot{2}(:,icp(j)),Psampleplot{2}(:,icp(i)),'filled',"MarkerFaceColor",scatcolor(2));
            m3 = scatter(Psampleplot{3}(:,icp(j)),Psampleplot{3}(:,icp(i)),'filled',"MarkerFaceColor",scatcolor(3));
            m4 = scatter(Psampleplot{4}(:,icp(j)),Psampleplot{4}(:,icp(i)),'filled',"MarkerFaceColor",scatcolor(4));
            if figmode == 0
                m5 = scatter(Psampleplot{5}(:,icp(j)),Psampleplot{5}(:,icp(i)),'filled',"MarkerFaceColor",scatcolor(5));
            end
            hold off          
            
            if i ~= sdim
                xlim(mpr(icp(j),:))
                set(gca,'XTick',[])
            end
            if i == sdim
                xlabel(paramnames{icp(j)})
                xticks(mpr(icp(j),:))
                xticklabels(axislabels(icp(j),:))
                set(gca,'FontSize',fsl)
            end
            if j ~= 1
                ylim(mpr(icp(i),:))
                set(gca,'YTick',[])
            end
            if j == 1
                ylabel(paramnames{icp(i)})
                set(gca,'YTick',mpr(icp(i),:),'YTickLabel',axislabels(icp(i),:),'FontSize',fsl)
            end
        end
        if i == j
            [y1,x1] = ksdensity(Psampleplot{1}(:,icp(i)));
            [y2,x2] = ksdensity(Psampleplot{2}(:,icp(i)));
            [y3,x3] = ksdensity(Psampleplot{3}(:,icp(i)));
            [y4,x4] = ksdensity(Psampleplot{4}(:,icp(i)));
            if figmode == 0 
                [y5,x5] = ksdensity(Psampleplot{5}(:,icp(i)));
            end
           
            nexttile(sdim*(i-1)+i)
            plot(x1,y1,'Color',scatcolor(1),'LineWidth',1);
            hold on
            plot(x2,y2,'Color',scatcolor(2),'LineWidth',1);
            plot(x3,y3,'Color',scatcolor(3),'LineWidth',1);
            plot(x4,y4,'Color',scatcolor(4),'LineWidth',1);
            if figmode == 0
                plot(x5,y5,'Color',scatcolor(5),'LineWidth',1);
            end
            hold off
     
            if i==sdim
                xticks(mpr(icp(i),:))
                xticklabels(axislabels(icp(i),:))
                set(gca,'YTick',[],'FontSize',fsl)
            else
                set(gca,'XTick',[],'YTick',[])
            end
            set(gcf,'Position',[0 300 1000 800],'PaperPositionMode','auto');
            if i == sdim
                xlabel(paramnames{icp(j)})
                xticks(mpr(icp(j),:))
                xticklabels(axislabels(icp(j),:))
            end
        end
    end
end 
if figmode == 0
    leg = legend(ax2,[m5,m4,m3,m2,m1],legt,'FontSize',fsl-3);
    leg.Layout.Tile = sdim-2;
else
    leg = legend(ax2,[m2,m1,m3,m4],legt,'FontSize',fsl);
    leg.Layout.Tile = sdim-2;
end

%% Check day of death / day of disease-free survival for Figure 6b
dayofdeath = zeros(n,1);
daydiseasefree = zeros(n,endptcheck);

for i = 1:n
    [i]
    
    % Day of death
    deaddays = find(tumoralldays(i,:) >= Cmortal);
    if isempty(deaddays) == 1
        dayofdeath(i) = endptcheck;
    else
        dayofdeath(i) = deaddays(1);
    end
    
    % Day of disease-free survival
    ddfi = find(tumoralldays(i,:) < 1); % index of days which are disease-free
    daydiseasefree(i,ddfi) = 1;
end


%% Figure 6b: survival / disease-free curves
fsl = 16; % font size

% Set up
percentsurviveday = zeros(endptcheck,2);
percentdiseasefree = zeros(endptcheck,2);
for i = 1:endptcheck
    percentsurviveday(i,1) = i;
    percentsurviveday(i,2) = length(find(dayofdeath > i))/n;
    percentdiseasefree(i,1) = i;
    percentdiseasefree(i,2) = length(nonzeros(daydiseasefree(:,i)))/n;
end


% Plot
figure(40)
percentfig = tiledlayout(1,2,'TileSpacing','compact');
% Percent survival (Kaplan-Meyer Survival Curve)
nexttile(1);
plot(percentsurviveday(1:endptcheck,1), 100*percentsurviveday(1:endptcheck,2),'LineWidth',2); 
ylim([0 100])
xlabel('Days after implantation')
ylabel('% Survival')
set(gca,'LineWidth',1,'FontSize',fsl,'FontName','Arial')
% Percent Disease-free
nexttile(2);
plot(percentdiseasefree(1:endptcheck,1), 100*percentdiseasefree(1:endptcheck,2),'LineWidth',2); 
ylim([0 100])
xlabel('Days after implantation')
ylabel('% Disease-free')
set(gca,'LineWidth',1,'FontSize',fsl,'FontName','Arial')


%% Calculate total time above threshold
function g = totaltime(input1)

if length(input1) >20000
    skip = ceil(length(input1)/15000);
    n = floor(length(input1)/skip);
    input = zeros(n,1);
    for i = 1:n
        input(i) = input1(skip*i);
    end
else 
    input = input1;
end



index = find(input == 0);
n = length(input);
k = length(index);
firstmatrix = zeros(n,k);
secondmatrix = zeros(n,k);
thirdvector = zeros(n,1);


if isempty(index) == 1
    g = max(input) - min(input);
else

for i = 1:k
    firstmatrix(:,i) = input(:) - [zeros(index(i)-1,1);input(index(i):end)];
end

for i = 1:k
    if i ==1 || sum(firstmatrix(:,i)) ==0
        secondmatrix(:,i) = firstmatrix(:,i);
    else
        secondmatrix(:,i) = firstmatrix(:,i) - firstmatrix(:,i-1);
    end
    if isempty(nonzeros(secondmatrix(:,i))) == 1
        thirdvector(i) = 0;
    else
        thirdvector(i) = max(secondmatrix(:,i))-min(nonzeros(secondmatrix(:,i)));
    end
end


g = sum(thirdvector);
end

end

%% Determine points to plot for different categories
function f = pp(data,ptstokeep)

v1 = data.*ptstokeep;
for i = 1:length(data)
    if data(i) == 0 && ptstokeep(i) == 1
        v1(i) = -1; % ensures we keep data points that are actually 0
    end
end

v2 = nonzeros(v1);
for i = 1:length(v2)
    if v2(i) == -1
        v2(i) = 0; % change them back to 0 after the nonzeros function
    end
end


f = v2;
end
