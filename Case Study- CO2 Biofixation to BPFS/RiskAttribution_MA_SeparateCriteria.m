%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% STAKEHOLDER'S RISK ATTRIBUTION & UNCERTAINTY ANALYSIS %%%%%%%%%%
%%%%%%%%%%%%%%%%         CO2 BIOFIXATION TO BPFS           %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Jeehwan Lee, KAIST, stevelee@kaist.ac.kr %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clc;
tic;
format longG;
progress = waitbar(0, 'Running...', 'Name', 'Running Risk Analysis (Individual Criteria)...');
total_steps = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [1] Initialize Parameters and Spaces
% (a) Define # of Model Evaluation (# of MC sampling of Parameter Sets)
n_sim = 5000;
% (b) Automatically generate Bayesian target values, or manually enter them
bayes_auto = 1;
bayes_target = questdlg('Automatically generate Stakeholders Sustainability Criteria?',...
        'USER INPUT',...
        'No','Yes','No');
switch bayes_target
    case 'Yes'
        bayes_auto = 1;
    case 'No'
        bayes_auto = 0;
end
% (c) Initialize Parameter Values for Key Parameters with kp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];        
alpha = 0.05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
% (d) Scalar # of Key Parameters. Count columns of kp array.
n_kp = size(kp, 2);
% (d) Define two identical parameter spaces for Sobol Analysis
%     If n_kp params are simulated n_sim times, space is n_sim x n_kp 
ParSpace = [];                 % Space for the Key Parameter of Interest
c_ParSpace = [];               % Space for Complementary Parameters to kp
% (e) Define Resolution of Kernel Distributions
np_1 = 20000;               % Resolution for Parameter KDs for MC sampling
np_2 = 20000;               % Resolution for Bayesian KDs for computing area difference
%%Progress Bar%%
waitbar(20/total_steps, progress, 'Generating Parameter Spaces...');




%% [2] Populate Parameter Space via MC Simulation
% Open bootstrap results
[FileName, PathName] = uigetfile('.mat', 'Select PE_Bootstrap.mat');
if isequal(FileName, 0)
    disp('User has selected Cancel')
else
    disp([fullfile(FileName), ' has been loaded'])
end
Bootstrap = load(FileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Define Parametric Distributions for Key Params]
% [List Kernel Data Points]
data_x1 = Bootstrap.Param_Array(1,:);                                      % Bootstrap KD for Mu_max
data_x2 = Bootstrap.Param_Array(2,:);                                      % Bootstrap KD for KI
data_x3 = Bootstrap.Param_Array(3,:);                                      % Bootstrap KD for T_opt
data_x4 = Bootstrap.Param_Array(4,:);                                      % Bootstrap KD for Theta
data_x5 = Bootstrap.Param_Array(5,:);                                      % Bootstrap KD for k_D
data_x6 = Bootstrap.Param_Array(6,:);                                      % Bootstrap KD for Pc                                              % [KOREA] GWI of Process Water, tonCO2eq/ton
%%%%%%%%%%%%%%% Remove Outliers
data_x1 = data_x1(data_x1(:) >= 0.0035 & data_x1(:) < 0.0038);             % Bootstrap KD for Mu_max
data_x2 = data_x2(data_x2(:) >= 385 & data_x2(:) < 390);                   % Bootstrap KD for KI
data_x3 = data_x3(data_x3(:) >= 25.5 & data_x3(:) < 28.5);                 % Bootstrap KD for T_opt
data_x4 = data_x4(data_x4(:) >= 8.5 & data_x4(:) < 9.5);                   % Bootstrap KD for Theta
data_x5 = data_x5(data_x5(:) >= 0.00103 & data_x5(:) < 0.00118);           % Bootstrap KD for k_D
data_x6 = data_x6(data_x6(:) >= 187.618 & data_x6(:) < 187.627);           % Bootstrap KD for Pc
%%%%%%%%%%%%%%%
% Endogenous: Separation
data_x7 = [0.9, 0.9, 0.91, 0.93];                                          % KD for DAF Recovery
data_x8 = [0.1, 0.175, 0.2, 0.225, 0.25];                                  % KD for slurry solids
% Exogenous: Price of Raw Materials
data_x9 = [90, 120, 140];                                                  % KD for NH4Cl Price
data_x10 = [1200, 1425, 1650];                                             % KD for KH2PO4 Price
data_x11 = [1200, 1350, 1400];                                             % KD for K2HPO4 Price
data_x12 = [75, 100, 128];                                                 % KD for MgSO4-7H2O Price
data_x13 = [100, 110, 140, 190];                                           % KD for CaCl2-2H2O Price
data_x14 = [1100, 1520, 1900];                                             % KD for Na2-EDTA Price
data_x15 = [1222, 1285, 1189, 1134, 1050];                                 % KD for LDPE Film
data_x16 = [0.681, 0.749, 0.619];                                          % KD for Water Price
% Exogenous: Price of Utilities
data_x17 = [13.41, 13.41, 15.08, 25.89, 18.73, 25.89, 45.27, 25.96, 39.3, 25.49, 18.35, 25.49];  % KD for Grid Mix Electricity Price                                            % KD for NH4Cl Price
% Exogenous: GWP of Raw Materials
data_x18 = [1.18, 0.91, 1.867464];                              % KD for Emission Factor NH4Cl Production
data_x19 = [2.8636, 2.522826];                                  % KD for Emission Factor KH2PO4, K2HPO4 Production
data_x20 = [0.3, 0.32, 0.15, 0.45];                             % KD for Emission Factor, MgSO4-7H2O
data_x21 = [0.89, 0.445, 1.335];                                % KD for Emission Factor, CaCl2-2H2O
data_x22 = [1.845, 2.06, 1.72, 2.472];                          % KD for Emission Factor, LDPE
data_x23 = [0.000173, 0.000123];                                % KD for Emission Factor, Process Water
% Exogenous: GWP of Utilities
data_x24 = [0.1852779, 0.153889, 0.1294445, 0.1375];            % KD for Emission Factor, Grid Mix

% [Define Kernel Distributions]
f_x1 = fitdist(data_x1', 'Kernel', 'Support', [0.0035, 0.0038]);   % Mu_max
f_x2 = fitdist(data_x2', 'Kernel', 'Support', [385, 390]);         % KI
f_x3 = fitdist(data_x3', 'Kernel', 'Support', [25.5, 28.5]);       % T_opt
f_x4 = fitdist(data_x4', 'Kernel', 'Support', [8.5, 9.5]);         % Theta
f_x5 = fitdist(data_x5', 'Kernel', 'Support', [0.00103, 0.00118]); % Death Rate
f_x6 = fitdist(data_x6', 'Kernel', 'Support', [187.618, 187.627]); % Periodicity Coefficient
f_x7 = fitdist(data_x7', 'Kernel', 'Support', [0.5, 1]);           % Biomass Recovery
f_x8 = fitdist(data_x8', 'Kernel', 'Support', [0.05, 1]);          % Slurry Solid %
f_x9 = fitdist(data_x9', 'Kernel', 'Support', 'Positive');         % Price, NH4Cl
f_x10 = fitdist(data_x10', 'Kernel', 'Support', 'Positive');       % Price, KH2PO4
f_x11 = fitdist(data_x11', 'Kernel', 'Support', 'Positive');       % Price, K2HPO4
f_x12 = fitdist(data_x12', 'Kernel', 'Support', 'Positive');       % Price, MgSO4-7H2O
f_x13 = fitdist(data_x13', 'Kernel', 'Support', 'Positive');       % Price, CaCL2-2H2O
f_x14 = fitdist(data_x14', 'Kernel', 'Support', 'Positive');       % Price, Na2-EDTA
f_x15 = fitdist(data_x15', 'Kernel', 'Support', 'Positive');       % Price, LDPE
f_x16 = fitdist(data_x16', 'Kernel', 'Support', [0.01, 5]);        % Price, Water
f_x17 = fitdist(data_x17', 'Kernel', 'Support', 'Positive');       % Price, Grid Mix Elec
f_x18 = fitdist(data_x18', 'Kernel', 'Support', 'Positive');       % Emission, NH4Cl
f_x19 = fitdist(data_x19', 'Kernel', 'Support', 'Positive');       % Emission, Monopotassium Phosphate
f_x20 = fitdist(data_x20', 'Kernel', 'Support', 'Positive');       % Emission, MgSO4-7H2O
f_x21 = fitdist(data_x21', 'Kernel', 'Support', 'Positive');       % Emission, CaCl2-2H2O
f_x22 = fitdist(data_x22', 'Kernel', 'Support', 'Positive');       % Emission, LDPE
f_x23 = fitdist(data_x23', 'Kernel', 'Support', [0, 1]);           % Emission, Process Water
f_x24 = fitdist(data_x24', 'Kernel', 'Support', [0, 1]);           % Emission, Grid Mix Elec
% [Sample from Kernel Distributions]
% Random Sample to Populate Parameter Space
ParSpace(:,1) = random(f_x1, n_sim, 1);
ParSpace(:,2) = random(f_x2, n_sim, 1);
ParSpace(:,3) = random(f_x3, n_sim, 1);
ParSpace(:,4) = random(f_x4, n_sim, 1);
ParSpace(:,5) = random(f_x5, n_sim, 1);
ParSpace(:,6) = random(f_x6, n_sim, 1);
ParSpace(:,7) = random(f_x7, n_sim, 1);
ParSpace(:,8) = random(f_x8, n_sim, 1);
ParSpace(:,9) = random(f_x9, n_sim, 1);
ParSpace(:,10) = random(f_x10, n_sim, 1);
ParSpace(:,11) = random(f_x11, n_sim, 1);
ParSpace(:,12) = random(f_x12, n_sim, 1);
ParSpace(:,13) = random(f_x13, n_sim, 1);
ParSpace(:,14) = random(f_x14, n_sim, 1);
ParSpace(:,15) = random(f_x15, n_sim, 1);
ParSpace(:,16) = random(f_x16, n_sim, 1);
ParSpace(:,17) = random(f_x17, n_sim, 1);
ParSpace(:,18) = random(f_x18, n_sim, 1);
ParSpace(:,19) = random(f_x19, n_sim, 1);
ParSpace(:,20) = random(f_x20, n_sim, 1);
ParSpace(:,21) = random(f_x21, n_sim, 1);
ParSpace(:,22) = random(f_x22, n_sim, 1);
ParSpace(:,23) = random(f_x23, n_sim, 1);
ParSpace(:,24) = random(f_x24, n_sim, 1);
% Random Sample to Populate Complementary Parameter Space
c_ParSpace(:,1) = random(f_x1, n_sim, 1);
c_ParSpace(:,2) = random(f_x2, n_sim, 1);
c_ParSpace(:,3) = random(f_x3, n_sim, 1);
c_ParSpace(:,4) = random(f_x4, n_sim, 1);
c_ParSpace(:,5) = random(f_x5, n_sim, 1);
c_ParSpace(:,6) = random(f_x6, n_sim, 1);
c_ParSpace(:,7) = random(f_x7, n_sim, 1);
c_ParSpace(:,8) = random(f_x8, n_sim, 1);
c_ParSpace(:,9) = random(f_x9, n_sim, 1);
c_ParSpace(:,10) = random(f_x10, n_sim, 1);
c_ParSpace(:,11) = random(f_x11, n_sim, 1);
c_ParSpace(:,12) = random(f_x12, n_sim, 1);
c_ParSpace(:,13) = random(f_x13, n_sim, 1);
c_ParSpace(:,14) = random(f_x14, n_sim, 1);
c_ParSpace(:,15) = random(f_x15, n_sim, 1);
c_ParSpace(:,16) = random(f_x16, n_sim, 1);
c_ParSpace(:,17) = random(f_x17, n_sim, 1);
c_ParSpace(:,18) = random(f_x18, n_sim, 1);
c_ParSpace(:,19) = random(f_x19, n_sim, 1);
c_ParSpace(:,20) = random(f_x20, n_sim, 1);
c_ParSpace(:,21) = random(f_x21, n_sim, 1);
c_ParSpace(:,22) = random(f_x22, n_sim, 1);
c_ParSpace(:,23) = random(f_x22, n_sim, 1);
c_ParSpace(:,24) = random(f_x22, n_sim, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Progress Bar%%
waitbar(100/total_steps, progress, 'Initializing Monte Carlo Simulations');
%% Prior Distribution Plots
Plot_String_Array = ["$\mathrm{{\mu}_{m}}$ maximum specific growth rate (1/min)",...
"Half light saturation constant ($\mathrm{{\mu}E/m^{2}s}$)",...
"Optimal growth temperature ($\mathrm{^{o}C}$)",...
"Temperature distribution parameter ($\mathrm{^{o}C}$)",...
"Specific cell death rate (1/min)",...
"Cell periodicity coefficient",...
"Dissolved air flotation biomass recovery",...
"Target biomass concentration in algal slurry (wt percent)",...
"Raw material price of $\mathrm{NH_{4}Cl}$, (USD/ton)",...
"Raw material price of $\mathrm{KH_{2}PO_{4}}$, (USD/ton)",...
"Raw material price of $\mathrm{K_{2}HPO_{4}}$, (USD/ton)",...
"Raw material price of $\mathrm{MgSO_{4}-7H_{2}O}$, (USD/ton)",...
"Raw material price of $\mathrm{CaCl_{2}-2H_{2}O}$, (USD/ton)",...
"Raw material price of $\mathrm{Na_{2}-EDTA}$, (USD/ton)",...
"Raw material price of LDPE, (USD/ton)",...
"Utility price of process water, (USD/$\mathrm{m^{3}}$)",...
"Utility price of electricity, (USD/$\mathrm{GJ_{elec}}$)",...
"Emission factor for $\mathrm{NH_{4}Cl}$ (ton-$\mathrm{NH_{4}Cl}$/ton)",...
"Emission factor for $\mathrm{KH_{2}PO_{4}}$ (ton-$\mathrm{KH_{2}PO_{4}}$/ton)",...
"Emission factor for $\mathrm{MgSO_{4}-7H_{2}O}$ (ton-$\mathrm{MgSO_{4}-7H_{2}O}$/ton)",...
"Emission factor for $\mathrm{CaCl_{2}-2H_{2}O}$ (ton-$\mathrm{CaCl_{2}-2H_{2}O}$/ton)",...
"Emission factor for LDPE (ton-LDPE/ton)",...
"Emission factor for process water production, Korea (ton-$\mathrm{CO_{2}eq}$/$\mathrm{m^{3}})$",...
"Emission factor for electricity production, Korea (ton-$\mathrm{CO_{2}eq}$/$\mathrm{GJ_{elec}}$)"];

% PROCESS: Maximum specific growth rate, Mu_m
figure(1)
[Y_1, X_1, BW1] = ksdensity(data_x1, 'npoints', np_2);
plot(X_1,Y_1,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(1), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param1.png')
% PROCESS: Half light saturation constant
figure(2)
[Y_2, X_2, BW2] = ksdensity(data_x2, 'npoints', np_2);
plot(X_2,Y_2,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(2), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param2.png')
% PROCESS: Optimal growth temperature
figure(3)
[Y_3, X_3, BW3] = ksdensity(data_x3, 'npoints', np_2);
plot(X_3,Y_3,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(3), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param3.png')
% PROCESS: Temperature distribution parameter
figure(4)
[Y_4, X_4, BW4] = ksdensity(data_x4, 'npoints', np_2);
plot(X_4,Y_4,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(4), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param4.png')
% PROCESS: Specific cell death rate
figure(5)
[Y_5, X_5, BW5] = ksdensity(data_x5, 'npoints', np_2);
plot(X_5,Y_5,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(5), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param5.png')
% PROCESS: Growth periodicity coefficient
figure(6)
[Y_6, X_6, BW6] = ksdensity(data_x6, 'npoints', np_2);
plot(X_6,Y_6,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(6), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param6.png')
% PROCESS: DAF Biomass Recovery %
figure(7)
[Y_7, X_7, BW7] = ksdensity(data_x7, 'npoints', np_2);
plot(X_7,Y_7,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(7), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param7.png')
% PROCESS: Target algal slurry biomass concentration, wt%
figure(8)
[Y_8, X_8, BW8] = ksdensity(data_x8, 'npoints', np_2);
plot(X_8,Y_8,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(8), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param8.png')
% COST: NH4Cl Price
figure(9)
[Y_9, X_9, BW9] = ksdensity(data_x9, 'npoints', np_2);
plot(X_9,Y_9,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(9), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param9.png')
% COST: KH2PO4 Price
figure(10)
[Y_10, X_10, BW10] = ksdensity(data_x10, 'npoints', np_2);
plot(X_10,Y_10,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(10), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param10.png')
% COST: K2HPO4 Price
figure(11)
[Y_11, X_11, BW11] = ksdensity(data_x11, 'npoints', np_2);
plot(X_11,Y_11,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(11), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param11.png')
% COST: MgSO4-7H2O
figure(12)
[Y_12, X_12, BW12] = ksdensity(data_x12, 'npoints', np_2);
plot(X_12,Y_12,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(12), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param12.png')
% COST: CaCl2-2H2O
figure(13)
[Y_13, X_13, BW13] = ksdensity(data_x13, 'npoints', np_2);
plot(X_13,Y_13,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(13), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param13.png')
% COST: Na2-EDTA
figure(14)
[Y_14, X_14, BW14] = ksdensity(data_x14, 'npoints', np_2);
plot(X_14,Y_14,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(14), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param14.png')
% COST: LDPE
figure(15)
[Y_15, X_15, BW15] = ksdensity(data_x15, 'npoints', np_2);
plot(X_15,Y_15,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(15), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param15.png')
% COST: Process water
figure(16)
[Y_16, X_16, BW16] = ksdensity(data_x16, 'npoints', np_2);
plot(X_16,Y_16,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(16), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param16.png')
% COST: Grid Mix Electricity
figure(17)
[Y_17, X_17, BW17] = ksdensity(data_x17, 'npoints', np_2);
plot(X_17,Y_17,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(17), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param17.png')
% EMISSION: NH4Cl
figure(18)
[Y_18, X_18, BW18] = ksdensity(data_x18, 'npoints', np_2);
plot(X_18,Y_18,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(18), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param18.png')
% EMISSION: K2HPO4-KH2PO4
figure(19)
[Y_19, X_19, BW19] = ksdensity(data_x19, 'npoints', np_2);
plot(X_19,Y_19,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(19), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param19.png')
% EMISSION: MgSO4-7H2O
figure(20)
[Y_20, X_20, BW20] = ksdensity(data_x20, 'npoints', np_2);
plot(X_20,Y_20,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(20), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param20.png')
% EMISSION: CaCl2-2H2O
figure(21)
[Y_21, X_21, BW21] = ksdensity(data_x21, 'npoints', np_2);
plot(X_21,Y_21,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(21), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param21.png')
% EMISSION: LDPE
figure(22)
[Y_22, X_22, BW22] = ksdensity(data_x22, 'npoints', np_2);
plot(X_22,Y_22,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(22), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param22.png')
% EMISSION: Proces water
figure(23)
[Y_23, X_23, BW23] = ksdensity(data_x23, 'npoints', np_2);
plot(X_23,Y_23,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(23), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param23.png')
% EMISSION: Grid Mix electricity
figure(24)
[Y_24, X_24, BW24] = ksdensity(data_x24, 'npoints', np_2);
plot(X_22,Y_22,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(24), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_MA_param24.png')




%% [3] Process Parameter Space for System Evaluation
% (a) Determine the # of System Model Outputs (# of Evaluation Metrics)
n_out = size(MODEL_MA(kp),2);
% (b) Create containers for Evaluated Model Outputs
% NOTE: eval_P and eval_C are n_sim x n_out x n_kp matrix. 
fx = zeros(n_sim, n_out);  % Evaluated outputs with all Params from ParSpace
fx_P = zeros(n_sim, n_out);% Evaluated outputs with i from ParSpace, ~i from c_ParSpace
fx_C = zeros(n_sim, n_out);% Evaluated outputs with i from c_ParSpacek, ~i from ParSpace
% (c) Evaluate Model from Monte Carlo Sampled Inputs (ParSpace)
parfor i = 1:n_sim
    % Each Parameter Set is a Row in the ParSpace matrix
    Parameter_Set = ParSpace(i,:); 
    fx(i,:) = MODEL_MA(Parameter_Set);
end
% (d) Generate Function Output Space based on i and ~i
for i = 1:n_sim
    for j = 1:n_kp
        %%%% fx_P = f(x_ik, x'_~ik)
        kp = [c_ParSpace(i,1:j-1), ParSpace(i,j), c_ParSpace(i,j+1:n_kp)];
        fx_P(i,:,j) = MODEL_MA(kp);
        %%%% fx_C = f(x'_ik, x_~ik)
        kp = [ParSpace(i,1:j-1), c_ParSpace(i,j), ParSpace(i,j+1:n_kp)];
        fx_C(i,:,j) = MODEL_MA(kp);
    end
    %%Progress Bar%%
    waitbar((100 + ((i/n_sim)*520))/total_steps, progress,...
    sprintf('Processing %d of %d Simulations', i, n_sim));
end
%%Progress Bar%%
waitbar(580/total_steps, progress, 'Computing Sobol Variances');




%% [4-SOBOL] Computing Variances
% (a) Initialize Integrals and Total Variances
% Initialize the Integral of Model Outputs (f0^2) 
f0 = zeros(1, n_out);            
% Initialize Total Variance of Model Outputs 
D = zeros(1, n_out);            
% (b) Compute the Average of Model Outputs, f0
parfor i = 1:n_out 
    for j = 1:n_sim
        f0(i) = f0(i) + fx(j,i);
    end
    f0(i) = f0(i)/n_sim;
end
% (c) Estimating Total Output Variance, D, using MC Model Outputs
parfor i = 1:n_out
    for j = 1:n_sim
        D(i) = D(i) + ((fx(j,i)^2)-(f0(i)^2));
    end
    D(i) = D(i)/n_sim;
end
% (d) Compute Partial Variances (1st Order Sobol)
%D_1st = zeros(n_kp, n_out);
%F_Par = zeros(n_kp, n_out);           % Initialize Partial Variance Factors
%parfor i = 1:n_kp
%    for j = 1:n_out
%        for k = 1:n_sim
%            F_Par(i,j) = F_Par(i,j) + (fx(k,j) - fx_P(k,j,i))^2;
%        end
%        D_1st(i,j) = D(j) - (F_Par(i,j)/(n_sim*2));
%    end
%end
%(e) Compute Total Variances (Total Sobol)
D_Tot = zeros(n_kp, n_out);
F_cPar = zeros(n_kp, n_out);          % Initialize Total Sobol Factors
parfor i = 1:n_kp
    for j = 1:n_out
        for k = 1:n_sim
            F_cPar(i,j) = F_cPar(i,j) + (fx(k,j) - fx_C(k,j,i))^2;
        end
        D_Tot(i,j) = F_cPar(i,j)/(n_sim*2);
    end
end
%%Progress Bar%%
waitbar(630/total_steps, progress, 'Calculating Sobol Indices for each System Output');




%% [4-SOBOL] Determine Sobol Indices
% (a) Compute Sum of Partial and Total Variances
%Sum_Partial = zeros(n_out);
Sum_Total = zeros(n_out);
for i = 1:n_out
    %Sum_Partial(i) = sum(D_1st(:,i));
    Sum_Total(i) = sum(D_Tot(:,i));
end
% (b) Compute Rank Scores based on Partial/Total Variances
%Sobol_1st = zeros(n_kp, n_out);
Sobol_Total = zeros(n_kp, n_out);
for i = 1:n_out
    for j = 1:n_kp
        %Sobol_1st(j,i) = D_1st(j,i)/D(i);
        Sobol_Total(j,i) = D_Tot(j,i)/D(i);
    end
end
% (c) Rank Input Parameters via Sobol Indices (Total)
for i = 1:n_out
    [score2, rank2] = sort(Sobol_Total(:,i), 'descend');
    fprintf('Total Sobol Ranks for Model Evaluation Output %.0f \n', i)
    sprintf('The Total Sobol Indices for the Above Order are: ')
    Sobol_Total(:,i)
    fprintf('Sum of Total Sobol Indices for Model Evaluation Output %.0f \n', i)
    sum(Sobol_Total(:,i))
end
waitbar(640/total_steps, progress, 'Registering Stakeholders Criteria');




%% [4-SOBOL] Normalized Total Sobol Indices
% (a) Define a vector of weights for each model output
%for i = 1:n_out
%    weight_prompt = sprintf('Enter Stakeholder Weight for Sustainability Criteria %.0f \n', i);
%    weight_prompt_title = 'Enter Stakeholder Weight';
%    dims1 = [1 35];
%    weight_answer = inputdlg(weight_prompt, weight_prompt_title, dims1);
%    weight(i) = str2double(weight_answer{1});
%end
% (b) Check if the sum of the weight vector = 1
%if sum(weight)~= 1
%    error('The sum of the metric weights must equal to 1')
%end
%if size(weight,2)~= n_out
%    error('A numerical weight must be assigned to each system evaluation output. For no weights, enter 0')
%end
% (c) Compute Normalized Total Sobol Indices
%NTS = zeros(1, n_kp);
%for i = 1:n_kp
%    NTS(i) = sum(weight.*Sobol_Total(i,:));
%end
%NTS_Sum = sum(NTS);
%waitbar(650/total_steps, progress, 'Classifying System Model Outputs');




%% [5-GNB-k] Stakeholder's Criteria for Sustainability per Indicator
% (a) Define Stakeholder's Target Criteria for Sustainability
targets = zeros(1, n_out);  % Array of Target Bayesian Output Metric Values
if bayes_auto == 1
    for i = 1:n_out
        targets(i) = mean(fx(:,i));
    end
else
    for i = 1:n_out
        prompt = sprintf('Enter Stakeholder Decision Criteria for Output Metric %.0f \n', i);
        prompt_title = 'Enter Stakeholder Decision Criteria Value';
        dims = [1 70];
        answer = inputdlg(prompt, prompt_title, dims);
        targets(i) = str2double(answer{1});
    end
end             
% (b) Generate Bayesian Classification Matrix
classmat = zeros(n_sim, n_out);
Success_Out = zeros(n_sim, n_out);
Failure_Out = zeros(n_sim, n_out);
parfor i = 1:n_sim
    for j = 1:n_out
        if fx(i,j) <= targets(j)
            classmat(i,j) = 1;
            Success_Out(i,j) = fx(i,j);
        else
            classmat(i,j) = 0;
            Failure_Out(i,j) = fx(i,j);
        end
    end
end
% (c) Plot the Overall Output Uncertainty Distributions 
Output_Prob = zeros(n_out,1);
for i = 1:n_out
    % Plot the Overall Posterior Distribution for Output
    Out_Success_K = nonzeros(Success_Out(:,i))';
    Out_Failure_K = nonzeros(Failure_Out(:,i))';
    [Out_F_s, Out_x_s] = ksdensity(Out_Success_K, 'npoints', np_2);
    [Out_F_f, Out_x_f] = ksdensity(Out_Failure_K, 'npoints', np_2);
    figure(i)
    plot(Out_x_s, Out_F_s, 'b','LineWidth',2);
    hold on
    plot(Out_x_f, Out_F_f,'r','LineWidth',2);
    title(['Posterior distribution for output ', num2str(i)])
    legend('Success Posterior', 'Failure Posterior')
    hold off
    saveas(gcf, ['OUTPUT_POSTERIOR_criteria_',num2str(i),'.png'])
    % Compute overall probability of achieving sustainability criteria
    Output_Prob(i) = length(Out_Success_K)/n_sim;
    fprintf('Probability of Achieving Stakeholders Criteria %.0f \n',i)
    Output_Prob(i)
end
%%Progress Bar%%
waitbar(700/total_steps, progress, 'Classifying Parameter Inputs based on Criteria');




%% [5-GNB-k] Categorization of Input Values based on Classification
% (a) Generate Empty Factorized Containers. The maximum dimensions of each
%     Label matrix is n_sim x n_kp x n_out
Success_Mat = [];
Failure_Mat = [];
% (b) Loop through "classmat" then append values from ParamSpace to 
%     Success_Mat or Failure_Mat accordingly
for i = 1:n_sim
    for j = 1:n_out
        if classmat(i,j) == 1
            Success_Mat(end+1,:,j) = ParSpace(i,:);
        else
            Failure_Mat(end+1,:,j) = ParSpace(i,:);
        end
    end
end
%%Progress Bar%%
waitbar(750/total_steps, progress, 'Generating Likelihoods based on Kernel Data');




%% [5-GNB-k] Likelihood Distributions using Kernel Density Estimation
% (a) Initialize Containers
S_KernelData = zeros(n_sim, n_kp, n_out); % Sampled matrix for "success" 
F_KernelData = zeros(n_sim, n_kp, n_out); % Sampled matrix for "failure" 
Sus_Risk_Score = zeros(n_kp, n_out);    % A n_kp x n_out matrix of ranks
Posterior_Diff = zeros(1,np_2);           
X_LIST = zeros(1,np_2);                   % Normalized X-Axis array
% (b) Populate Stakeholders Risk Score by evaluating the differences
for j = 1:n_kp
    figure()
    set(gcf, 'Position',  [100, 100, 1600, 400])
    for i = 1:n_out
        % Generate array of sample data for Bayesian Kernel Distributions
        S_KernelData = nonzeros(Success_Mat(:,j,i))' ;   
        F_KernelData = nonzeros(Failure_Mat(:,j,i))';
        % Generate the Kernel Distribution for each Bayesian Outcome
        [F_s, x_s] = ksdensity(S_KernelData, 'npoints', np_2);
        [F_f, x_f] = ksdensity(F_KernelData, 'npoints', np_2);
        % Because Success/Failure are 2 data sets, need to interpolate on a
        % consistent X-Axis. Generate the consistent X-Axis array
        X_APPEND = [x_s, x_f];
        X_MIN = min(X_APPEND);
        X_MAX = max(X_APPEND);
        X_LIST = linspace(X_MIN,X_MAX,np_2);
        % Interpolate the Success/Failure probabilities with above array
        Y_s = interp1(x_s, F_s, X_LIST);
        Y_f = interp1(x_f, F_f, X_LIST);
        % Calculate the Posterior Probability differences as Risk Scores
        Posterior_Diff = abs(Y_s-Y_f);
        Posterior_Diff(isnan(Posterior_Diff))=0;
        Sus_Risk_Score(j,i)  = trapz(X_LIST, Posterior_Diff);
        % Plot the Posterior Probability Differences
        subplot(1,n_out,i)
        plot(X_LIST,Y_s,'k-','LineWidth',2.5);
        hold on
        plot(X_LIST,Y_f,'k:','LineWidth',2.5);
        % Draw and color the posterior probability difference
        ProbDiff = area(X_LIST, abs(Y_s-Y_f));
        ProbDiff.FaceColor = [0.8 0.15 0.1];
        ProbDiff.EdgeColor = [0.8 0.15 0.1];
        ProbDiff.FaceAlpha = 0.38;
        ProbDiff.EdgeAlpha = 0.38;
        ProbDiff.LineWidth = 1;
        % Legend and Axes
        ax = gca;
        ax.FontSize = 13;
        legend(['$P(x_{' num2str(j) '}|L_{s})$'], ['$P(x_{' num2str(j) '}|L_{f})$'], ['$|P(x_{' num2str(j) '}|L_{s}) - P(x_{' num2str(j) '}|L_{f})|$'], 'FontUnits', 'points', 'interpreter', 'latex', 'FontSize', 12, 'FontName', 'Times', 'Location', 'NorthEast', 'box', 'off')
        xlabel(Plot_String_Array(j), 'interpreter', 'latex', 'FontSize', 14)
        ylabel(["Density"], 'interpreter', 'latex', 'FontSize', 14)
        if i == 1
            title('Unit COGM criteria', 'interpreter', 'latex', 'FontSize', 14);
        elseif i == 2
            title('Specific carbon footprint criteria', 'interpreter', 'latex', 'FontSize', 14);
        end
        hold off
    end
    saveas(gcf, ['LIKELIHOOD_MA_output_',num2str(i),'_parameter_',num2str(j),'.png'])
end
for i = 1:n_out
    fprintf('The Sustainability Risk Score for Sustainability Criteria %.0f \n',i)
    Sus_Risk_Score(:,i)
end
%%Progress Bar%%
waitbar(990/total_steps, progress, 'Plotting Distributions of System Outputs');




%% [6-RESULTS] Plot and Export Results
% (a) Plot Distributions of Model Outputs
%nbins = 50;             % Define resolution of output metric's historgram
% (b) Generate plots (Example for via For loop below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nbins = 100;
    figure('Name', 'Variance of Cost of Goods Manufactured')
    histogram(fx(:,1), nbins, 'facecolor', [0, 0, 0])
    ax = gca;
    ax.FontSize = 13;
    xlabel('Unit COGM, USD/ton Formic Acid', 'interpreter', 'latex', 'FontSize', 14)
    ylabel('Frequency', 'interpreter', 'latex', 'FontSize', 14)
    saveas(gcf, 'COGM_Distribution.png')
    
    figure('Name', 'Variance of Carbon Footprint')
    histogram(fx(:,2), nbins, 'facecolor', [0, 0, 0])
    ax = gca;
    ax.FontSize = 13;
    xlabel('Specific carbon footprint, ton$\mathrm{CO_{2-eq}}$/ton Formic Acid', 'interpreter', 'latex', 'FontSize', 14)
    ylabel('Frequency', 'interpreter', 'latex', 'FontSize', 14)
    saveas(gcf, 'CarbFoot_Distribution.png')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Save Results
% Define Vars
Sort_COGM = sort(fx(:,1));
Sort_GWI = sort(fx(:,2));
INTQUART = zeros(2,n_out);
TAILS = zeros(2,n_out);
VAR = zeros(1,n_out);
CVAR = zeros(1,n_out);
% Inter-quartile Range
INTQUART(1,1) = Sort_COGM(n_sim*0.25);
INTQUART(2,1) = Sort_GWI(n_sim*0.25);
INTQUART(1,2) = Sort_COGM(n_sim*0.75);
INTQUART(2,2) = Sort_GWI(n_sim*0.75);
% 5% - 95% range
TAILS(1,1) = Sort_COGM(n_sim*0.05);
TAILS(2,1) = Sort_GWI(n_sim*0.05);
TAILS(1,2) = Sort_COGM(n_sim*0.95);
TAILS(2,2) = Sort_GWI(n_sim*0.95);
% Value-at-Risk for alpha=0.05
VAR(1) = Sort_COGM((1-alpha)*n_sim);
VAR(2) = Sort_GWI((1-alpha)*n_sim);
% Conditional Value-at-Risk for alpha = 0.05
CVAR(1) = mean(Sort_COGM((1-alpha)*n_sim:n_sim));
CVAR(2) = mean(Sort_GWI((1-alpha)*n_sim:n_sim));
%%Progress Bar%%
waitbar(1000/total_steps, progress, 'Complete');
delete(progress)
toc




%% EXPORT RESULTS
RESULTS = struct();
RESULTS.Input = ParSpace;
RESULTS.SobolTot = Sobol_Total;
RESULTS.RiskScore = Sus_Risk_Score;
RESULTS.Output = fx;
RESULTS.IQ = INTQUART;
RESULTS.TAILS = TAILS;
RESULTS.VAR = VAR;
RESULTS.CVAR = CVAR;
% Export to CSV
MA_Res = struct2cell(RESULTS);
writecell( MA_Res, 'MA_Results_Separate.csv');