%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% STAKEHOLDER'S RISK ATTRIBUTION & UNCERTAINTY ANALYSIS %%%%%%%%%%
%%%%%%%%%%%%%         CO2-BASED FORMIC ACID PROCESS           %%%%%%%%%%%%%
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
n_sim = 30000;
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
kp = [1 1 1 1 1 1 1 1 1 1 1 1 1 1];      
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Define Parametric Distributions for Key Params]
% [List Kernel Data Points]
data_x1 = [0.638, 0.669, 0.411, 0.656, 0.333, 0.583, 0.545, 0.448, 0.569, 0.546, 0.501, 0.436, 0.480, 0.457, 0.429, 0.368, 0.479];% [GLOBAL] CO2 hyrdrogenation conversion, %
data_x2 = [0.08, 0.1, 0.125];                                                                % [GLOBAL] TREA_loss_rate, %
data_x3 = [0.08, 0.1, 0.125];                                                                 % [GLOBAL] NBIM_loss_rate, %
data_x4 = [13.41, 13.41, 15.08, 25.89, 18.73, 25.89, 45.27, 25.96, 39.3, 25.49, 18.35, 25.49]; % [KOREA] Elec. cost, USD/GJ
data_x5 = [17.81, 16.63, 16.69, 19.03, 17.85, 17.91, 18.92, 17.75, 17.81];                     % [KOREA] MP steam cost, USD/GJ
data_x6 = [268.63, 250.84, 251.87, 287.04, 269.24, 270.27, 285.52, 267.73, 268.76];            % [KOREA] NG cost, USD/ton
data_x7 = [40, 60, 80];                                                                        % [GLOBAL] Captured CO2 (Power Plant) cost, USD/ton
data_x8 = [2544.3, 1500, 2200, 2150, 2390];                                                    % [GLOBAL] TREA cost, USD/ton
data_x9 = [8000, 4000, 10000];                                                                 % [GLOBAL] NBIM cost, USD/ton
data_x10 = [-0.8, -0.86, -0.8, -0.67];                                                         % [KOREA] GWI of Captured CO2 from Power Plant, tonCO2eq/ton
data_x11 = [0.185, 0.154, 0.129, 0.1375];                                                              % [KOREA] GWI of Grid Mix Electricity, tonCO2eq/GJ
data_x12 = [0.0691, 0.0889, 0.0597];                                                           % [KOREA] GWI of Industrial Steam, tonCO2eq/ton
data_x13 = [0.301, 0.614, 0.102];                                                              % [KOREA] GWI of NG, tonCO2eq/ton
data_x14 = [0.000173, 0.00019895, 0.000123];                                                 % [KOREA] GWI of Process Water, tonCO2eq/ton


% [Define Kernel Distributions]
f_x1 = fitdist(data_x1', 'Kernel', 'Support', [0.01,1]);
f_x2 = fitdist(data_x2', 'Kernel', 'Support', [0, 1]);
f_x3 = fitdist(data_x3', 'Kernel', 'Support', [0, 1]);
f_x4 = fitdist(data_x4', 'Kernel', 'Support', [0, 100]);
f_x5 = fitdist(data_x5', 'Kernel', 'Support', [0, 100]);
f_x6 = fitdist(data_x6', 'Kernel', 'Support', [0, 1000]);
f_x7 = fitdist(data_x7', 'Kernel', 'Support', [0, 500]);
f_x8 = fitdist(data_x8', 'Kernel', 'Support', [100, 50000]);
f_x9 = fitdist(data_x9', 'Kernel', 'Support', [100, 100000]);
f_x10 = fitdist(data_x10', 'Kernel', 'Support', [-1, 0]);
f_x11 = fitdist(data_x11', 'Kernel', 'Support', [0, 10]);
f_x12 = fitdist(data_x12', 'Kernel', 'Support', [0, 10]);
f_x13 = fitdist(data_x13', 'Kernel', 'Support', [0, 10]);
f_x14 = fitdist(data_x14', 'Kernel', 'Support', [0, 1]);


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Progress Bar%%
waitbar(100/total_steps, progress, 'Initializing Monte Carlo Simulations');
%% Prior Distribution Plots
Plot_String_Array = ["$\mathrm{CO_{2}}$ conversion in hydrogenation",...
"Triethylamine fractional degradation",...
"$n$-butyl imidazole fractional degradation",...
"Price of Korean grid mix electricity (USD/$\mathrm{GJ_{elec}}$)",...
"Price of Korean industrial steam (USD/GJ)",...
"Price of Korean natural gas (USD/ton)",...
"Price of captured $\mathrm{CO_{2}}$ from power plants (USD/ton)",...
"Price of Triethylamine (USD/ton)",...
"Price of $n$-butyl imidazole (USD/ton)",...
"Emission factor for captured $\mathrm{CO_{2}}$ from power plants (ton-$\mathrm{CO_{2}eq}$/ton)",...
"Emission factor for electricity production, Korea (ton-$\mathrm{CO_{2}eq}$/$\mathrm{GJ_{elec}}$)",...
"Emission factor for industrial steam production, Korea (ton-$\mathrm{CO_{2}eq}$/GJ)",...
"Emission factor for Korean natural gas (ton-$\mathrm{CO_{2}eq}$/ton)",...
"Emission factor for process water production, Korea (ton-$\mathrm{CO_{2}eq}$/$\mathrm{m^{3}})$"];

% PROCESS: Hydrogenation Reactor Conversion
figure(1)
[Y_1, X_1, BW1] = ksdensity(data_x1, 'npoints', np_2);
plot(X_1,Y_1,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(1), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param1.png')
% PROCESS: TREA Recovery Rate
figure(2)
[Y_2, X_2, BW2] = ksdensity(data_x2, 'npoints', np_2);
plot(X_2,Y_2,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(2), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param2.png')
% PROCESS: NBIM Recovery Rate
figure(3)
[Y_3, X_3, BW3] = ksdensity(data_x3, 'npoints', np_2);
plot(X_3,Y_3,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(3), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param3.png')
% COST: Electricity Price
figure(4)
[Y_4, X_4, BW4] = ksdensity(data_x4, 'npoints', np_2);
plot(X_4,Y_4,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(4), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param4.png')
% COST: Steam (MP/HP)
figure(5)
[Y_5, X_5, BW5] = ksdensity(data_x5, 'npoints', np_2);
plot(X_5,Y_5,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(5), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param5.png')
% COST: Natural Gas Heat
figure(6)
[Y_6, X_6, BW6] = ksdensity(data_x6, 'npoints', np_2);
plot(X_6,Y_6,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(6), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param6.png')
% COST: Captured CO2
figure(7)
[Y_7, X_7, BW7] = ksdensity(data_x7, 'npoints', np_2);
plot(X_7,Y_7,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(7), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param7.png')
% COST: TREA
figure(8)
[Y_8, X_8, BW8] = ksdensity(data_x8, 'npoints', np_2);
plot(X_8,Y_8,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(8), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param8.png')
% COST: NBIM
figure(9)
[Y_9, X_9, BW9] = ksdensity(data_x9, 'npoints', np_2);
plot(X_9,Y_9,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(9), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param9.png')
% GWI: Captured CO2 Consumption
figure(10)
[Y_10, X_10, BW10] = ksdensity(data_x10, 'npoints', np_2);
plot(X_10,Y_10,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(10), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param10.png')
% GWI: Electricity
figure(11)
[Y_11, X_11, BW11] = ksdensity(data_x11, 'npoints', np_2);
plot(X_11,Y_11,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(11), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param11.png')
% GWI: Steam (MP/HP)
figure(12)
[Y_12, X_12, BW12] = ksdensity(data_x12, 'npoints', np_2);
plot(X_12,Y_12,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(12), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param12.png')
% GWI: Natural Gas Heat
figure(13)
[Y_13, X_13, BW13] = ksdensity(data_x13, 'npoints', np_2);
plot(X_13,Y_13,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
xlabel(Plot_String_Array(13), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param13.png')
% GWI: Water
figure(14)
[Y_14, X_14, BW14] = ksdensity(data_x14, 'npoints', np_2);
plot(X_14,Y_14,'k-','LineWidth',2);
ax = gca;
ax.FontSize = 13;
ax.XAxis.Exponent = 0;
xlabel(Plot_String_Array(14), 'interpreter', 'latex', 'FontSize', 14)
set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf, 'PRIOR_FA_param14.png')




%% [3] Process Parameter Space for System Evaluation
% (a) Determine the # of System Model Outputs (# of Evaluation Metrics)
n_out = size(MODEL_FA(kp),2);
% (b) Create containers for Evaluated Model Outputs
% NOTE: eval_P and eval_C are n_sim x n_out x n_kp matrix. 
fx = zeros(n_sim, n_out);  % Evaluated outputs with all Params from ParSpace
fx_P = zeros(n_sim, n_out);% Evaluated outputs with i from ParSpace, ~i from c_ParSpace
fx_C = zeros(n_sim, n_out);% Evaluated outputs with i from c_ParSpacek, ~i from ParSpace
% (c) Evaluate Model from Monte Carlo Sampled Inputs (ParSpace)
parfor i = 1:n_sim
    % Each Parameter Set is a Row in the ParSpace matrix
    Parameter_Set = ParSpace(i,:); 
    fx(i,:) = MODEL_FA(Parameter_Set);
end
% (d) Generate Function Output Space based on i and ~i
for i = 1:n_sim
    for j = 1:n_kp
        %%%% fx_P = f(x_ik, x'_~ik)
        kp = [c_ParSpace(i,1:j-1), ParSpace(i,j), c_ParSpace(i,j+1:n_kp)];
        fx_P(i,:,j) = MODEL_FA(kp);
        %%%% fx_C = f(x'_ik, x_~ik)
        kp = [ParSpace(i,1:j-1), c_ParSpace(i,j), ParSpace(i,j+1:n_kp)];
        fx_C(i,:,j) = MODEL_FA(kp);
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



min(fx(:,1))
max(fx(:,1))
min(fx(:,2))
max(fx(:,2))
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
    saveas(gcf, ['OUTPUT_LIKELIHOOD_criteria_',num2str(i),'.png'])
    % Compute overall probability of achieving sustainability criteria
    Output_Prob(i) = length(Out_Success_K)/n_sim;
    fprintf('Probability of Achieving Stakeholders Criteria %.0f \n',i)
    Output_Prob(i)
end
%%Progress Bar%%
waitbar(700/total_steps, progress, 'Classifying Parameter Inputs based on Criteria');




%% [5-GNB-k] Categorization of Input Values based on Classification
% (a) Generate Empty Factorized Containers
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
        legend(['$P(x_{' num2str(j) '}|L_{s})$'], ['$P(x_{' num2str(j) '}|L_{f})$'], ['$|P(x_{' num2str(j) '}|L_{s}) - P(x_{' num2str(j) '}|L_{f})|$'], 'FontUnits', 'points', 'interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times', 'Location', 'NorthEast', 'box', 'off')
        xlabel(Plot_String_Array(j), 'interpreter', 'latex', 'FontSize', 14)
        ylabel(["Density"], 'interpreter', 'latex', 'FontSize', 14)
        if i == 1
            title('Unit COGM criteria', 'interpreter', 'latex', 'FontSize', 14);
        elseif i == 2
            title('Specific carbon footprint criteria', 'interpreter', 'latex', 'FontSize', 14);
        end
        hold off
    end
    saveas(gcf, ['POSTERIOR_FA_output_',num2str(i),'_parameter_',num2str(j),'.png'])
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
FA_Res = struct2cell(RESULTS);
writecell( FA_Res, 'FA_Results_Separate.csv');