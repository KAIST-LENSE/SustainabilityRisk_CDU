%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Process + Evaluation Model for Microalgae Production as BPFS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EvalMetrics = MODEL_MA(kp)             
format long g
%% SCENARIOS
% GG = DYCOPS_BioAdditive([0.0035, 384.5, 27.9, 9.83, 0.0011, 187.62, 0.95, 0.25, 115, 1150, 1150, 140, 115, 2575, 1200, 0.68, 25, 1.319, 2.69, 0.3, 0.89, 1.85, 0.000173, 0.154])

%% KEY PARAMETERS WITH UNCERTAINTY
% [1] Cultivation Process
Mu_max = kp(1);       % Maximum specific growth rate, 0.0035 min^-1
K_I = kp(2);          % Half-light saturation constant, 384.5 ?E/m2s
T_opt = kp(3);        % Optimal growth temperature, 27.9 degC
Theta = kp(4);        % Temperature distribution parameter, 9.83 degC
k_D = kp(5);          % Cell death, 0.0011 min^-1
C_D = kp(6);          % Light dependence Coefficient, 187.62

% [2] Harvest+Dewater Process
rec_CENTR = kp(7);    % Biomass recovery during centrifuge harvest, 95%
conc_CENTR = kp(8);   % Exit conc. of biomass @ Centrifuge, 25 dwt%

% [3] Blowdown + Media Recycle
% [4] Convective Flue Gas Drying

% [TEA] Techno-Economic Analysis Parameters
%%% [EQUIPMENT] Equipment Sizing Costs
%%% [OPERATING] Plant Costs
%%% [OPERATING] Raw Material Costs
OP_RM_NH4Cl = kp(9);   % $/ton, NH4Cl, 115
OP_RM_KH2PO4 = kp(10); % $/ton, KH2PO4, 1150
OP_RM_K2HPO4 = kp(11); % $/ton, K2HPO4, same as above
OP_RM_MgSO4 = kp(12);  % $/ton, MgSO4, 140
OP_RM_CaCl2 = kp(13);  % $/ton, CaCl2, 115
OP_RM_NaEDTA = kp(14); % $/ton, Na2-EDTA, 2575
OP_RM_LDPE = kp(15);   % $/ton, LDPE Film, 1200
OP_RM_WATER = kp(16);  % $/ton, 0.68 from IPON Quote
%%% [OPERATING] Utility Costs
UT_ELEC = kp(17);      % $/GJ, 25 from IPON Quote

% [LCA] CO2 Life Cycle Assessment Parameters
%%% [MATERIAL] Global Warming Potentials for Material Consumption
GWI_NH4Cl = kp(18);    % ton CO2-eq/ton NH4Cl, 1.32
GWI_KH2PO4 = kp(19);   % ton CO2-eq/ton KH2PO4/K2HPO4, 2.69
GWI_MgSO4 = kp(20);    % ton CO2-eq/ton MgSO4-7H2O, 0.3
GWI_CaCl2 = kp(21);    % ton CO2-eq/ton CaCl2-2H2O, 0.89
GWI_LDPE = kp(22);     % ton CO2-eq/ton LDPE, 1.85
GWI_WATER = kp(23);    % ton CO2-eq/ton Water, 0.000173
%%% [UTILITY] Global Warming Potentials for Utility Consumption
GWI_ELEC_GRID = kp(24);% ton CO2-eq/GJ Elec, 0.154




%% NON VARYING FIXED-PARAMETERS & DESIGN SPECS
% [1] Cultivation Process
Cx0 = 0.0262;          % Average innoculum concentration in g/L
HAR_T = 150;           % Cultivation time in hours, 150 hrs
STOCK_T = 0.5;         % Time to transfer stock soln. to GM Tank, hrs
GM_T = 1.0;            % Time to transfer growth media to Cult., hrs
REC_T = 1.0;           % Time to transfer cult. broth back to GM Tank, hrs
V_Cult_T = 15000000;   % Total cultivation volume, liters  (Approx 100 ton/week)
V_Cult_Unit = 1000;    % Liters per unit Cultivation Module
Land_Cult_Unit = 5;    % 3.3 m^2 per 1000L PBR + 1.7 m^2 Settling tank 
VVM = 0.05;            % Volume fluegas per volume
Daylight = 0.5;        % Percent of day photosynthesis occurs
FG_MM = [44.01, 18.02, 32, 28.01];           % NGCC, W. Zhang et al., International J. Greenhouse Gas Control, 2016
% FG_Moles = [0.1383, 0.0832, 0.0006, 0.036, 0.7419];  % For COAL [CO2, Water, H2S, NO2, O2, N2]
FG_Moles = [0.041, 0.079, 0.121, 0.759];

% [2] Harvest+Dewater Process
SEP_T = 10;            % Operation time/batch for culture separation, including DAF, hrs

% [3] Blowdown + Media Recycle
BD_LOSS = [0.01, 0.01, 0, 0, 0, 0, 0, 0];         % [Water, Biomass, NH4Cl, KH2PO4, K2HPO4, MgSO4, CaCl2, Na2-EDTA]

% [4] Convective Flue Gas Dryer
Target_Moisture = 0.1;                   % 10% moisture content [A. Giostri paper]
FG_Comp = [0.079, 0.759, 0.041, 0.121];    % [H2O, N2, CO2, O2] NGCC
FG_MW = [18.02, 28.014, 44.01, 32];        % [H2O, N2, CO2, O2, H2S] NGCC
FG_Temp_in = 106;                        % degC [A. Giostri paper]
FG_Temp_out = 45;                        % degC [A. Giostri paper]
FG_Pressure = 825.068;                   % mmHg
Antoine_Coeff_Water = [8.14, 1810.94, 244.485];     % Water from 99 - 374 degC
CP_FG = 1.054;                           % kJ/kgC
CP_H2O = 4.187;                          % kJ/kgC
H_lat = 2260;                            % kJ/kg
DRY_T_m = 242.7;                         % y = mx+b where y = Drying time (mins) and x = moisture content (%) [H. Hosseinizand et al., 2017]
DRY_T_b = 83.015;                        % y = mx+b where y = Drying time (mins) and x = moisture content (%) [H. Hosseinizand et al., 2017]
AREA_LOAD = 30;                          % kg/m2. Kiranoudia and Markatos (2000) suggested a maximum unit loading of 50 kg/m2 (wet)

% [TEA] Techno-Economic Analysis Parameters
%%% Chemical Engineering Plant Cost Indices (CEPCI)
CEPCI = 619.2;                        % Based on CEPCI for 2019
CEPCI_PBR = 619.2;                    % 2019 CEPCI based on [21] publish yr
CEPCI_PumpsBlowers = 541.3;           % 2016 CEPCI based on [15] publish yr
CEPCI_Separator = 541.3;              % 2016 CEPCI based on [15] publish yr
CEPCI_Upstr_Tanks = 603.1;            % 2018 CEPCI from [22]
%%% Equipment Efficiencies
NU_Pump = 0.73;                       % Efficiency for Water & Solutions
NU_Motor = 0.90;                      % Motor efficiency for 50-250kW Range
NU_HydFan = 0.70;                     % Hydraulic fan efficiency for blower, [A. Giostri et al., 2016]
NU_MechFan = 0.94;                    % Mechanical-electric fan efficiency for blower, [A. Giostri et al., 2016] 
%%% Capital Cost Lang Factors
OSBL_OS = 0.25;                       % Offsite Operations, from [23]
OSBL_DE = 0.2;                        % Design and Engineering, from [23]
OSBL_CN = 0.3;                        % Contingency, from [23]
LIFET = 30;                           % Plant Lifetime, years
DISCO = 0.05;                         % Annual Discount Rate % of TCI
%%% Equipment Costs
CAP_SS_T = 500;                       % m^3 Storage Capacity for SS Tank
CAP_GM_T = 10000;                     % m^3 Storage Capacity for GM Tank
C_PBR = 81867;                        % USD/acre Disregarding greenhouse
C_FG_Blower = 5803;                   % 7.5kW 3-lobe Blower, [15], 5803 USD/unit
C_Stock_Pump = 1826;                  % 5.5kW Centrifugal Pump, [15], 1826 USD/unit
C_Rec_Pump = C_Stock_Pump;            % 5.5kW Centrifugal Pump, [15], 1826 USD/unit
C_GM_Pump = 20234;                    % 10kW Submersible Pump, [15], 20234 USD/unit
C_SEP = 62591;                        % 7.5kW Centrifugal Separator, 62591 USD/unit
C_SS_T = 54159;                       % Cost per Tank, Carbon Steel, 54159 USD/unit
C_GM_T = 478113;                      % Cost per Tank, Carbon Steel, 478113 USD/unit
%%% [OPERATING COST] Raw Material Costs
OP_RM_FG = 0;                         % $/kg of Flue Gas
OP_PE_MAINT = 1090;                   % $/Acre from [21]. Land maintainence + Clean-in-Place costs, 340+750

% [LCA] CO2 Life Cycle Assessment Parameters
%%% [BOUNDARY] Binary Parameters to Set LCA Boundary for Plant
B_DIR = 1;    % Include direct plant emissions? Default=1
B_ENR = 1;    % '' indirect effects from energy consumption? Default=1
B_MAT = 1;    % '' indirect effects from RawMat. consumption? Default=1
B_CaS = 0;    % '' indirect effects from plant construction & salvage?
B_PC = 0;     % '' indirect effects from product consumption?
%%% [MATERIAL] Global Warming Potentials for Material Consumption
LDPE_per_Module = 1.5;                % kg LDPE per Cultivation Module (3kg per 2 years per 1000L bag)
%%% [UTILITY] Global Warming Potentials for Utility Consumption
GWI_WASTEWATER = 0.000553;            % tonCO2eq/m3 Wastewater
%%% [CONSUMPTION] Product Biocrude Consumption
GWI_Biofiller_Cons = 0;               % For biodegradable plastic additive




%% DERIVED PARAMETERS AND DESIGN SPECS
% [1] Cultivation Process
% [2] Harvest+Dewater Process
% [3] Blowdown Recycle Treatment
% [4] Convective Flue Gas Drying
% [TEA] Techno-Economic Analysis Parameters
OP_Eff = HAR_T/(HAR_T+STOCK_T+SEP_T+max([GM_T,REC_T])); % Operating eff, %
Num_Har = floor((OP_Eff*8760)/HAR_T);                   % # of Harvests/yr
% [LCA] CO2 Life Cycle Assessment Parameters




%% PROCESS MASS BALANCE MODEL
% [1] Cultivation Process (Kinetic)
% (a) Load Solar Intensity and Temperature Data
Nat_Data = load('PAR_Temp_Data-Full.mat');
Temp = Nat_Data.Temp';
PAR = Nat_Data.L_int';
Time = Nat_Data.Exp_t';

%Nat_Data = load('PAR_Temp_08.mat');
%Temp = Nat_Data.Temp;
%PAR = Nat_Data.L_int;
%Time = Nat_Data.Exp_t;

% (b) Calculate Growth Rate at each Timestep
Model_Mu = Mu_max*(PAR./(K_I+PAR)).*(exp(-((Temp-T_opt)./Theta).^2));
Model_D = k_D.*(exp(-C_D*PAR));

% (c) Calculate Concentration at Harvest using Growth Model
Num_Steps = size(Time,2);
%Num_Steps = size(Time,1);
Model_Cx = zeros(1,Num_Steps);
Model_Cx(1) = Cx0;
for i = 2:Num_Steps
    Model_Cx(i) = Model_Cx(i-1)*exp((Model_Mu(i-1)-Model_D(i-1))*(Time(i)-Time(i-1)));
end
%Model_Cx;
Final_Conc = Model_Cx(end);

% figure()
% plot(Time, PAR)
% hold on
% plot(Time2, PAR2)
% hold off
% figure()
% plot(Time, Temp)
% hold on
% plot(Time2, Temp2)
% hold off
% figure()
% plot(Time, Model_Cx)
% xlabel('Time, mins')
% ylabel('Concentration, g/L')
 



% (d) Calculate Harvest Fraction based on Final Concentration
Har_Frac = 1 - (Cx0/Final_Conc);    
prod_bm = (Har_Frac * Final_Conc * V_Cult_T)/1000;      % g/L * L * kg/1000g
% (e) Complete Cultivation Mass Balance
% Initial Nutrient loading of [NH4Cl, KH2PO4, K2HPO4, MgSO4, CaCl2, Na2-EDTA], (kg/L)
nutri_conc = [0.000375, 0.000288, 0.000144, 0.0001, 0.00005, 0.00005];      
% Specific (kg/kg) Biomass Consumption of [Water, CO2, NH4Cl, KH2PO4, K2HPO4, MgSO4, CaCl2, Na2-EDTA]
cult_stoic = [0.71, 1.54, 0.28, 0.018, 0.009, 0.011, 0.0014, 0.0014];   
cons_cult = cult_stoic.*prod_bm;                                  
% (f) Exit Flow Composition: [Water, Biomass, NH4Cl, KH2PO4, K2HPO4, MgSO4, CaCl2, Na2-EDTA], kg/hr 
exit_cult = [];
exit_cult(1) = (V_Cult_T - cons_cult(1))*Har_Frac/HAR_T;
exit_cult(2) = prod_bm*Har_Frac/HAR_T;
exit_cult(3:8) = ((V_Cult_T*nutri_conc)-cons_cult(3:8)).*(Har_Frac/HAR_T);

% [2] Harvest+Dewater Process  (Black Box)
% (a) Dissolved Air Flotation (Assume imperfect recovery is used as
% innoculum, which is calculated by harvest frac. So just take into account
% capital and energy expense of DAF
% (b) Centrifugation
bm_CENTR = rec_CENTR*exit_cult(2);             % Recovered biomass, kg/hr
water_CENTR = (bm_CENTR/conc_CENTR)-bm_CENTR;  % Exit water flowrate, kg/hr
exit_CENTR = [water_CENTR, bm_CENTR, water_CENTR*nutri_conc]; 
% (b) Recycle Culture to Blowdown
recycle_CENTR = exit_cult - exit_CENTR;        % kg/hr

% [3] Blowdown + Media Recycle (Black Box)
% (a) Exit Streams
blowdown_loss = recycle_CENTR.*BD_LOSS;                   % kg/hr
wet_bm = exit_CENTR(1:2);                                 % kg/hr 
% (b) Make-Up Stream, kg
Cult_Initial = nutri_conc.*V_Cult_T;
Cult_Consumption = cult_stoic(3:8).*prod_bm;
Cult_Final = Cult_Initial - Cult_Consumption;
Cult_Leftover = Cult_Final.*(1-Har_Frac);
Cult_Recycle = (Cult_Final.*Har_Frac) - blowdown_loss(3:8);
MU_Nutri = Cult_Initial - Cult_Leftover - Cult_Recycle;
MU_Water = cons_cult(1) + (BD_LOSS(1)*V_Cult_T*Har_Frac);

% [4] Convective Flue Gas Dryer (Black Box)
%%% Assume an adiabatic dryer where the Biomass is dried from direct
%%% convective contact with the flue gas and not via evaporation from heat
%%% conduction. The direct contact of flue gas to the Biomass does transfer
%%% heat, but it is still adiabatic because the evaporation is released with
%%% the flue gas at the outlet.
% (a) Water Evaporation Rate
Exit_Water_Flow = Target_Moisture * (wet_bm(2)/(1 - Target_Moisture));% kg/hr
Evaporation_Rate = wet_bm(1) - Exit_Water_Flow;                       % kg/hr
% (b) Humidity of Inlet
FG_Comp_Mass = (FG_MW.*FG_Comp)/sum(FG_MW.*FG_Comp);    
Mois_FG_in = FG_Comp_Mass(1);
Psat_FG_in = 10^(Antoine_Coeff_Water(1) - (Antoine_Coeff_Water(2)/(Antoine_Coeff_Water(3) + FG_Temp_in)));
Hum_in = (FG_MW(1)/mean(FG_MW.*FG_Comp)) * (Mois_FG_in*Psat_FG_in/(FG_Pressure - (Mois_FG_in*Psat_FG_in)));
% (c) Inlet FG Enthalpy (assume adiabatic process Hf_in = Hf_out)
Hf_in = ((CP_FG + CP_H2O*Hum_in)*(FG_Temp_in+273)) + (H_lat*Hum_in);  % kJ/kg, equals Hf_out
% (d) Calculate Exit Humidity from Adiabatic Assumption
Hum_out = (Hf_in - CP_FG*FG_Temp_out)/((CP_H2O*FG_Temp_out) + H_lat);
% (e) Calculate Flue Gas flowrate from Humidity Difference
FG_Flow = Evaporation_Rate/(Hum_out - Hum_in);                        % kg/hr flue gas flowrate
Prod_BM = wet_bm - [Evaporation_Rate, 0];                             % kg/hr

% ORGANIZE PRODUCT STREAMS
PROD_BIOFILLER = sum(Prod_BM);   % kg/hr




%% PROCESS ENERGY BALANCE MODEL
% [1] ENERGY REQUIREMENTS: Cultivation Process
% (a) Flue Gas Centrifugal Blower Energy (DAYTIME) (in MW) [15]
RATING_FG = 7.5;                                  % kW per Blower, 0.12 bar
FG_UNIT_FWRATE = 520;                             % Nm^3/hr capacity
NonNormal_FWRATE = FG_UNIT_FWRATE*(1.01325/(FG_Pressure/750.062))*((273.15+FG_Temp_in)/273.15);
FWRATE_FG = VVM*V_Cult_T*0.06;                    % m3/hr Flue Gas Feedrate
NUM_FG = ceil(FWRATE_FG/NonNormal_FWRATE);        % Number of Blower Modules
PWR_FG = (FWRATE_FG/NonNormal_FWRATE)*RATING_FG;  % kW Power Consumption
ENR_FG = PWR_FG*Daylight*HAR_T*Num_Har*0.0036;    % Daylight FG Bubbling, GJ
% (b) Air Centrifugal Blower Energy (NIGHTTIME) (in MW) [15]
PWR_Air = PWR_FG*(31/73);                         % kW Power Consumption
ENR_Air = PWR_Air*(1-Daylight)*HAR_T*Num_Har*0.0036;  % Nightime FG Bubbling, GJ
% (c) Stock Solution Pump Energy (in MW) [15]
RATING_Stock = 5.5;                       % kW per Pump, 2.7 bar diffP
FWRATE_Stock = 40;                        % m^3/hr flowrate per pump
NUM_Stock = ceil((MU_Water/(1000*HAR_T))/FWRATE_Stock/STOCK_T);  
PWR_Stock = (((MU_Water/(1000*HAR_T))/FWRATE_Stock/STOCK_T)*RATING_Stock)/(NU_Pump*NU_Motor); % kW Power Cons.
ENR_Stock = PWR_Stock*STOCK_T*Num_Har*0.0036; % Stock Transf. Energy, GJ
% (d) Growth Media Pump Energy (in MW) [15]
RATING_GM = 10;                           % kW per Pump, 2.7 bar diffP
FWRATE_GM = 300;                          % m^3/hr flowrate per pump
NUM_GM = ceil((V_Cult_T*Har_Frac/1000)/FWRATE_GM/GM_T);      
PWR_GM = (((V_Cult_T*Har_Frac/1000)/FWRATE_GM/GM_T)*RATING_GM)/(NU_Pump*NU_Motor); % kW Power Cons.
ENR_GM = PWR_GM*GM_T*Num_Har*0.0036;      % GM Transf. Energy, GJ 

% [2] ENERGY REQUIREMENTS: Harvest+Dewater Process
% (a) Centrifugal Separator [15]
RATING_SEP = 7.5;                 % kW per Centrifugal Separator
CAPACITY_SEP = 4;                 % Culture processing capacity, m^3/hr
NUM_SEP = ceil((V_Cult_T*Har_Frac/(1000*HAR_T))/CAPACITY_SEP);
PWR_SEP = ((V_Cult_T*Har_Frac/(1000*HAR_T))/CAPACITY_SEP)*RATING_SEP;             % kW Separation Power Cons.
ENR_SEP = PWR_SEP*SEP_T*Num_Har*0.0036;   % Centrif. separation energy, GJ

% [3] ENERGY REQUIREMENTS: Blowdown Recycle Treatment
% (a) Recycle Growth Media Pump (w/ Filter) [15]
RATING_Rec = 5.5;                 % kW per Recycle Pump
FWRATE_Rec = 40;                  % m^3/hr flowrate processing capacity
NUM_Rec = ceil((recycle_CENTR(1)/1000)/FWRATE_Rec); 
PWR_Rec = (((recycle_CENTR(1)/1000)/FWRATE_Rec)*RATING_Rec)/(NU_Pump*NU_Motor); % kW Power Consumption
ENR_Rec = PWR_Rec*REC_T*Num_Har*0.0036;   % Broth Recycle pump energy, GJ

% [4] ENERGY REQUIREMENTS: Convective Flue Gas Dryer
% (a) Belt Dryer Conveyor Mechanical Belt
RATING_MechBelt = 90;                     % kJ/kg-H2O Evap [A. Giostri et al, 2016] 
FWRATE_MechBelt = Evaporation_Rate/3600;  % kg-H2O/s
PWR_MechBelt = RATING_MechBelt * FWRATE_MechBelt;  % kW mech belt power consumption
DRY_T = DRY_T_m*(1-conc_CENTR) + DRY_T_b; % Drying time in minutes
ENR_MechBelt = PWR_MechBelt*DRY_T*60/1000000;      % Belt mechanical energy, GJ
% (b) Flue Gas Centrifugal Blower 
P_DROP_FG = 0.5;                          % kPa FG pressure drop to drive blower
DENSITY_FG = 1.224;                       % kg/m^3 
PWR_FGBlower = ((FG_Flow*P_DROP_FG)/(DENSITY_FG*NU_HydFan*NU_MechFan))/3600;  % kW blower duty
ENR_FGBlower = PWR_FGBlower*DRY_T*60/1000000;

% ORGANIZE ENERGY DUTIES FOR UTILITY CALCULATIONS
DUTY_Elec = ENR_FG+ENR_Air+ENR_Stock+ENR_GM+ENR_SEP+ENR_Rec+ENR_MechBelt+ENR_FGBlower;




%% TECHNOECONOMIC EVALUATION [TEA] MODEL
%============================ CAPITAL EXPENSES ===========================%
% [1] Cultivation Process
% (a) Vertical Airlift Plastic Bag Cultivators [21]
Num_Cult_Unit = ceil(V_Cult_T/V_Cult_Unit);           % Number of 1m3 PBR Modules
T_Land = Num_Cult_Unit*Land_Cult_Unit                % Total Land Occupied
EQ_PBR_Cult = C_PBR*(T_Land/4046.86)*(CEPCI/CEPCI_PBR);
% (b) Flue Gas Centrifugal Blowers
EQ_FG_Blower = C_FG_Blower*NUM_FG*(CEPCI/CEPCI_PumpsBlowers);  
% (c) Stock Solution Centrifugal Pumps
EQ_Stock_Pump = C_Stock_Pump*NUM_Stock*(CEPCI/CEPCI_PumpsBlowers);
% (d) Growth Media Centrifugal Pumps
EQ_GM_Pump = C_GM_Pump*NUM_GM*(CEPCI/CEPCI_PumpsBlowers);     

% [2] Harvest+Dewater Process
% (a) Centrifugal Heavy Duty Separator 
EQ_SEP = C_SEP*NUM_SEP*(CEPCI/CEPCI_Separator)/SEP_T;

% [3] Blowdown + Media Recycle
% (a) Recycle Growth Media Centrifugal Pump
EQ_Rec_Pump = C_Rec_Pump*NUM_Rec*(CEPCI/CEPCI_PumpsBlowers);
% (b) Stock Solution Holding Tanks (Small Field Erected)
NUM_SS_Tank = ceil((MU_Water/1000)/CAP_SS_T);
EQ_SS_Tank = C_SS_T*NUM_SS_Tank*(CEPCI/CEPCI_Upstr_Tanks);
% (c) Growth Media Holding Tanks (Large Field Erected)
NUM_GM_Tank = ceil((V_Cult_T*Har_Frac/1000)/CAP_GM_T);
EQ_GM_Tank = C_GM_T*NUM_GM_Tank*(CEPCI/CEPCI_Upstr_Tanks);

% [4] Convective Flue Gas Dryer
% (a) Cross Sectional Area Calculation
AREA_DRY = sum(wet_bm)*(1+(wet_bm(2)/sum(wet_bm)))*(DRY_T/60)/AREA_LOAD;  %m^2 area
EQ_DRYER = 2700*AREA_DRY;          % From H. Li et al.

% CAPITAL COST SUMMARY
ISBL_Cult = EQ_PBR_Cult+EQ_FG_Blower+EQ_Stock_Pump+EQ_GM_Pump;
ISBL_HarDew = EQ_SEP;
ISBL_Media = EQ_Rec_Pump+EQ_SS_Tank+EQ_GM_Tank;
ISBL_Drying = EQ_DRYER;
ISBL_Total = ISBL_Cult + ISBL_HarDew + ISBL_Media + ISBL_Drying;
FIXED_CAPEX = ISBL_Total*(1+OSBL_OS)*(1+OSBL_DE+OSBL_CN);
CRF = (DISCO*(DISCO+1)^LIFET)/(((DISCO+1)^LIFET)-1);
ANNUALIZED_CAPEX = CRF*FIXED_CAPEX;

%EQ_PBR_Cult*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%EQ_FG_Blower*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%EQ_Stock_Pump*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%EQ_GM_Pump*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%EQ_SEP*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%EQ_Rec_Pump*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%EQ_SS_Tank*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%EQ_GM_Tank*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%EQ_DRYER*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%ISBL_Total*OSBL_OS*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%ISBL_Total*(1+OSBL_OS)*OSBL_DE*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
%ISBL_Total*(1+OSBL_OS)*OSBL_CN*1000/(PROD_BIOFILLER*HAR_T*Num_Har)

%=========================== OPERATING EXPENSES ==========================%
% [1] Plant Expenses
% (a) Land, Clean-In-Place & PBR Plastic Replacement Costs
COST_PE_LCP = OP_PE_MAINT*(T_Land/4046.86);
COST_PE_LDPE = (LDPE_per_Module*Num_Cult_Unit/1000)*OP_RM_LDPE;
% (b) Innoculum Costs
COST_PE_INOC = sum([OP_RM_NH4Cl, OP_RM_KH2PO4, OP_RM_K2HPO4, OP_RM_MgSO4, OP_RM_CaCl2, OP_RM_NaEDTA].*(nutri_conc*V_Cult_T/1000))/LIFET;

% [2] Raw Material Costs
% (a) Flue Gas Utilization Costs
MOLES_FG = (prod_bm*cult_stoic(2)*Num_Har*(1/44.01))/FG_Moles(1);
COST_RM_FG = OP_RM_FG*MOLES_FG*sum(FG_MM.*FG_Moles);
% (b) Nutrient Stock Make-Up Costs
MU_Nutri = MU_Nutri*Num_Har;  % kg/yr Makeup Nutrients
COST_RM_NU = sum([OP_RM_NH4Cl, OP_RM_KH2PO4, OP_RM_K2HPO4, OP_RM_MgSO4, OP_RM_CaCl2, OP_RM_NaEDTA].*(MU_Nutri/1000));
% (c) Water Make-up Costs 
MU_Water = MU_Water*Num_Har;  % kg/yr Makeup Water
COST_RM_H2O = OP_RM_WATER*(MU_Water/1000);      

% [3] Utility Costs
COST_UT_ELEC = UT_ELEC*DUTY_Elec;          % Yearly Electrical Costs

%======================== FIXED OPERATING EXPENSES =======================%
% Labor
LAND_ACRE = T_Land/4046.86;
INF_2014 = 1.098;
SAL_PLANT_MANAGER = 155617*INF_2014*LAND_ACRE/5000;
SAL_PLANT_ENGINEER = 82050*INF_2014*LAND_ACRE/5000;
SAL_MAINT_SUPERVISOR = 60341*INF_2014*LAND_ACRE/5000;
SAL_MODULE_OPERATOR = 38590*INF_2014*LAND_ACRE/5000;
SAL_CLERK = 38110*INF_2014*LAND_ACRE/5000;
SAL_FIELD_EMPLOYEE = 3500*LAND_ACRE*INF_2014;
LABOR = SAL_FIELD_EMPLOYEE+SAL_PLANT_MANAGER+SAL_PLANT_ENGINEER+SAL_MAINT_SUPERVISOR+SAL_MODULE_OPERATOR+SAL_CLERK;
% Maintenance
%MAINT_PBR = 0.05*ISBL_Cult;
MAINT_Else = 0.03*(ISBL_HarDew + ISBL_Media + ISBL_Drying);
MAINTAINENCE = MAINT_Else;
% Administration + Overhead
ADMIN_OVERHEAD = 0.90*LABOR;
LABORATORY = 0.01*(ISBL_Cult + ISBL_HarDew + ISBL_Media + ISBL_Drying);
TOTAL_FIXED = LABOR + MAINTAINENCE + ADMIN_OVERHEAD + LABORATORY;
% OPERATING COST SUMMARY
TOTAL_OPEX_PE = COST_PE_LCP + COST_PE_LDPE + COST_PE_INOC;
TOTAL_OPEX_RM = COST_RM_FG + COST_RM_NU + COST_RM_H2O;
TOTAL_OPEX_UT = COST_UT_ELEC;
%============================== COST SUMMARY =============================%
ANNUAL_PROD_COST = TOTAL_OPEX_PE + TOTAL_OPEX_RM + TOTAL_OPEX_UT + ANNUALIZED_CAPEX + TOTAL_FIXED;
% Cost of Goods Manufactured, befitting TRL 3&4
COGM = ANNUAL_PROD_COST*1000/(PROD_BIOFILLER*HAR_T*Num_Har);  % USD/ton



%% CO2 LIFE CYCLE ASSESSMENT MODEL
% [1] Direct Plant Emissions of GHG
% (a) Cultivation Off Gases **ASSUME FG IS FED AS BYPASS**
CULT_CO2_In_Moles = (VVM*V_Cult_T*0.06/FG_Moles(1))*(FG_Moles(1)/1000)*HAR_T*Num_Har;
CULT_CO2_In = CULT_CO2_In_Moles*FG_MM(1);      % kg CO2 in Total
CULT_CO2_Consumed = cons_cult(2)*Num_Har;
% (b) Sequestration as TIC in Media
O2_Production = prod_bm*1.2567409*Num_Har;     % kg O2 in Total
Bicarb_Production = prod_bm*0.3091724*Num_Har; % kg Bicarb in Total
CULT_CO2_Out = CULT_CO2_In - CULT_CO2_Consumed - Bicarb_Production;
% Sum all Direct Emissions
TOT_EM_DIR = -CULT_CO2_Consumed/1000;     % BYPASS APPROACH, tons CO2eq

% [2] Indirect Plant Emissions of GHG
% (a) Indirect Emissions from Energy Consumption
EM_ENR_ELEC = GWI_ELEC_GRID*DUTY_Elec;
% (b) Indirect Emissions from Wastewater Treatment
EM_WASTEWATER = GWI_WASTEWATER*(sum(blowdown_loss)*HAR_T*Num_Har/1000); % ton CO2
TOT_EM_ENR = EM_ENR_ELEC + EM_WASTEWATER;
% (c) Indirect Emissions from Raw Material Consumption
EM_RM_NH4Cl = ((((nutri_conc(1)*V_Cult_T)/LIFET)+MU_Nutri(1))/1000)*GWI_NH4Cl;
EM_RM_KH2PO4 = ((((nutri_conc(2)*V_Cult_T)/LIFET)+sum(MU_Nutri(2:3)))/1000)*GWI_KH2PO4;
EM_RM_MGSO4 = ((((nutri_conc(3)*V_Cult_T)/LIFET)+MU_Nutri(4))/1000)*GWI_MgSO4;
EM_RM_CACL2 = ((((nutri_conc(4)*V_Cult_T)/LIFET)+MU_Nutri(5))/1000)*GWI_CaCl2;
EM_RM_NAEDTA = ((((nutri_conc(5)*V_Cult_T)/LIFET)+MU_Nutri(6))/1000)*0.053;
EM_RM_WATER = (MU_Water/1000)*GWI_WATER;
EM_PBR_REPLACE = GWI_LDPE*(LDPE_per_Module*Num_Cult_Unit/1000);
TOT_EM_RM = EM_RM_NH4Cl + EM_RM_KH2PO4 + EM_RM_MGSO4 + EM_RM_CACL2 + EM_RM_NAEDTA + EM_RM_WATER + EM_PBR_REPLACE;
% (d) Indirect Emissions from Plant Construction and Salvage
TOT_EM_CaS = 0;             % Not Studied
% (e) Indirect Emissions from Product Consumption
TOT_EM_PC = GWI_Biofiller_Cons*(PROD_BIOFILLER*HAR_T*Num_Har);  

% CO2 GWI SUMMARY (in kg CO2-eq)
TOTAL_GWI = B_DIR*TOT_EM_DIR + B_ENR*TOT_EM_ENR + B_MAT*TOT_EM_RM + B_CaS*TOT_EM_CaS + B_PC*TOT_EM_PC;
SPECIFIC_GWI = TOTAL_GWI*1000/(PROD_BIOFILLER*HAR_T*Num_Har);




% %% Baseline: Table
% % Raw Materials-
% TOTAL_PRODUCT = (PROD_BIOFILLER*HAR_T*Num_Har)/1000
% TOTAL_LAND = T_Land
% R_NH4Cl = MU_Nutri(1)/(PROD_BIOFILLER*HAR_T*Num_Har)
% R_KH2PO4 = MU_Nutri(2)/(PROD_BIOFILLER*HAR_T*Num_Har)
% R_K2HPO4 = MU_Nutri(3)/(PROD_BIOFILLER*HAR_T*Num_Har)
% R_MGSO4 = MU_Nutri(4)/(PROD_BIOFILLER*HAR_T*Num_Har)
% R_CaCl2 = MU_Nutri(5)/(PROD_BIOFILLER*HAR_T*Num_Har)
% R_NaEDTA = MU_Nutri(6)/(PROD_BIOFILLER*HAR_T*Num_Har)
% R_WATER = MU_Water/(PROD_BIOFILLER*HAR_T*Num_Har)
% R_FLUE_GAS = CULT_CO2_In/(PROD_BIOFILLER*HAR_T*Num_Har)
% % Utilities
% U_ELEC = DUTY_Elec*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
ENR_FG*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
ENR_Air*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% %ENR_Stock*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% %ENR_GM*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% %ENR_SEP*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% %ENR_Rec*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% %ENR_MechBelt*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% %ENR_FGBlower*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% %FG_Feed_Total = CULT_CO2_In/(PROD_BIOFILLER*HAR_T*Num_Har)
% %U_WASTE_HEAT = FG_Flow*(DRY_T/60)*Num_Har/(PROD_BIOFILLER*HAR_T*Num_Har)
% U_WASTE_HEAT_GJ = (Hf_in*FG_Flow/1000000)*(DRY_T/60)*Num_Har/(PROD_BIOFILLER*HAR_T*Num_Har)
% % Direct Emissions
% CO2_Bypass = CULT_CO2_Out/(PROD_BIOFILLER*HAR_T*Num_Har)
% O2_Byproduct = O2_Production/(PROD_BIOFILLER*HAR_T*Num_Har)
% Bicarb_Waste = Bicarb_Production/(PROD_BIOFILLER*HAR_T*Num_Har)
% Wastewater = sum(blowdown_loss)/PROD_BIOFILLER
% 
% %% Baseline: COGM Breakdown
% OPEX_Plant_Expenses = TOTAL_OPEX_PE*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% OPEX_Utilities = TOTAL_OPEX_UT*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% OPEX_Raw_Materials = TOTAL_OPEX_RM*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% Annual_CAPEX = ANNUALIZED_CAPEX*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% OPEX_Fixed = TOTAL_FIXED*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% 
% %% Baseline: Carbon Footprint Breakdown
% Direct_Consumption = TOT_EM_DIR*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% Indirect_Energy = TOT_EM_ENR*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% LDPE_Consumption = (LDPE_per_Module*Num_Cult_Unit/1000)*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% Indirect_RawMaterials = TOT_EM_RM*1000/(PROD_BIOFILLER*HAR_T*Num_Har)
% 
% %% Mass Balancce Info
% % Flue Gas Composition Stream 1a, 2a (day)
% Stream1_FG = ((VVM*V_Cult_T)/1000)*60*DENSITY_FG*FG_Comp;   % [H2O, N2, CO2, O2]
% % Flue Gas Composition Stream 1b, 2b (night)
% Stream2_FG = ((VVM*(31/73)*V_Cult_T)/1000)*60*DENSITY_FG*FG_Comp;
% % Water makeup per batch Stream 3
% Stream3_Water = MU_Water/Num_Har;
% % Nutrient makeup per batch Stream 4
% Stream4_Nutri = MU_Nutri/Num_Har;
% % Recycle media per batch Stream 11
% Stream11_Nutri = Cult_Recycle;    % NH4Cl, KH2PO4, K2HPO4, MgSO4, CaCl2, Na2-EDTA
% Stream11_Water = (exit_cult(1)*HAR_T)-(MU_Water/Num_Har)-(blowdown_loss(1)*HAR_T); % Water
% % Bypass to Stack Stream 7
% Sequestration_Rate = (cons_cult(2) + (prod_bm*0.3091724))/HAR_T; %kg/hr
% Stream7 = Stream1_FG;
% Stream7(3) = Stream7(3) - Sequestration_Rate;
% Stream7(4) = Stream7(4) + (O2_Production/(HAR_T*Num_Har));
% Stream7;
% % Blowdown Loss Stream 10
% Stream10_Blowdown = (blowdown_loss*HAR_T);
% % PBR Exit Stream 8
% Stream8_Exit = exit_cult*HAR_T; % Water, Biomass, NH4Cl, KH2PO4, K2HPO4, MgSO4, CaCl2, Na2-EDTA
% Stream8_Bicarb = Bicarb_Production/Num_Har;
% % Sludge to Dryer, Stream 12
% Stream12_Sludge = wet_bm*HAR_T;
% % Stream 11 Biomass Recycle to Cult
% Stream11_Biomass = Stream8_Exit(2) - Stream12_Sludge(2);
% % Stream 13 Flue Gas Flow 
% Stream13 = FG_Flow*HAR_T*FG_Comp;
% % Stream 15 Produced Biomass
% Stream15 = Prod_BM*HAR_T;
% 
% %% Volumetric Productivity
% Productivity = prod_bm/(V_Cult_T/1000)/(149/24)   % kg/m3/day




%% EXPORT SUSTAINABILITY CRITERIA METRICS
% Criteria 1: Cost of Goods Manufactured (Technoeconomic Criteria)
% Criteria 2: CO2 Global Warming Impact (Environmental Emissions Criteria)
EvalMetrics = [COGM, SPECIFIC_GWI];    
end