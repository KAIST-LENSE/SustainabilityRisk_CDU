%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Process+Evaluation Model for CO2-based Formic Acid Production %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Revision_version_in_terms_of_Formic Acid Purity 99wt% %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EvalMetrics = MODEL_FA(kp)   
%% Declare Global Variables
format long g
%% KEY PARAMETERS WITH UNCERTAINTY
% Process Parameters
CO2_Conv = kp(1);                       % Conversion, 0.35 (%)
TREA_loss_rate = kp(2);                 % Loss Rate, 0.35 (%)
NBIM_losss_rate = kp(3);                % Loss Rate, 0.35 (%)
% [TEA] Techno-Economic Analysis Parameters
REFCOST_electricity = kp(4);            % 25.9 USD/GJ
REFCOST_MP_steam = kp(5);               % 17.8 USD/GJ
REFCOST_NG = kp(6);                     % 250 USD/ton
REFCOST_CO2 = kp(7);                    % 60 USD/ton
REFCOST_TREA = kp(8);                   % 2544.3 USD/ton
REFCOST_NBIM = kp(9);                   % 8000 USD/ton
% [LCA] CO2 Life Cycle Assessment Parameters
GWI_CO2_captured = kp(10);              % -0.86 tonCO2eq/ton
GWI_Elec = kp(11);                      % 0.185 tonCO2eq/GJ
GWI_MP_Steam = kp(12);                  % 0.069 tonCO2eq/GJ
GWI_NG = kp(13);                        % 0.301 tonCO2eq/ton
GWI_Water = kp(14);                     % 0.000173 tonCO2eq/ton


% Process Parameters
% CO2_Conv = 0.35;                       % Conversion, 0.35 (%)
% TREA_loss_rate = 0.35;                 % Loss Rate, 0.35 (%)
% NBIM_losss_rate = 0.35;                % Loss Rate, 0.35 (%)
% % [TEA] Techno-Economic Analysis Parameters
% REFCOST_electricity = 25.9;            % 25.9 USD/GJ
% REFCOST_MP_steam = 17.8;               % 17.8 USD/GJ
% REFCOST_NG = 250;                     % 250 USD/ton
% REFCOST_CO2 = 60;                    % 60 USD/ton
% REFCOST_TREA = 2544.3;                   % 2544.3 USD/ton
% REFCOST_NBIM = 8000;                   % 8000 USD/ton
% % [LCA] CO2 Life Cycle Assessment Parameters
% GWI_CO2_captured = -0.86;              % -0.86 tonCO2eq/ton
% GWI_Elec = 0.185;                      % 0.185 tonCO2eq/GJ
% GWI_MP_Steam =0.069;                  % 0.069 tonCO2eq/GJ
% GWI_NG =0.301;                        % 0.301 tonCO2eq/ton
% GWI_Water = 0.000173;                     % 0.000173 tonCO2eq/ton


%% NON VARYING FIXED-PARAMETERS & DESIGN SPECS
% Molar Weights
MW_FA = 46.026;              %kg/kmol = g/mol
MW_NG = 16;                  %kg/kmol = g/mol
MW_CO2 = 44;
MW_TREA = 101.192;
MW_NBIM = 124.18;
MW_H2O = 18;
% Process Params (SMR + Hydrogenation Reactor)
L_SMR=0.00069 ;              % SMR Reactor_length (m)
D_SMR=0.1016;                % SMR Reactor_diameter(m)
R=8.314472;                  % J/mol K
p_SMR=2900;                  % kPa, constant pressure
n_SMR=1000;                  % Number of for loop
num_SMR=210;                 % reactor_divide
TREA_concentration = 4;      % 3M
T_SMR = 1223;                % T_SMR=950+273; % K
production_scale = 1;        % FA production scale(ton/yr)
% [TEA] Techno-Economic Analysis Parameters
REFCOST_comp = 180000;              %USD
REFCOST_pump = 10000;               %USD
REFCOST_reactor = 110000;           %USD, Gas catalytic reactor
REFCOST_FUR = 250000;               %USD
REFCOST_HEX = 70000;                %USD
REFCOST_column = 100000;            %USD, Vertical column
REFCOST_Col_1 = 8142400;
REFCOST_Reb_1 = 201500;
REFCOST_Col_2 = 265600;
REFCOST_Reb_2 = 205900;
REFCOST_Cond_2 = 96400;
REFCOST_Col_3 = 134800;
REFCOST_Reb_3 = 35800;
REFCOST_Cond_3 = 14400;
REF_electricity_comp = 100;         %kW
REF_electricity_pump = 23;          %kW
REF_gas_cat_reactor = 20;           %m^2.5 = length*diameter^1.5
REF_HEX = 100;                      %heat transfer area (m^2)
REF_FUR = 1;                        %MW
REF_COL = 20;                       %m^2.5 = length*diameter^1.5
% [LANG] Capital Cost Scaling Factors
NTH_PMP = 0.79;
NTH_COMP = 0.79;
NTH_GRXTOR = 0.81;
NTH_HEX = 0.71;
NTH_FUR = 0.74;
NTH_COL = 0.81;
NTH_Tray = 0.78;
NTH_Reb = 0.78;
f_IEC = 1.1;                        % installation factor
f_DCC = 1.44;                       % Direct capital cost factor
f_ICC = 0.6;                        % Indirect capital cost factor
f_WC = 0.05;                        % working capital
DISCO = 0.07;                       % annual discount rate
CEPCI = 619.2;                      %2019 based
CEPCI_ref = 1000;
Lifetime_plant = 30;
Operating_hour = 8100;
REFCOST_Labor = 970797;             % USD/year, 200,000 tonne/yr production based
f_Maintenance = 0.03;
f_admin = 0.2;
f_overhead = 0.6;
f_Laboratory = 0.01;
% [LCA]
GWI_CO2 = -1;                       %ton/tonFA = kg/kgFA




%% PROCESS MASS BALANCE MODEL
% [1] SMR Process
% (a) Feed Preparation
NG_flowrate = 21.0377*production_scale;             % kmol/hr - 30000tone FA 
H2O_flowrate = 62.5075*production_scale;          % kmol/hr - 30000tone FA 
% (b) Steam Methane Reforming Rate Constants
SMR_inlet_flowrate = 86.4127*production_scale;      % kmol/hr
unit_change = 1000/60;                              % kmol/hr -> mol/min
CO2_frac=0.0164121;                                 % SMR_inlet_molfrac 
H2_frac=0.0626852;
CH4_frac=0.230199;
H2O_frac=0.690357;
CO_frac=0.000179439;
N2_frac= 0.0002435;
p_CH4=p_SMR*CH4_frac/num_SMR;                       % partial pressure kPa
p_H2O=p_SMR*H2O_frac/num_SMR;
p_CO=p_SMR*CO_frac/num_SMR;
p_H2=p_SMR*H2_frac/num_SMR;
p_CO2=p_SMR*CO2_frac/num_SMR;
k_1=7.382*10^14*exp(-261.9*10^3/(R*T_SMR));         % molkPa^0.5/gcat min
K_1=1.198*10^17*exp(-26830/T_SMR);                  % kPa^2
k_2=3.783*10^4*exp(-38.9*10^3/(R*T_SMR));           % mol/gcat min Pa
K_2=1.767*10^-2*exp(4400/T_SMR);                    % kPa^0
k_3=1.636*10^14*exp(-243.7*10^3/(R*T_SMR));         % molkPa^0.5/gcat min
K_3=2.117*10^15*exp(-22430/T_SMR);                  % kPa^2
K_CO=8.825*10^-7*exp(85900/(R*T_SMR));              % kPa^-1
K_H2=6.051*10^-11*exp(79700/(R*T_SMR));
K_CH4=7.015*10^-6*exp(42300/(R*T_SMR));
K_H2O=2.206*10^5*exp(-82900/(R*T_SMR));
DEN=1+K_CO*p_CO+K_H2*p_H2+K_CH4*p_CH4+K_H2O*p_H2O/p_H2;
r1=k_1*(p_CH4*p_H2O/p_H2^2.5-p_H2^0.5*p_CO/K_1)/DEN^2;
r2=k_2*(p_CO*p_H2O/p_H2-p_CO2/K_2)/DEN^2;
r3=k_3*(p_CH4*p_H2O^2/p_H2^3.5-p_H2^0.5*p_CO2/K_3)/DEN^2;
r_CO=r1-r2;
r_CO2=r2+r3;
r_CH4=-r1-r3;
r_H2=3*r1+r2+4*r3;
r_H2O=-r1-r2-2*r3;
w_cat=pi*((D_SMR/2)^2)*L_SMR*2355.5*0.6*1000;         % g_cat
F_COin=0.25843/num_SMR;
F_CO2in=23.6369/num_SMR;
F_CH4in=331.42/num_SMR;
F_H2in=90.2799/num_SMR;
F_H2Oin=994.26/num_SMR;
% (c) SMR Mass Balance
for i = 1:n_SMR
    F_COout=r_CO*w_cat/n_SMR+F_COin;
    F_COin=F_COout;
    L_SMR(i)=i;
    F_CO(i)=F_COout;
    F_CO2out=r_CO2*w_cat/n_SMR+F_CO2in;
    F_CO2in=F_CO2out;
    F_CO2(i)=F_CO2out;
    F_CH4out=r_CH4*w_cat/n_SMR+F_CH4in;
    F_CH4in=F_CH4out;
    F_CH4(i)=F_CH4out;
    F_H2out=r_H2*w_cat/n_SMR+F_H2in;
    F_H2in=F_H2out;
    F_H2(i)=F_H2out;
    F_H2Oout=r_H2O*w_cat/n_SMR+F_H2Oin;
    F_H2Oin=F_H2Oout;
    F_H2O(i)=F_H2Oout;
    F_tot=F_COin+F_CO2in+F_CH4in+F_H2in+F_H2Oin;
    p_CO=p_SMR*F_COin/F_tot/num_SMR;
    p_CO2=p_SMR*F_CO2in/F_tot/num_SMR;
    p_CH4=p_SMR*F_CH4in/F_tot/num_SMR;
    p_H2=p_SMR*F_H2in/F_tot/num_SMR;
    p_H2O=p_SMR*F_H2Oin/F_tot/num_SMR;
    DEN=1+K_CO*p_CO+K_H2*p_H2+K_CH4*p_CH4+K_H2O*p_H2O/p_H2;
    r1=k_1*(p_CH4*p_H2O/p_H2^2.5-1/K_1*p_H2^0.5*p_CO)/DEN^2;
    r2=k_2*(p_CO*p_H2O/p_H2-p_CO2/K_2)/DEN^2;
    r3=k_3*(p_CH4*p_H2O^2/p_H2^3.5-p_H2^0.5*p_CO2/K_3)/DEN^2;
    r_CO=r1-r2;
    r_CO2=r2+r3;
    r_CH4=-r1-r3;
    r_H2=3*r1+r2+4*r3;
    r_H2O=-r1-r2-2*r3;
end
F_CH4out=F_CH4out*num_SMR;
F_CO2out=F_CO2out*num_SMR;
F_H2Oout=F_H2Oout*num_SMR;
F_H2out=F_H2out*num_SMR;
F_COout=F_COout*num_SMR;
% Component Indices : CO2(1), H2(2), O2(3), N2(4), CH4(5), TREA(6), H2O(7), HCOOH(8), NBIM(9) TEA-FA(10), NBIM-FA(11), CO(12)
SMR_reactor_out = zeros(1,12);
SMR_reactor_out(1) = F_CO2out/unit_change;
SMR_reactor_out(2) = F_H2out/unit_change;
SMR_reactor_out(4) = N2_frac*SMR_inlet_flowrate;
SMR_reactor_out(5) = F_CH4out/unit_change;
SMR_reactor_out(7) = F_H2Oout/unit_change;
SMR_reactor_out(12) = F_COout/unit_change;
% (d) High Temperature Water-Gas-Shift (HT-WGS)
HT_WGS_out = SMR_reactor_out;   % equilibrium based reaction
HT_WGS_Keq = 31.0836;           % at 320 celius
HT_eqn= [SMR_reactor_out(12),SMR_reactor_out(7),SMR_reactor_out(1),SMR_reactor_out(2),HT_WGS_Keq];
a = HT_eqn(5)-1;
b = (HT_eqn(1)+HT_eqn(2))*HT_eqn(5)+HT_eqn(3)+HT_eqn(4);
c = HT_eqn(1)*HT_eqn(2)*HT_eqn(5)-HT_eqn(3)*HT_eqn(4);
WGS_S = (b-sqrt(b^2-4*a*c))/(2*a);
if WGS_S > min(SMR_reactor_out(12),SMR_reactor_out(7))
    WGS_S = (b+sqrt(b^2-4*a*c))/(2*a);
end
%%%DEBUGGING CODE%%%
if 4*a*c > b^2
    WGS_S = b/(2*a);
end
%%%%%%%%%%%%%%%%%%%%%
HT_WGS_out(1) = SMR_reactor_out(1) + WGS_S;
HT_WGS_out(2) = SMR_reactor_out(2) + WGS_S;
HT_WGS_out(12) = SMR_reactor_out(12) - WGS_S;
HT_WGS_out(7) = SMR_reactor_out(7) - WGS_S;
% (e) Low Temperature Water-Gas-Shift (LT_WGS)
LT_WGS_out = HT_WGS_out;    % equilibrium based reaction
LT_WGS_Keq = 295.436;       % at 190 celius
LT_eqn= [HT_WGS_out(12),HT_WGS_out(7),HT_WGS_out(1),HT_WGS_out(2),LT_WGS_Keq];
a = LT_eqn(5)-1;
b = (LT_eqn(1)+LT_eqn(2))*LT_eqn(5)+LT_eqn(3)+LT_eqn(4);
c = LT_eqn(1)*LT_eqn(2)*LT_eqn(5)-LT_eqn(3)*LT_eqn(4);
LT_WGS_S = (b-sqrt(b^2-4*a*c))/(2*a);
if LT_WGS_S > min(HT_WGS_out(12),HT_WGS_out(7))
    LT_WGS_S = (b+sqrt(b^2-4*a*c))/(2*a);
end
%%%DEBUGGING CODE%%%
if 4*a*c > b^2
    LT_WGS_S = b/(2*a);
end
%%%%%%%%%%%%%%%%%%%%%
LT_WGS_out(1) = HT_WGS_out(1) + LT_WGS_S;
LT_WGS_out(2) = HT_WGS_out(2) + LT_WGS_S;
LT_WGS_out(12) = HT_WGS_out(12) - LT_WGS_S;
LT_WGS_out(7) = HT_WGS_out(7) - LT_WGS_S;
% (f) Gas Separation Flash Drum
FL_flowrate = sum(LT_WGS_out);
FL_top_molfrac = [0.194947,0.783896,0,0.000210377,0.0158805,0,0,0,0,0,0,0.00231169];
FL_top=FL_flowrate*FL_top_molfrac*100/123.0143;

% [2] CO2 Hydrogenation Reactor
% (a) Conversion Reaction: CO2+H2+TREA -> TEA-HCOOH
conversion_factor = (CO2_Conv/0.35)^(-0.721);
SMR_FUR = 27.0381*production_scale*[1,0,0,0,0,0,0,0,0,0,0,0];
SMR_EFFL = FL_top;
CO2_feed_flowrate = FL_top(2)-(FL_top(1)+SMR_FUR(1));
% Component Indices : CO2(1), H2(2), O2(3), N2(4), CH4(5), TREA(6), H2O(7), HCOOH(8), NBIM(9) TEA-FA(10), NBIM-FA(11), CO(12)
CO2 = [CO2_feed_flowrate,0,0,0,0,0,0,0,0,0,0,0];
TEA_H2O_initial = [0,0,0,0,0,90,75,0,0,0,0,0];
TEA_H2O = TEA_H2O_initial;
TEA_H2O(7) = TEA_H2O(6)*(500/(9*TREA_concentration)-7.74334);
Recycle_G_initial = 300*[0.427445,0.427083,0,0,0.12,0,0,0,0,0,0,0];
Recycle_G = Recycle_G_initial;
for i = 1:500
    % (b) Hydrogenation Mass Balanace (Stoichiometric)
    CO2_hydrogenation_in = CO2 + TEA_H2O + SMR_FUR + SMR_EFFL + Recycle_G;
    TEA_H2O(6) = CO2_hydrogenation_in(2)*(CO2_Conv+0.01);
    TEA_H2O(7) = TEA_H2O(6)*(500/(9*TREA_concentration)-7.74334);
    CO2_hydrogenation_in = CO2 + TEA_H2O + SMR_FUR + SMR_EFFL + Recycle_G;
    stoichiometric = [-1,-1,0,0,0,-1,0,0,0,1,0,0];
    CO2_hydrogenation_out= CO2_hydrogenation_in + (CO2_hydrogenation_in(1)*CO2_Conv).*stoichiometric;
    % (c) Purge and Separation of Unwanted
    Unreact_G = CO2_hydrogenation_out.*[1,1,1,1,1,0,0,0,0,0,0,1];
    Purge=0.04*Unreact_G;
    Recycle_G = Unreact_G*0.96;
    reaction_product_out = CO2_hydrogenation_out.*[0,0,0,0,0,1,1,1,1,1,1,0];
end

% [3] COLUMN 1: TEA-Solvent Recovery
% (a) Separation Fractions
column_1_bottom = reaction_product_out.*[0,0,0,0,0,0,0,0,0,1,0,0];
column_1_top = reaction_product_out.*[0,0,0,0,0,1,1,0,0,0,0,0];

% [4] COLUMN 2: NBIM-replacement
% (a) Column Mass Balance: TEA-FA + NBIM -> TREA + NBIM-FA
NBIM_initial = zeros(1,12);
NBIM_initial(9) = column_1_bottom(10);
column_2_feed = zeros(1, 12);
column_2_bottom = zeros(1, 12);
column_2_top = zeros(1, 12);
column_2_feed(6) = NBIM_initial(9);
column_2_feed(11) = column_1_bottom(10);
column_2_bottom(11) = column_2_feed(11);    % NBIM-FA
column_2_top(6) = column_2_feed(6);         % TREA

% [5] COLUMN 3: FA Purification 
% (a) Column MAss Balance: NBIM-FA -> FA + NBIM
column_3_feed = zeros(1, 12);
column_3_bottom = zeros(1, 12);
column_3_top = zeros(1, 12);
column_3_feed(8) = column_2_bottom(11);
column_3_feed(9) = column_2_bottom(11);
column_3_bottom(9) = column_3_feed(9);   % NBIM
column_3_top(8) = column_3_feed(8);     % FA
column_2_reboiler_heat = 5396.7/72.7704*sum(column_3_bottom);   %MJ/kmol*kmol

% ORGANIZE PRODUCT STREAMS
C_FA = column_3_top(8);     %kmol/hr

W_FA = C_FA*MW_FA/1000;      %tone/hr  
% C_Water = (15*C_FA*MW_FA/85)/MW_H2O;
% ORGANIZE RECYCLE STREAMS
TEA_H2O = column_1_top + column_2_top;




%% PROCESS ENERGY BALANCE MODEL
% [1] SMR Process
% (a) Heat Duty
NG_heater_duty = 745.76 * production_scale;             % MJ/hr
NG_combustion_duty = 6167.01 * production_scale;        % kg/hr
% (a) Electricity Duty
NG_compressing_duty = 41.9172 * production_scale;       % MJ/hr
H2O_pump_duty = 13.0285 * production_scale;             % MJ/hr
SMR_electricity = NG_compressing_duty+H2O_pump_duty;

% [2] CO2 Hydrogenation Reactor
% (a) Heat Duty
CO2_heater_enthalpy = 288.390206457/32.0152;            % MJ/kmol
CO2_heater_duty = CO2_heater_enthalpy*CO2(1)*conversion_factor;       % MJ
TEA_H2O_heater_enthalpy = 1760.3409090830385/150;       % MJ/kmol
TEA_H2O_heater_duty = TEA_H2O_heater_enthalpy*sum(TEA_H2O)*conversion_factor;
CO2_hydrogenation_heat = CO2_heater_duty+TEA_H2O_heater_duty;       % MJ/hr enthalpy change
% (b) Electricity Duty 
SMR_EFFL_comp1_electric = 321.893/100;                  % MJ/kmol
SMR_EFFL_comp2_electric = 344.822/100;
SMR_EFFL_comp1_duty = SMR_EFFL_comp1_electric*sum(SMR_EFFL);
SMR_EFFL_comp2_duty = SMR_EFFL_comp2_electric*sum(SMR_EFFL);
TEA_H2O_pump_electric = 260.036/150*conversion_factor;
TEA_H2O_pump_duty = TEA_H2O_pump_electric*sum(TEA_H2O)*conversion_factor;
SMR_FUR_multicomp1_electric = 197.766/27.0381;
SMR_FUR_multicomp2_electric = 116.28/27.0381;
SMR_FUR_multicomp3_electric = 116.28/27.0381;
SMR_FUR_multicomp4_electric = 116.28/27.0381;
SMR_FUR_multicomp5_electric = 116.28/27.0381;
SMR_FUR_duty = (SMR_FUR_multicomp1_electric+SMR_FUR_multicomp2_electric+SMR_FUR_multicomp3_electric+SMR_FUR_multicomp4_electric+SMR_FUR_multicomp5_electric)*sum(SMR_FUR);
CO2_hydrogenation_electricity = SMR_EFFL_comp1_duty+SMR_EFFL_comp2_duty+TEA_H2O_pump_duty+SMR_FUR_duty;

% [3] COLUMN 1: TEA-Solvent Recovery
% (a) Heat Duty
column_1_reboiler_heat = 12645.7/150;       % MJ/kmol
column_1_reboiler_heat_duty = column_1_reboiler_heat*sum(reaction_product_out)/2;
TEA_solvent_recovery_heat = column_1_reboiler_heat_duty;
%electricity
column_1_vacuum_pump_elec = 48.561/150;     % MJ/kmol
column_1_vacuum_pump_elec_duty = column_1_vacuum_pump_elec*sum(reaction_product_out);
TEA_solvent_recovery_electricity = column_1_vacuum_pump_elec_duty;

% [4] COLUMN 2: NBIM-replacement
% (a) Heat Duty
column_2_reaction_heat = 33422.6/73.0955;   % MJ/kmol
column_2_reboiler_heat = 4963.35/72.8362;   % MJ/kmol
column_2_reaction_heat_duty = column_2_reaction_heat*column_1_bottom(10);
column_2_reboiler_heat_duty = column_2_reboiler_heat*sum(column_2_bottom);
NBIM_replacement_heat = column_2_reaction_heat_duty + column_2_reboiler_heat_duty;

% [5] COLUMN 3: FA Purification 
% (a) Heat Duty
column_3_reboiler_heat = 5396.7/72.7704;    % MJ/kmol
column_3_reboiler_heat_duty = column_3_reboiler_heat*sum(column_3_bottom);
FA_purification_heat = column_3_reboiler_heat_duty;

% ORGANIZE ENERGY DUTIES FOR UTILITY CALCULATIONS
Elec_duty = SMR_electricity+CO2_hydrogenation_electricity+TEA_solvent_recovery_electricity; %MJ
Heat_duty = 0.75*(CO2_hydrogenation_heat+TEA_solvent_recovery_heat+NBIM_replacement_heat+FA_purification_heat); %MJ
NG_combustion_usage = 91.77 * production_scale; % NG usage - kg/hr



%% TECHNOECONOMIC EVALUATION [TEA] MODEL
%============================ CAPITAL EXPENSES ===========================%
% [1] SMR Process
% (a) Feed Preparation
FOB_NG_COMP = REFCOST_comp * ((NG_compressing_duty/3.6)/REF_electricity_comp)^NTH_COMP * CEPCI/CEPCI_ref;
FOB_H2O_PUMP = REFCOST_pump * ((H2O_pump_duty/3.6)/REF_electricity_pump)^NTH_PMP * CEPCI/CEPCI_ref;
% (b) SMR Mass Balance
FOB_SMR_reactor= REFCOST_reactor * ((12*(0.1016^1.5)*(production_scale^(2.5/3)))/REF_gas_cat_reactor)^NTH_GRXTOR*4.68*CEPCI/CEPCI_ref; % alloy factor : 3.6 , pressure factor : 1.3
FOB_FUR = REFCOST_FUR *(1.71306*production_scale/REF_FUR)^NTH_FUR*0.35*CEPCI/CEPCI_ref; %pressure factor : 0.35
% (c) Water Gas Shift Reactions
% Heat Exchangers
FOB_HEX_1 = REFCOST_HEX * (3.03419*(production_scale^(2/3))/REF_HEX)^NTH_HEX*2*CEPCI/CEPCI_ref; % 304 s/s factor : 2
FOB_HEX_2 = REFCOST_HEX * (0.79402*(production_scale^(2/3))/REF_HEX)^NTH_HEX*2*CEPCI/CEPCI_ref;
FOB_HEX_3 = REFCOST_HEX * (3.52969*(production_scale^(2/3))/REF_HEX)^NTH_HEX*2*CEPCI/CEPCI_ref;
% WGS Reactor
FOB_HT_WGSR = REFCOST_reactor * ((12*(0.1016^1.5)*(production_scale^(2.5/3)))/REF_gas_cat_reactor)^NTH_GRXTOR*4.68*CEPCI/CEPCI_ref; % alloy factor : 3.6 , pressure factor : 1.3
FOB_LT_WGSR = REFCOST_reactor * ((12*(0.1016^1.5)*(production_scale^(2.5/3)))/REF_gas_cat_reactor)^NTH_GRXTOR*4.68*CEPCI/CEPCI_ref; % alloy factor : 3.6 , pressure factor : 1.3
% (d) Gas Separation Flash 
FOB_FLASH = REFCOST_column*((3.6576*(0.9144^1.5)*(production_scale^(2.5/3)))/REF_COL)^NTH_COL*CEPCI/CEPCI_ref;
EQ_SMR = FOB_NG_COMP+FOB_H2O_PUMP+FOB_SMR_reactor+FOB_FUR+FOB_HEX_1+FOB_HEX_2+FOB_HEX_3+FOB_HT_WGSR+FOB_LT_WGSR+FOB_FLASH;

% [2] CO2 Hydrogenation Reactor
% (a) Reactor Cost
adj_pressure = 2.3; % 100bar=10MPa - 2.3 

FOB_CO2_RXTOR = REFCOST_column*(adj_pressure*(3.6576*(0.9144^1.5)*conversion_factor*(production_scale^(2.5/3)))/REF_COL)^NTH_COL*CEPCI/CEPCI_ref;
% (b) Compressor Costs
FOB_CO2_COMP_1 = REFCOST_comp * ((SMR_EFFL_comp1_duty/3.6)/REF_electricity_comp)^NTH_COMP*2.5*CEPCI/CEPCI_ref;
FOB_CO2_COMP_2 = REFCOST_comp * ((SMR_EFFL_comp2_duty/3.6)/REF_electricity_comp)^NTH_COMP*2.5*CEPCI/CEPCI_ref;
FOB_SMREFFL_COOL_1 = REFCOST_HEX * (16.533*(production_scale^(2/3))/REF_HEX)^NTH_HEX*2*CEPCI/CEPCI_ref;
FOB_SMFEFFL_COOL_2 = REFCOST_HEX * (20.3185*(production_scale^(2/3))/REF_HEX)^NTH_HEX*2*CEPCI/CEPCI_ref;
FOB_CO2_Heat = REFCOST_HEX * (1.22008*(production_scale^(2/3))/REF_HEX)^NTH_HEX*2*CEPCI/CEPCI_ref;
FOB_TEAH2O_Heat = REFCOST_HEX * (14.7981*(production_scale^(2/3))/REF_HEX)^NTH_HEX*2*CEPCI/CEPCI_ref;
FOB_TEAH2O_PUMP = REFCOST_pump * ((TEA_H2O_pump_duty/3.6)/REF_electricity_pump)^NTH_PMP * CEPCI/CEPCI_ref;
% (c) MultiCompressor Costs
FOB_Multi_comp_1 = REFCOST_comp * ((SMR_FUR_multicomp1_electric*sum(SMR_FUR)/3.6)/REF_electricity_comp)^NTH_COMP*CEPCI/CEPCI_ref;
FOB_Multi_comp_2 = REFCOST_comp * ((SMR_FUR_multicomp2_electric*sum(SMR_FUR)/3.6)/REF_electricity_comp)^NTH_COMP*CEPCI/CEPCI_ref;
FOB_Multi_comp_3 = REFCOST_comp * ((SMR_FUR_multicomp3_electric*sum(SMR_FUR)/3.6)/REF_electricity_comp)^NTH_COMP*CEPCI/CEPCI_ref;
FOB_Multi_comp_4 = REFCOST_comp * ((SMR_FUR_multicomp4_electric*sum(SMR_FUR)/3.6)/REF_electricity_comp)^NTH_COMP*CEPCI/CEPCI_ref;
FOB_Multi_comp_5 = REFCOST_comp * ((SMR_FUR_multicomp5_electric*sum(SMR_FUR)/3.6)/REF_electricity_comp)^NTH_COMP*CEPCI/CEPCI_ref;
FOB_Multi_comp = FOB_Multi_comp_1+FOB_Multi_comp_2+FOB_Multi_comp_3+FOB_Multi_comp_4+FOB_Multi_comp_5;
EQ_CO2_Hydrogenation = FOB_CO2_RXTOR+FOB_CO2_COMP_1+FOB_CO2_COMP_2+FOB_SMREFFL_COOL_1+FOB_SMFEFFL_COOL_2+FOB_CO2_Heat+FOB_TEAH2O_Heat+FOB_TEAH2O_PUMP+FOB_Multi_comp;

% [3] COLUMN 1: TEA-Solvent Recovery
% (a) Distillation Column Estimation
FOB_Col1_PUMP = REFCOST_pump * ((column_1_vacuum_pump_elec_duty/3.6)/REF_electricity_pump)^NTH_PMP * CEPCI/CEPCI_ref;
FOB_Col_1 = REFCOST_Col_1*production_scale^NTH_Tray*CEPCI/CEPCI_ref;
FOB_Reb_1 = REFCOST_Reb_1*production_scale^NTH_Reb*CEPCI/CEPCI_ref;
EQ_TEA_solvent_recovery = FOB_Col1_PUMP + FOB_Col_1 +FOB_Reb_1;

% [4] COLUMN 2: NBIM-replacement
% (a) Distillation Column Estimation
FOB_Col_2 = REFCOST_Col_2*production_scale^NTH_Tray*CEPCI/CEPCI_ref;
FOB_Reb_2 = REFCOST_Reb_2*production_scale^NTH_Reb*CEPCI/CEPCI_ref;
FOB_Cond_2 = REFCOST_Cond_2*production_scale^NTH_HEX*CEPCI/CEPCI_ref;
EQ_NBIM_replacement=FOB_Col_2+FOB_Reb_2+FOB_Cond_2;

% [5] COLUMN 3: FA Purification 
% (a) Distillation Column Estimation
FOB_Col_3 = REFCOST_Col_3*production_scale^NTH_Tray*CEPCI/CEPCI_ref;
FOB_Reb_3 = REFCOST_Reb_3*production_scale^NTH_Reb*CEPCI/CEPCI_ref;
FOB_Cond_3 = REFCOST_Cond_3*production_scale^NTH_HEX*CEPCI/CEPCI_ref;
FOB_Col3_PUMP = REFCOST_pump * ((32.4889*production_scale)/REF_electricity_pump)^NTH_PMP * CEPCI/CEPCI_ref;
EQ_FA_purification = FOB_Col_3+FOB_Reb_3+FOB_Cond_3+FOB_Col3_PUMP;

%%% CAPITAL COST SUMMARY
IEC = (EQ_SMR + EQ_CO2_Hydrogenation + EQ_TEA_solvent_recovery + EQ_NBIM_replacement + EQ_FA_purification) * f_IEC;
DCC = f_DCC * IEC;
ICC = f_ICC * IEC;
FCI = DCC + ICC; %fixed capital investment
WC = f_WC * FCI;
Total_CAPEX = FCI + WC;
CRF = (DISCO*(DISCO+1)^Lifetime_plant)/(((DISCO+1)^Lifetime_plant)-1);
ANNUALIZED_CAPEX = CRF*Total_CAPEX;

%=========================== OPERATING EXPENSES ==========================%
% [1] Raw Material Costs
% (a) Natural Gas for SMR
COST_NG = MW_NG*NG_flowrate*Operating_hour/1000*REFCOST_NG;
% (b) Captured CO2 Cost
COST_CO2 = MW_CO2*CO2(1)*Operating_hour/1000*REFCOST_CO2;
% (c) Amine Cost (Make-Up
COST_TREA = MW_TREA*TEA_H2O(6)*(TREA_loss_rate/675)*Operating_hour/1000*REFCOST_TREA;
COST_NBIM = MW_NBIM*column_3_feed(9)*(NBIM_losss_rate/675)*Operating_hour/1000*REFCOST_NBIM;

% [2] Utility Costs
COST_elec = Elec_duty*Operating_hour/1000*REFCOST_electricity;   %(MJ * # of Hours/1000) >> GJ >> USD/GJ
COST_heat = Heat_duty*Operating_hour/1000*REFCOST_MP_steam;      %(MJ * # of Hours/1000) >> GJ >> USD/GJ
COST_NG_combustion = NG_combustion_usage*Operating_hour/1000*REFCOST_NG;

% [3] Fixed OPEX
COST_FO_Labor = REFCOST_Labor * W_FA * Operating_hour / 200000;
COST_FO_Maintenance = IEC * f_Maintenance;
COST_FO_admin = COST_FO_Labor * f_admin;
COST_FO_overhead = COST_FO_Labor * f_overhead;
COST_FO_Laboratory = IEC * f_Laboratory;

% OPERATING COST SUMMARY
TOTAL_OPEX_FO = COST_FO_Labor + COST_FO_Maintenance + COST_FO_admin + COST_FO_overhead + COST_FO_Laboratory;
TOTAL_OPEX_Variable = COST_NG + COST_CO2 + COST_TREA + COST_NBIM + COST_elec + COST_heat + COST_NG_combustion;

%============================== COST SUMMARY =============================%
ANNUAL_PROD_COST = ANNUALIZED_CAPEX + TOTAL_OPEX_FO + TOTAL_OPEX_Variable; % USD/yr
COGM = ANNUAL_PROD_COST/(W_FA * Operating_hour); % $/ton




%% CO2 LIFE CYCLE ASSESSMENT MODEL
% [1] Direct Emissions from Process
EM_CO2_direct = -Purge(1)*MW_CO2/(W_FA*1000)*GWI_CO2;

% [2] Indirect Emissions from Energy Consumption
% (a) Electricity
EM_Elec = Elec_duty/(W_FA*1000) * GWI_Elec; 
% (b) Heat Consumption
EM_MP = Heat_duty/(W_FA*1000) * GWI_MP_Steam;
% (c) NG Combustion
EM_NG_combustion = NG_combustion_usage/(W_FA*1000) * GWI_NG;

% [3] Indirect Emissions from Raw Materials
% (a) Raw Material Production
EM_NG = NG_flowrate*MW_NG/(W_FA*1000) * GWI_NG;
EM_H2O = (H2O_flowrate)*MW_H2O/(W_FA*1000) * GWI_Water;
% EM_H2O = (H2O_flowrate+C_Water)*MW_H2O/(W_FA*1000) * GWI_Water;

% [4] Direct Consumption
EM_CO2_con = CO2(1)*MW_CO2/(W_FA*1000)*GWI_CO2_captured;

% CO2 GWI SUMMARY (in ton CO2-eq)
SPECIFIC_GWI = EM_Elec+EM_MP+EM_NG_combustion+EM_NG+EM_CO2_con+EM_CO2_direct+EM_H2O;




%% EXPORT SUSTAINABILITY CRITERIA METRICS
% Criteria 1: Technoeconomics = Cost of Goods Manufactured
% Criteria 2: Global Warming Emissions = Specific Global Warming Impact
EvalMetrics = [COGM, SPECIFIC_GWI];  
end