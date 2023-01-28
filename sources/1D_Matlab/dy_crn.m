%===============================================================================
% CellML file:   Z:\home\DATA\tmp\courtemanche_ramirez_nattel_1998-43c29fd00a0d\courtemanche_ramirez_nattel_1998.cellml
% CellML model:  courtemanche_1998
% Date and time: 2/1/2012 at 9:27:33 at AM
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2012 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = dy_crn(time, Y)

% global Istim;
global total_current;

gk1mult=1.0-0.34;
gtomult=1.0-0.41;

%**********************************************
% Stimulus
%**********************************************
% % This uses previous datarow in stim-matrix for the value of stimuluscurrent at time t. In other words: interpolating to previous time-point.   
% st = stim( find(stim(:,1) <= time & stim(:,1) >= 0) , 2 ); % Find all (time, current)-data which t is smaller or equal to t
% if isempty(st)
%     i_st = 0;
% else
%     i_st = st(end); % use the last value (highest t)
% end

% i_st = Istim; % this Istim is now only in the driver/pde solver.

i_st = 0.0;

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [2.35e-112, 1.0, 0.9992, 1.367e-4, 7.755e-1, 9.996e-1, 9.649e-1, 9.775e-1, 2.908e-3, 1.013e-4, 1.488, 1.488, 1.39e2, 1.117e1, -81.18, 3.296e-5, 1.869e-2, 3.043e-2, 9.992e-1, 4.966e-3, 9.986e-1];

% YNames = {'u', 'v', 'w', 'd', 'f_Ca', 'f', 'h', 'j', 'm', 'Ca_i', 'Ca_rel', 'Ca_up', 'K_i', 'Na_i', 'V', 'xr', 'xs', 'oa', 'oi', 'ua', 'ui'};
% YUnits = {'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millivolt', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless'};
% YComponents = {'Ca_release_current_from_JSR_u_gate', 'Ca_release_current_from_JSR_v_gate', 'Ca_release_current_from_JSR_w_gate', 'L_type_Ca_channel_d_gate', 'L_type_Ca_channel_f_Ca_gate', 'L_type_Ca_channel_f_gate', 'fast_sodium_current_h_gate', 'fast_sodium_current_j_gate', 'fast_sodium_current_m_gate', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'membrane', 'rapid_delayed_rectifier_K_current_xr_gate', 'slow_delayed_rectifier_K_current_xs_gate', 'transient_outward_K_current_oa_gate', 'transient_outward_K_current_oi_gate', 'ultrarapid_delayed_rectifier_K_current_ua_gate', 'ultrarapid_delayed_rectifier_K_current_ui_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: u (dimensionless) (in Ca_release_current_from_JSR_u_gate)
% 2: v (dimensionless) (in Ca_release_current_from_JSR_v_gate)
% 3: w (dimensionless) (in Ca_release_current_from_JSR_w_gate)
% 4: d (dimensionless) (in L_type_Ca_channel_d_gate)
% 5: f_Ca (dimensionless) (in L_type_Ca_channel_f_Ca_gate)
% 6: f (dimensionless) (in L_type_Ca_channel_f_gate)
% 7: h (dimensionless) (in fast_sodium_current_h_gate)
% 8: j (dimensionless) (in fast_sodium_current_j_gate)
% 9: m (dimensionless) (in fast_sodium_current_m_gate)
% 10: Ca_i (millimolar) (in intracellular_ion_concentrations)
% 11: Ca_rel (millimolar) (in intracellular_ion_concentrations)
% 12: Ca_up (millimolar) (in intracellular_ion_concentrations)
% 13: K_i (millimolar) (in intracellular_ion_concentrations)
% 14: Na_i (millimolar) (in intracellular_ion_concentrations)
% 15: V (millivolt) (in membrane)
% 16: xr (dimensionless) (in rapid_delayed_rectifier_K_current_xr_gate)
% 17: xs (dimensionless) (in slow_delayed_rectifier_K_current_xs_gate)
% 18: oa (dimensionless) (in transient_outward_K_current_oa_gate)
% 19: oi (dimensionless) (in transient_outward_K_current_oi_gate)
% 20: ua (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
% 21: ui (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ui_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

CMDN_max = 0.05;   % millimolar (in Ca_buffers)
CSQN_max = 10.0;   % millimolar (in Ca_buffers)
Km_CMDN = 0.00238;   % millimolar (in Ca_buffers)
Km_CSQN = 0.8;   % millimolar (in Ca_buffers)
Km_TRPN = 0.0005;   % millimolar (in Ca_buffers)
TRPN_max = 0.07;   % millimolar (in Ca_buffers)
Ca_up_max = 15.0;   % millimolar (in Ca_leak_current_by_the_NSR)
K_rel = 30.0;   % per_millisecond (in Ca_release_current_from_JSR)
I_up_max = 0.005;   % millimolar_per_millisecond (in Ca_uptake_current_by_the_NSR)
K_up = 0.00092;   % millimolar (in Ca_uptake_current_by_the_NSR)
g_Ca_L = 0.12375;   % nanoS_per_picoF (in L_type_Ca_channel)
I_NaCa_max = 1600.0;   % picoA_per_picoF (in Na_Ca_exchanger_current)
K_mCa = 1.38;   % millimolar (in Na_Ca_exchanger_current)
K_mNa = 87.5;   % millimolar (in Na_Ca_exchanger_current)
K_sat = 0.1;   % dimensionless (in Na_Ca_exchanger_current)
gamma = 0.35;   % dimensionless (in Na_Ca_exchanger_current)
g_B_Ca = 0.001131;   % nanoS_per_picoF (in background_currents)
g_B_K = 0.0;   % nanoS_per_picoF (in background_currents)
g_B_Na = 0.0006744375;   % nanoS_per_picoF (in background_currents)
g_Na = 7.8;   % nanoS_per_picoF (in fast_sodium_current)
V_cell = 20100.0;   % micrometre_3 (in intracellular_ion_concentrations)
Cm = 100.0;   % picoF (in membrane)
F = 96.4867;   % coulomb_per_millimole (in membrane)
R = 8.3143;   % joule_per_mole_kelvin (in membrane)
T = 310.0;   % kelvin (in membrane)
% stim_amplitude = -2000.0;   % picoA (in membrane)
% stim_duration = 2.0;   % millisecond (in membrane)
% stim_end = 50000.0;   % millisecond (in membrane)
% stim_period = 1000.0;   % millisecond (in membrane)
% stim_start = 50.0;   % millisecond (in membrane)
g_Kr = 0.029411765;   % nanoS_per_picoF (in rapid_delayed_rectifier_K_current)
i_CaP_max = 0.275;   % picoA_per_picoF (in sarcolemmal_calcium_pump_current)
g_Ks = 0.12941176;   % nanoS_per_picoF (in slow_delayed_rectifier_K_current)
Km_K_o = 1.5;   % millimolar (in sodium_potassium_pump)
Km_Na_i = 10.0;   % millimolar (in sodium_potassium_pump)
i_NaK_max = 0.59933874;   % picoA_per_picoF (in sodium_potassium_pump)
Ca_o = 1.8;   % millimolar (in standard_ionic_concentrations)
K_o = 5.4;   % millimolar (in standard_ionic_concentrations)
Na_o = 140.0;   % millimolar (in standard_ionic_concentrations)
g_K1 = gk1mult*0.09;   % nanoS_per_picoF (in time_independent_potassium_current)
tau_tr = 180.0;   % millisecond (in transfer_current_from_NSR_to_JSR)
K_Q10 = 3.0;   % dimensionless (in transient_outward_K_current)
g_to = gtomult*0.1652;   % nanoS_per_picoF (in transient_outward_K_current)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% Ca_CMDN (millimolar) (in Ca_buffers)
% Ca_CSQN (millimolar) (in Ca_buffers)
% Ca_TRPN (millimolar) (in Ca_buffers)
% i_up_leak (millimolar_per_millisecond) (in Ca_leak_current_by_the_NSR)
% tau_u (millisecond) (in Ca_release_current_from_JSR_u_gate)
% u_infinity (dimensionless) (in Ca_release_current_from_JSR_u_gate)
% tau_v (millisecond) (in Ca_release_current_from_JSR_v_gate)
% v_infinity (dimensionless) (in Ca_release_current_from_JSR_v_gate)
% tau_w (millisecond) (in Ca_release_current_from_JSR_w_gate)
% w_infinity (dimensionless) (in Ca_release_current_from_JSR_w_gate)
% Fn (dimensionless) (in Ca_release_current_from_JSR)
% i_rel (millimolar_per_millisecond) (in Ca_release_current_from_JSR)
% i_up (millimolar_per_millisecond) (in Ca_uptake_current_by_the_NSR)
% d_infinity (dimensionless) (in L_type_Ca_channel_d_gate)
% tau_d (millisecond) (in L_type_Ca_channel_d_gate)
% f_Ca_infinity (dimensionless) (in L_type_Ca_channel_f_Ca_gate)
% tau_f_Ca (millisecond) (in L_type_Ca_channel_f_Ca_gate)
% f_infinity (dimensionless) (in L_type_Ca_channel_f_gate)
% tau_f (millisecond) (in L_type_Ca_channel_f_gate)
% i_Ca_L (picoA) (in L_type_Ca_channel)
% i_NaCa (picoA) (in Na_Ca_exchanger_current)
% E_Ca (millivolt) (in background_currents)
% i_B_Ca (picoA) (in background_currents)
% i_B_K (picoA) (in background_currents)
% i_B_Na (picoA) (in background_currents)
% alpha_h (per_millisecond) (in fast_sodium_current_h_gate)
% beta_h (per_millisecond) (in fast_sodium_current_h_gate)
% h_inf (dimensionless) (in fast_sodium_current_h_gate)
% tau_h (millisecond) (in fast_sodium_current_h_gate)
% alpha_j (per_millisecond) (in fast_sodium_current_j_gate)
% beta_j (per_millisecond) (in fast_sodium_current_j_gate)
% j_inf (dimensionless) (in fast_sodium_current_j_gate)
% tau_j (millisecond) (in fast_sodium_current_j_gate)
% alpha_m (per_millisecond) (in fast_sodium_current_m_gate)
% beta_m (per_millisecond) (in fast_sodium_current_m_gate)
% m_inf (dimensionless) (in fast_sodium_current_m_gate)
% tau_m (millisecond) (in fast_sodium_current_m_gate)
% E_Na (millivolt) (in fast_sodium_current)
% i_Na (picoA) (in fast_sodium_current)
% B1 (millimolar_per_millisecond) (in intracellular_ion_concentrations)
% B2 (dimensionless) (in intracellular_ion_concentrations)
% V_i (micrometre_3) (in intracellular_ion_concentrations)
% V_rel (micrometre_3) (in intracellular_ion_concentrations)
% V_up (micrometre_3) (in intracellular_ion_concentrations)
% i_st (picoA) (in membrane)
% alpha_xr (per_millisecond) (in rapid_delayed_rectifier_K_current_xr_gate)
% beta_xr (per_millisecond) (in rapid_delayed_rectifier_K_current_xr_gate)
% tau_xr (millisecond) (in rapid_delayed_rectifier_K_current_xr_gate)
% xr_infinity (dimensionless) (in rapid_delayed_rectifier_K_current_xr_gate)
% i_Kr (picoA) (in rapid_delayed_rectifier_K_current)
% i_CaP (picoA) (in sarcolemmal_calcium_pump_current)
% alpha_xs (per_millisecond) (in slow_delayed_rectifier_K_current_xs_gate)
% beta_xs (per_millisecond) (in slow_delayed_rectifier_K_current_xs_gate)
% tau_xs (millisecond) (in slow_delayed_rectifier_K_current_xs_gate)
% xs_infinity (dimensionless) (in slow_delayed_rectifier_K_current_xs_gate)
% i_Ks (picoA) (in slow_delayed_rectifier_K_current)
% f_NaK (dimensionless) (in sodium_potassium_pump)
% i_NaK (picoA) (in sodium_potassium_pump)
% sigma (dimensionless) (in sodium_potassium_pump)
% E_K (millivolt) (in time_independent_potassium_current)
% i_K1 (picoA) (in time_independent_potassium_current)
% i_tr (millimolar_per_millisecond) (in transfer_current_from_NSR_to_JSR)
% alpha_oa (per_millisecond) (in transient_outward_K_current_oa_gate)
% beta_oa (per_millisecond) (in transient_outward_K_current_oa_gate)
% oa_infinity (dimensionless) (in transient_outward_K_current_oa_gate)
% tau_oa (millisecond) (in transient_outward_K_current_oa_gate)
% alpha_oi (per_millisecond) (in transient_outward_K_current_oi_gate)
% beta_oi (per_millisecond) (in transient_outward_K_current_oi_gate)
% oi_infinity (dimensionless) (in transient_outward_K_current_oi_gate)
% tau_oi (millisecond) (in transient_outward_K_current_oi_gate)
% i_to (picoA) (in transient_outward_K_current)
% alpha_ua (per_millisecond) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
% beta_ua (per_millisecond) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
% tau_ua (millisecond) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
% ua_infinity (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ua_gate)
% alpha_ui (per_millisecond) (in ultrarapid_delayed_rectifier_K_current_ui_gate)
% beta_ui (per_millisecond) (in ultrarapid_delayed_rectifier_K_current_ui_gate)
% tau_ui (millisecond) (in ultrarapid_delayed_rectifier_K_current_ui_gate)
% ui_infinity (dimensionless) (in ultrarapid_delayed_rectifier_K_current_ui_gate)
% g_Kur (nanoS_per_picoF) (in ultrarapid_delayed_rectifier_K_current)
% i_Kur (picoA) (in ultrarapid_delayed_rectifier_K_current)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (millisecond)

Ca_CMDN = CMDN_max*Y(10)/(Y(10)+Km_CMDN);
Ca_TRPN = TRPN_max*Y(10)/(Y(10)+Km_TRPN);
Ca_CSQN = CSQN_max*Y(11)/(Y(11)+Km_CSQN);
i_up_leak = I_up_max*Y(12)/Ca_up_max;
V_rel = 0.0048*V_cell;
i_rel = K_rel*Y(1)^2.0*Y(2)*Y(3)*(Y(11)-Y(10));
i_Ca_L = Cm*g_Ca_L*Y(4)*Y(6)*Y(5)*(Y(15)-65.0);
i_NaCa = Cm*I_NaCa_max*(exp(gamma*F*Y(15)/(R*T))*Y(14)^3.0*Ca_o-exp((gamma-1.0)*F*Y(15)/(R*T))*Na_o^3.0*Y(10))/((K_mNa^3.0+Na_o^3.0)*(K_mCa+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*Y(15)*F/(R*T))));
Fn = 1.0e3*(1.0e-15*V_rel*i_rel-1.0e-15/(2.0*F)*(0.5*i_Ca_L-0.2*i_NaCa));
tau_u = 8.0;
u_infinity = (1.0+exp(-(Fn-3.4175e-13)/13.67e-16))^-1.0;
dY(1, 1) = (u_infinity-Y(1))/tau_u;
tau_v = 1.91+2.09*(1.0+exp(-(Fn-3.4175e-13)/13.67e-16))^-1.0;
v_infinity = 1.0-(1.0+exp(-(Fn-6.835e-14)/13.67e-16))^-1.0;
dY(2, 1) = (v_infinity-Y(2))/tau_v;

if (abs(Y(15)-7.9) < 1.0e-10)
   tau_w = 6.0*0.2/1.3;
else
   tau_w = 6.0*(1.0-exp(-(Y(15)-7.9)/5.0))/((1.0+0.3*exp(-(Y(15)-7.9)/5.0))*1.0*(Y(15)-7.9));
end;

w_infinity = 1.0-(1.0+exp(-(Y(15)-40.0)/17.0))^-1.0;
dY(3, 1) = (w_infinity-Y(3))/tau_w;
i_up = I_up_max/(1.0+K_up/Y(10));
d_infinity = (1.0+exp((Y(15)+10.0)/-8.0))^-1.0;

if (abs(Y(15)+10.0) < 1.0e-10)
   tau_d = 4.579/(1.0+exp((Y(15)+10.0)/-6.24));
else
   tau_d = (1.0-exp((Y(15)+10.0)/-6.24))/(0.035*(Y(15)+10.0)*(1.0+exp((Y(15)+10.0)/-6.24)));
end;

dY(4, 1) = (d_infinity-Y(4))/tau_d;
f_Ca_infinity = (1.0+Y(10)/0.00035)^-1.0;
tau_f_Ca = 2.0;
dY(5, 1) = (f_Ca_infinity-Y(5))/tau_f_Ca;
f_infinity = exp(-(Y(15)+28.0)/6.9)/(1.0+exp(-(Y(15)+28.0)/6.9));
tau_f = 9.0*(0.0197*exp(-0.0337^2.0*(Y(15)+10.0)^2.0)+0.02)^-1.0;
dY(6, 1) = (f_infinity-Y(6))/tau_f;
E_Ca = R*T/(2.0*F)*log(Ca_o/Y(10));
E_Na = R*T/F*log(Na_o/Y(14));
i_B_Na = Cm*g_B_Na*(Y(15)-E_Na);
i_B_Ca = Cm*g_B_Ca*(Y(15)-E_Ca);
E_K = R*T/F*log(K_o/Y(13));
i_B_K = Cm*g_B_K*(Y(15)-E_K);
i_Na = Cm*g_Na*Y(9)^3.0*Y(7)*Y(8)*(Y(15)-E_Na);

if (Y(15) < -40.0)
   alpha_h = 0.135*exp((Y(15)+80.0)/-6.8);
else
   alpha_h = 0.0;
end;

if (Y(15) < -40.0)
   beta_h = 3.56*exp(0.079*Y(15))+3.1e5*exp(0.35*Y(15));
else
   beta_h = 1.0/(0.13*(1.0+exp((Y(15)+10.66)/-11.1)));
end;

h_inf = alpha_h/(alpha_h+beta_h);
tau_h = 1.0/(alpha_h+beta_h);
dY(7, 1) = (h_inf-Y(7))/tau_h;

if (Y(15) < -40.0)
   alpha_j = (-1.2714e5*exp(0.2444*Y(15))-3.474e-5*exp(-0.04391*Y(15)))*(Y(15)+37.78)/(1.0+exp(0.311*(Y(15)+79.23)));
else
   alpha_j = 0.0;
end;

if (Y(15) < -40.0)
   beta_j = 0.1212*exp(-0.01052*Y(15))/(1.0+exp(-0.1378*(Y(15)+40.14)));
else
   beta_j = 0.3*exp(-2.535e-7*Y(15))/(1.0+exp(-0.1*(Y(15)+32.0)));
end;

j_inf = alpha_j/(alpha_j+beta_j);
tau_j = 1.0/(alpha_j+beta_j);
dY(8, 1) = (j_inf-Y(8))/tau_j;

if (Y(15) == -47.13)
   alpha_m = 3.2;
else
   alpha_m = 0.32*(Y(15)+47.13)/(1.0-exp(-0.1*(Y(15)+47.13)));
end;

beta_m = 0.08*exp(-Y(15)/11.0);
m_inf = alpha_m/(alpha_m+beta_m);
tau_m = 1.0/(alpha_m+beta_m);
dY(9, 1) = (m_inf-Y(9))/tau_m;
V_i = V_cell*0.68;
V_up = 0.0552*V_cell;
sigma = 1.0/7.0*(exp(Na_o/67.3)-1.0);
f_NaK = (1.0+0.1245*exp(-0.1*F*Y(15)/(R*T))+0.0365*sigma*exp(-F*Y(15)/(R*T)))^-1.0;
i_NaK = Cm*i_NaK_max*f_NaK*1.0/(1.0+(Km_Na_i/Y(14))^1.5)*K_o/(K_o+Km_K_o);
dY(14, 1) = (-3.0*i_NaK-(3.0*i_NaCa+i_B_Na+i_Na))/(V_i*F);
i_K1 = Cm*g_K1*(Y(15)-E_K)/(1.0+exp(0.07*(Y(15)+80.0)));
i_to = Cm*g_to*Y(18)^3.0*Y(19)*(Y(15)-E_K);
g_Kur = 0.005+0.05/(1.0+exp((Y(15)-15.0)/-13.0));
i_Kur = Cm*g_Kur*Y(20)^3.0*Y(21)*(Y(15)-E_K);
i_Kr = Cm*g_Kr*Y(16)*(Y(15)-E_K)/(1.0+exp((Y(15)+15.0)/22.4));
i_Ks = Cm*g_Ks*Y(17)^2.0*(Y(15)-E_K);
dY(13, 1) = (2.0*i_NaK-(i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))/(V_i*F);
i_CaP = Cm*i_CaP_max*Y(10)/(0.0005+Y(10));
B1 = (2.0*i_NaCa-(i_CaP+i_Ca_L+i_B_Ca))/(2.0*V_i*F)+(V_up*(i_up_leak-i_up)+i_rel*V_rel)/V_i;
B2 = 1.0+TRPN_max*Km_TRPN/(Y(10)+Km_TRPN)^2.0+CMDN_max*Km_CMDN/(Y(10)+Km_CMDN)^2.0;
dY(10, 1) = B1/B2;
i_tr = (Y(12)-Y(11))/tau_tr;
dY(12, 1) = i_up-(i_up_leak+i_tr*V_rel/V_up);
dY(11, 1) = (i_tr-i_rel)*(1.0+CSQN_max*Km_CSQN/(Y(11)+Km_CSQN)^2.0)^-1.0;

% if ((time >= stim_start) && (time <= stim_end) && (time-stim_start-floor((time-stim_start)/stim_period)*stim_period <= stim_duration))
%    i_st = stim_amplitude;
% else
%    i_st = 0.0;
% end;

% dY(15, 1) = -(i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+i_st)/Cm;

dY(15,1) = 0.0; % I am doing this calculation in the PDE.
total_current = -(i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+i_st)/Cm;

if (abs(Y(15)+14.1) < 1.0e-10)
   alpha_xr = 0.0015;
else
   alpha_xr = 0.0003*(Y(15)+14.1)/(1.0-exp((Y(15)+14.1)/-5.0));
end;

if (abs(Y(15)-3.3328) < 1.0e-10)
   beta_xr = 3.7836118e-4;
else
   beta_xr = 0.000073898*(Y(15)-3.3328)/(exp((Y(15)-3.3328)/5.1237)-1.0);
end;

tau_xr = (alpha_xr+beta_xr)^-1.0;
xr_infinity = (1.0+exp((Y(15)+14.1)/-6.5))^-1.0;
dY(16, 1) = (xr_infinity-Y(16))/tau_xr;

if (abs(Y(15)-19.9) < 1.0e-10)
   alpha_xs = 0.00068;
else
   alpha_xs = 0.00004*(Y(15)-19.9)/(1.0-exp((Y(15)-19.9)/-17.0));
end;

if (abs(Y(15)-19.9) < 1.0e-10)
   beta_xs = 0.000315;
else
   beta_xs = 0.000035*(Y(15)-19.9)/(exp((Y(15)-19.9)/9.0)-1.0);
end;

tau_xs = 0.5*(alpha_xs+beta_xs)^-1.0;
xs_infinity = (1.0+exp((Y(15)-19.9)/-12.7))^-0.5;
dY(17, 1) = (xs_infinity-Y(17))/tau_xs;
alpha_oa = 0.65*(exp((Y(15)-(-10.0))/-8.5)+exp((Y(15)-(-10.0)-40.0)/-59.0))^-1.0;
beta_oa = 0.65*(2.5+exp((Y(15)-(-10.0)+72.0)/17.0))^-1.0;
tau_oa = (alpha_oa+beta_oa)^-1.0/K_Q10;
oa_infinity = (1.0+exp((Y(15)-(-10.0)+10.47)/-17.54))^-1.0;
dY(18, 1) = (oa_infinity-Y(18))/tau_oa;
alpha_oi = (18.53+1.0*exp((Y(15)-(-10.0)+103.7)/10.95))^-1.0;
beta_oi = (35.56+1.0*exp((Y(15)-(-10.0)-8.74)/-7.44))^-1.0;
tau_oi = (alpha_oi+beta_oi)^-1.0/K_Q10;
oi_infinity = (1.0+exp((Y(15)-(-10.0)+33.1)/5.3))^-1.0;
dY(19, 1) = (oi_infinity-Y(19))/tau_oi;
alpha_ua = 0.65*(exp((Y(15)-(-10.0))/-8.5)+exp((Y(15)-(-10.0)-40.0)/-59.0))^-1.0;
beta_ua = 0.65*(2.5+exp((Y(15)-(-10.0)+72.0)/17.0))^-1.0;
tau_ua = (alpha_ua+beta_ua)^-1.0/K_Q10;
ua_infinity = (1.0+exp((Y(15)-(-10.0)+20.3)/-9.6))^-1.0;
dY(20, 1) = (ua_infinity-Y(20))/tau_ua;
alpha_ui = (21.0+1.0*exp((Y(15)-(-10.0)-195.0)/-28.0))^-1.0;
beta_ui = 1.0/exp((Y(15)-(-10.0)-168.0)/-16.0);
tau_ui = (alpha_ui+beta_ui)^-1.0/K_Q10;
ui_infinity = (1.0+exp((Y(15)-(-10.0)-109.45)/27.48))^-1.0;
dY(21, 1) = (ui_infinity-Y(21))/tau_ui;

%===============================================================================
% End of file
%===============================================================================
