clc; clear; 
set(0,'defaultfigurecolor',[1 1 1])
set(0, 'defaultfigureposition', [100, 100, 600, 600])
%% UPDATE HERE PLEASE

%% =========================================================
%  =================== USER INPUTS =========================
%  =========================================================

% --- Operating Conditions --- %
m_flowrate      = 0.0013;
reactor_length  = 10;
P0              = 1.6;
T0              = 666;

% --- Inlet Mole Fractions --- %
y0_CH3OH = 0.11;
y0_O2    = 0.06;
y0_H2O   = 0.02;
y0_HCHO  = 0.00;
y0_CO    = 0.00;
y0_N2    = 1 - y0_CH3OH - y0_O2 - y0_H2O - y0_HCHO - y0_CO;

% --- Catalyst/Reactor Parameters --- %
alpha         = 0.44;
rho_cat       = 1500;
bed_voidage   = 0.4;
D_p           = 7.04e-3;
tube_diameter = 20e-3;

% --- Reaction Enthalpies --- %
Hr1 = -159e3;
Hr2 = -200e3;

% --- Solver Options --- %
Kinetic_type    = "LHHW"; % "LHHW", "Simple Power"
Pressure_drop   = "Yes"; % "Yes", "No"
Temp_dependence = "Yes"; % "Yes", "No"

% --- Event Limits --- %
T_max        = 673;
T_max_active = "Yes"; % "Yes", "No"
P_min        = 1.1/1.01325;
P_min_active = "Yes"; % "Yes", "No"

% --- Cooling Strategy --- %
cooling_type = "Isothermal"; % "Isothermal", "Cold Shot", "Interstage"

% === Cold Shot ===
split_ratio  = 0.5;
cs_n_stages  = 4;

% === Interstage ===
isc_no_stages = 8;
T_subcool     = 610;

% --- Scale-up --- %
goal_fr = 38.6;

%  =========================================================
%  ================== END OF USER INPUTS ===================
%  =========================================================

inlet         = [y0_CH3OH, y0_O2, y0_H2O, y0_HCHO, y0_CO, y0_N2];
reactor_param = [alpha, rho_cat, bed_voidage, D_p, tube_diameter, Hr1, Hr2];
option        = [Kinetic_type, Pressure_drop, Temp_dependence];
option_limits = {T_max, T_max_active, P_min, P_min_active};

% --- Build feed stream and solve ---
[feed_stream, ~, F0_CH3OH, F0_HCHO] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);

if cooling_type == "Cold Shot"
    Temp_dependence = "Yes";
    options = [Kinetic_type Pressure_drop Temp_dependence];
    [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, P, T, extent1, extent2, ...
     F0_CH3OH, F0_HCHO, stage_boundaries] = ...
        coldshot_reactor_multistage(m_flowrate, P0, T0, reactor_length, split_ratio, cs_n_stages, ...
                                    option_limits, options, inlet, reactor_param);
    Q_dot_per_length = [];
    Q_total = 0;

elseif cooling_type == "Interstage"
    Temp_dependence = "Yes";
    options = [Kinetic_type Pressure_drop Temp_dependence];
    [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, P, T, extent1, extent2, F0_CH3OH, F0_HCHO, stage_boundaries, Q_total] = ...
          reactor_interstage_cooling(m_flowrate, isc_no_stages, T_subcool, P0, T0, option_limits, options, inlet, reactor_param);
    Q_dot_per_length = [];

elseif cooling_type == "None"
    Temp_dependence = "Yes";
    options = [Kinetic_type Pressure_drop Temp_dependence];
    [feed_stream, ~, F0_CH3OH, F0_HCHO] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);
    [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
    solve_reactor(feed_stream, reactor_length, P0, T0, option_limits, options, reactor_param);
    stage_boundaries = [];
    Q_dot_per_length = [];
    Q_total = 0;

elseif cooling_type == "Isothermal"
    Temp_dependence = "No";
    options = [Kinetic_type Pressure_drop Temp_dependence];
    [feed_stream, ~, F0_CH3OH, F0_HCHO] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);
    [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
    solve_reactor(feed_stream, reactor_length, P0, T0, option_limits, options, reactor_param);
    [Q_dot_per_length, Q_total] = calc_isothermal_duty(z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, P, T, reactor_param);
    stage_boundaries = [];
end
plot_reactor_results(z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, ...
                     P, T, extent1, extent2, F0_CH3OH, F0_HCHO, ...
                     [], Q_dot_per_length, reactor_param);

% --- Catalyst weight axis ---
dW_dz = rho_cat * (1 - bed_voidage) * pi*(tube_diameter/2)^2;
W     = z_points * dW_dz;

num_tubes = calc_num_tubes(F_HCHO, goal_fr);
fprintf("Total number of tubes required: %.0f \n", num_tubes)
fprintf('Total isothermal cooling duty: %.2f W (%.2f MW)\n', Q_total*num_tubes, Q_total*num_tubes/1000000);


%% Calling for Graphs (Step 2 onwards)
plot_r1_simplepower()
plot_r2_simplepower()
step3_graphs(reactor_param)
step4_graphs(reactor_param)
step5_graphs(reactor_param)
plot_extent_vs_catalyst_weight(reactor_param, inlet)
plot_rate_profiles_lhhw_430(reactor_param, inlet)
plot_temperature_profile_6cases(reactor_param, inlet)
plot_rate_vs_temperature_comparison(reactor_param, inlet, P0)

sensitivity_analysis()
heatmap_conversion_tubes(reactor_param, inlet, P0, reactor_length, goal_fr)
pareto_conversion_tubes(reactor_param, inlet, P0, reactor_length, goal_fr)
reactor_comparision_graph(m_flowrate,P0,inlet,reactor_param,reactor_length)
fprintf('\n=== Outlet Summary ===\n')
fprintf('X_CH3OH       = %.4f\n',   (F_CH3OH(1)-F_CH3OH(end))/F_CH3OH(1))
fprintf('Y_HCHO        = %.4f\n', F_HCHO(end)/F_CH3OH(1))
fprintf('F_HCHO out    = %.4e mol/s\n', F_HCHO(end))
fprintf('F_CO out      = %.4e mol/s\n', F_CO(end))
fprintf('P out         = %.4f atm\n',   P(end))
fprintf('T out         = %.2f K\n',     T(end))
fprintf('W total       = %.4f kg\n',    W(end))



%% =========================================================
%% FUNCTIONS
%% =========================================================

%% STEP 2 - Reaction 1 Graphs Function
function plot_r1_simplepower()
    R = 8.314; % J/mol/K
    fprintf("\n")
    fprintf("=============================\n")
    fprintf("Reaction 1 - Simple Power Law Determination \n")
    fprintf("=============================\n")
    % --- Activation Energy parameters (converted kJ/mol -> J/mol) --- %
    A_kCH3OH = 1.50e7;   Ea_kCH3OH = 86.00e3;
    A_KCH3OH = 2.60e-4;  Ea_KCH3OH = -56.78e3;
    A_KO2    = 1.42e-5;  Ea_KO2    = -60.32e3;
    A_KH2O   = 5.50e-7;  Ea_KH2O   = -86.45e3;
    A_kCO    = 3.50e2;   Ea_kCO    = 46.00e3;
    
    % --- Rate constants at 250C --- %
    T1 = 250 + 273.15;
    kCH3OH  = A_kCH3OH * exp(-Ea_kCH3OH/(R*T1));
    KCH3OH  = A_KCH3OH * exp(-Ea_KCH3OH/(R*T1));
    KO2     = A_KO2    * exp(-Ea_KO2/(R*T1));
    KH2O    = A_KH2O   * exp(-Ea_KH2O/(R*T1));
    kCO     = A_kCO    * exp(-Ea_kCO/(R*T1));
    
    % --- Rate constants at 270C --- %
    T2 = 270 + 273.15;
    kCH3OH_270 = A_kCH3OH * exp(-Ea_kCH3OH/(R*T2));
    KCH3OH_270 = A_KCH3OH * exp(-Ea_KCH3OH/(R*T2));
    KO2_270    = A_KO2    * exp(-Ea_KO2/(R*T2));
    KH2O_270   = A_KH2O   * exp(-Ea_KH2O/(R*T2));
    kCO_270    = A_kCO    * exp(-Ea_kCO/(R*T2));
    
    alpha = 0.44;
    
    %==========================================================================
    % DATA
    %==========================================================================
    
    % --- 220C datasets --- 
    H2O_c1_220 = [0.002 0.005 0.010] * 100;        
    R_HCHO_c1_220 = [2.025 1.396 0.920] * 10^-3;   % convert mmol to mol
    
    O2_c2_220 = [0.005 0.010 0.020] * 100;
    R_HCHO_c2_220 = [0.842 0.920 0.984] * 10^-3;
    
    CH3OH_c3_220 = [0.002 0.005 0.010] * 100;
    R_HCHO_c3_220 = [0.226 0.520 0.920] * 10^-3;
    
    % --- 260C datasets --- 
    H2O_d1_260 = [0.002 0.005 0.010] * 100;
    R_HCHO_d1_260 = [5.535 4.557 3.520] * 10^-3;
    
    O2_d2_260 = [0.005 0.010 0.020] * 100;
    R_HCHO_d2_260 = [2.954 3.520 4.073] * 10^-3;
    
    CH3OH_d3_260 = [0.002 0.005 0.010] * 100;
    R_HCHO_d3_260 = [0.895 2.030 3.520] * 10^-3;

    % --- 250C datasets ---
    H2O_a1_250  = [1.5961 2.2488 2.4685 2.7851 3.3150];
    R_HCHO_a1_250 = 10^-3.*[9.6658 7.1342 6.4438 5.9836 7.1342];
    
    O2_a2_250   = [1.0363 2.1111 3.2218 5.1195 7.0852 10.2314];
    R_HCHO_a2_250 = 10^-3.*[5.2358 6.9431 7.2846 7.0569 7.2087 6.4499];
    
    CH3OH_a3_250 = [1.0554 2.3443 3.6676 5.6237 5.6607 6.3680 6.9861 7.8738 9.6124];
    R_HCHO_a3_250 = 10^-3.*[3.5956 5.6612 6.3115 5.9672 6.5792 7.1913 8.8743 8.8743 8.4536];
    
    % --- 270C datasets ---
    CH3OH_b1_270 = [1.0287 2.3103 3.6594 5.6661 6.3575 6.9646 7.8752 9.5953];
    R_HCHO_b1_270 = 10^-3.*[5.4221 10.4288 13.0400 13.1025 15.3457 17.7447 18.6869 17.7584];
    
    O2_b2_270   = [1.0745 2.1317 3.2409 5.1300 7.0884 10.2600];
    R_HCHO_b2_270 = 10^-3.*[8.0262 12.3003 12.5638 12.4426 13.3567 13.2410];
    
    H2O_b3_270  = [0.8894 1.5865 2.2436 2.3317 2.4679 2.7724 3.3093 3.8301];
    R_HCHO_b3_270 = 10^-3.*[28.3815 17.8622 15.3651 16.4867 13.1957 16.0785 13.4228 14.0562];
    
    % --- Combined 220C and 260C arrays ---
    CH3OH_220 = [0.002 0.005 0.010 0.010 0.010 0.010 0.010] * 100;
    O2_220    = [0.010 0.010 0.010 0.005 0.020 0.010 0.010] * 100;
    H2O_220   = [0.010 0.010 0.010 0.010 0.010 0.005 0.002] * 100;
    rate_220  = [0.226 0.520 0.920 0.842 0.984 1.396 2.025] * 10^-3;
    
    CH3OH_260 = [0.002 0.005 0.010 0.010 0.010 0.010 0.010] * 100;
    O2_260    = [0.010 0.010 0.010 0.005 0.020 0.010 0.010] * 100;
    H2O_260   = [0.010 0.010 0.010 0.010 0.010 0.005 0.002] * 100;
    rate_260  = [0.895 2.030 3.520 2.954 4.073 4.557 5.535] * 10^-3;

    CH3OH_250 = [6 6 6 6 6 6 1.0554 2.3443 3.6676 5.6237 5.6607 6.3680 6.9861 7.8738 9.6124 6 6 6 6 6];
    O2_250    = [1.0363 2.1111 3.2218 5.1195 7.0852 10.2314 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
    H2O_250   = [2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 1.5961 2.2488 2.4685 2.7851 3.3150];
    rate_250  = 10^-3.*[5.2358 6.9431 7.2846 7.0569 7.2087 6.4499 3.5956 5.6612 6.3115 5.9672 6.5792 7.1913 8.8743 8.8743 8.4536 9.6658 7.1342 6.4438 5.9836 7.1342];
    
    CH3OH_270 = [6 6 6 6 6 6 1.0287 2.3103 3.6594 5.6661 6.3575 6.9646 7.8752 9.5953 6 6 6 6 6 6 6 6];
    O2_270    = [1.0745 2.1317 3.2409 5.1300 7.0884 10.2600 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
    H2O_270   = [2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 0.8894 1.5865 2.2436 2.3317 2.4679 2.7724 3.3093 3.8301];
    rate_270  = 10^-3.*[8.0262 12.3003 12.5638 12.4426 13.3567 13.2410 5.4221 10.4288 13.0400 13.1025 15.3457 17.7447 18.6869 17.7584 28.3815 17.8622 15.3651 16.4867 13.1957 16.0785 13.4228 14.0562];

    %==========================================================================
    % FIGURE 1: Reaction Order Determination at 220C (1x3)
    %==========================================================================
    figure("Name","Rxn 1 - Reaction Order Determination at 220C")
    
    subplot(1,3,1)
    ln_x = log(H2O_c1_220);
    ln_y = log(R_HCHO_c1_220);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    H2O_Trendline = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{H_2O}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("H_2O Order | 0.01atm CH3OH, 0.01atm O_2, 220°C")
    hold off
    
    subplot(1,3,2)
    ln_x = log(O2_c2_220);
    ln_y = log(R_HCHO_c2_220);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    O2_Trendline = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{O_2}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("O_2 Order | 0.01atm CH3OH, 0.01atm H_2O, 250°C")
    hold off
    
    subplot(1,3,3)
    ln_x = log(CH3OH_c3_220);
    ln_y = log(R_HCHO_c3_220);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    CH3OH_Trendline = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{CH3OH}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("CH3OH Order | 0.01atm O_2, 0.01atm H_2O, 250°C")
    hold off
    
    sgtitle("Rxn 1 (CH_3OH \rightarrow HCHO): Reaction Order Determination at 220°C")
    
    rxn1_H2O_order_220  = H2O_Trendline(1);
    rxn1_O2_order_220   = O2_Trendline(1);
    rxn1_CH3OH_order_220 = CH3OH_Trendline(1);
    
    fprintf('Rxn1 220C Orders — CH3OH: %.4f, O2: %.4f, H2O: %.4f\n', ...
        rxn1_CH3OH_order_220, rxn1_O2_order_220, rxn1_H2O_order_220)
    
    %==========================================================================
    % FIGURE 2: Reaction Order Determination at 250C (1x3)
    %==========================================================================
    figure("Name","Rxn 1 - Reaction Order Determination at 250C")
    
    subplot(1,3,1)
    ln_x = log(H2O_a1_250);
    ln_y = log(R_HCHO_a1_250);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    H2O_Trendline = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{H_2O}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("H_2O Order | 6% CH3OH, 10% O_2, 250°C")
    hold off
    
    subplot(1,3,2)
    ln_x = log(O2_a2_250);
    ln_y = log(R_HCHO_a2_250);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    O2_Trendline = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{O_2}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("O_2 Order | 6% CH3OH, 2.5% H_2O, 250°C")
    hold off
    
    subplot(1,3,3)
    ln_x = log(CH3OH_a3_250);
    ln_y = log(R_HCHO_a3_250);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    CH3OH_Trendline = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{CH3OH}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("CH3OH Order | 10% O_2, 2.5% H_2O, 250°C")
    hold off
    
    sgtitle("Rxn 1 (CH_3OH \rightarrow HCHO): Reaction Order Determination at 250°C")
    
    rxn1_H2O_order_250  = H2O_Trendline(1);
    rxn1_O2_order_250   = O2_Trendline(1);
    rxn1_CH3OH_order_250 = CH3OH_Trendline(1);
    
    fprintf('Rxn1 250C Orders — CH3OH: %.4f, O2: %.4f, H2O: %.4f\n', ...
        rxn1_CH3OH_order_250, rxn1_O2_order_250, rxn1_H2O_order_250)
    
    %==========================================================================
    % FIGURE 3: Reaction Order Determination at 260C (1x3)
    %==========================================================================
    figure("Name","Rxn 1 - Reaction Order Determination at 260C")
    
    subplot(1,3,1)
    ln_x = log(CH3OH_d3_260);
    ln_y = log(R_HCHO_d3_260);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    CH3OH_Trendline_b = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{CH_3OH}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("CH3OH Order | 0.010atm O_2, 0.010atm H_2O, 260°C")
    hold off
    
    subplot(1,3,2)
    ln_x = log(O2_d2_260);
    ln_y = log(R_HCHO_d2_260);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    O2_Trendline_b = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{O_2}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("O_2 Order | 0.010atm CH_3OH, 0.010 H_2O, 260°C")
    hold off
    
    subplot(1,3,3)
    ln_x = log(H2O_d1_260);
    ln_y = log(R_HCHO_d1_260);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    H2O_Trendline_b = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{H_2O}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("H_2O Order | 0.010atm CH_3OH, 0.010atm O_2, 260°C")
    hold off
    
    sgtitle("Rxn 1 (CH_3OH \rightarrow HCHO): Reaction Order Determination at 260°C")
    
    rxn1_CH3OH_order_260 = CH3OH_Trendline_b(1);
    rxn1_O2_order_260    = O2_Trendline_b(1);
    rxn1_H2O_order_260   = H2O_Trendline_b(1);
    
    fprintf('Rxn1 260C Orders — CH3OH: %.4f, O2: %.4f, H2O: %.4f\n', ...
        rxn1_CH3OH_order_260, rxn1_O2_order_260, rxn1_H2O_order_260)
    

    %==========================================================================
    % FIGURE 4: Reaction Order Determination at 270C (1x3)
    %==========================================================================
    figure("Name","Rxn 1 - Reaction Order Determination at 270C")
    
    subplot(1,3,1)
    ln_x = log(CH3OH_b1_270);
    ln_y = log(R_HCHO_b1_270);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    CH3OH_Trendline_b = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{CH_3OH}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("CH3OH Order | 10% O_2, 2.5% H_2O, 270°C")
    hold off
    
    subplot(1,3,2)
    ln_x = log(O2_b2_270);
    ln_y = log(R_HCHO_b2_270);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    O2_Trendline_b = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{O_2}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("O_2 Order | 6% CH_3OH, 2.5% H_2O, 270°C")
    hold off
    
    subplot(1,3,3)
    ln_x = log(H2O_b3_270);
    ln_y = log(R_HCHO_b3_270);
    coeff = polyfit(ln_x, ln_y, 1);
    y_fit = polyval(coeff, ln_x);
    mdl = fitlm(ln_x, ln_y);
    r2 = mdl.Rsquared.Ordinary;
    H2O_Trendline_b = coeff;
    scatter(ln_x, ln_y, 'x', 'SizeData', 100, 'LineWidth', 2)
    hold on
    plot(ln_x, y_fit, 'r-', 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 = %.4f', r2), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    legend('Data Points','Linear Fit','Location','southeast')
    xlabel("ln(C_{H_2O}) (%)"); ylabel("ln(r_{HCHO}) (mol/g_{cat}/s)")
    title("H_2O Order | 6% CH_3OH, 10% O_2, 270°C")
    hold off
    
    sgtitle("Rxn 1 (CH_3OH \rightarrow HCHO): Reaction Order Determination at 270°C")
    
    rxn1_CH3OH_order_270 = CH3OH_Trendline_b(1);
    rxn1_O2_order_270    = O2_Trendline_b(1);
    rxn1_H2O_order_270   = H2O_Trendline_b(1);
    
    fprintf('Rxn1 270C Orders — CH3OH: %.4f, O2: %.4f, H2O: %.4f\n', ...
        rxn1_CH3OH_order_270, rxn1_O2_order_270, rxn1_H2O_order_270)
    
   
    %==========================================================================
    % Calculating Average Order 
    %==========================================================================
    rxn1_CH3OH_order_avg = (rxn1_CH3OH_order_260+rxn1_CH3OH_order_270+rxn1_CH3OH_order_220+rxn1_CH3OH_order_250)/4;
    rxn1_O2_order_avg = (rxn1_O2_order_260+rxn1_O2_order_270+rxn1_O2_order_220+rxn1_O2_order_250)/4;
    rxn1_H2O_order_avg = (rxn1_H2O_order_260+rxn1_H2O_order_270+rxn1_H2O_order_220+rxn1_H2O_order_250)/4;
    fprintf("Average Order of CH3OH for Reaction 1: %.3f \n", rxn1_CH3OH_order_avg)
    fprintf("Average Order of O2 for Reaction 1: %.3f \n", rxn1_O2_order_avg)
    fprintf("Average Order of H2O for Reaction 1: %.3f \n", rxn1_H2O_order_avg)

    %==========================================================================
    % Calculating k for each temperature
    %==========================================================================
    k = rate_220 ./ (CH3OH_220.^rxn1_CH3OH_order_220 .* O2_220.^rxn1_O2_order_220 .* H2O_220.^rxn1_H2O_order_220);
    k1_220 = mean(k);
    fprintf('Rxn1 k1 at 220C = %.8f\n', k1_220)
    k = rate_250 ./ (CH3OH_250.^rxn1_CH3OH_order_250 .* O2_250.^rxn1_O2_order_250 .* H2O_250.^rxn1_H2O_order_250);
    k1_250 = mean(k);
    fprintf('Rxn1 k1 at 250C = %.8f\n', k1_250)
    k = rate_260 ./ (CH3OH_260.^rxn1_CH3OH_order_260 .* O2_260.^rxn1_O2_order_260 .* H2O_260.^rxn1_H2O_order_260);
    k1_260 = mean(k);
    fprintf('Rxn1 k1 at 260C = %.8f\n', k1_260)
    k = rate_270 ./ (CH3OH_270.^rxn1_CH3OH_order_270 .* O2_270.^rxn1_O2_order_270 .* H2O_270.^rxn1_H2O_order_270);
    k1_270 = mean(k);
    fprintf('Rxn1 k1 at 270C = %.8f\n', k1_270)
    %==========================================================================
    % FIGURE 3: Kinetics Model Comparison at 250C (1x3)
    %==========================================================================
    figure("Name","Rxn 1 - Kinetics Model Comparison at 250C")
    
    % --- Subplot 1: O2 vs R(HCHO) at 250C ---
    subplot(1,3,1)
    o2_conc = linspace(0,15,100);
    R_pow  = k1_250 .* (6.^rxn1_CH3OH_order_avg) .* (o2_conc.^rxn1_O2_order_avg) .* (2.5.^rxn1_H2O_order_avg);
    nom    = alpha*kCH3OH*KCH3OH*(6/100)*KO2.*((o2_conc./100).^0.5);
    den    = (1+KCH3OH*(6/100)+KH2O*(2.5/100)).*(1+KO2.*(o2_conc./100).^0.5);
    R_lhhw = nom./den;
    r2_pow_O2_250  = compute_r2_pow(O2_a2_250, R_HCHO_a2_250, k1_250, rxn1_CH3OH_order_avg, rxn1_O2_order_avg, rxn1_H2O_order_avg, 6, 2.5, 'O2');
    r2_lhhw_O2_250 = compute_r2_lhhw1(O2_a2_250, R_HCHO_a2_250, alpha, kCH3OH, KCH3OH, KO2, KH2O, 6/100, 2.5/100, 'O2');
    plot(o2_conc, R_pow, 'LineWidth', 1.5); hold on
    scatter(O2_a2_250, R_HCHO_a2_250, 'x', 'SizeData', 100, 'LineWidth', 2)
    plot(o2_conc, R_lhhw, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow_O2_250, r2_lhhw_O2_250), 'Units','normalized','VerticalAlignment','top','FontSize',8)
    legend("Power Law","Experimental","LHHW",'Location','southeast')
    xlabel("C_{O_2} (%)"); ylabel("r_{HCHO} (mol/g_{cat}/s)")
    title("O_2 vs r_{HCHO} | 6% CH3OH, 2.5% H_2O, 250°C")
    hold off
    
    % --- Subplot 2: H2O vs R(HCHO) at 250C ---
    subplot(1,3,2)
    H2O_conc = linspace(0.5,5,100);
    R_pow  = k1_250 .* (6.^rxn1_CH3OH_order_avg) .* (10.^rxn1_O2_order_avg) .* (H2O_conc.^rxn1_H2O_order_avg);
    nom    = alpha*kCH3OH*KCH3OH*(6/100)*KO2.*((10/100).^0.5);
    den    = (1+KCH3OH*(6/100)+KH2O*(H2O_conc./100)).*(1+KO2.*(10/100).^0.5);
    R_lhhw = nom./den;
    r2_pow_H2O_250  = compute_r2_pow(H2O_a1_250, R_HCHO_a1_250, k1_250, rxn1_CH3OH_order_avg, rxn1_O2_order_avg, rxn1_H2O_order_avg, 6, 10, 'H2O');
    r2_lhhw_H2O_250 = compute_r2_lhhw1(H2O_a1_250, R_HCHO_a1_250, alpha, kCH3OH, KCH3OH, KO2, KH2O, 6/100, 10/100, 'H2O');
    plot(H2O_conc, R_pow, 'LineWidth', 1.5); hold on
    scatter(H2O_a1_250, R_HCHO_a1_250, 'x', 'SizeData', 100, 'LineWidth', 2)
    plot(H2O_conc, R_lhhw, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow_H2O_250, r2_lhhw_H2O_250), 'Units','normalized','VerticalAlignment','top','FontSize',8)
    legend("Power Law","Experimental","LHHW",'Location','northeast')
    xlabel("C_{H_2O} (%)"); ylabel("r_{HCHO} (mol/g_{cat}/s)")
    title("H_2O vs r_{HCHO} | 6% CH3OH, 10% O_2, 250°C")
    hold off
    
    % --- Subplot 3: CH3OH vs R(HCHO) at 250C ---
    subplot(1,3,3)
    CH3OH_conc = linspace(0,15,100);
    R_pow  = k1_250 .* (CH3OH_conc.^rxn1_CH3OH_order_avg) .* (10.^rxn1_O2_order_avg) .* (2.5.^rxn1_H2O_order_avg);
    nom    = alpha*kCH3OH*KCH3OH*(CH3OH_conc./100).*KO2.*((10/100).^0.5);
    den    = (1+KCH3OH*(CH3OH_conc./100)+KH2O*(2.5/100)).*(1+KO2.*(10/100).^0.5);
    R_lhhw = nom./den;
    r2_pow_CH3OH_250  = compute_r2_pow(CH3OH_a3_250, R_HCHO_a3_250, k1_250, rxn1_CH3OH_order_avg, rxn1_O2_order_avg, rxn1_H2O_order_avg, 10, 2.5, 'CH3OH');
    r2_lhhw_CH3OH_250 = compute_r2_lhhw1(CH3OH_a3_250, R_HCHO_a3_250, alpha, kCH3OH, KCH3OH, KO2, KH2O, 10/100, 2.5/100, 'CH3OH');
    plot(CH3OH_conc, R_pow, 'LineWidth', 1.5); hold on
    scatter(CH3OH_a3_250, R_HCHO_a3_250, 'x', 'SizeData', 100, 'LineWidth', 2)
    plot(CH3OH_conc, R_lhhw, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow_CH3OH_250, r2_lhhw_CH3OH_250), 'Units','normalized','VerticalAlignment','top','FontSize',8)
    legend("Power Law","Experimental","LHHW",'Location','southeast')
    xlabel("C_{CH3OH} (%)"); ylabel("r_{HCHO} (mol/g_{cat}/s)")
    title("CH3OH vs r_{HCHO} | 10% O_2, 2.5% H_2O, 250°C")
    hold off
    
    sgtitle("Rxn 1 (CH_3OH \rightarrow HCHO): Kinetics Model Comparison at 250°C")
    
    %==========================================================================
    % FIGURE 4: Kinetics Model Comparison at 270C (1x3)
    %==========================================================================
    figure("Name","Rxn 1 - Kinetics Model Comparison at 270C")
    
    % --- Subplot 1: CH3OH vs R(HCHO) at 270C ---
    subplot(1,3,1)
    CH3OH_conc = linspace(0,15,100);
    R_pow  = k1_270 .* (CH3OH_conc.^rxn1_CH3OH_order_avg) .* (10.^rxn1_O2_order_avg) .* (2.5.^rxn1_H2O_order_avg);
    nom    = alpha*kCH3OH_270*KCH3OH_270*(CH3OH_conc./100).*KO2_270.*((10/100).^0.5);
    den    = (1+KCH3OH_270*(CH3OH_conc./100)+KH2O_270*(2.5/100)).*(1+KO2_270.*(10/100).^0.5);
    R_lhhw = nom./den;
    r2_pow_CH3OH_270  = compute_r2_pow(CH3OH_b1_270, R_HCHO_b1_270, k1_270, rxn1_CH3OH_order_avg, rxn1_O2_order_avg, rxn1_H2O_order_avg, 10, 2.5, 'CH3OH');
    r2_lhhw_CH3OH_270 = compute_r2_lhhw1(CH3OH_b1_270, R_HCHO_b1_270, alpha, kCH3OH_270, KCH3OH_270, KO2_270, KH2O_270, 10/100, 2.5/100, 'CH3OH');
    plot(CH3OH_conc, R_pow, 'LineWidth', 1.5); hold on
    scatter(CH3OH_b1_270, R_HCHO_b1_270, 'x', 'SizeData', 100, 'LineWidth', 2)
    plot(CH3OH_conc, R_lhhw, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow_CH3OH_270, r2_lhhw_CH3OH_270), 'Units','normalized','VerticalAlignment','top','FontSize',8)
    legend("Power Law","Experimental","LHHW",'Location','southeast')
    xlabel("C_{CH_3OH} (%)"); ylabel("r_{HCHO} (mol/g_{cat}/s)")
    title("CH3OH vs r_{HCHO} | 10% O_2, 2.5% H_2O, 270°C")
    hold off
    
    % --- Subplot 2: O2 vs R(HCHO) at 270C ---
    subplot(1,3,2)
    o2_conc = linspace(0,15,100);
    R_pow  = k1_270 .* (6.^rxn1_CH3OH_order_avg) .* (o2_conc.^rxn1_O2_order_avg) .* (2.5.^rxn1_H2O_order_avg);
    nom    = alpha*kCH3OH_270*KCH3OH_270*(6/100)*KO2_270.*((o2_conc./100).^0.5);
    den    = (1+KCH3OH_270*(6/100)+KH2O_270*(2.5/100)).*(1+KO2_270.*(o2_conc./100).^0.5);
    R_lhhw = nom./den;
    r2_pow_O2_270  = compute_r2_pow(O2_b2_270, R_HCHO_b2_270, k1_270, rxn1_CH3OH_order_avg, rxn1_O2_order_avg, rxn1_H2O_order_avg, 6, 2.5, 'O2');
    r2_lhhw_O2_270 = compute_r2_lhhw1(O2_b2_270, R_HCHO_b2_270, alpha, kCH3OH_270, KCH3OH_270, KO2_270, KH2O_270, 6/100, 2.5/100, 'O2');
    plot(o2_conc, R_pow, 'LineWidth', 1.5); hold on
    scatter(O2_b2_270, R_HCHO_b2_270, 'x', 'SizeData', 100, 'LineWidth', 2)
    plot(o2_conc, R_lhhw, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow_O2_270, r2_lhhw_O2_270), 'Units','normalized','VerticalAlignment','top','FontSize',8)
    legend("Power Law","Experimental","LHHW",'Location','southeast')
    xlabel("C_{O_2} (%)"); ylabel("r_{HCHO} (mol/g_{cat}/s)")
    title("O_2 vs r_{HCHO} | 6% CH_3OH, 2.5% H_2O, 270°C")
    hold off
    
    % --- Subplot 3: H2O vs R(HCHO) at 270C ---
    subplot(1,3,3)
    H2O_conc = linspace(0.5,5,100);
    R_pow  = k1_270 .* (6.^rxn1_CH3OH_order_avg) .* (10.^rxn1_O2_order_avg) .* (H2O_conc.^rxn1_H2O_order_avg);
    nom    = alpha*kCH3OH_270*KCH3OH_270*(6/100)*KO2_270.*(10/100).^0.5;
    den    = (1+KCH3OH_270*(6/100)+KH2O_270*(H2O_conc./100)).*(1+KO2_270*(10/100).^0.5);
    R_lhhw = nom./den;
    r2_pow_H2O_270  = compute_r2_pow(H2O_b3_270, R_HCHO_b3_270, k1_270, rxn1_CH3OH_order_avg, rxn1_O2_order_avg, rxn1_H2O_order_avg, 6, 10, 'H2O');
    r2_lhhw_H2O_270 = compute_r2_lhhw1(H2O_b3_270, R_HCHO_b3_270, alpha, kCH3OH_270, KCH3OH_270, KO2_270, KH2O_270, 6/100, 10/100, 'H2O');
    plot(H2O_conc, R_pow, 'LineWidth', 1.5); hold on
    scatter(H2O_b3_270, R_HCHO_b3_270, 'x', 'SizeData', 100, 'LineWidth', 2)
    plot(H2O_conc, R_lhhw, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow_H2O_270, r2_lhhw_H2O_270), 'Units','normalized','VerticalAlignment','top','FontSize',8)
    legend("Power Law","Experimental","LHHW",'Location','northeast')
    xlabel("C_{H_2O} (%)"); ylabel("r_{HCHO} (mol/g_{cat}/s)")
    title("H_2O vs r_{HCHO} | 6% CH_3OH, 10% O_2, 270°C")
    hold off
    
    sgtitle("Rxn 1 (CH_3OH \rightarrow HCHO): Kinetics Model Comparison at 270°C")
    
    %==========================================================================
    % Arrhenius plot
    %==========================================================================
    T_arr = [T1 T2 220+273.15 260+273.15];
    inverseT = 1./T_arr;
    lnk_powerlaw = log([k1_250 k1_270 k1_220 k1_260]);
    
    figure("Name","Rxn 1 - Arrhenius Plot - Power Law")
    scatter(inverseT, lnk_powerlaw, 'x', 'SizeData', 100, 'LineWidth', 1.5)
    hold on
    Arrhenius_coeff = polyfit(inverseT, lnk_powerlaw, 1);
    plot(inverseT, polyval(Arrhenius_coeff, inverseT), 'r-', 'LineWidth', 1.5)
    hold off
    xlabel("1/T (K^{-1})"); ylabel("ln(k_1)")
    title("Rxn 1 (CH_3OH \rightarrow HCHO): Arrhenius Plot - Simple Power Law")
    legend("ln(k) data","Linear Fit",'Location','best')
    Ea_rxn1 = -Arrhenius_coeff(1) * 8.314 / 1000;
    k0_rxn1 = exp(Arrhenius_coeff(2));
    fprintf('Rxn1 Ea (Power Law) = %.4f kJ/mol\n', Ea_rxn1)
    fprintf('Rxn1 k0 (Power Law) = %.6e\n', k0_rxn1)
    
    
   %==========================================================================
    % Parity Plot - Reaction 1 (All temperatures, average orders)
    %==========================================================================
    figure("Name","Rxn 1 - Parity Plot")
    
    R = 8.314;
    T_220 = 220 + 273.15;
    T_250 = 250 + 273.15;
    T_260 = 260 + 273.15;
    T_270 = 270 + 273.15;
    
    % --- LHHW rate constants at each temperature ---
    A_kCH3OH = 1.50e7;   Ea_kCH3OH = 86.00e3;
    A_KCH3OH = 2.60e-4;  Ea_KCH3OH = -56.78e3;
    A_KO2    = 1.42e-5;  Ea_KO2    = -60.32e3;
    A_KH2O   = 5.50e-7;  Ea_KH2O   = -86.45e3;
    
    kCH3OH_220 = A_kCH3OH*exp(-Ea_kCH3OH/(R*T_220)); KCH3OH_220 = A_KCH3OH*exp(-Ea_KCH3OH/(R*T_220));
    KO2_220    = A_KO2   *exp(-Ea_KO2   /(R*T_220)); KH2O_220   = A_KH2O  *exp(-Ea_KH2O  /(R*T_220));
    
    kCH3OH_250 = A_kCH3OH*exp(-Ea_kCH3OH/(R*T_250)); KCH3OH_250 = A_KCH3OH*exp(-Ea_KCH3OH/(R*T_250));
    KO2_250    = A_KO2   *exp(-Ea_KO2   /(R*T_250)); KH2O_250   = A_KH2O  *exp(-Ea_KH2O  /(R*T_250));
    
    kCH3OH_260 = A_kCH3OH*exp(-Ea_kCH3OH/(R*T_260)); KCH3OH_260 = A_KCH3OH*exp(-Ea_KCH3OH/(R*T_260));
    KO2_260    = A_KO2   *exp(-Ea_KO2   /(R*T_260)); KH2O_260   = A_KH2O  *exp(-Ea_KH2O  /(R*T_260));
    
    kCH3OH_270 = A_kCH3OH*exp(-Ea_kCH3OH/(R*T_270)); KCH3OH_270 = A_KCH3OH*exp(-Ea_KCH3OH/(R*T_270));
    KO2_270    = A_KO2   *exp(-Ea_KO2   /(R*T_270)); KH2O_270   = A_KH2O  *exp(-Ea_KH2O  /(R*T_270));
    
    % --- Power law predictions using average orders ---
    % 220C
    rate_220_power = k1_220 .* (CH3OH_220.^rxn1_CH3OH_order_avg) .* (O2_220.^rxn1_O2_order_avg) .* (H2O_220.^rxn1_H2O_order_avg);
    % 250C
    rate_250_power = k1_250 .* (CH3OH_250.^rxn1_CH3OH_order_avg) .* (O2_250.^rxn1_O2_order_avg) .* (H2O_250.^rxn1_H2O_order_avg);
    % 260C
    rate_260_power = k1_260 .* (CH3OH_260.^rxn1_CH3OH_order_avg) .* (O2_260.^rxn1_O2_order_avg) .* (H2O_260.^rxn1_H2O_order_avg);
    % 270C
    rate_270_power = k1_270 .* (CH3OH_270.^rxn1_CH3OH_order_avg) .* (O2_270.^rxn1_O2_order_avg) .* (H2O_270.^rxn1_H2O_order_avg);
    
    % --- LHHW predictions ---
    % 220C — concentrations already in atm
    nom_220_LHHW  = alpha.*kCH3OH_220.*KCH3OH_220.*(CH3OH_220).*KO2_220.*((O2_220).^0.5);
    den_220_LHHW  = (1+KCH3OH_220.*(CH3OH_220)+KH2O_220.*(H2O_220)).*(1+KO2_220.*(O2_220).^0.5);
    rate_220_LHHW = nom_220_LHHW./den_220_LHHW;
    
    % 250C — concentrations in percent, divide by 100
    nom_250_LHHW  = alpha.*kCH3OH_250.*KCH3OH_250.*(CH3OH_250./100).*KO2_250.*((O2_250./100).^0.5);
    den_250_LHHW  = (1+KCH3OH_250.*(CH3OH_250./100)+KH2O_250.*(H2O_250./100)).*(1+KO2_250.*(O2_250./100).^0.5);
    rate_250_LHHW = nom_250_LHHW./den_250_LHHW;
    
    % 260C — concentrations already in atm
    nom_260_LHHW  = alpha.*kCH3OH_260.*KCH3OH_260.*(CH3OH_260).*KO2_260.*((O2_260).^0.5);
    den_260_LHHW  = (1+KCH3OH_260.*(CH3OH_260)+KH2O_260.*(H2O_260)).*(1+KO2_260.*(O2_260).^0.5);
    rate_260_LHHW = nom_260_LHHW./den_260_LHHW;
    
    % 270C — concentrations in percent, divide by 100
    nom_270_LHHW  = alpha.*kCH3OH_270.*KCH3OH_270.*(CH3OH_270./100).*KO2_270.*((O2_270./100).^0.5);
    den_270_LHHW  = (1+KCH3OH_270.*(CH3OH_270./100)+KH2O_270.*(H2O_270./100)).*(1+KO2_270.*(O2_270./100).^0.5);
    rate_270_LHHW = nom_270_LHHW./den_270_LHHW;
    % --- Combine all experimental and predicted for axis limits ---
    all_exp  = [rate_220, rate_250, rate_260, rate_270];
    all_pred = [rate_220_power, rate_250_power, rate_260_power, rate_270_power, ...
                rate_220_LHHW,  rate_250_LHHW,  rate_260_LHHW,  rate_270_LHHW];
    ax_max   = max([all_exp, all_pred]) * 1.1;
    x_line   = linspace(0, ax_max, 200);
    
    % --- Scatter: Power Law (all temps, same marker/colour) ---
    scatter([rate_220, rate_250, rate_260, rate_270], ...
            [rate_220_power, rate_250_power, rate_260_power, rate_270_power], ...
            60, 'x', 'r', 'LineWidth', 1.5, 'DisplayName', 'Simple Power Law')
    hold on
    
    % --- Scatter: LHHW (all temps, same marker/colour) ---
    scatter([rate_220, rate_250, rate_260, rate_270], ...
            [rate_220_LHHW,  rate_250_LHHW,  rate_260_LHHW,  rate_270_LHHW], ...
            60, 'x', 'b', 'LineWidth', 1.5, 'DisplayName', 'LHHW')
    
    % --- y = x line ---
    plot(x_line, x_line, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
    hold off
    
    xlabel('r_{HCHO} experimental (mol/kg_{cat}/s)')
    ylabel('r_{HCHO} predicted (mol/kg_{cat}/s)')
    title('Rxn 1 (CH_3OH \rightarrow HCHO): Parity Plot (220–270°C, Average Orders)')
    legend('Location', 'best')
    axis equal
    xlim([0 0.03]); ylim([0 0.03])
    grid on
    

end

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

function r2 = compute_r2_pow(x_data, y_data, k1, ch3oh_ord, o2_ord, h2o_ord, const1, const2, varying)
    switch varying
        case 'CH3OH'
            y_pred = k1 .* (x_data.^ch3oh_ord) .* (const1.^o2_ord)  .* (const2.^h2o_ord);
        case 'O2'
            y_pred = k1 .* (const1.^ch3oh_ord) .* (x_data.^o2_ord)  .* (const2.^h2o_ord);
        case 'H2O'
            y_pred = k1 .* (const1.^ch3oh_ord) .* (const2.^o2_ord)  .* (x_data.^h2o_ord);
    end
    ss_res = sum((y_data - y_pred).^2);
    ss_tot = sum((y_data - mean(y_data)).^2);
    r2 = 1 - ss_res/ss_tot;
end

function r2 = compute_r2_lhhw1(x_data, y_data, alpha, kCH3OH, KCH3OH, KO2, KH2O, const1, const2, varying)
    switch varying
        case 'CH3OH'
            nom = alpha.*kCH3OH.*KCH3OH.*(x_data./100).*KO2.*(const1.^0.5);
            den = (1+KCH3OH.*(x_data./100)+KH2O.*const2).*(1+KO2.*const1.^0.5);
        case 'O2'
            nom = alpha.*kCH3OH.*KCH3OH.*const1.*KO2.*((x_data./100).^0.5);
            den = (1+KCH3OH.*const1+KH2O.*const2).*(1+KO2.*(x_data./100).^0.5);
        case 'H2O'
            nom = alpha.*kCH3OH.*KCH3OH.*const1.*KO2.*(const2.^0.5);
            den = (1+KCH3OH.*const1+KH2O.*(x_data./100)).*(1+KO2.*const2.^0.5);
    end
    y_pred = nom./den;
    ss_res = sum((y_data - y_pred).^2);
    ss_tot = sum((y_data - mean(y_data)).^2);
    r2 = 1 - ss_res/ss_tot;
end
%% STEP 2 - Reaction 2 Graphs Function
function plot_r2_simplepower()
    fprintf("\n")
    fprintf("=============================\n")
    fprintf("Reaction 2 - Simple Power Law Determination \n")
    fprintf("=============================\n")
    R = 8.314; % J/mol/K
    T = 250 + 273.15;
    A_kCH3OH = 1.50e7;   Ea_kCH3OH = 86.00e3;
    A_KCH3OH = 2.60e-4;  Ea_KCH3OH = -56.78e3;
    A_KO2    = 1.42e-5;  Ea_KO2    = -60.32e3;
    A_KH2O   = 5.50e-7;  Ea_KH2O   = -86.45e3;
    A_kCO    = 3.50e2;   Ea_kCO    = 46.00e3;
    
    kCH3OH = A_kCH3OH * exp(-Ea_kCH3OH/(R*T));
    KCH3OH = A_KCH3OH * exp(-Ea_KCH3OH/(R*T));
    KO2    = A_KO2    * exp(-Ea_KO2/(R*T));
    KH2O   = A_KH2O   * exp(-Ea_KH2O/(R*T));
    kCO    = A_kCO    * exp(-Ea_kCO/(R*T));
    alpha  = 0.44;

    % --- Data ---
    HCHO_a_250 = [0.734917733 1.524680073 2.413162706 2.654478976 3.411334552 4.091407678 5.188299817 6.175502742];
    R_CO_a_250 = 10^-3.*[0.021367521 0.051994302 0.061253561 0.06980057 0.095441595 0.113960114 0.134615385 0.153846154];
    ln_HCHO_a_250 = log(HCHO_a_250);
    ln_R_CO_a_250 = log(R_CO_a_250);
    
    O2_b_250 = [1.085502878 1.98787339 2.888599616 3.886544259 4.885996163 10.02959715 14.71245547];
    R_CO_b_250 = 10^-3.*[0.069213483 0.069606742 0.079438202 0.091235955 0.094382022 0.11011236 0.110505618];
    ln_O2_b_250 = log(O2_b_250);
    ln_R_CO_b_250 = log(R_CO_b_250);
    
    HCHO_a_290 = [0.764044944 1.382022472 2.303370787 2.820224719 3.04494382 3.06741573 4.078651685 5.191011236 6.168539326];
    R_CO_a_290 = 10^-3.*[0.06282318 0.089151381 0.131491618 0.154381809 0.161262098 0.184058389 0.245788918 0.289305035 0.31912033];
    ln_HCHO_a_290 = log(HCHO_a_290);
    ln_R_CO_a_290 = log(R_CO_a_290);
    
    O2_b_290 = [1.12012987 1.899350649 3.043831169 3.993506494 5.089285714 10.08116883 14.75649351];
    R_CO_b_290 = 10^-3.*[0.120923913 0.151494565 0.161005435 0.171195652 0.180706522 0.203125 0.222826087];
    ln_O2_b_290 = log(O2_b_290);
    ln_R_CO_b_290 = log(R_CO_b_290);
    
    %==========================================================================
    % FIGURE 1: Reaction Order Subplots (2x2)
    % Reaction 2: HCHO + 1/2 O2 -> CO + H2O
    %==========================================================================
    figure("Name","Rxn 2 - Reaction Order Determination")
    
    % --- Subplot 1: HCHO order at 250C ---
    subplot(2,2,1)
    coeff_a = polyfit(ln_HCHO_a_250, ln_R_CO_a_250, 1);
    y_fit_a_250 = polyval(coeff_a, ln_HCHO_a_250);
    mdl = fitlm(ln_HCHO_a_250, ln_R_CO_a_250);
    eqn2_a_250_rsquared = mdl.Rsquared.Ordinary;
    scatter(ln_HCHO_a_250, ln_R_CO_a_250, 'x', "SizeData", 100, "LineWidth", 2);
    hold on
    plot(ln_HCHO_a_250, y_fit_a_250, 'r-', 'LineWidth', 1.5);
    text(0.05, 0.95, sprintf('R^2 = %.4f', eqn2_a_250_rsquared), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9);
    legend('Data Points', 'Linear Fit', 'Location', 'southeast');
    xlabel("ln(C_{HCHO}) (%)")
    ylabel("ln(r_{CO}) (mol/g_{cat}/s)")
    title("HCHO Order | 5% O_2, 250°C")
    hold off

    % --- Subplot 2: O2 order at 250C ---
    subplot(2,2,2)
    coeff_b = polyfit(ln_O2_b_250, ln_R_CO_b_250, 1);
    y_fit_b_250 = polyval(coeff_b, ln_O2_b_250);
    mdl = fitlm(ln_O2_b_250, ln_R_CO_b_250);
    eqn2_b_250_rsquared = mdl.Rsquared.Ordinary;
    scatter(ln_O2_b_250, ln_R_CO_b_250, 'x', "SizeData", 100, "LineWidth", 2);
    hold on
    plot(ln_O2_b_250, y_fit_b_250, 'r-', 'LineWidth', 1.5);
    text(0.05, 0.95, sprintf('R^2 = %.4f', eqn2_b_250_rsquared), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9);
    legend('Data Points', 'Linear Fit', 'Location', 'southeast');
    xlabel("ln(C_{O_2}) (%)")
    ylabel("ln(r_{CO}) (mol/g_{cat}/s)")
    title("O_2 Order | 3.5% HCHO, 250°C")
    hold off

    % --- Subplot 3: HCHO order at 290C ---
    subplot(2,2,3)
    coeff_a_290 = polyfit(ln_HCHO_a_290, ln_R_CO_a_290, 1);
    y_fit_a_290 = polyval(coeff_a_290, ln_HCHO_a_290);
    mdl = fitlm(ln_HCHO_a_290, ln_R_CO_a_290);
    eqn2_a_290_rsquared = mdl.Rsquared.Ordinary;
    scatter(ln_HCHO_a_290, ln_R_CO_a_290, 'x', "SizeData", 100, "LineWidth", 2);
    hold on
    plot(ln_HCHO_a_290, y_fit_a_290, 'r-', 'LineWidth', 1.5);
    text(0.05, 0.95, sprintf('R^2 = %.4f', eqn2_a_290_rsquared), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9);
    legend('Data Points', 'Linear Fit', 'Location', 'southeast');
    xlabel("ln(C_{HCHO}) (%)")
    ylabel("ln(r_{CO}) (mol/g_{cat}/s)")
    title("HCHO Order | 5% O_2, 290°C")
    hold off

    % --- Subplot 4: O2 order at 290C ---
    subplot(2,2,4)
    coeff_b_290 = polyfit(ln_O2_b_290, ln_R_CO_b_290, 1);
    y_fit_b_290 = polyval(coeff_b_290, ln_O2_b_290);
    mdl = fitlm(ln_O2_b_290, ln_R_CO_b_290);
    eqn2_b_290_rsquared = mdl.Rsquared.Ordinary;
    scatter(ln_O2_b_290, ln_R_CO_b_290, 'x', "SizeData", 100, "LineWidth", 2);
    hold on
    plot(ln_O2_b_290, y_fit_b_290, 'r-', 'LineWidth', 1.5);
    text(0.05, 0.95, sprintf('R^2 = %.4f', eqn2_b_290_rsquared), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9);
    legend('Data Points', 'Linear Fit', 'Location', 'southeast');
    xlabel("ln(C_{O_2}) (%)")
    ylabel("ln(r_{CO}) (mol/g_{cat}/s)")
    title("O_2 Order | 3.5% HCHO, 290°C")
    hold off
    
    sgtitle("Rxn 2 (HCHO \rightarrow CO): Reaction Order Determination")
    
    % --- Reaction orders ---
    eqn2_HCHO_order_250 = coeff_a(1);
    eqn2_O2_order_250   = coeff_b(1);
    eqn2_HCHO_order_290 = coeff_a_290(1);
    eqn2_O2_order_290   = coeff_b_290(1);
    eqn2_HCHO_order = (eqn2_HCHO_order_290 + eqn2_HCHO_order_250) / 2;
    eqn2_O2_order   = (eqn2_O2_order_250   + eqn2_O2_order_290)   / 2;
    fprintf("Average Order of HCHO for Reaction 2: %.3f\n",eqn2_HCHO_order)
    fprintf("Average Order of O2 for Reaction 2: %.3f\n",eqn2_O2_order)
    
    % --- Back-calculate k2 at 250C from data ---
    k_vals_250 = [R_CO_a_250 ./ (HCHO_a_250.^eqn2_HCHO_order .* 5.^eqn2_O2_order), ...
                  R_CO_b_250 ./ (3.5.^eqn2_HCHO_order .* O2_b_250.^eqn2_O2_order)];
    k2_250c = mean(k_vals_250);
    fprintf('Rxn2 k2 at 250C (back-calculated) = %.6e\n', k2_250c)
    
    % --- Back-calculate k2 at 290C from data ---
    k_vals_290 = [R_CO_a_290 ./ (HCHO_a_290.^eqn2_HCHO_order .* 5.^eqn2_O2_order), ...
                  R_CO_b_290 ./ (3.5.^eqn2_HCHO_order .* O2_b_290.^eqn2_O2_order)];
    k2_290c = mean(k_vals_290);
    fprintf('Rxn2 k2 at 290C (back-calculated) = %.6e\n', k2_290c)
    %==========================================================================
    % FIGURE 2: Kinetics Model Comparison (2x2)
    % Reaction 2: HCHO + 1/2 O2 -> CO + H2O
    %==========================================================================
    figure("Name","Rxn 2 - Kinetics Model Comparison")
    
    % --- Subplot 1: HCHO vs rate of CO at 250C ---
    T = 250 + 273.15;
    kCO_250  = A_kCO * exp(-Ea_kCO/(R*T));
    KO2_250  = A_KO2 * exp(-Ea_KO2/(R*T));
    
    subplot(2,2,1)
    HCHO_values = linspace(0,7,100);
    R_CO_a_predicted = k2_250c .* HCHO_values.^eqn2_HCHO_order .* 5.^eqn2_O2_order;
    nom2 = alpha * kCO_250 .* HCHO_values/100 * KO2_250 * sqrt(0.05);
    den2 = 1 + KO2_250 * sqrt(0.05);
    R_lhhw_a = nom2 ./ den2;
    r2_pow     = compute_r2(HCHO_a_250, R_CO_a_250, k2_250c, eqn2_HCHO_order, eqn2_O2_order, 5, 'HCHO');
    r2_lhhw_val = compute_r2_lhhw(HCHO_a_250, R_CO_a_250, alpha, kCO_250, KO2_250, 0.05, 'HCHO');
    plot(HCHO_values, R_CO_a_predicted, 'LineWidth', 1.5)
    hold on
    scatter(HCHO_a_250, R_CO_a_250, 'x', "SizeData", 100, "LineWidth", 2)
    plot(HCHO_values, R_lhhw_a, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow, r2_lhhw_val), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    legend("Power Law","Experimental","LHHW",'Location','southeast')
    xlabel("C_{HCHO} (%)")
    ylabel("r_{CO} (mol/g_{cat}/s)")
    title("HCHO vs r_{CO} | 5% O_2, 250°C")
    hold off
    
    % --- Subplot 2: O2 vs rate of CO at 250C ---
    subplot(2,2,2)
    o2_values = linspace(0,15,100);
    R_CO_b_predicted = k2_250c .* 3.5.^eqn2_HCHO_order .* o2_values.^eqn2_O2_order;
    nom2 = alpha * kCO_250 * 3.5/100 * KO2_250 .* (o2_values/100).^0.5;
    den2 = 1 + KO2_250 .* (o2_values/100).^0.5;
    R_lhhw_b = nom2 ./ den2;
    r2_pow     = compute_r2(O2_b_250, R_CO_b_250, k2_250c, eqn2_HCHO_order, eqn2_O2_order, 3.5, 'O2');
    r2_lhhw_val = compute_r2_lhhw(O2_b_250, R_CO_b_250, alpha, kCO_250, KO2_250, 3.5/100, 'O2');
    plot(o2_values, R_CO_b_predicted, 'LineWidth', 1.5)
    hold on
    scatter(O2_b_250, R_CO_b_250, 'x', "SizeData", 100, "LineWidth", 2)
    plot(o2_values, R_lhhw_b, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow, r2_lhhw_val), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    legend("Power Law","Experimental","LHHW",'Location','southeast')
    xlabel("C_{O_2} (%)")
    ylabel("r_{CO} (mol/g_{cat}/s)")
    title("O_2 vs r_{CO} | 3.5% HCHO, 250°C")
    hold off
    
    % --- Subplot 3: HCHO vs rate of CO at 290C ---
    T = 290 + 273.15;
    kCO_290  = A_kCO * exp(-Ea_kCO/(R*T));
    KO2_290  = A_KO2 * exp(-Ea_KO2/(R*T));
    
    subplot(2,2,3)
    HCHO_values = linspace(0,7,100);
    R_CO_a_predicted = k2_290c .* HCHO_values.^eqn2_HCHO_order .* 5.^eqn2_O2_order;
    nom2 = alpha * kCO_290 .* HCHO_values/100 * KO2_290 * sqrt(0.05);
    den2 = 1 + KO2_290 * sqrt(0.05);
    R_lhhw_a = nom2 ./ den2;
    r2_pow     = compute_r2(HCHO_a_290, R_CO_a_290, k2_290c, eqn2_HCHO_order, eqn2_O2_order, 5, 'HCHO');
    r2_lhhw_val = compute_r2_lhhw(HCHO_a_290, R_CO_a_290, alpha, kCO_290, KO2_290, 0.05, 'HCHO');
    plot(HCHO_values, R_CO_a_predicted, 'LineWidth', 1.5)
    hold on
    scatter(HCHO_a_290, R_CO_a_290, 'x', "SizeData", 100, "LineWidth", 2)
    plot(HCHO_values, R_lhhw_a, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow, r2_lhhw_val), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    legend("Power Law","Experimental","LHHW",'Location','southeast')
    xlabel("C_{HCHO} (%)")
    ylabel("r_{CO} (mol/g_{cat}/s)")
    title("HCHO vs r_{CO} | 5% O_2, 290°C")
    hold off
    
    % --- Subplot 4: O2 vs rate of CO at 290C ---
    subplot(2,2,4)
    o2_values = linspace(0,15,100);
    R_CO_b_predicted = k2_290c .* 3.5.^eqn2_HCHO_order .* o2_values.^eqn2_O2_order;
    nom2 = alpha * kCO_290 * 3.5/100 * KO2_290 .* (o2_values/100).^0.5;
    den2 = 1 + KO2_290 .* (o2_values/100).^0.5;
    R_lhhw_b = nom2 ./ den2;
    r2_pow     = compute_r2(O2_b_290, R_CO_b_290, k2_290c, eqn2_HCHO_order, eqn2_O2_order, 3.5, 'O2');
    r2_lhhw_val = compute_r2_lhhw(O2_b_290, R_CO_b_290, alpha, kCO_290, KO2_290, 3.5/100, 'O2');
    plot(o2_values, R_CO_b_predicted, 'LineWidth', 1.5)
    hold on
    scatter(O2_b_290, R_CO_b_290, 'x', "SizeData", 100, "LineWidth", 2)
    plot(o2_values, R_lhhw_b, 'LineWidth', 1.5)
    text(0.05, 0.95, sprintf('R^2 PL = %.4f\nR^2 LHHW = %.4f', r2_pow, r2_lhhw_val), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    legend("Power Law","Experimental","LHHW",'Location','southeast')
    xlabel("C_{O_2} (%)")
    ylabel("r_{CO} (mol/g_{cat}/s)")
    title("O_2 vs r_{CO} | 3.5% HCHO, 290°C")
    hold off
    
    sgtitle("Rxn 2 (HCHO \rightarrow CO): Kinetics Model Comparison")

    %==========================================================================
    % Parity Plot - Reaction 2 (All temperatures, 2 legend entries)
    %==========================================================================
    figure("Name","Rxn 2 - Parity Plot")
    
    % --- 250C predictions ---
    % Power Law
    R_CO_pow_HCHO_250 = k2_250c .* HCHO_a_250.^eqn2_HCHO_order .* 5.^eqn2_O2_order;
    R_CO_pow_O2_250   = k2_250c .* 3.5.^eqn2_HCHO_order .* O2_b_250.^eqn2_O2_order;
    % LHHW
    R_CO_lhhw_HCHO_250 = alpha.*kCO_250.*(HCHO_a_250./100).*KO2_250.*sqrt(0.05) ...
                         ./ (1 + KO2_250.*sqrt(0.05));
    R_CO_lhhw_O2_250   = alpha.*kCO_250.*(3.5/100).*KO2_250.*sqrt(O2_b_250./100) ...
                         ./ (1 + KO2_250.*sqrt(O2_b_250./100));
    
    % --- 290C predictions ---
    % Power Law
    R_CO_pow_HCHO_290 = k2_290c .* HCHO_a_290.^eqn2_HCHO_order .* 5.^eqn2_O2_order;
    R_CO_pow_O2_290   = k2_290c .* 3.5.^eqn2_HCHO_order .* O2_b_290.^eqn2_O2_order;
    % LHHW
    R_CO_lhhw_HCHO_290 = alpha.*kCO_290.*(HCHO_a_290./100).*KO2_290.*sqrt(0.05) ...
                         ./ (1 + KO2_290.*sqrt(0.05));
    R_CO_lhhw_O2_290   = alpha.*kCO_290.*(3.5/100).*KO2_290.*sqrt(O2_b_290./100) ...
                         ./ (1 + KO2_290.*sqrt(O2_b_290./100));
    
    % --- Combine all experimental and predicted for axis limits ---
    all_exp  = [R_CO_a_250, R_CO_b_250, R_CO_a_290, R_CO_b_290];
    all_pred = [R_CO_pow_HCHO_250, R_CO_pow_O2_250, R_CO_lhhw_HCHO_250, R_CO_lhhw_O2_250, ...
                R_CO_pow_HCHO_290, R_CO_pow_O2_290, R_CO_lhhw_HCHO_290, R_CO_lhhw_O2_290];
    ax_max = max([all_exp, all_pred]) * 1.1;
    x_line = linspace(0, ax_max, 200);
    
    % --- Scatter: Power Law (all temps, same marker/colour) ---
    scatter([R_CO_a_250, R_CO_b_250, R_CO_a_290, R_CO_b_290], ...
            [R_CO_pow_HCHO_250, R_CO_pow_O2_250, R_CO_pow_HCHO_290, R_CO_pow_O2_290], ...
            60, 'x', 'r', 'LineWidth', 1.5, 'DisplayName', 'Simple Power Law')
    hold on
    
    % --- Scatter: LHHW (all temps, same marker/colour) ---
    scatter([R_CO_a_250, R_CO_b_250, R_CO_a_290, R_CO_b_290], ...
            [R_CO_lhhw_HCHO_250, R_CO_lhhw_O2_250, R_CO_lhhw_HCHO_290, R_CO_lhhw_O2_290], ...
            60, 'x', 'b', 'LineWidth', 1.5, 'DisplayName', 'LHHW')
    
    % --- y = x line ---
    plot(x_line, x_line, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
    hold off
    
    xlabel('r_{CO} experimental (mol/kg_{cat}/s)')
    ylabel('r_{CO} predicted (mol/kg_{cat}/s)')
    title('Rxn 2 (HCHO \rightarrow CO): Parity Plot (250–290°C)')
    legend('Location', 'best')
    axis equal
    xlim([0 ax_max]); ylim([0 ax_max])
    grid on

    %==========================================================================
    % Arrhenius Plot - Reaction 2
    %==========================================================================
    T_250K = 250 + 273.15;
    T_290K = 290 + 273.15;
    inverseT_r2  = 1./[T_250K, T_290K];
    lnk_r2       = log([k2_250c, k2_290c]);
    
    figure("Name","Rxn 2 - Arrhenius Plot - Power Law")
    scatter(inverseT_r2, lnk_r2, 'x', 'SizeData', 100, 'LineWidth', 1.5)
    hold on
    Arrhenius_coeff_r2 = polyfit(inverseT_r2, lnk_r2, 1);
    plot(inverseT_r2, polyval(Arrhenius_coeff_r2, inverseT_r2), 'r-', 'LineWidth', 1.5)
    hold off
    xlabel("1/T (K^{-1})"); ylabel("ln(k_2)")
    title("Rxn 2 (HCHO \rightarrow CO): Arrhenius Plot - Simple Power Law")
    legend("ln(k) data","Linear Fit",'Location','best')
    
    Ea_rxn2 = -Arrhenius_coeff_r2(1) * 8.314 / 1000;
    k0_rxn2 = exp(Arrhenius_coeff_r2(2));
    fprintf('Rxn2 Ea (Power Law) = %.4f kJ/mol\n', Ea_rxn2)
    fprintf('Rxn2 k0 (Power Law) = %.6e\n', k0_rxn2)
end

% --- Helper: R² for Power Law (Reaction 2) ---
function r2 = compute_r2(x_data, y_data, k2, hcho_order, o2_order, const_conc, varying)
    if strcmp(varying, 'HCHO')
        y_pred = k2 .* x_data.^hcho_order .* const_conc.^o2_order;
    else
        y_pred = k2 .* const_conc.^hcho_order .* x_data.^o2_order;
    end
    ss_res = sum((y_data - y_pred).^2);
    ss_tot = sum((y_data - mean(y_data)).^2);
    r2 = 1 - ss_res / ss_tot;
end

% --- Helper: R² for LHHW (Reaction 2) ---
function r2 = compute_r2_lhhw(x_data, y_data, alpha, kCO, KO2, const_conc, varying)
    if strcmp(varying, 'HCHO')
        y_pred = alpha * kCO .* (x_data/100) * KO2 * sqrt(const_conc) ...
                 ./ (1 + KO2 * sqrt(const_conc));
    else
        y_pred = alpha * kCO * const_conc * KO2 .* sqrt(x_data/100) ...
                 ./ (1 + KO2 .* sqrt(x_data/100));
    end
    ss_res = sum((y_data - y_pred).^2);
    ss_tot = sum((y_data - mean(y_data)).^2);
    r2 = 1 - ss_res / ss_tot;
end

%% Getting Initial Conditions
function [feed_stream, coldshot_stream, F0_CH3OH, F0_HCHO] = get_initial_conditions(m_flowrate, P0, T0, split_ratio, inlet)
    Mr_CH3OH = 32.042e-3;
    Mr_O2    = 32.00e-3;
    Mr_H2O   = 18.00e-3;
    Mr_N2    = 28.014e-3;

    y0_CH3OH = inlet(1); y0_O2   = inlet(2);
    y0_H2O   = inlet(3); y0_HCHO = inlet(4);
    y0_CO    = inlet(5); y0_N2   = inlet(6);

    Mr_mix = y0_CH3OH*Mr_CH3OH + y0_O2*Mr_O2 + y0_H2O*Mr_H2O + y0_N2*Mr_N2;
    nT0    = m_flowrate / Mr_mix;

    nT0_cs = nT0 * split_ratio;
    nT0_f  = nT0 - nT0_cs;

    n0_CH3OH_f = y0_CH3OH * nT0_f;  n0_O2_f   = y0_O2   * nT0_f;
    n0_H2O_f   = y0_H2O   * nT0_f;  n0_HCHO_f = y0_HCHO * nT0_f;
    n0_CO_f    = y0_CO    * nT0_f;  n0_N2_f   = y0_N2   * nT0_f;

    n0_CH3OH_cs = y0_CH3OH * nT0_cs; n0_O2_cs   = y0_O2   * nT0_cs;
    n0_H2O_cs   = y0_H2O   * nT0_cs; n0_HCHO_cs = y0_HCHO * nT0_cs;
    n0_CO_cs    = y0_CO    * nT0_cs; n0_N2_cs   = y0_N2   * nT0_cs;

    F0_CH3OH = n0_CH3OH_f + n0_CH3OH_cs;
    F0_HCHO  = n0_HCHO_f  + n0_HCHO_cs;

    feed_stream = [n0_CH3OH_f, n0_O2_f, n0_HCHO_f, n0_H2O_f, n0_CO_f, n0_N2_f, 0, 0, P0, T0];
    coldshot_stream = [y0_CH3OH, y0_O2, y0_H2O, y0_HCHO, y0_CO, y0_N2, ...
                       n0_CH3OH_cs, n0_O2_cs, n0_HCHO_cs, n0_H2O_cs, n0_CO_cs, n0_N2_cs, ...
                       0, 0, P0, T0, nT0_cs];
end

%% Getting Outlet Conditions
function outlet_conditions = get_outlet_conditions(z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T)
    outlet_conditions = [z_points(end), F_CH3OH(end), F_O2(end), F_HCHO(end), ...
                         F_H2O(end), F_CO(end), F_N2(end), extent1(end), extent2(end), P(end), T(end)];
end

%% Cp calculation
function y_Cp = Cp(T, y_CH3OH, y_O2, y_H2O, y_N2, y_HCHO, y_CO)
    y = [y_CH3OH, y_O2, y_H2O, y_N2, y_HCHO, y_CO];
    a = [40.0460; 29.5260; 33.9330; 29.3420; 29.6750; 29.5560];
    b = [-3.8287e-02; -8.8999e-03; -8.4186e-03; -3.5395e-03;  1.8937e-02; -6.5807e-03];
    c = [ 2.4529e-04;  3.8083e-05;  2.9906e-05;  1.0076e-05;  2.8739e-05;  2.0130e-05];
    d = [-2.1679e-07; -3.2629e-08; -1.7825e-08; -4.3116e-09; -2.0092e-08; -1.2227e-08];
    e = [ 5.9909e-11;  8.8607e-12;  3.6934e-12;  2.5935e-13;  3.4333e-12;  2.2617e-12];
    Cp_ls = zeros(6,1);
    for i = 1:6
        Cp_ls(i) = y(i)*(a(i) + b(i)*T + c(i)*T^2 + d(i)*T^3 + e(i)*T^4);
    end
    y_Cp = Cp_ls;
end

function y_Cp_mean = Cp_mean(T1, T2, y_CH3OH, y_O2, y_H2O, y_N2, y_HCHO, y_CO)
    y = [y_CH3OH, y_O2, y_H2O, y_N2, y_HCHO, y_CO];
    a = [40.0460; 29.5260; 33.9330; 29.3420; 29.6750; 29.5560];
    b = [-3.8287e-02; -8.8999e-03; -8.4186e-03; -3.5395e-03;  1.8937e-02; -6.5807e-03];
    c = [ 2.4529e-04;  3.8083e-05;  2.9906e-05;  1.0076e-05;  2.8739e-05;  2.0130e-05];
    d = [-2.1679e-07; -3.2629e-08; -1.7825e-08; -4.3116e-09; -2.0092e-08; -1.2227e-08];
    e = [ 5.9909e-11;  8.8607e-12;  3.6934e-12;  2.5935e-13;  3.4333e-12;  2.2617e-12];
    Cp_mean_ls = zeros(6,1);
    for i = 1:6
        integral_val = a(i)*(T2-T1) + b(i)/2*(T2^2-T1^2) + c(i)/3*(T2^3-T1^3) ...
                     + d(i)/4*(T2^4-T1^4) + e(i)/5*(T2^5-T1^5);
        Cp_mean_ls(i) = y(i) * integral_val / (T2 - T1);
    end
    y_Cp_mean = Cp_mean_ls;
end

%% Event functions
function [value, isterminal, direction] = T_max_event(z, F, T_max)
    value = F(10) - T_max; isterminal = 1; direction = 1;
end
function [value, isterminal, direction] = P_min_event(z, F, P_min)
    value = F(9) - P_min; isterminal = 1; direction = -1;
end
function [value, isterminal, direction] = combine_events(z, F, event_list)
    value = []; isterminal = []; direction = [];
    for k = 1:length(event_list)
        [v, it, d] = event_list{k}(z, F);
        value = [value; v(:)]; isterminal = [isterminal; it(:)]; direction = [direction; d(:)];
    end
end

%% Rate functions
function [r1, r2] = rates_LHHW(T, alpha, P_CH3OH, P_O2, P_H2O, P_HCHO)
    R = 8.314;
    A_kCH3OH = 1.50e7;  Ea_kCH3OH = 86.00e3;
    A_KCH3OH = 2.60e-4; Ea_KCH3OH = -56.78e3;
    A_KO2    = 1.42e-5; Ea_KO2    = -60.32e3;
    A_KH2O   = 5.50e-7; Ea_KH2O   = -86.45e3;
    A_kCO    = 3.50e2;  Ea_kCO    = 46.00e3;

    kCH3OH = A_kCH3OH * exp(-Ea_kCH3OH/(R*T));
    KCH3OH = A_KCH3OH * exp(-Ea_KCH3OH/(R*T));
    KO2    = A_KO2    * exp(-Ea_KO2/(R*T));
    KH2O   = A_KH2O   * exp(-Ea_KH2O/(R*T));
    kCO    = A_kCO    * exp(-Ea_kCO/(R*T));

    P_O2_sqrt = sqrt(max(P_O2, 0));
    nom1 = alpha*kCH3OH*KCH3OH*P_CH3OH*KO2*P_O2_sqrt;
    den1 = (1 + KCH3OH*P_CH3OH + KH2O*P_H2O) * (1 + KO2*P_O2_sqrt);
    r1   = nom1 / den1;

    nom2 = alpha*kCO*P_HCHO*KO2*P_O2_sqrt;
    den2 = 1 + KO2*P_O2_sqrt;
    r2   = nom2 / den2;
end

function [r1, r2] = rates_powerlaw(T, y_CH3OH, y_O2, y_HCHO, y_H2O)
    R   = 8.314;
    Ea1 = 80.0173e3; Ea2 = 44.1522e4;
    k01 = 3.117346e+05;    k02 = 5.707655e-01;
    k1  = k01 * exp(-Ea1/(R*T));
    k2  = k02 * exp(-Ea2/(R*T));
    a=0.660; b=0.156; c=0.8635; d=0.2149; e=-0.442;
    r1 = k1*(y_CH3OH^a)*(y_O2^b)*(y_H2O^e);
    r2 = k2*(y_HCHO^c)*(y_O2^d);
end

%% Energy balance
function dT_dz = energy_bal(r1, r2, y, n_T, T, Hr1, Hr2)
    a = [40.0460; 29.5260; 33.9330; 29.3420; 29.6750; 29.5560];
    b = [-3.8287e-02; -8.8999e-03; -8.4186e-03; -3.5395e-03;  1.8937e-02; -6.5807e-03];
    c = [ 2.4529e-04;  3.8083e-05;  2.9906e-05;  1.0076e-05;  2.8739e-05;  2.0130e-05];
    d = [-2.1679e-07; -3.2629e-08; -1.7825e-08; -4.3116e-09; -2.0092e-08; -1.2227e-08];
    e = [ 5.9909e-11;  8.8607e-12;  3.6934e-12;  2.5935e-13;  3.4333e-12;  2.2617e-12];
    y_Cp = zeros(6,1);
    for i = 1:6
        y_Cp(i) = y(i)*(a(i) + b(i)*T + c(i)*T^2 + d(i)*T^3 + e(i)*T^4);
    end
    dT_dz = -(r1*Hr1 + r2*Hr2) / (n_T * sum(y_Cp));
end

%% Reactor ODE 
function dni_dz = reactor(t, F, options, reactor_param)
    if any(F(1:6) < 0)
        fprintf('WARNING: Negative flows at z=%.2f\n', t);
        dni_dz = zeros(10,1); return;
    end

    alpha         = reactor_param(1); rho_cat       = reactor_param(2);
    bed_voidage   = reactor_param(3); D_p           = reactor_param(4);
    tube_diameter = reactor_param(5); Hr1           = reactor_param(6);
    Hr2           = reactor_param(7);

    R = 8.314;
    Mr_CH3OH=32.042e-3; Mr_O2=32.00e-3; Mr_HCHO=30.026e-3;
    Mr_H2O=18.00e-3;    Mr_CO=28e-3;    Mr_N2=28.014e-3;
    A = pi*(tube_diameter/2)^2;

    n_CH3OH=F(1); n_O2=F(2); n_HCHO=F(3); n_H2O=F(4);
    n_CO=F(5);    n_N2=F(6);
    P=F(9); P_pascal=P*101325; T=F(10);

    n_T = n_CH3OH+n_O2+n_HCHO+n_H2O+n_CO+n_N2;
    y_CH3OH=n_CH3OH/n_T; y_O2=n_O2/n_T; y_H2O=n_H2O/n_T;
    y_N2=n_N2/n_T;       y_HCHO=n_HCHO/n_T; y_CO=n_CO/n_T;
    y = [y_CH3OH; y_O2; y_H2O; y_N2; y_HCHO; y_CO];

    P_CH3OH=y_CH3OH*P; P_O2=y_O2*P; P_H2O=y_H2O*P; P_HCHO=y_HCHO*P;
    Mr_total = y_CH3OH*Mr_CH3OH + y_O2*Mr_O2 + y_H2O*Mr_H2O + ...
               y_HCHO*Mr_HCHO + y_CO*Mr_CO + y_N2*Mr_N2;
    G       = (Mr_total*n_T) / A;
    rho_mix = (P_pascal*Mr_total) / (R*T);
    dw_dz   = rho_cat*(1-bed_voidage)*A;

    if options(2) == "Yes"
        dP_dz = -(1.75*(G^2)*(1-bed_voidage)) / (D_p*rho_mix*(bed_voidage)^3) / 101325;
    else
        dP_dz = 0;
    end

    if options(1) == "LHHW"
        [r1, r2] = rates_LHHW(T, alpha, P_CH3OH, P_O2, P_H2O, P_HCHO);
    elseif options(1) == "Simple Power"
        [r1, r2] = rates_powerlaw(T, y_CH3OH, y_O2, y_HCHO, y_H2O);
    else
        error('Unknown kinetic model.')
    end

    if options(3) == "Yes"
        dT_dz = energy_bal(r1, r2, y, n_T, T, Hr1, Hr2) * dw_dz;
    else
        dT_dz = 0;
    end

    dni_dz = [-r1*dw_dz;
              -(r1+r2)/2*dw_dz;
              (r1-r2)*dw_dz;
              (r1+r2)*dw_dz;
              r2*dw_dz;
              0;
              A*rho_cat*(1-bed_voidage)*r1;
              A*rho_cat*(1-bed_voidage)*r2;
              dP_dz;
              dT_dz];
end

%% solve_reactor 
function [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
         solve_reactor(initial_conditions, reactor_length, P0, T0, option_limits, options, reactor_param)

    T_max = option_limits{1}; T_max_active = option_limits{2};
    P_min = option_limits{3}; P_min_active = option_limits{4};

    active_events = {};
    if T_max_active == "Yes"
        active_events{end+1} = @(z,F) T_max_event(z, F, T_max);
    end
    if P_min_active == "Yes"
        active_events{end+1} = @(z,F) P_min_event(z, F, P_min);
    end

    if ~isempty(active_events)
        combined_event = @(z,F) combine_events(z, F, active_events);
        opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'Events',combined_event);
    else
        opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
    end

    [z_points, solution, z_stop] = ode45( ...
        @(z,F) reactor(z, F, options, reactor_param), ...
        [0 reactor_length], initial_conditions, opts);

    if ~isempty(z_stop)
        if T_max_active == "Yes" && solution(end,10) >= T_max
            fprintf('Reactor stopped at z = %.4f m (T_max = %g K)\n', z_stop(end), T_max);
        end
        if P_min_active == "Yes" && solution(end,9) <= P_min
            fprintf('Reactor stopped at z = %.4f m (P_min = %g atm)\n', z_stop(end), P_min);
        end
    end

    solution = solution';
    F_CH3OH=solution(1,:); F_O2=solution(2,:);    F_HCHO=solution(3,:);
    F_H2O  =solution(4,:); F_CO=solution(5,:);    F_N2  =solution(6,:);
    extent1=solution(7,:); extent2=solution(8,:);
    P      =solution(9,:); T    =solution(10,:);
end

%% Isothermal heat duty
function [Q_dot_per_length, Q_total] = calc_isothermal_duty(z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, P, T, reactor_param)
    alpha=reactor_param(1); rho_cat=reactor_param(2); bed_voidage=reactor_param(3);
    tube_diameter=reactor_param(5); Hr1=reactor_param(6); Hr2=reactor_param(7);
    A     = pi*(tube_diameter/2)^2;
    dW_dz = rho_cat*(1-bed_voidage)*A;
    Q_dot_per_length = zeros(length(z_points),1);
    for i = 1:length(z_points)
        n_T_i = F_CH3OH(i)+F_O2(i)+F_HCHO(i)+F_H2O(i)+F_CO(i)+F_N2(i);
        P_CH3OH_i=(F_CH3OH(i)/n_T_i)*P(i); P_O2_i=(F_O2(i)/n_T_i)*P(i);
        P_H2O_i =(F_H2O(i)/n_T_i)*P(i);   P_HCHO_i=(F_HCHO(i)/n_T_i)*P(i);
        [r1_i,r2_i] = rates_LHHW(T(i), alpha, P_CH3OH_i, P_O2_i, P_H2O_i, P_HCHO_i);
        Q_dot_per_length(i) = (r1_i*Hr1 + r2_i*Hr2)*dW_dz;
    end
    Q_total = trapz(z_points, Q_dot_per_length);
end

%% Plot function for manual user input
function plot_reactor_results(z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2, ...
                               P, T, extent1, extent2, F0_CH3OH, F0_HCHO, ...
                               stage_boundaries, Q_dot_per_length, reactor_param)

    alpha = reactor_param(1);

    figure("Name","Molar Flowrate of Species Across the Reactor");
    plot(z_points, F_CH3OH, 'b-', 'LineWidth', 2); hold on;
    plot(z_points, F_O2,    'r-', 'LineWidth', 2);
    plot(z_points, F_HCHO,  'g-', 'LineWidth', 2);
    plot(z_points, F_H2O,   'c-', 'LineWidth', 2);
    plot(z_points, F_CO,    'm-', 'LineWidth', 2);
    xlabel('Reactor Length (m)'); ylabel('Molar Flow Rate (mol/s)');
    xlim([0 z_points(end)]); legend('CH_3OH','O_2','HCHO','H_2O','CO');
    title("Molar Flowrate of Species Across the Reactor");
    set(gcf,'color','w'); grid on;

    figure("Name","Pressure Drop Across Reactor");
    plot(z_points, P, 'LineWidth', 2);
    xlabel('Reactor Length (m)'); ylabel('Pressure (atm)');
    xlim([0 z_points(end)]); title("Pressure Drop Across Reactor");
    set(gcf,'color','w'); grid on;

    figure("Name","Extent Of Reaction");
    plot(z_points, extent1, 'b-', 'LineWidth', 2); hold on;
    plot(z_points, extent2, 'r-', 'LineWidth', 2);
    title('Extent of Reaction'); xlim([0 z_points(end)]);
    legend("Reaction 1","Reaction 2"); set(gcf,'color','w'); grid on;

    figure("Name","Temperature Profile");
    plot(z_points, T, 'LineWidth', 2);
    title('Temperature Profile'); legend('T');
    xlabel('Reactor Length (m)'); ylabel('Temperature (K)');
    xlim([0 z_points(end)]); set(gcf,'color','w'); grid on;

    % --- Reaction rates ---
    r1_profile = zeros(length(z_points), 1);
    r2_profile = zeros(length(z_points), 1);
    for i = 1:length(z_points)
        n_T_i     = F_CH3OH(i)+F_O2(i)+F_HCHO(i)+F_H2O(i)+F_CO(i)+F_N2(i);
        y_CH3OH_i = F_CH3OH(i)/n_T_i;  y_O2_i   = F_O2(i)/n_T_i;
        y_H2O_i   = F_H2O(i)/n_T_i;    y_HCHO_i = F_HCHO(i)/n_T_i;
        P_CH3OH_i = y_CH3OH_i*P(i);    P_O2_i   = y_O2_i*P(i);
        P_H2O_i   = y_H2O_i*P(i);      P_HCHO_i = y_HCHO_i*P(i);
        [r1_profile(i), r2_profile(i)] = rates_LHHW(T(i), alpha, P_CH3OH_i, P_O2_i, P_H2O_i, P_HCHO_i);
    end
    figure("Name","Reaction Rates Across Reactor");
    plot(z_points, r1_profile, 'b-', 'LineWidth', 2); hold on;
    plot(z_points, r2_profile, 'r-', 'LineWidth', 2);
    xlabel('Reactor Length (m)'); ylabel('Rate (mol/kg_{cat}/s)');
    legend('r_1 (CH_3OH \rightarrow HCHO)','r_2 (HCHO \rightarrow CO)');
    xlim([0 z_points(end)]); title('Reaction Rates Across the Reactor');
    grid on; set(gcf,'color','w');

    % --- Conversion ---
    X_CH3OH = (F0_CH3OH - F_CH3OH) / F0_CH3OH;
    figure("Name","Conversion of CH3OH Across Reactor");
    plot(z_points, X_CH3OH, 'b-', 'LineWidth', 2);
    for k = 1:length(stage_boundaries)
        xline(stage_boundaries(k), "LineStyle", "--", "LineWidth", 2, "Color", "r");
    end
    xlabel('Reactor Length (m)'); ylabel('Conversion X_{CH_3OH}');
    title('Conversion of CH_3OH Across the Reactor');
    xlim([0 z_points(end)]); ylim([0 1]); set(gcf,'color','w'); grid on;

    % --- HCHO yield ---
    X_HCHO = (F_HCHO - F0_HCHO) / F0_CH3OH;
    figure("Name","HCHO Yield Across Reactor");
    plot(z_points, X_HCHO, 'b-', 'LineWidth', 2);
    for k = 1:length(stage_boundaries)
        xline(stage_boundaries(k), "LineStyle", "--", "LineWidth", 2, "Color", "r");
    end
    xlabel('Reactor Length (m)'); ylabel('Yield Y_{HCHO}');
    title('HCHO Yield Across the Reactor');
    xlim([0 z_points(end)]); ylim([0 1]); set(gcf,'color','w'); grid on;

    % --- Isothermal cooling duty (optional) ---
    if ~isempty(Q_dot_per_length)
        figure("Name","Isothermal Cooling Duty Profile");
        plot(z_points, abs(Q_dot_per_length)/1000, 'r-', 'LineWidth', 2);
        xlabel('Reactor Length (m)'); ylabel('Cooling Duty (kW/m)');
        title('Isothermal Cooling Duty per Unit Length');
        xlim([0 z_points(end)]); set(gcf,'color','w'); grid on;
    end

end
%% ========================================================================
%  GRAPHS FOR STEP 3 ONWARDS
% Independent of user input
%% =========================================================================
%% Step 3 - All Graphs
function step3_graphs(reactor_param)

    % Isothermal, Isobaric Reactor, with allocated inlet conditions
    %
    % --- Operating Conditions --- %
    m_flowrate     = 0.001;
    reactor_length = 10;
    P0             = 1.6;
    T0             = 430;

    % --- Inlet Mole Fractions --- %
    y0_CH3OH = 0.11;
    y0_O2    = 0.06;
    y0_H2O   = 0.02;
    y0_HCHO  = 0.00;
    y0_CO    = 0.00;
    y0_N2    = 1 - y0_CH3OH - y0_O2 - y0_H2O - y0_HCHO - y0_CO;
    inlet    = [y0_CH3OH, y0_O2, y0_H2O, y0_HCHO, y0_CO, y0_N2];

    % --- Option Limits --- %
    T_max        = 673;   T_max_active = "Yes";
    P_min        = 1.1/1.01325;   P_min_active = "Yes";
    option_limits = {T_max, T_max_active, P_min, P_min_active};

    % --- Reactor/catalyst parameters --- %
    alpha       = reactor_param(1);
    rho_cat     = reactor_param(2);
    bed_voidage = reactor_param(3);
    tube_diam   = reactor_param(5);
    A_tube      = pi*(tube_diam/2)^2;
    dW_dz       = rho_cat * (1 - bed_voidage) * A_tube;  % kg_cat/m

    %==========================================================
    % Run LHHW solver
    %==========================================================
    options_LHHW = ["LHHW", "No", "No"];
    [feed_LHHW, ~, F0_CH3OH, ~] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);

    [z_L, F_CH3OH_L, F_O2_L, F_HCHO_L, F_H2O_L, F_CO_L, F_N2_L, ...
     extent1_L, extent2_L, P_L, T_L, ~] = ...
        solve_reactor(feed_LHHW, reactor_length, P0, T0, option_limits, options_LHHW, reactor_param);

    %==========================================================
    % Run Simple Power Law solver
    %==========================================================
    options_PL = ["Simple Power", "No", "No"];
    [feed_PL, ~, ~, ~] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);

    [z_P, F_CH3OH_P, F_O2_P, F_HCHO_P, F_H2O_P, F_CO_P, F_N2_P, ...
     extent1_P, extent2_P, P_P, T_P, ~] = ...
        solve_reactor(feed_PL, reactor_length, P0, T0, option_limits, options_PL, reactor_param);

    %==========================================================
    % Catalyst weight axis: W = dW_dz * z
    %==========================================================
    W_L = z_L * dW_dz;   % [kg_cat]
    W_P = z_P * dW_dz;   % [kg_cat]

    %==========================================================
    % FIGURE 1: Molar Flowrates - LHHW vs Power Law
    % One subplot per species (essential components only)
    %==========================================================
    species_names  = {'CH_3OH', 'O_2', 'HCHO', 'H_2O', 'CO'};
    F_LHHW = {F_CH3OH_L, F_O2_L, F_HCHO_L, F_H2O_L, F_CO_L};
    F_PL   = {F_CH3OH_P, F_O2_P, F_HCHO_P, F_H2O_P, F_CO_P};
    colors = {'b','r','g','c','m'};

    figure("Name","Step 3 - Molar Flowrates: LHHW vs Power Law")
    hold on
    for i = 1:5
        plot(z_L, F_LHHW{i}, '-',  'Color', colors{i}, 'LineWidth', 2, ...
             'DisplayName', [species_names{i} ' (LHHW)'])
        plot(z_P, F_PL{i},   '--', 'Color', colors{i}, 'LineWidth', 2, ...
             'DisplayName', [species_names{i} ' (Power Law)'])
    end
    hold off
    xlabel('Reactor Length (m)')
    ylabel('Molar Flow Rate (mol/s)')
    title(sprintf('Molar Flowrates (Isothermal, Isobaric, T_0 = %d K)', T0))
    legend('Location','bestoutside')
    grid on
    set(gcf,'color','w')

    %==========================================================
    % FIGURE 2: Catalyst Weight vs Extent of Reaction
    % Both reactions, both models
    %==========================================================
    figure("Name","Step 3 - Catalyst Weight vs Extent of Reaction")

    % --- Subplot 1: Reaction 1 extent ---
    subplot(1,2,1)
    plot(W_L, extent1_L, 'b-',  'LineWidth', 2)
    hold on
    plot(W_P, extent1_P, 'b--', 'LineWidth', 2)
    hold off
    xlabel('Catalyst Weight W (kg_{cat})')
    ylabel('\xi_1 (mol/s)')
    title('Rxn 1 (CH_3OH \rightarrow HCHO): Extent of Reaction')
    legend('LHHW','Power Law','Location','best')
    grid on

    % --- Subplot 2: Reaction 2 extent ---
    subplot(1,2,2)
    plot(W_L, extent2_L, 'r-',  'LineWidth', 2)
    hold on
    plot(W_P, extent2_P, 'r--', 'LineWidth', 2)
    hold off
    xlabel('Catalyst Weight W (kg_{cat})')
    ylabel('\xi_2 (mol/s)')
    title('Rxn 2 (HCHO \rightarrow CO): Extent of Reaction')
    legend('LHHW','Power Law','Location','best')
    grid on

    sgtitle(sprintf('Step 3: Catalyst Weight vs Extent of Reaction (T_0 = %d K)', T0))
    
    figure("Name","Step 3 - Catalyst Weight vs Extent of Reaction (R1 vs R2")
    plot(W_L, extent1_L, 'b-',  'LineWidth', 2)
    hold on
    plot(W_L, extent2_L, 'r-',  'LineWidth', 2)
    hold off 
    grid on 
    legend('Reaction 1', 'Reaction 2','Location','best')
    xlabel('Catalyst Weight W (kg_{cat})')
    ylabel('Extent of Reaction (mol s^{-1})')
    title('Extent of Reaction vs Catalyst Weight (LHHW, T_0 = 430 K)')
    %==========================================================
    % Print summary to console
    %==========================================================
    W_total_L = W_L(end);
    W_total_P = W_P(end);
    X_CH3OH_L = (F0_CH3OH - F_CH3OH_L(end)) / F0_CH3OH;
    X_CH3OH_P = (F0_CH3OH - F_CH3OH_P(end)) / F0_CH3OH;

    fprintf('\n=== Step 3 Summary ===\n')
    fprintf('%-25s %10s %10s\n', 'Quantity', 'LHHW', 'Power Law')
    fprintf('%-25s %10.4f %10.4f\n', 'Total catalyst W (kg)',  W_total_L, W_total_P)
    fprintf('%-25s %10.4f %10.4f\n', 'X_CH3OH at outlet',      X_CH3OH_L, X_CH3OH_P)
    fprintf('%-25s %10.6f %10.6f\n', 'F_HCHO out (mol/s)',     F_HCHO_L(end), F_HCHO_P(end))
    fprintf('%-25s %10.6f %10.6f\n', 'F_CO out (mol/s)',       F_CO_L(end),   F_CO_P(end))
    fprintf('%-25s %10.4f %10.4f\n', 'xi_1 at outlet (mol/s)', extent1_L(end), extent1_P(end))
    fprintf('%-25s %10.4f %10.4f\n', 'xi_2 at outlet (mol/s)', extent2_L(end), extent2_P(end))
end
%% Step 4 - All Graphs
function step4_graphs(reactor_param)
    inlet_sup = [0.11, 0.06, 0.02, 0.00, 0.00, 0.81];
    rl_sup    = 5;
    P0_sup    = 1.6;

    %% --- Pressure profiles for different mass flowrates ---
    opt_P    = ["LHHW", "Yes", "Yes"];
    olim_P   = {800, "No", 1.1/1.01325, "Yes"};
    T0_P     = 430;
    m_list   = [0.001, 0.0008, 0.0006];
    labels_P = {'0.001 kg/s', '0.0008 kg/s', '0.0006 kg/s'};

    figure("Name","Pressure Profiles for Different Mass Flowrates");
    hold on;
    for k = 1:3
        [fs_k, ~, ~, ~] = get_initial_conditions(m_list(k), P0_sup, T0_P, 0, inlet_sup);
        [zk, ~, ~, ~, ~, ~, ~, ~, ~, Pk, ~, ~] = ...
            solve_reactor(fs_k, rl_sup, P0_sup, T0_P, olim_P, opt_P, reactor_param);
        plot(zk, Pk, 'LineWidth', 2, 'DisplayName', labels_P{k});
    end
    hold off;
    xlabel('Reactor Length (m)'); ylabel('Pressure (atm)');
    title('Pressure Profile at Different Mass Flowrates');
    legend('Location', 'best'); grid on; set(gcf,'color','w');

    %% --- Reaction Rate vs Total Pressure ---
    alpha   = reactor_param(1);
    T_rate  = 430;                      % temperature [K]
    y_CH3OH = inlet_sup(1);
    y_O2    = inlet_sup(2);
    y_H2O   = inlet_sup(3);
    y_HCHO  = inlet_sup(4);
    P_range = linspace(1.1,1.6 , 200);  % [atm]

    r1_vec = zeros(size(P_range));
    r2_vec = zeros(size(P_range));
    for i = 1:length(P_range)
        P_total     = P_range(i);
        P_CH3OH_i   = y_CH3OH * P_total;
        P_O2_i      = y_O2    * P_total;
        P_H2O_i     = y_H2O   * P_total;
        P_HCHO_i    = y_HCHO  * P_total;
        [r1_vec(i), r2_vec(i)] = rates_LHHW(T_rate, alpha, P_CH3OH_i, P_O2_i, P_H2O_i, P_HCHO_i);
    end

    figure("Name","Reaction Rate vs Total Pressure");
    plot(P_range, r1_vec, 'b-', 'LineWidth', 2); hold on;
    plot(P_range, r2_vec, 'r-', 'LineWidth', 2);
    hold off;
    xlabel('Total Pressure (atm)');
    ylabel('Reaction Rate (mol kg^{-1} s^{-1})');
    title(sprintf('Reaction Rate vs Total Pressure at T = %.0f K', T_rate));
    legend('r_1 (CH_3OH \rightarrow HCHO)','r_2 (HCHO \rightarrow CO)');
    grid on; set(gcf,'color','w');

    %% --- Molar Flowrate: Isobaric vs Non-Isobaric ---
    T0_molar = 430;
    m_molar  = 0.001;
    rl_molar = 2.5;
    olim_molar = {800, "No", 1.1/1.01325, "Yes"};  % P_min off so both run full length
    
    % Isobaric (no pressure drop)
    opt_isobaric    = ["LHHW", "No", "No"];
    [fs_iso, ~, ~, ~] = get_initial_conditions(m_molar, P0_sup, T0_molar, 0, inlet_sup);
    [z_iso, F_CH3OH_iso, F_O2_iso, F_HCHO_iso, F_H2O_iso, F_CO_iso, F_N2_iso, ...
     ~,~,~,~,~] = solve_reactor(fs_iso, rl_molar, P0_sup, T0_molar, olim_molar, opt_isobaric, reactor_param);
    
    % Non-isobaric (with pressure drop)
    opt_nonisobaric = ["LHHW", "Yes", "No"];
    [fs_noniso, ~, ~, ~] = get_initial_conditions(m_molar, P0_sup, T0_molar, 0, inlet_sup);
    [z_noniso, F_CH3OH_noniso, F_O2_noniso, F_HCHO_noniso, F_H2O_noniso, F_CO_noniso, F_N2_noniso, ...
     ~,~,~,~,~] = solve_reactor(fs_noniso, rl_molar, P0_sup, T0_molar, olim_molar, opt_nonisobaric, reactor_param);
    
    figure("Name", "Molar Flowrate Profiles (LHHW Kinetics)");
    hold on
    
    % CH3OH
    plot(z_iso,    F_CH3OH_iso,    '-',  'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'CH3OH isobaric')
    plot(z_noniso, F_CH3OH_noniso, '--', 'Color', 'b', 'LineWidth', 3, 'DisplayName', 'CH3OH non-isobaric')
    
    % O2
    plot(z_iso,    F_O2_iso,    '-',  'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'O2 isobaric')
    plot(z_noniso, F_O2_noniso, '--', 'Color', 'r', 'LineWidth', 3, 'DisplayName', 'O2 non-isobaric')
    
    % HCHO
    plot(z_iso,    F_HCHO_iso,    '-',  'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'HCHO isobaric')
    plot(z_noniso, F_HCHO_noniso, '--', 'Color', 'g', 'LineWidth', 3, 'DisplayName', 'HCHO non-isobaric')
    
    % H2O
    plot(z_iso,    F_H2O_iso,    '-',  'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'H2O isobaric')
    plot(z_noniso, F_H2O_noniso, '--', 'Color', 'c', 'LineWidth', 3, 'DisplayName', 'H2O non-isobaric')
    
    % CO
    plot(z_iso,    F_CO_iso,    '-',  'Color', 'm', 'LineWidth', 1.5, 'DisplayName', 'CO isobaric')
    plot(z_noniso, F_CO_noniso, '--', 'Color', 'm', 'LineWidth', 3, 'DisplayName', 'CO non-isobaric')
    
    hold off
    xlabel('Reactor Length (m)')
    ylabel('Molar Flowrate (mol/s)')
    title('Molar Flowrate Profiles (LHHW Kinetics)')
    legend('Location', 'bestoutside')
    grid on
    xlim([0 rl_molar])
end
%% Step 5 - All Graphs
function step5_graphs(reactor_param)

    %% --- Temperature profiles for 5 LHHW cases ---
    inlet_sup = [0.11, 0.06, 0.02, 0.00, 0.00, 0.81];
    rl_sup    = 4;
    P0_sup    = 1.6;
    opt_T  = ["LHHW", "Yes", "Yes"];
    olim_T = {673, "Yes", 1.1/1.01325, "Yes"};  % T_max = 673 K active

    cases = {0.0007, 430, '430 K, 0.0007 kg/s';
             0.0006, 430, '430 K, 0.0006 kg/s';
             0.0005, 430, '430 K, 0.0005 kg/s';
             0.0006, 460, '460 K, 0.0006 kg/s';
             0.0006, 490, '490 K, 0.0006 kg/s'};

    figure("Name","Temperature Profiles: 5 LHHW Cases");
    hold on;
    for k = 1:size(cases, 1)
        m_k  = cases{k, 1};
        T0_k = cases{k, 2};
        [fs_k, ~, ~, ~] = get_initial_conditions(m_k, P0_sup, T0_k, 0, inlet_sup);
        [zk, ~, ~, ~, ~, ~, ~, ~, ~, ~, Tk, ~] = ...
            solve_reactor(fs_k, rl_sup, P0_sup, T0_k, olim_T, opt_T, reactor_param);
        plot(zk, Tk, 'LineWidth', 1.5, 'DisplayName', cases{k, 3});
    end
    hold off;
    xlabel('Reactor Length (m)'); ylabel('Temperature (K)');
    title('Temperature Profiles for Different Inlet Conditions');
    legend('Location', 'best'); grid on; set(gcf,'color','w');
    xlim([0 rl_sup]);
    %% --- Concentration Profile for one representative case ---
    % Using case: 430 K, 0.0006 kg/s (case 2 — moderate conditions)
    m_conc  = 0.0006;
    T0_conc = 430;
    olim_T = {673, "No", 1.1/1.01325, "Yes"};
    [fs_conc, ~, ~, ~] = get_initial_conditions(m_conc, P0_sup, T0_conc, 0, inlet_sup);
    [z_points, F_CH3OH, F_O2, F_HCHO, F_H2O, F_CO, F_N2,~,~, P_conc, T_conc, ~] = ...
        solve_reactor(fs_conc, rl_sup, P0_sup, T0_conc, olim_T, opt_T, reactor_param);
    P = P_conc; T = T_conc;

    R = 8.314;
    npts = length(z_points);
    C_CH3OH = zeros(npts,1); C_O2 = zeros(npts,1);
    C_HCHO  = zeros(npts,1); C_H2O = zeros(npts,1);
    C_CO    = zeros(npts,1); C_N2  = zeros(npts,1);
    for i = 1:npts
        F_T_i     = F_CH3OH(i)+F_O2(i)+F_HCHO(i)+F_H2O(i)+F_CO(i)+F_N2(i);
        C_total   = (P(i)*101325) / (R*T(i));
        C_CH3OH(i) = (F_CH3OH(i)/F_T_i) * C_total;
        C_O2(i)    = (F_O2(i)   /F_T_i) * C_total;
        C_HCHO(i)  = (F_HCHO(i) /F_T_i) * C_total;
        C_H2O(i)   = (F_H2O(i)  /F_T_i) * C_total;
        C_CO(i)    = (F_CO(i)   /F_T_i) * C_total;
        C_N2(i)    = (F_N2(i)   /F_T_i) * C_total;
    end

    figure('Name', 'Concentration Profiles: 430 K, 0.0006 kg/s');
    plot(z_points, C_CH3OH, 'b-', 'LineWidth', 2); hold on
    plot(z_points, C_O2,    'r-', 'LineWidth', 2)
    plot(z_points, C_HCHO,  'g-', 'LineWidth', 2)
    plot(z_points, C_H2O,   'c-', 'LineWidth', 2)
    plot(z_points, C_CO,    'm-', 'LineWidth', 2)
    hold off
    xlabel('Reactor Length (m)')
    ylabel('Concentration (mol m^{-3})')
    title('Concentration Profiles Across the Reactor (430 K, 0.0006 kg/s, LHHW)')
    legend('CH_3OH','O_2','HCHO','H_2O','CO','Location','best')
    xlim([0 z_points(end)])
    grid on
    set(gcf,'color','w')
end

%% Interstage cooling
function [z_combined, F_CH3OH_combined, F_O2_combined, F_HCHO_combined, ...
          F_H2O_combined, F_CO_combined, F_N2_combined, ...
          P_combined, T_combined, extent1_combined, extent2_combined, ...
          F0_CH3OH, F0_HCHO, stage_boundaries, Q_interstage_total] = ...
          reactor_interstage_cooling(m_flowrate, n_stages, T_subcool, P0, T0, option_limits, options, inlet, reactor_param)

    max_stage_length = 4;
    [IC, ~, F0_CH3OH, F0_HCHO] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);

    z_combined=[]; F_CH3OH_combined=[]; F_O2_combined=[]; F_HCHO_combined=[];
    F_H2O_combined=[]; F_CO_combined=[]; F_N2_combined=[];
    extent1_combined=[]; extent2_combined=[]; P_combined=[]; T_combined=[];
    stage_boundaries=[]; z_offset=0; Q_interstage_total=0;
    s=0;
    P_min        = option_limits{3};
    P_min_active = option_limits{4};

    while true
        s = s+1;
        fprintf('--- Stage %d | T_in=%.1f K | P_in=%.4f atm ---\n', s, IC(10), IC(9));
        [z_s,F_CH3OH_s,F_O2_s,F_HCHO_s,F_H2O_s,F_CO_s,F_N2_s, ...
         extent1_s,extent2_s,P_s,T_s,~] = ...
            solve_reactor(IC, max_stage_length, P0, T0, option_limits, options, reactor_param);
        fprintf('    Stage %d ended | z=%.4f m | T_out=%.1f K | X=%.3f\n', ...
            s, z_s(end), T_s(end), (F0_CH3OH-F_CH3OH_s(end))/F0_CH3OH);

        z_combined=[z_combined; z_s(:)+z_offset];
        F_CH3OH_combined=[F_CH3OH_combined; F_CH3OH_s(:)];
        F_O2_combined   =[F_O2_combined;    F_O2_s(:)];
        F_HCHO_combined =[F_HCHO_combined;  F_HCHO_s(:)];
        F_H2O_combined  =[F_H2O_combined;   F_H2O_s(:)];
        F_CO_combined   =[F_CO_combined;     F_CO_s(:)];
        F_N2_combined   =[F_N2_combined;     F_N2_s(:)];
        extent1_combined=[extent1_combined;  extent1_s(:)];
        extent2_combined=[extent2_combined;  extent2_s(:)];
        P_combined=[P_combined; P_s(:)]; T_combined=[T_combined; T_s(:)];
        z_offset = z_offset + z_s(end);

        if (P_min_active=="Yes" && P_s(end)<=P_min)
            fprintf('    P_min reached at stage %d.\n', s); break
        end
        if s >= n_stages
            fprintf('    Max stages (%d) reached.\n', n_stages); break
        end

        stage_boundaries=[stage_boundaries; z_offset];
        n_T_s=F_CH3OH_s(end)+F_O2_s(end)+F_HCHO_s(end)+F_H2O_s(end)+F_CO_s(end)+F_N2_s(end);
        y_CH3OH_s=F_CH3OH_s(end)/n_T_s; y_O2_s=F_O2_s(end)/n_T_s;
        y_H2O_s =F_H2O_s(end)/n_T_s;   y_N2_s=F_N2_s(end)/n_T_s;
        y_HCHO_s=F_HCHO_s(end)/n_T_s;  y_CO_s=F_CO_s(end)/n_T_s;
        Cp_s = sum(Cp_mean(T_subcool, T_s(end), y_CH3OH_s, y_O2_s, y_H2O_s, y_N2_s, y_HCHO_s, y_CO_s));
        Q_stage_s = n_T_s*Cp_s*(T_subcool - T_s(end));
        Q_interstage_total = Q_interstage_total + Q_stage_s;
        fprintf('    Cooling duty stage %d: %.2f W\n', s, Q_stage_s);
        IC = [F_CH3OH_s(end),F_O2_s(end),F_HCHO_s(end),F_H2O_s(end), ...
              F_CO_s(end),F_N2_s(end),extent1_s(end),extent2_s(end),P_s(end),T_subcool];
    end
    fprintf('\n=== Interstage complete: %d stages, total length=%.4f m ===\n', s, z_offset);
end

%% Cold-shot mixing
function mix_out = mix_coldshot(cs_opts, out_opts)
    ycs_CH3OH=cs_opts(1); ycs_O2=cs_opts(2);  ycs_H2O=cs_opts(3);
    ycs_HCHO =cs_opts(4); ycs_CO=cs_opts(5);  ycs_N2=cs_opts(6);
    n0_CH3OH_cs=cs_opts(7); n0_O2_cs=cs_opts(8);  n0_HCHO_cs=cs_opts(9);
    n0_H2O_cs  =cs_opts(10); n0_CO_cs=cs_opts(11); n0_N2_cs =cs_opts(12);
    nT0_cs=cs_opts(17); P_cs=cs_opts(15); T_cs=cs_opts(16);

    F_CH3OH_out=out_opts(2); F_O2_out=out_opts(3); F_HCHO_out=out_opts(4);
    F_H2O_out  =out_opts(5); F_CO_out=out_opts(6); F_N2_out  =out_opts(7);
    P_out=out_opts(10); T_out=out_opts(11);

    Ftot_out  = F_CH3OH_out+F_O2_out+F_HCHO_out+F_H2O_out+F_CO_out+F_N2_out;
    yout_CH3OH=F_CH3OH_out/Ftot_out; yout_O2  =F_O2_out/Ftot_out;
    yout_HCHO =F_HCHO_out/Ftot_out;  yout_H2O =F_H2O_out/Ftot_out;
    yout_CO   =F_CO_out/Ftot_out;    yout_N2  =F_N2_out/Ftot_out;

    Cp_cs  = sum(Cp_mean(T_cs,T_out,ycs_CH3OH, ycs_O2, ycs_H2O, ycs_N2, ycs_HCHO, ycs_CO));
    Cp_out = sum(Cp_mean(T_cs,T_out,yout_CH3OH,yout_O2,yout_H2O,yout_N2,yout_HCHO,yout_CO));
    T_mix  = (Ftot_out*Cp_out*T_out + nT0_cs*Cp_cs*T_cs) / (Ftot_out*Cp_out + nT0_cs*Cp_cs);

    mix_out = [F_CH3OH_out+n0_CH3OH_cs, F_O2_out+n0_O2_cs, F_HCHO_out+n0_HCHO_cs, ...
               F_H2O_out+n0_H2O_cs,     F_CO_out+n0_CO_cs, F_N2_out+n0_N2_cs, ...
               0, 0, P_out, T_mix];
end

%% Cold-shot multi-stage
function [z_combined, F_CH3OH_combined, F_O2_combined, F_HCHO_combined, ...
          F_H2O_combined, F_CO_combined, F_N2_combined, ...
          P_combined, T_combined, extent1_combined, extent2_combined, ...
          F0_CH3OH, F0_HCHO, stage_boundaries] = ...
    coldshot_reactor_multistage(m_flowrate, P0, T0, reactor_length, split_ratio, n_stages, ...
                                option_limits, options, inlet, reactor_param)

    split_per_stage = split_ratio / n_stages;
    [feed_stream, ~, F0_CH3OH, F0_HCHO] = get_initial_conditions(m_flowrate, P0, T0, split_ratio, inlet);
    [~, coldshot_one, ~, ~] = get_initial_conditions(m_flowrate, P0, T0, split_per_stage, inlet);

    z_combined=[]; F_CH3OH_combined=[]; F_O2_combined=[]; F_HCHO_combined=[];
    F_H2O_combined=[]; F_CO_combined=[]; F_N2_combined=[];
    extent1_combined=[]; extent2_combined=[]; P_combined=[]; T_combined=[];
    stage_boundaries=[]; z_offset=0;
    P_min=option_limits{3}; P_min_active=option_limits{4};
    IC = feed_stream;

    for s = 1:n_stages
        fprintf('--- Cold-Shot Stage %d | T_in=%.1f K | P_in=%.4f atm ---\n',s,IC(10),IC(9));
        [z_s,F_CH3OH_s,F_O2_s,F_HCHO_s,F_H2O_s,F_CO_s,F_N2_s, ...
         extent1_s,extent2_s,P_s,T_s,~] = ...
            solve_reactor(IC, reactor_length, P0, T0, option_limits, options, reactor_param);
        fprintf('    Stage %d | z=%.4f m | T_out=%.1f K | X=%.3f\n', ...
            s, z_s(end), T_s(end), (F0_CH3OH-F_CH3OH_s(end))/F0_CH3OH);

        z_combined=[z_combined; z_s(:)+z_offset];
        F_CH3OH_combined=[F_CH3OH_combined; F_CH3OH_s(:)];
        F_O2_combined   =[F_O2_combined;    F_O2_s(:)];
        F_HCHO_combined =[F_HCHO_combined;  F_HCHO_s(:)];
        F_H2O_combined  =[F_H2O_combined;   F_H2O_s(:)];
        F_CO_combined   =[F_CO_combined;     F_CO_s(:)];
        F_N2_combined   =[F_N2_combined;     F_N2_s(:)];
        extent1_combined=[extent1_combined;  extent1_s(:)];
        extent2_combined=[extent2_combined;  extent2_s(:)];
        P_combined=[P_combined; P_s(:)]; T_combined=[T_combined; T_s(:)];
        z_offset = z_offset + z_s(end);

        if (P_min_active=="Yes" && P_s(end)<=P_min)
            fprintf('    P_min reached at stage %d.\n',s); break
        end
        if s == n_stages
            fprintf('    All %d cold-shot stages complete.\n',n_stages); break
        end

        stage_boundaries=[stage_boundaries; z_offset];
        outlet = get_outlet_conditions(z_s,F_CH3OH_s,F_O2_s,F_HCHO_s,F_H2O_s, ...
                                       F_CO_s,F_N2_s,extent1_s,extent2_s,P_s,T_s);
        mix_out = mix_coldshot(coldshot_one, outlet);
        fprintf('    Cold-shot injected | T_mix=%.1f K\n', mix_out(10));
        IC = mix_out;
    end
    fprintf('\n=== Multi-stage cold-shot: %d stages, total length=%.4f m ===\n',n_stages,z_offset);
end

%% Number of tubes
function num_tubes = calc_num_tubes(F_HCHO, goal_fr)

    num_tubes = goal_fr / F_HCHO(end);
end
%% Sensitivity Analysis (function)
function sensitivity_analysis()
    %% Operating Conditions of optimal isothermal reactor
    % --- Operating Conditions --- %
    m_flowrate      = 0.0013;
    reactor_length  = 0.8734;
    P0              = 1.6;
    T0              = 668;
    
    % --- Inlet Mole Fractions --- %
    y0_CH3OH = 0.11;
    y0_O2    = 0.06;
    y0_H2O   = 0.02;
    y0_HCHO  = 0.00;
    y0_CO    = 0.00;
    y0_N2    = 1 - y0_CH3OH - y0_O2 - y0_H2O - y0_HCHO - y0_CO;
    
    % --- Catalyst/Reactor Parameters --- %
    alpha         = 0.44;
    rho_cat       = 1500;
    bed_voidage   = 0.4;
    D_p           = 7.04e-3;
    tube_diameter = 20e-3;
    
    % --- Reaction Enthalpies --- %
    Hr1 = -159e3;
    Hr2 = -200e3;
    
    % --- Solver Options --- %
    Kinetic_type    = "LHHW";
    Pressure_drop   = "Yes";
    Temp_dependence = "No";
    
    % --- Event Limits --- %
    T_max        = 673;
    T_max_active = "No";
    P_min        = 1.1/1.01325;
    P_min_active = "Yes";

    inlet         = [y0_CH3OH, y0_O2, y0_H2O, y0_HCHO, y0_CO, y0_N2];
    reactor_param = [alpha, rho_cat, bed_voidage, D_p, tube_diameter, Hr1, Hr2];
    option_limits = {T_max, T_max_active, P_min, P_min_active};
    %% Isothermal Reactor Sensitivity Analysis : Temperature 
    T0_range = T0 * [0.8, 0.9, 1.0, 1.1, 1.2];
    percent_inc = [-20,-10,0,+10,+20];
    figure("Name","Sensitivity Analysis : Inlet Temperature")
    for j = 1:length(T0_range)
        Temp_dependence_j = "No";
        options_j = [Kinetic_type Pressure_drop Temp_dependence_j];
        [feed_stream, ~, F0_CH3OH_j, ~] = get_initial_conditions(m_flowrate, P0, T0_range(j), 0, inlet);
        [z_j, F_CH3OH, F_O2, F_HCHO_j, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
            solve_reactor(feed_stream, reactor_length, P0, T0_range(j), option_limits, options_j, reactor_param);
        plot(z_j, F_HCHO_j/F0_CH3OH_j, 'LineWidth', 2, ...
             'DisplayName', sprintf('T_0 = %.0f K (%+.0f%%)', T0_range(j), percent_inc(j)))
        hold on 
    end 
    hold off
    xlabel('Reactor Length (m)')
    ylabel('HCHO Yield  F_{HCHO}/F_{0,CH_3OH}')
    title('Sensitivity Analysis: Inlet Temperature')
    legend('show', 'Location', 'best')
    xlim([0, reactor_length])
    set(gcf, 'color', 'w')
    grid on
    fprintf("Sensitivity Analysis : Inlet Temperature Completed\n")
    
    %% Isothermal Reactor Sensitivity Analysis : Mass Flowrate 
    mflow_range = m_flowrate * [0.8, 0.9, 1.0, 1.1, 1.2];
    percent_inc = [-20,-10,0,+10,+20];
    figure("Name","Sensitivity Analysis : Inlet Mass Flowrate")
    for j = 1:length(mflow_range)
        Temp_dependence_j = "No";
        options_j = [Kinetic_type Pressure_drop Temp_dependence_j];
        [feed_stream, ~, F0_CH3OH_j, ~] = get_initial_conditions(mflow_range(j), P0, T0, 0, inlet);
        [z_j, F_CH3OH, F_O2, F_HCHO_j, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
            solve_reactor(feed_stream, reactor_length, P0, T0, option_limits, options_j, reactor_param);
        plot(z_j, F_HCHO_j/F0_CH3OH_j, 'LineWidth', 2, ...
             'DisplayName', sprintf('ṁ = %.5f kg/s (%+.0f%%)', mflow_range(j), percent_inc(j)))
        hold on 
    end 
    hold off
    xlabel('Reactor Length (m)')
    ylabel('HCHO Yield  F_{HCHO}/F_{0,CH_3OH}')
    title('Sensitivity Analysis: Inlet Mass Flowrate')
    legend('show', 'Location', 'best')
    xlim([0, reactor_length])
    set(gcf, 'color', 'w')
    grid on
    fprintf("Sensitivity Analysis : Mass Flowrate Completed\n")
    %% Isothermal Reactor Sensitivity Analysis : Alpha
    alpha_range = alpha * [0.8, 0.9, 1.0, 1.1, 1.2];
    percent_inc = [-20,-10,0,+10,+20];
    figure("Name","Sensitivity Analysis : Alpha")
    for j = 1:length(alpha_range)
        reactor_param_j    = reactor_param;
        reactor_param_j(1) = alpha_range(j);
        Temp_dependence_j = "No";
        options_j = [Kinetic_type Pressure_drop Temp_dependence_j];
        [feed_stream, ~, F0_CH3OH_j, ~] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);
        [z_j, F_CH3OH, F_O2, F_HCHO_j, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
            solve_reactor(feed_stream, reactor_length, P0, T0, option_limits, options_j, reactor_param_j);
        plot(z_j, F_HCHO_j/F0_CH3OH_j, 'LineWidth', 2, ...
             'DisplayName', sprintf('\\alpha = %.3f (%+.0f%%)', alpha_range(j), percent_inc(j)))
        hold on
    end
    hold off
    xlabel('Reactor Length (m)')
    ylabel('HCHO Yield  F_{HCHO}/F_{0,CH_3OH}')
    title('Sensitivity Analysis: \alpha')
    legend('show', 'Location', 'best')
    xlim([0, reactor_length])
    set(gcf, 'color', 'w')
    grid on
    fprintf("Sensitivity Analysis : Alpha Completed\n")
    %% Isothermal Reactor Sensitivity Analysis : Inlet O2 Composition
    y0_O2_range = y0_O2 * [0.8, 0.9, 1.0, 1.1, 1.2];
    percent_inc = [-20,-10,0,+10,+20];
    figure("Name","Sensitivity Analysis : Inlet O2 Composition")
    for j = 1:length(y0_O2_range)
        inlet_j    = inlet;
        inlet_j(2) = y0_O2_range(j);
        inlet_j(6) = 1 - y0_CH3OH - y0_O2_range(j) - y0_H2O - y0_HCHO - y0_CO;
        Temp_dependence_j = "No";
        options_j = [Kinetic_type Pressure_drop Temp_dependence_j];
        [feed_stream, ~, F0_CH3OH_j, ~] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet_j);
        [z_j, F_CH3OH, F_O2, F_HCHO_j, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
            solve_reactor(feed_stream, reactor_length, P0, T0, option_limits, options_j, reactor_param);
        plot(z_j, F_HCHO_j/F0_CH3OH_j, 'LineWidth', 2, ...
             'DisplayName', sprintf('y_{O2} = %.4f (%+.0f%%)', y0_O2_range(j), percent_inc(j)))
        hold on
    end
    hold off
    xlabel('Reactor Length (m)')
    ylabel('HCHO Yield  F_{HCHO}/F_{0,CH_3OH}')
    title('Sensitivity Analysis: Inlet O_2 Composition')
    legend('show', 'Location', 'best')
    xlim([0, reactor_length])
    set(gcf, 'color', 'w')
    grid on
    
    fprintf("Sensitivity Analysis : Inlet O_2 Composition Completed\n")
    %% Isothermal Reactor Sensitivity Analysis : Tube Diameter
    d_range = tube_diameter * [0.8, 0.9, 1.0, 1.1, 1.2];
    percent_inc = [-20,-10,0,+10,+20];
    figure("Name","Sensitivity Analysis : Tube Diameter")
    for j = 1:length(d_range)
        Temp_dependence_j = "No";
        options_j = [Kinetic_type Pressure_drop Temp_dependence_j];
        reactor_param_j = [alpha, rho_cat, bed_voidage, D_p, d_range(j), Hr1, Hr2];
        [feed_stream, ~, F0_CH3OH_j, ~] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);
        [z_j, F_CH3OH, F_O2, F_HCHO_j, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
            solve_reactor(feed_stream, reactor_length, P0, T0, option_limits, options_j, reactor_param_j);
        plot(z_j, F_HCHO_j/F0_CH3OH_j, 'LineWidth', 2, ...
             'DisplayName', sprintf('d = %.1fmm (%+.0f%%)', d_range(j)*1000, percent_inc(j)))
        hold on 
    end 
    hold off
    xlabel('Reactor Length (m)')
    ylabel('HCHO Yield  F_{HCHO}/F_{0,CH_3OH}')
    title('Sensitivity Analysis: Tube Diameter')
    legend('show', 'Location', 'best')
    xlim([0, reactor_length])
    set(gcf, 'color', 'w')
    grid on
    
    fprintf("Sensitivity Analysis : Tube Diameter Completed\n")
    %% Isothermal Reactor Sensitivity Analysis : Bed Voidage
    bed_voidage_range = bed_voidage * [0.8, 0.9, 1.0, 1.1, 1.2];
    percent_inc = [-20,-10,0,+10,+20];
    figure("Name","Sensitivity Analysis : Bed Voidage")
    for j = 1:length(bed_voidage_range)
        Temp_dependence_j = "No";
        options_j = [Kinetic_type Pressure_drop Temp_dependence_j];
        reactor_param_j = [alpha, rho_cat, bed_voidage_range(j), D_p, tube_diameter, Hr1, Hr2];
        [feed_stream, ~, F0_CH3OH_j, ~] = get_initial_conditions(m_flowrate, P0, T0, 0, inlet);
        [z_j, F_CH3OH, F_O2, F_HCHO_j, F_H2O, F_CO, F_N2, extent1, extent2, P, T, z_stop] = ...
            solve_reactor(feed_stream, reactor_length, P0, T0, option_limits, options_j, reactor_param_j);
        plot(z_j, F_HCHO_j/F0_CH3OH_j, 'LineWidth', 2, ...
             'DisplayName', sprintf('bed voidage = %.2f (%+.0f%%)', bed_voidage_range(j), percent_inc(j)))
        hold on 
    end 
    hold off
    xlabel('Reactor Length (m)')
    ylabel('HCHO Yield  F_{HCHO}/F_{0,CH_3OH}')
    title('Sensitivity Analysis: Bed Voidage')
    legend('show', 'Location', 'best')
    xlim([0, reactor_length])
    set(gcf, 'color', 'w')
    grid on
    
    fprintf("Sensitivity Analysis : Bed Voidage Completed\n")
end
%% Heatmap for conversion and no. of tubes
function heatmap_conversion_tubes(reactor_param, inlet, P0, reactor_length, goal_fr)
    % heatmap_conversion_tubes - Heatmap of CH3OH conversion and number of tubes
    %
    % INPUTS:
    %   reactor_param   - [alpha, rho_cat, bed_voidage, D_p, tube_diameter, Hr1, Hr2]
    %   inlet           - [y0_CH3OH, y0_O2, y0_H2O, y0_HCHO, y0_CO, y0_N2]
    %   P0              - inlet pressure [atm]
    %   reactor_length  - reactor length [m]
    %   goal_fr         - target HCHO molar flowrate for tube count [mol/s]

    % --- Grid --- %
    m_flow_range = linspace(0.001, 0.0025, 100);
    T0_range     = linspace(600, 750, 100);

    % --- Solver options: isothermal, with pressure drop --- %
    option_heatmap        = ["LHHW", "Yes", "No"];
    option_limits_heatmap = {1000, "No", 1.1, "Yes"};

    % --- Preallocate --- %
    X_CH3OH_heatmap  = zeros(length(m_flow_range), length(T0_range));
    no_tubes_heatmap = zeros(length(m_flow_range), length(T0_range));

    total_iters  = length(m_flow_range) * length(T0_range);
    current_iter = 0;
    wb = waitbar(0, 'Running heatmap...', 'Name', 'Heatmap Progress');

    for i = 1:length(m_flow_range)
        for j = 1:length(T0_range)
            try
                [feed_stream_hm, ~, F0_CH3OH_hm, ~] = ...
                    get_initial_conditions(m_flow_range(i), P0, T0_range(j), 0, inlet);

                [~, F_CH3OH_hm, ~, F_HCHO_hm, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                    solve_reactor(feed_stream_hm, reactor_length, P0, T0_range(j), ...
                                  option_limits_heatmap, option_heatmap, reactor_param);

                X_CH3OH_heatmap(i,j) = (F0_CH3OH_hm - F_CH3OH_hm(end)) / F0_CH3OH_hm;

                if F_HCHO_hm(end) > 0
                    no_tubes_heatmap(i,j) = goal_fr / F_HCHO_hm(end);
                else
                    no_tubes_heatmap(i,j) = NaN;
                end

            catch
                X_CH3OH_heatmap(i,j)  = NaN;
                no_tubes_heatmap(i,j) = NaN;
            end

            current_iter = current_iter + 1;
            waitbar(current_iter / total_iters, wb, ...
                sprintf('%.1f%% complete', 100*current_iter/total_iters));
        end
    end
    close(wb);

    % ----------------------------------------------------------------
    % FIGURE 1: Conversion Heatmap
    % ----------------------------------------------------------------
    figure("Name", "Conversion Heatmap: Mass Flowrate vs Inlet Temperature")
    imagesc(T0_range, m_flow_range * 1000, X_CH3OH_heatmap)
    set(gca, 'YDir', 'normal')
    xlabel('Inlet Temperature T_0 (K)')
    ylabel('Mass Flowrate (g s^{-1})')
    title(sprintf('HCHO Yield Y_{HCHO}\n(Isothermal, P_0 = %.1f atm)', ...
          P0))
    clim([0 1])
    c = colorbar;
    c.Label.String = 'Y_{HCHO}';
    colormap(gca, "turbo")
    set(gcf, 'color', 'w')

    % ----------------------------------------------------------------
    % FIGURE 2: Number of Tubes Heatmap
    % ----------------------------------------------------------------
    figure("Name", "Number of Tubes Heatmap: Mass Flowrate vs Inlet Temperature")
    imagesc(T0_range, m_flow_range * 1000, no_tubes_heatmap)
    set(gca, 'YDir', 'normal')
    xlabel('Inlet Temperature T_0 (K)')
    ylabel('Mass Flowrate (g s^{-1})')
    title(sprintf('Number of Tubes to Achieve %.1f mol/s HCHO\n(Isothermal, P_0 = %.1f atm)', ...
          goal_fr, P0))
    c = colorbar;
    c.Label.String = 'Number of Tubes';
    colormap(gca, "turbo")
    set(gcf, 'color', 'w')

    valid_tubes = no_tubes_heatmap(isfinite(no_tubes_heatmap));
    if ~isempty(valid_tubes)
        clim([min(valid_tubes), prctile(valid_tubes, 95)])
    end

    % ----------------------------------------------------------------
    % Console summary
    % ----------------------------------------------------------------
    [max_X, idx_X] = max(X_CH3OH_heatmap(:));
    [row_X, col_X] = ind2sub(size(X_CH3OH_heatmap), idx_X);
    [min_T, idx_T] = min(no_tubes_heatmap(:));
    [row_T, col_T] = ind2sub(size(no_tubes_heatmap), idx_T);

    fprintf('\n=== Heatmap Summary ===\n')
    fprintf('Max conversion    = %.4f\n',      max_X)
    fprintf('  at T0           = %.1f K\n',    T0_range(col_X))
    fprintf('  at m_flow       = %.5f kg/s\n', m_flow_range(row_X))
    fprintf('Min tubes needed  = %.0f\n',      min_T)
    fprintf('  at T0           = %.1f K\n',    T0_range(col_T))
    fprintf('  at m_flow       = %.5f kg/s\n', m_flow_range(row_T))
end

%% Pareto Front plot for optimisation

function pareto_conversion_tubes(reactor_param, inlet, P0, reactor_length, goal_fr)
    % pareto_conversion_tubes - Pareto front: maximise conversion, minimise tubes
    %
    % INPUTS:
    %   reactor_param   - [alpha, rho_cat, bed_voidage, D_p, tube_diameter, Hr1, Hr2]
    %   inlet           - [y0_CH3OH, y0_O2, y0_H2O, y0_HCHO, y0_CO, y0_N2]
    %   P0              - inlet pressure [atm]
    %   reactor_length  - reactor length [m]
    %   goal_fr         - target HCHO molar flowrate [mol/s]

    option_pareto        = ["LHHW", "Yes", "No"];
    option_limits_pareto = {1000, "No", 1.1, "Yes"};

    lb = [0.001,  550];
    ub = [0.0025, 673];

    obj_fun_2d = @(x) reactor_objectives_2d(x, reactor_length, P0, ...
                                             option_limits_pareto, option_pareto, ...
                                             inlet, reactor_param, goal_fr);

    ga_opts = optimoptions('gamultiobj', ...
        'PopulationSize',   150,   ...
        'MaxGenerations',   150,   ...
        'ParetoFraction',   0.5,   ...
        'FunctionTolerance', 1e-4, ...
        'Display',          'iter');

    [x_pareto_2d, f_pareto_2d] = gamultiobj(obj_fun_2d, 2, [], [], [], [], lb, ub, ga_opts);

    % --- Unpack and sort by conversion --- %
    X_pareto_2d     = -f_pareto_2d(:, 1);
    tubes_pareto_2d =  f_pareto_2d(:, 2);

    [~, sort_idx] = sort(X_pareto_2d);
    X_sorted      = X_pareto_2d(sort_idx);
    tubes_sorted  = tubes_pareto_2d(sort_idx);
    mflow_sorted  = x_pareto_2d(sort_idx, 1);
    T0_sorted     = x_pareto_2d(sort_idx, 2);

    % ----------------------------------------------------------------
    % FIGURE: Pareto Front
    % ----------------------------------------------------------------
    figure("Name", "2D Pareto Front: Yield vs Tubes")
    plot(X_sorted, tubes_sorted, 'k-o', ...
        'LineWidth',       1.5, ...
        'MarkerSize',      5,   ...
        'MarkerFaceColor', 'k')
    xlabel('HCHO Yield  Y_{HCHO}', 'FontSize', 12)
    ylabel('Number of Tubes',               'FontSize', 12)
    title(sprintf('Pareto Front: Minimise Tubes vs Maximise Yield\n(Isothermal, P_0 = %.1f atm)', ...
          P0), 'FontSize', 13)
    grid on
    set(gcf, 'color', 'w')

    % ----------------------------------------------------------------
    % Console table
    % ----------------------------------------------------------------
    fprintf('\n=== Pareto Front Summary ===\n')
    fprintf('%-6s %-12s %-12s %-14s %-10s\n', ...
        'Point', 'X_CH3OH', 'Num Tubes', 'm_flow (g/s)', 'T0 (K)')
    fprintf('%s\n', repmat('-', 1, 58))
    for k = 1:length(X_sorted)
        fprintf('%-6d %-12.4f %-12.0f %-14.4f %-10.1f\n', ...
            k, X_sorted(k), tubes_sorted(k), mflow_sorted(k)*1000, T0_sorted(k));
    end
end

% ----------------------------------------------------------------
% Objective function (local helper)
% ----------------------------------------------------------------
function objectives = reactor_objectives_2d(x, reactor_length, P0, ...
                                             option_limits, options, ...
                                             inlet, reactor_param, goal_fr)
    m_flow = x(1);
    T0     = x(2);

    try
        [fs, ~, F0_CH3OH, ~] = get_initial_conditions(m_flow, P0, T0, 0, inlet);

        [~, ~, ~, F_HCHO, ~, ~, ~, ~, ~, ~, ~, ~] = ...
            solve_reactor(fs, reactor_length, P0, T0, option_limits, options, reactor_param);

        X       = F_HCHO(end) / F0_CH3OH;
        n_tubes = goal_fr / F_HCHO(end);

        objectives = [-X, n_tubes];

    catch
        objectives = [0, 1e6];
    end
end

function reactor_comparision_graph(m_flowrate,P0,inlet,reactor_param,reactor_length)
    %% =========================================================
    %% Comparison : Isothermal vs Interstage vs Cold Shot
    %% =========================================================
    
    % --- Shared fixed parameters --- %
    cmp_m_flowrate    = m_flowrate;
    cmp_P0            = P0;
    cmp_reactor_param = reactor_param;
    cmp_inlet         = inlet;
    cmp_P_min         = 1.1;
    cmp_P_min_active  = "Yes";
    cmp_Kinetic_type  = "LHHW";
    cmp_Pressure_drop = "Yes";
    
    % --- Cold Shot parameters --- %
    cmp_cs_T0      = 430;
    cmp_cs_split   = 0.5;
    cmp_cs_stages  = 4;
    cmp_cs_Tmax    = 673;
    cmp_cs_limits  = {cmp_cs_Tmax, "Yes", cmp_P_min, cmp_P_min_active};
    cmp_cs_options = [cmp_Kinetic_type, cmp_Pressure_drop, "Yes"];
    
    % --- Interstage parameters --- %
    cmp_isc_T0     = 560;
    cmp_isc_Tmax   = 673;
    cmp_isc_stages = 8;
    cmp_isc_Tsub   = 610;
    cmp_isc_limits  = {cmp_isc_Tmax, "Yes", cmp_P_min, cmp_P_min_active};
    cmp_isc_options = [cmp_Kinetic_type, cmp_Pressure_drop, "Yes"];
    
    % --- Isothermal parameters --- %
    cmp_iso_T0      = 666;
    cmp_iso_limits  = {720, "No", cmp_P_min, cmp_P_min_active};
    cmp_iso_options = [cmp_Kinetic_type, cmp_Pressure_drop, "No"];
    
    %% --- Run Cold Shot --- %
    fprintf('\n--- Running Cold Shot comparison ---\n')
    [cmp_cs_z, cmp_cs_FCH3OH, cmp_cs_FO2, cmp_cs_FHCHO, cmp_cs_FH2O, cmp_cs_FCO, cmp_cs_FN2, ...
     cmp_cs_P, cmp_cs_T, ~, ~, cmp_cs_F0CH3OH, ~, cmp_cs_boundaries] = ...
        coldshot_reactor_multistage( ...
            cmp_m_flowrate, cmp_P0, cmp_cs_T0, reactor_length, ...
            cmp_cs_split, cmp_cs_stages, ...
            cmp_cs_limits, cmp_cs_options, cmp_inlet, cmp_reactor_param);
    
    %% --- Run Interstage --- %
    fprintf('\n--- Running Interstage comparison ---\n')
    [cmp_isc_z, cmp_isc_FCH3OH, cmp_isc_FO2, cmp_isc_FHCHO, cmp_isc_FH2O, cmp_isc_FCO, cmp_isc_FN2, ...
     cmp_isc_P, cmp_isc_T, ~, ~, cmp_isc_F0CH3OH, ~, cmp_isc_boundaries, ~] = ...
        reactor_interstage_cooling( ...
            cmp_m_flowrate, cmp_isc_stages, cmp_isc_Tsub, ...
            cmp_P0, cmp_isc_T0, ...
            cmp_isc_limits, cmp_isc_options, cmp_inlet, cmp_reactor_param);
    
    %% --- Run Isothermal --- %
    fprintf('\n--- Running Isothermal comparison ---\n')
    [cmp_iso_feed, ~, cmp_iso_F0CH3OH, ~] = ...
        get_initial_conditions(cmp_m_flowrate, cmp_P0, cmp_iso_T0, 0, cmp_inlet);
    [cmp_iso_z, cmp_iso_FCH3OH, cmp_iso_FO2, cmp_iso_FHCHO, cmp_iso_FH2O, cmp_iso_FCO, cmp_iso_FN2, ...
     ~, ~, cmp_iso_P, cmp_iso_T, ~] = ...
        solve_reactor(cmp_iso_feed, reactor_length, cmp_P0, cmp_iso_T0, ...
                      cmp_iso_limits, cmp_iso_options, cmp_reactor_param);
    
    %% --- Compute reaction rate profiles (mol/m_reactor/s) --- %
    A_cmp    = pi * (cmp_reactor_param(5)/2)^2;
    dWdz_cmp = cmp_reactor_param(2) * (1 - cmp_reactor_param(3)) * A_cmp;
    
    % Cold Shot rates
    cmp_cs_r1 = zeros(length(cmp_cs_z), 1);
    cmp_cs_r2 = zeros(length(cmp_cs_z), 1);
    for i = 1:length(cmp_cs_z)
        n_T_i     = cmp_cs_FCH3OH(i)+cmp_cs_FO2(i)+cmp_cs_FHCHO(i)+cmp_cs_FH2O(i)+cmp_cs_FCO(i)+cmp_cs_FN2(i);
        P_CH3OH_i = (cmp_cs_FCH3OH(i)/n_T_i) * cmp_cs_P(i);
        P_O2_i    = (cmp_cs_FO2(i)/n_T_i)    * cmp_cs_P(i);
        P_H2O_i   = (cmp_cs_FH2O(i)/n_T_i)   * cmp_cs_P(i);
        P_HCHO_i  = (cmp_cs_FHCHO(i)/n_T_i)  * cmp_cs_P(i);
        [r1_i, r2_i] = rates_LHHW(cmp_cs_T(i), cmp_reactor_param(1), P_CH3OH_i, P_O2_i, P_H2O_i, P_HCHO_i);
        cmp_cs_r1(i) = r1_i * dWdz_cmp;
        cmp_cs_r2(i) = r2_i * dWdz_cmp;
    end
    
    % Interstage rates
    cmp_isc_r1 = zeros(length(cmp_isc_z), 1);
    cmp_isc_r2 = zeros(length(cmp_isc_z), 1);
    for i = 1:length(cmp_isc_z)
        n_T_i     = cmp_isc_FCH3OH(i)+cmp_isc_FO2(i)+cmp_isc_FHCHO(i)+cmp_isc_FH2O(i)+cmp_isc_FCO(i)+cmp_isc_FN2(i);
        P_CH3OH_i = (cmp_isc_FCH3OH(i)/n_T_i) * cmp_isc_P(i);
        P_O2_i    = (cmp_isc_FO2(i)/n_T_i)    * cmp_isc_P(i);
        P_H2O_i   = (cmp_isc_FH2O(i)/n_T_i)   * cmp_isc_P(i);
        P_HCHO_i  = (cmp_isc_FHCHO(i)/n_T_i)  * cmp_isc_P(i);
        [r1_i, r2_i] = rates_LHHW(cmp_isc_T(i), cmp_reactor_param(1), P_CH3OH_i, P_O2_i, P_H2O_i, P_HCHO_i);
        cmp_isc_r1(i) = r1_i * dWdz_cmp;
        cmp_isc_r2(i) = r2_i * dWdz_cmp;
    end
    
    % Isothermal rates
    cmp_iso_r1 = zeros(length(cmp_iso_z), 1);
    cmp_iso_r2 = zeros(length(cmp_iso_z), 1);
    for i = 1:length(cmp_iso_z)
        n_T_i     = cmp_iso_FCH3OH(i)+cmp_iso_FO2(i)+cmp_iso_FHCHO(i)+cmp_iso_FH2O(i)+cmp_iso_FCO(i)+cmp_iso_FN2(i);
        P_CH3OH_i = (cmp_iso_FCH3OH(i)/n_T_i) * cmp_iso_P(i);
        P_O2_i    = (cmp_iso_FO2(i)/n_T_i)    * cmp_iso_P(i);
        P_H2O_i   = (cmp_iso_FH2O(i)/n_T_i)   * cmp_iso_P(i);
        P_HCHO_i  = (cmp_iso_FHCHO(i)/n_T_i)  * cmp_iso_P(i);
        [r1_i, r2_i] = rates_LHHW(cmp_iso_T(i), cmp_reactor_param(1), P_CH3OH_i, P_O2_i, P_H2O_i, P_HCHO_i);
        cmp_iso_r1(i) = r1_i * dWdz_cmp;
        cmp_iso_r2(i) = r2_i * dWdz_cmp;
    end
    
    %% --- Compute yields --- %
    cmp_cs_Y  = cmp_cs_FHCHO  ./ cmp_cs_F0CH3OH;
    cmp_isc_Y = cmp_isc_FHCHO ./ cmp_isc_F0CH3OH;
    cmp_iso_Y = cmp_iso_FHCHO ./ cmp_iso_F0CH3OH;
    
    %% --- Colours --- %
    col_iso = [0.00, 0.45, 0.70];   % blue
    col_isc = [0.85, 0.33, 0.10];   % orange
    col_cs  = [0.47, 0.67, 0.19];   % green
    
    %% --- Figure 1: Temperature Profile --- %
    figure('Name', 'Comparison: Temperature Profile', 'Color', 'w');
    hold on
    
    plot(cmp_cs_z,  cmp_cs_T,  '-', 'Color', col_cs,  'LineWidth', 2, ...
         'DisplayName', sprintf('Cold Shot (T_0=%dK, split=%.1f, %d stages)', ...
         cmp_cs_T0, cmp_cs_split, cmp_cs_stages));
    
    plot(cmp_isc_z, cmp_isc_T, '-', 'Color', col_isc, 'LineWidth', 2, ...
         'DisplayName', sprintf('Interstage (T_0=%dK, T_{sub}=%dK, %d stages)', ...
         cmp_isc_T0, cmp_isc_Tsub, cmp_isc_stages));
    
    plot(cmp_iso_z, cmp_iso_T, '-', 'Color', col_iso, 'LineWidth', 2, ...
         'DisplayName', sprintf('Isothermal (T_0=%dK)', cmp_iso_T0));
    
    yline(673, 'k--', 'LineWidth', 1.2, 'DisplayName', 'T_{max} = 673 K');
    
    xlabel('Reactor Length (m)'); ylabel('Temperature (K)');
    title('Temperature Profile Comparison');
    legend('Location', 'best'); grid on; hold off
    
    %% --- Figure 2: HCHO Yield --- %
    figure('Name', 'Comparison: HCHO Yield', 'Color', 'w');
    hold on
    
    plot(cmp_cs_z,  cmp_cs_Y,  '-', 'Color', col_cs,  'LineWidth', 2, ...
         'DisplayName', sprintf('Cold Shot (T_0=%dK, split=%.1f, %d stages)', ...
         cmp_cs_T0, cmp_cs_split, cmp_cs_stages));
    
    plot(cmp_isc_z, cmp_isc_Y, '-', 'Color', col_isc, 'LineWidth', 2, ...
         'DisplayName', sprintf('Interstage (T_0=%dK, T_{sub}=%dK, %d stages)', ...
         cmp_isc_T0, cmp_isc_Tsub, cmp_isc_stages));
    
    plot(cmp_iso_z, cmp_iso_Y, '-', 'Color', col_iso, 'LineWidth', 2, ...
         'DisplayName', sprintf('Isothermal (T_0=%dK)', cmp_iso_T0));
    
    % Final yield value annotations
    text(cmp_cs_z(end),  cmp_cs_Y(end),  sprintf('  %.3f', cmp_cs_Y(end)),  ...
         'Color', col_cs,  'FontSize', 9, 'FontWeight', 'bold');
    text(cmp_isc_z(end), cmp_isc_Y(end), sprintf('  %.3f', cmp_isc_Y(end)), ...
         'Color', col_isc, 'FontSize', 9, 'FontWeight', 'bold');
    text(cmp_iso_z(end), cmp_iso_Y(end), sprintf('  %.3f', cmp_iso_Y(end)), ...
         'Color', col_iso, 'FontSize', 9, 'FontWeight', 'bold');
    
    xlabel('Reactor Length (m)'); ylabel('HCHO Yield  F_{HCHO} / F_{0,CH_3OH}');
    title('HCHO Yield Comparison');
    legend('Location', 'best'); grid on; hold off
    
    %% --- Figure 3: Reaction Rate r1 --- %
    figure('Name', 'Comparison: Reaction Rate r1', 'Color', 'w');
    hold on
    
    plot(cmp_cs_z,  cmp_cs_r1,  '-', 'Color', col_cs,  'LineWidth', 2, ...
         'DisplayName', sprintf('Cold Shot (T_0=%dK, split=%.1f, %d stages)', ...
         cmp_cs_T0, cmp_cs_split, cmp_cs_stages));
    
    plot(cmp_isc_z, cmp_isc_r1, '-', 'Color', col_isc, 'LineWidth', 2, ...
         'DisplayName', sprintf('Interstage (T_0=%dK, T_{sub}=%dK, %d stages)', ...
         cmp_isc_T0, cmp_isc_Tsub, cmp_isc_stages));
    
    
    plot(cmp_iso_z, cmp_iso_r1, '-', 'Color', col_iso, 'LineWidth', 2, ...
         'DisplayName', sprintf('Isothermal (T_0=%dK)', cmp_iso_T0));
    
    xlabel('Reactor Length (m)');
    ylabel('r_1 (mol m^{-1} s^{-1})');
    title('Reaction Rate r_1 (CH_3OH \rightarrow HCHO) Comparison');
    legend('Location', 'best'); grid on; hold off
    
    % --- Print final summary to command window --- %
    fprintf('\n=== Comparison Summary ===\n')
    fprintf('  Cold Shot  : Final Y_HCHO = %.4f | T_out = %.1f K | z_total = %.3f m\n', ...
        cmp_cs_Y(end),  cmp_cs_T(end),  cmp_cs_z(end))
    fprintf('  Interstage : Final Y_HCHO = %.4f | T_out = %.1f K | z_total = %.3f m\n', ...
        cmp_isc_Y(end), cmp_isc_T(end), cmp_isc_z(end))
    fprintf('  Isothermal : Final Y_HCHO = %.4f | T_out = %.1f K | z_total = %.3f m\n', ...
        cmp_iso_Y(end), cmp_iso_T(end), cmp_iso_z(end))
end

%% ---------------------------------------------------------
%% Temperature Profile Comparison — 6 Cases
function plot_temperature_profile_6cases(reactor_param, inlet)

    reactor_length  = 20;
    P0              = 1.6/1.01325;
    Pressure_drop   = "Yes";
    Temp_dependence = "Yes";
    T_max           = 800;
    T_max_active    = "No";
    P_min           = 1.1/1.01325;
    P_min_active    = "Yes";
    option_limits   = {T_max, T_max_active, P_min, P_min_active};

    % ---- Case 1: 430 K, LHHW, 0.001 kg/s ----
    m1 = 0.001; T01 = 430;
    [fs1, ~, ~, ~] = get_initial_conditions(m1, P0, T01, 0, inlet);
    [z1,~,~,~,~,~,~,~,~,~,T1,~] = solve_reactor(fs1, reactor_length, P0, T01, option_limits, ["LHHW", Pressure_drop, Temp_dependence], reactor_param);

    % ---- Case 2: 430 K, Simple Power, 0.001 kg/s ----
    m2 = 0.001; T02 = 430;
    [fs2, ~, ~, ~] = get_initial_conditions(m2, P0, T02, 0, inlet);
    [z2,~,~,~,~,~,~,~,~,~,T2,~] = solve_reactor(fs2, reactor_length, P0, T02, option_limits, ["Simple Power", Pressure_drop, Temp_dependence], reactor_param);

    % ---- Case 3: 430 K, LHHW, 0.0006 kg/s ----
    m3 = 0.0006; T03 = 430;
    [fs3, ~, ~, ~] = get_initial_conditions(m3, P0, T03, 0, inlet);
    [z3,~,~,~,~,~,~,~,~,~,T3,~] = solve_reactor(fs3, reactor_length, P0, T03, option_limits, ["LHHW", Pressure_drop, Temp_dependence], reactor_param);

    % ---- Case 4: 430 K, Simple Power, 0.0006 kg/s ----
    m4 = 0.0006; T04 = 430;
    [fs4, ~, ~, ~] = get_initial_conditions(m4, P0, T04, 0, inlet);
    [z4,~,~,~,~,~,~,~,~,~,T4,~] = solve_reactor(fs4, reactor_length, P0, T04, option_limits, ["Simple Power", Pressure_drop, Temp_dependence], reactor_param);

    % ---- Case 5: 460 K, LHHW, 0.001 kg/s ----
    m5 = 0.001; T05 = 460;
    [fs5, ~, ~, ~] = get_initial_conditions(m5, P0, T05, 0, inlet);
    [z5,~,~,~,~,~,~,~,~,~,T5,~] = solve_reactor(fs5, reactor_length, P0, T05, option_limits, ["LHHW", Pressure_drop, Temp_dependence], reactor_param);

    % ---- Case 6: 460 K, Simple Power, 0.001 kg/s ----
    m6 = 0.001; T06 = 460;
    [fs6, ~, ~, ~] = get_initial_conditions(m6, P0, T06, 0, inlet);
    [z6,~,~,~,~,~,~,~,~,P6,T6,~] = solve_reactor(fs6, reactor_length, P0, T06, option_limits, ["Simple Power", Pressure_drop, Temp_dependence], reactor_param);
    fprintf('Case 6 final P = %.4f atm, P_min = %.4f atm\n', P6(end), P_min)
    % ---- Plot ----
    figure('Name', 'Temperature Profiles: 6 Cases')
    plot(z1, T1, 'LineWidth', 1.5,'Color','g','LineStyle','-'); hold on
    plot(z2, T2, 'LineWidth', 1.5,'Color','g','LineStyle','--')
    plot(z3, T3, 'LineWidth', 1.5,'Color','r')
    plot(z4, T4, 'LineWidth', 1.5,'Color','r','LineStyle','--')
    plot(z5, T5, 'LineWidth', 1.5,'Color','b')
    plot(z6, T6, 'LineWidth', 1.5','Color','b',LineStyle=':')
    xlabel('Reactor Length (m)')
    ylabel('Temperature (K)')
    title('Temperature Profiles for Different Kinetics, Inlet Temperatures, and Mass Flowrates')
    legend( ...
        '430 K, LHHW, 0.001 kg/s', ...
        '430 K, Simple Power, 0.001 kg/s', ...
        '430 K, LHHW, 0.0006 kg/s', ...
        '430 K, Simple Power, 0.0006 kg/s', ...
        '460 K, LHHW, 0.001 kg/s', ...
        '460 K, Simple Power, 0.001 kg/s', ...
        'Location', 'best')
    grid on
    set(gcf, 'color', 'w')
    xlim([0 7])

end


%% ---------------------------------------------------------
%% Reaction Rate Profiles — LHHW, 430 K
function plot_rate_profiles_lhhw_430(reactor_param, inlet)

    reactor_length = 5;
    P0             = 1.6/1.01325;
    T0             = 430;
    T_max          = 800;
    T_max_active   = "No";
    P_min          = 1.1/1.01325;
    P_min_active   = "No";
    option_limits  = {T_max, T_max_active, P_min, P_min_active};
    options        = ["LHHW", "Yes", "Yes"];
    alpha          = reactor_param(1);

    % ---- Case 1: 0.001 kg/s ----
    m1 = 0.001;
    [fs1, ~, ~, ~] = get_initial_conditions(m1, P0, T0, 0, inlet);
    [z1, F_CH3OH_1, F_O2_1, F_HCHO_1, F_H2O_1, F_CO_1, F_N2_1, ~,~, P1, T1, ~] = ...
        solve_reactor(fs1, reactor_length, P0, T0, option_limits, options, reactor_param);

    % ---- Case 2: 0.0006 kg/s ----
    m2 = 0.0006;
    [fs2, ~, ~, ~] = get_initial_conditions(m2, P0, T0, 0, inlet);
    [z2, F_CH3OH_2, F_O2_2, F_HCHO_2, F_H2O_2, F_CO_2, F_N2_2, ~,~, P2, T2, ~] = ...
        solve_reactor(fs2, reactor_length, P0, T0, option_limits, options, reactor_param);

    % ---- Compute rates for case 1 ----
    r1_case1 = zeros(length(z1), 1);
    r2_case1 = zeros(length(z1), 1);
    for i = 1:length(z1)
        n_T = F_CH3OH_1(i)+F_O2_1(i)+F_HCHO_1(i)+F_H2O_1(i)+F_CO_1(i)+F_N2_1(i);
        [r1_case1(i), r2_case1(i)] = rates_LHHW(T1(i), alpha, ...
            (F_CH3OH_1(i)/n_T)*P1(i), (F_O2_1(i)/n_T)*P1(i), ...
            (F_H2O_1(i)/n_T)*P1(i),   (F_HCHO_1(i)/n_T)*P1(i));
    end

    % ---- Compute rates for case 2 ----
    r1_case2 = zeros(length(z2), 1);
    r2_case2 = zeros(length(z2), 1);
    for i = 1:length(z2)
        n_T = F_CH3OH_2(i)+F_O2_2(i)+F_HCHO_2(i)+F_H2O_2(i)+F_CO_2(i)+F_N2_2(i);
        [r1_case2(i), r2_case2(i)] = rates_LHHW(T2(i), alpha, ...
            (F_CH3OH_2(i)/n_T)*P2(i), (F_O2_2(i)/n_T)*P2(i), ...
            (F_H2O_2(i)/n_T)*P2(i),   (F_HCHO_2(i)/n_T)*P2(i));
    end

    % ---- Plot ----
    figure('Name', 'Reaction Rate Profiles: LHHW, 430 K')
    plot(z1, r1_case1, 'LineWidth', 1.5); hold on
    plot(z1, r2_case1, 'LineWidth', 1.5)
    plot(z2, r1_case2, 'LineWidth', 1.5)
    plot(z2, r2_case2, 'LineWidth', 1.5)
    xlabel('Reactor Length (m)')
    ylabel('Reaction Rate (mol kg^{-1} s^{-1})')
    title('Reaction Rate Profiles for LHHW Kinetics at 430 K')
    legend('r_1, 0.001 kg/s', 'r_2, 0.001 kg/s', ...
           'r_1, 0.0006 kg/s', 'r_2, 0.0006 kg/s', 'Location', 'best')
    grid on
    set(gcf, 'color', 'w')
    xlim([1.5 4])
    ylim([0 0.055])

end

%% Extent of Reaction vs Catalyst Weight
% ---------------------------------------------------------
function plot_extent_vs_catalyst_weight(reactor_param, inlet)

    reactor_length = 5;
    P0             = 1.6/1.01325;
    T0             = 430;
    T_max          = 800;
    T_max_active   = "No";
    P_min          = 1.1/1.01325;
    P_min_active   = "No";
    option_limits  = {T_max, T_max_active, P_min, P_min_active};
    options        = ["LHHW", "Yes", "Yes"];

    rho_cat     = reactor_param(2);
    bed_voidage = reactor_param(3);
    tube_diam   = reactor_param(5);
    dW_dz       = rho_cat * (1 - bed_voidage) * pi*(tube_diam/2)^2;

    % ---- Case 1: 0.001 kg/s ----
    m1 = 0.001;
    [fs1, ~, ~, ~] = get_initial_conditions(m1, P0, T0, 0, inlet);
    [z1,~,~,~,~,~,~, extent1_1, extent2_1, ~,~,~] = ...
        solve_reactor(fs1, reactor_length, P0, T0, option_limits, options, reactor_param);
    W1 = z1 * dW_dz;

    % ---- Case 2: 0.0006 kg/s ----
    m2 = 0.0006;
    [fs2, ~, ~, ~] = get_initial_conditions(m2, P0, T0, 0, inlet);
    [z2,~,~,~,~,~,~, extent1_2, extent2_2, ~,~,~] = ...
        solve_reactor(fs2, reactor_length, P0, T0, option_limits, options, reactor_param);
    W2 = z2 * dW_dz;

    % ---- Plot ----
    figure('Name', 'Extent of Reaction vs Catalyst Weight')
    plot(W1, extent1_1, 'LineWidth', 1.5); hold on
    plot(W1, extent2_1, 'LineWidth', 1.5)
    plot(W2, extent1_2, 'LineWidth', 1.5)
    plot(W2, extent2_2, 'LineWidth', 1.5)
    xlabel('Catalyst Weight (kg)')
    ylabel('Extent of Reaction (mol/s)')
    title('Extent of Reaction vs Catalyst Weight at 430 K (LHHW)')
    legend('\xi_1, 0.001 kg/s', '\xi_2, 0.001 kg/s', ...
           '\xi_1, 0.0006 kg/s', '\xi_2, 0.0006 kg/s', 'Location', 'best')
    grid on
    set(gcf, 'color', 'w')
    xlim([0.4 1.1])

end
function plot_rate_vs_temperature_comparison(reactor_param, inlet, P0)
% Compares LHHW vs Simple Power Law rates (r1 and r2) over a temperature range
% at fixed inlet composition and pressure.
%
% INPUTS:
%   reactor_param - [alpha, rho_cat, bed_voidage, D_p, tube_diameter, Hr1, Hr2]
%   inlet         - [y0_CH3OH, y0_O2, y0_H2O, y0_HCHO, y0_CO, y0_N2]
%   P0            - total pressure [atm]

    alpha = reactor_param(1);

    % --- Fixed inlet mole fractions ---
    y_CH3OH = inlet(1);
    y_O2    = inlet(2);
    y_H2O   = inlet(3);
    y_HCHO  = inlet(4);  % 0 at inlet

    % --- Temperature sweep ---
    T_range = linspace(430, 1000, 300);  % [K]

    r1_LHHW  = zeros(size(T_range));
    r2_LHHW  = zeros(size(T_range));
    r1_PL    = zeros(size(T_range));
    r2_PL    = zeros(size(T_range));

    y_HCHO_r2 = 0.05;
    P_HCHO_r2 = y_HCHO_r2 * P0;
    
    r2_LHHW_combined = zeros(size(T_range));
    r2_PL_combined   = zeros(size(T_range));
    for i = 1:length(T_range)
        [r1_LHHW(i), r2_LHHW_combined(i)] = rates_LHHW(T_range(i), alpha, ...
            y_CH3OH*P0, y_O2*P0, y_H2O*P0, P_HCHO_r2);
        [r1_PL(i), r2_PL_combined(i)] = rates_powerlaw(T_range(i), y_CH3OH, ...
            y_O2, y_HCHO_r2, y_H2O);
    end
    
    figure('Name', 'Rate vs Temperature: All Reactions Superimposed', 'Color', 'w')
    plot(T_range, r1_LHHW,          'b-',  'LineWidth', 2, 'DisplayName', 'r_1 LHHW')
    hold on
    plot(T_range, r2_LHHW_combined, 'r-',  'LineWidth', 2, 'DisplayName', sprintf('r_2 LHHW [y_{HCHO}=%.2f]', y_HCHO_r2))
    hold off
    xlabel('Temperature (K)')
    ylabel('Reaction Rate (mol kg_{cat}^{-1} s^{-1})')
    title(sprintf('r_1 and r_2 vs Temperature\n(P=%.1f atm)', ...
        P0))
    legend('Location', 'best')
    grid on
end