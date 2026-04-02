
clear; clc; close all;

% 1. Define Constant Parameters
p.v_M2  = 5.7;      % uM/s (Modified from 15.0 for Problem 3)
p.v_M3  = 40.0;     % s^-1
p.v_in  = 0.05;     % uM/s
p.v_p   = 0.05;     % uM/s
p.k_2   = 0.1;      % uM
p.k_CaA = 0.15;     % uM
p.k_CaI = 0.15;     % uM
p.k_ip3 = 0.1;      % uM
p.k_p   = 0.3;      % uM
p.k_deg = 0.08;     % s^-1
p.k_out = 0.5;      % s^-1
p.k_f   = 0.5;      % s^-1
p.n     = 2.02;
p.m     = 2.2;

% 2. Initial Conditions [X; Y; Z]
% X = [Ca2+]cyt, Y = [Ca2+]ER, Z = [IP3]
y0 = [0.1; 1.5; 0.1]; 

% 3. Time Span for simulation (1000 seconds)
tspan = [0 1000];

% 4. Simulate the ODE system
[t, y] = ode45(@(t,y) lh_model(t, y, p), tspan, y0);

% Extract calcium concentrations
Ca_cyt = y(:,1);
Ca_ER  = y(:,2);

% 5. Phase Plane Plot
figure('Name', 'Problem 3: Phase Plane Analysis', 'Color', 'w');
plot(Ca_cyt, Ca_ER, 'b', 'LineWidth', 1.2);
xlabel('[Ca^{2+}]_{cyt}');
ylabel('[Ca^{2+}]_{ER}');
title('Phase Plane Analysis (v_{M2} = 5.7)');
grid on;
set(gca, 'FontSize', 10);

% =========================================================================
% ODE Function Definition
% =========================================================================
function dydt = lh_model(~, y, p)
    % State variables
    X = y(1); % Cytoplasmic Calcium
    Y = y(2); % ER Calcium
    Z = y(3); % IP3
    
    % Reaction Rates
    v_serca = p.v_M2 * (X^2 / (X^2 + p.k_2^2));
    v_PLC   = p.v_p * (X^2 / (X^2 + p.k_p^2));
    
    % CICR Rate 
    term1 = (p.k_CaA^p.n * X^p.n) / ((X^p.n + p.k_CaA^p.n) * (X^p.n + p.k_CaI^p.n));
    term2 = (Z^p.m) / (Z^p.m + p.k_ip3^p.m);
    v_CICR = 4 * p.v_M3 * term1 * term2 * (Y - X);
    
    % Differential Equations
    dXdt = p.v_in - p.k_out*X + v_CICR - v_serca + p.k_f*(Y - X);
    dYdt = v_serca - v_CICR - p.k_f*(Y - X);
    dZdt = v_PLC - p.k_deg*Z;
    
    % Return column vector
    dydt = [dXdt; dYdt; dZdt];
end
