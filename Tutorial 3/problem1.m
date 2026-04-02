clear; clc; close all;

% 1. Define Parameter Values
p.v_M2 = 15.0;     % uM/s
p.v_M3 = 40.0;     % s^-1
p.v_in = 0.05;     % uM/s
p.v_p = 0.05;      % uM/s
p.k_2 = 0.1;       % uM
p.k_CaA = 0.15;    % uM
p.k_CaI = 0.15;    % uM
p.k_ip3 = 0.1;     % uM
p.k_p = 0.3;       % uM
p.k_deg = 0.08;    % s^-1
p.k_out = 0.5;     % s^-1
p.k_f = 0.5;       % s^-1
p.n = 2.02;
p.m = 2.2;

% 2. Initial Conditions [X; Y; Z]
% X = [Ca2+]cyt, Y = [Ca2+]ER, Z = [IP3]
y0 = [0.1; 1.5; 0.1]; 

% 3. Time Span for Simulation (2000 seconds)
tspan = [0 2000];

% 4. Solve the ODE system using ode45
[t, y] = ode45(@(t,y) lh_model(t, y, p), tspan, y0);

% Extract variables for plotting
Ca_cyt = y(:,1);
Ca_ER  = y(:,2);
IP3    = y(:,3);

% 5. Plotting the results
figure('Name', 'Problem 1: Time Series', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Plot [Ca2+]cyt
subplot(3,1,1);
plot(t, Ca_cyt, 'b', 'LineWidth', 1.5);
ylabel('[Ca^{2+}]_{cyt}');
xlabel('Time, s');
xlim([0 2000]);
set(gca, 'FontSize', 10);

% Plot [Ca2+]ER
subplot(3,1,2);
plot(t, Ca_ER, 'r', 'LineWidth', 1.5);
ylabel('[Ca^{2+}]_{ER}');
xlabel('Time, s');
xlim([0 2000]);
set(gca, 'FontSize', 10);

% Plot [IP3]
subplot(3,1,3);
plot(t, IP3, 'm', 'LineWidth', 1.5);
ylabel('[IP_3]');
xlabel('Time, s');
xlim([0 2000]);
set(gca, 'FontSize', 10);

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
