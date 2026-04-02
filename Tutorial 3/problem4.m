
% Parametric Sensitivity Analysis (1D)
clear; clc; close all;

% 1. Define Constant Parameters
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
y0 = [0.1; 1.5; 0.1]; 

% 3. Time Span for simulation
% FIX: Force ode45 to evaluate exactly every 1 second
tspan = 0:1:1000; 

% 4. Initialize v_M2 range and Frequency array
v_M2_values = 1:20;
Freq = zeros(1, length(v_M2_values));

% 5. For loop to iterate through values of v_M2
for i = 1:length(v_M2_values)
    % Update v_M2 for the current iteration
    p.v_M2 = v_M2_values(i);
    
    % Simulate the ODE system
    [t, y] = ode45(@(t,y) lh_model(t, y, p), tspan, y0);
    
    % Extract cytoplasmic calcium concentration (1st column)
    Ca_cyt = y(:,1);
    
    % Find peaks (Since t is spaced by 1s, MinPeakDistance of 50 means 50 seconds)
    [pks, ~] = findpeaks(Ca_cyt, 'MinPeakDistance', 50, 'MinPeakHeight', 0.4);
    
    % Store the number of peaks found
    Freq(i) = length(pks);
end

% 6. Plot the Bar Chart
figure('Name', 'Problem 4: 1D Parametric Sensitivity', 'Color', 'w');
bar(v_M2_values, Freq, 0.8, 'FaceColor', '#0072BD'); 
xlabel('v_{M2}');
ylabel('Number of peaks');
xlim([0 21]); 
set(gca, 'FontSize', 10);
box on;

% =========================================================================
% ODE Function Definition
% =========================================================================
function dydt = lh_model(~, y, p)
    X = y(1); Y = y(2); Z = y(3); 
    
    v_serca = p.v_M2 * (X^2 / (X^2 + p.k_2^2));
    v_PLC   = p.v_p * (X^2 / (X^2 + p.k_p^2));
    
    term1 = (p.k_CaA^p.n * X^p.n) / ((X^p.n + p.k_CaA^p.n) * (X^p.n + p.k_CaI^p.n));
    term2 = (Z^p.m) / (Z^p.m + p.k_ip3^p.m);
    v_CICR = 4 * p.v_M3 * term1 * term2 * (Y - X);
    
    dXdt = p.v_in - p.k_out*X + v_CICR - v_serca + p.k_f*(Y - X);
    dYdt = v_serca - v_CICR - p.k_f*(Y - X);
    dZdt = v_PLC - p.k_deg*Z;
    
    dydt = [dXdt; dYdt; dZdt];
end
