
clear; clc; close all;

% 1. Define time span and initial conditions
tspan = [-3 22];    % Simulating from 3 days before symptom onset
V0 = 0.31;          % Initial viral concentration V(0)
T0_init = 10^6;     % Initial T cells T(0)
y0 = [V0; T0_init]; 

% 2. Solve the ODE system
[t, y] = ode45(@covidSystem, tspan, y0);

% Extract V and T
V = y(:, 1);
T = y(:, 2);

% 3. Plotting to replicate exactly
figure('Position', [100, 100, 600, 500], 'Color', [0.94 0.94 0.94]);

% Subplot 1: log10(V)
subplot(2,1,1);
plot(t, log10(V), 'b-', 'LineWidth', 1.5);
xlim([-3 22]);
ylim([-1 10]); % Matching the y-axis ticks of the desired plot
set(gca, 'YTick', 0:2:10, 'XTick', 0:5:20);
ylabel('log_{10} (V)', 'FontSize', 14);
xlabel('Days postonset of symptoms', 'FontSize', 14);
set(gca, 'FontSize', 12, 'LineWidth', 1, 'Color', 'w');
box on;

% Subplot 2: log10(T)
subplot(2,1,2);
plot(t, log10(T), 'r-', 'LineWidth', 1.5);
xlim([-3 22]);
ylim([6 6.6]); % Matching the y-axis ticks of the desired plot
set(gca, 'YTick', 6:0.2:6.6, 'XTick', 0:5:20);
ylabel('log_{10} (T)', 'FontSize', 14);
xlabel('Days postonset of symptoms', 'FontSize', 14);
set(gca, 'FontSize', 12, 'LineWidth', 1, 'Color', 'w');
box on;

% =========================================================================
% Local ODE Function
% =========================================================================
function dydt = covidSystem(~, y)
    V = y(1);
    T = y(2);
    
    % Parameters
    K = 7.11e9;
    p = 6.31;
    cT = 1.58e-6;
    r = 0.794;
    kT = 7.94e7;
    m = 2;
    c = 3;
    
   
    
    deltaT = 0.1;  
    
    T0 = 10^6;
    sT = deltaT * T0; 
    
    % Equations
    dVdt = p * V * (1 - V/K) - cT * V * T - c * V;
    dTdt = sT + r * T * (V^m / (V^m + kT^m)) - deltaT * T;
    
    dydt = [dVdt; dTdt];
end
