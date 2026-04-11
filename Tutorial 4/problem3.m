% --- Problem 3: Data Fitting and Simulation ---
clear
close all
clc

load('Data.mat')
tspan = [-3 22];
P0 = 3;             %% This is the initial guess (Do not change).

% Complete this line of code using the given hint. 
% We pass empty matrices for linear constraints and set a lower bound (lb) of 0 
% for parameter P to ensure the clearance rate doesn't go negative.
options = optimoptions('fmincon','Display','off');
[P, fval] = fmincon(@objcovid, P0, [], [], [], [], 0, [], [], options);

display(P);
display(fval);

y0 = [0.31; 10^6];

% Use a finer time grid for a smoother fitted curve
tspan_plot = -3:0.1:22;
[t, y] = ode45(@(t,y) covidODE(t,y,P), tspan_plot, y0);

%% Plot the simulation result and experimental data together.
figure('Name', 'Problem 3: Fitted Simulation', 'Position', [100, 100, 600, 450], 'Color', [0.94 0.94 0.94]);

% 1. Plot the fitted simulation
plot(t, log10(y(:,1)), 'b-', 'LineWidth', 1.5);
hold on;

% 2. Plot the patient data
plot(Data(:,1), log10(Data(:,2)), 'ro', 'MarkerSize', 6, 'LineWidth', 1);

% 3. Formatting to strictly replicate the image
xlim([-5 25]);
ylim([-1 9]);
set(gca, 'XTick', -5:5:25, 'YTick', -1:1:9);

xlabel('Days postonset of symptoms', 'FontSize', 14);
ylabel('log_{10} (V)', 'FontSize', 14);

% 4. Add Legend
legend('Fitted simulation', 'Patient data', 'Location', 'northeast', 'FontSize', 12, 'EdgeColor', 'k');

% 5. Aesthetics
set(gca, 'FontSize', 12, 'LineWidth', 1, 'Color', 'w');
box on;
hold off;
