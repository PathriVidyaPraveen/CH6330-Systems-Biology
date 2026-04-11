% --- Problem 2: Patient Data Scatter Plot ---
clear; clc; close all;

% 1. Fetch the data
load('Data.mat');

% First column is time, second is virus levels
t_data = Data(:, 1);
V_data = Data(:, 2);

% 2. Create the figure to match the desired result exactly
figure('Name', 'Problem 2: Patient Data', 'Position', [100, 100, 550, 450], 'Color', [0.94 0.94 0.94]);

% Plot the scatter plot using red open circles ('ro') and log10 scale for V
plot(t_data, log10(V_data), 'ro', 'MarkerSize', 6, 'LineWidth', 1);

% 3. Format the axes and labels
xlim([4 22]);
ylim([3 9]);

% Set exact ticks to match the provided image
set(gca, 'XTick', 4:2:22);
set(gca, 'YTick', 3:1:9);

% Set labels
xlabel('Days postonset of symptoms', 'FontSize', 14);
ylabel('log_{10} (V)', 'FontSize', 14);

% Match the box and background styling
set(gca, 'FontSize', 12, 'LineWidth', 1, 'Color', 'w');
box on;
