%% =========================================================
%% Li-Rinzel Model – Parameter Search via Monte Carlo
%% =========================================================
%  Two-variable ODE model for calcium oscillations:
%    dCa/dt = c1*v1*pinf^3*ninf^3*h^3*(CaER-Ca)
%             + c1*v2*(CaER-Ca)
%             - v3*Ca^2/(Ca^2+k3^2)
%    dh/dt  = alphah*(1-h) - betah*h
%
%  Auxiliary quantities:
%    CaER  = (c0 - Ca)/c1
%    pinf  = IP3/(IP3+d1)
%    ninf  = Ca/(Ca+d5)
%    Q2    = d2*(IP3+d1)/(IP3+d3)
%    tauh  = 1/(a2*(Q2+Ca))
%    hinf  = Q2/(Q2+Ca)
%    alphah= hinf/tauh
%    betah = (1-hinf)/tauh
% =========================================================

clear; clc; close all;

%% =========================================================
%% QUESTION 1: Euler Scheme Simulation (nominal parameters)
%% =========================================================

% --- Nominal parameters ---
IP3 = 0.5;  v1  = 6;      v2 = 0.11;
v3  = 0.9;  c0  = 2;      c1 = 0.185;
a2  = 0.2;  d1  = 0.13;   d2 = 1.049;
d3  = 0.9434; d5 = 0.08234; k3 = 0.1;

% --- Time vector & step ---
time = 0 : 0.1 : 100;
dt   = 0.1;

% --- Pre-allocate state variables ---
Ca = zeros(1, numel(time));
h  = zeros(1, numel(time));

% --- Initial conditions ---
Ca(1) = 0.1;
h(1)  = 0.1;

% --- Euler integration loop ---
for t = 2 : numel(time)

    CaER   = (c0 - Ca(t-1)) / c1;            % ER calcium
    pinf   = IP3 / (IP3 + d1);               % IP3-receptor activation
    ninf   = Ca(t-1) / (Ca(t-1) + d5);       % Ca-activation
    Q2     = d2 * (IP3 + d1) / (IP3 + d3);   % lumped parameter
    tauh   = 1 / (a2 * (Q2 + Ca(t-1)));      % time constant
    hinf   = Q2 / (Q2 + Ca(t-1));            % steady-state h
    alphah = hinf / tauh;                     % opening rate
    betah  = (1 - hinf) / tauh;              % closing rate

    dh  = alphah * (1 - h(t-1))  -  betah * h(t-1);

    dCa = c1*v1 * (pinf^3) * (ninf^3) * (h(t-1)^3) * (CaER - Ca(t-1)) ...
        + c1*v2 * (CaER - Ca(t-1)) ...
        - v3 * Ca(t-1)^2 / (Ca(t-1)^2 + k3^2);

    Ca(t) = Ca(t-1) + dCa * dt;
    h(t)  = h(t-1)  + dh  * dt;
end

% --- Plot nominal simulation ---
figure(1);
subplot(2,1,1);
plot(time, Ca, 'b', 'LineWidth', 2);
xlabel('Time (s)',        'FontSize', 14, 'FontWeight', 'bold');
ylabel('[Ca^{2+}] (\muM)', 'FontSize', 14, 'FontWeight', 'bold');
title('Li-Rinzel Model: Calcium Oscillations (Euler Scheme)', ...
      'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

subplot(2,1,2);
plot(time, h, 'r', 'LineWidth', 2);
xlabel('Time (s)',             'FontSize', 14, 'FontWeight', 'bold');
ylabel('h (gating variable)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

%% =========================================================
%% MONTE CARLO PARAMETER SEARCH
%% =========================================================

% --- Time vector for MC simulations ---
time1 = 0 : 0.1 : 240;   % extend to 240 s

% --- Generate "experimental" reference using nominal parameters ---
%     (run model and normalise to [0,1])
Ca_nom = lirinzel_ns(0.11, 0.2, 0.5, 0.08234, 0.1, time1);
normalised_exp = (Ca_nom - min(Ca_nom)) / (max(Ca_nom) - min(Ca_nom));

% --- Parameter bounds [lower, upper] ---
%     rows: v2 | a2 | IP3 | d5 | k3
param_bounds = [
    0.055, 0.07;   % v2
    0.18,  0.23;   % a2
    0.55,  0.63;   % IP3
    0.06,  0.12;   % d5
    0.06,  0.10;   % k3
];

% --- MC setup ---
num_simulations = 1000;
simulations     = zeros(num_simulations, 6);  % [v2 a2 IP3 d5 k3 RMSE]
rmse_values     = zeros(1, num_simulations);

fprintf('Running %d Monte Carlo simulations...\n', num_simulations);

for i = 1 : num_simulations

    % Sample 5 parameters from uniform distributions
    params = param_bounds(:,1) + rand(5,1) .* ...
             (param_bounds(:,2) - param_bounds(:,1));

    v2  = params(1);
    a2  = params(2);
    IP3 = params(3);
    d5  = params(4);
    k3  = params(5);

    % Simulate and normalise
    Ca_sim = lirinzel_ns(v2, a2, IP3, d5, k3, time1);

    % RMSE against "experimental" data
    rmse = sqrt(mean((normalised_exp - Ca_sim).^2));

    simulations(i, :) = [v2, a2, IP3, d5, k3, rmse];
    rmse_values(i)    = rmse;
end

fprintf('Done. Best RMSE = %.6f\n', min(rmse_values));

% --- Show 4 sample simulations (indices spread across 1000 runs) ---
sample_idx = randperm(num_simulations,4);
panel_titles = {'Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4'};

figure(2);
sgtitle('Monte Carlo – 4 Sample Simulations', 'FontSize', 14, 'FontWeight', 'bold');

for k = 1 : 4
    idx  = sample_idx(k);
    p    = simulations(idx, 1:5);
    Ca_sim = lirinzel_ns(p(1), p(2), p(3), p(4), p(5), time1);

    subplot(2, 2, k);
    plot(time1, normalised_exp, 'c--', 'LineWidth', 2); hold on;
    plot(time1, Ca_sim,         'LineWidth', 2);
    xlim([0 240]);
    xlabel('Time (s)',         'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Normalized Ca',   'FontSize', 12, 'FontWeight', 'bold');
    title(panel_titles{k},    'FontSize', 12, 'FontWeight', 'bold');
    legend('Experiment', 'Simulation', 'Location', 'best');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
end

% --- Find best match (minimum RMSE) ---
[best_rmse, idx] = min(rmse_values);
best_params = simulations(idx, :);

best_v2  = best_params(1);
best_a2  = best_params(2);
best_IP3 = best_params(3);
best_d5  = best_params(4);
best_k3  = best_params(5);

best_Ca_sim = lirinzel_ns(best_v2, best_a2, best_IP3, best_d5, best_k3, time1);

% --- Plot best match ---
figure(3000);
plot(time1, normalised_exp, 'k--', 'LineWidth', 2); hold on;
plot(time1, best_Ca_sim,    'r',   'LineWidth', 2);
xlabel('Time (s)',          'FontSize', 14, 'FontWeight', 'bold');
ylabel('Normalized Calcium','FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Best MC Match  (RMSE = %.5f)', best_rmse), ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Experimental', 'Best Simulation', 'Location', 'best');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 240]);

% --- Print best parameters to console ---
fprintf('\n--- Best Parameters (MC Search) ---\n');
fprintf('  v2   = %.6f\n', best_v2);
fprintf('  a2   = %.6f\n', best_a2);
fprintf('  IP3  = %.6f\n', best_IP3);
fprintf('  d5   = %.6f\n', best_d5);
fprintf('  k3   = %.6f\n', best_k3);
fprintf('  RMSE = %.6f\n', best_rmse);


%% =========================================================
%% LOCAL FUNCTION: lirinzel_ns
%% =========================================================
%  Simulates the Li-Rinzel model with Euler's scheme.
%  Returns NORMALISED Ca (scaled to [0,1]).
%
%  Inputs:
%    v2, a2, IP3, d5, k3 – sampled parameters
%    time                 – time vector (uniform spacing)
%
%  Output:
%    Ca  – normalised calcium trace (same length as time)
% =========================================================
function Ca = lirinzel_ns(v2, a2, IP3, d5, k3, time)

    % Fixed (non-sampled) parameters
    v1 = 6;  v3 = 0.9;  c0 = 2;  c1 = 0.185;
    d1 = 0.13; d2 = 1.049; d3 = 0.9434;

    dt = time(2) - time(1);
    N  = numel(time);

    Ca_raw = zeros(1, N);
    h      = zeros(1, N);
    Ca_raw(1) = 0.1;
    h(1)      = 0.1;

    for t = 2 : N
        CaER   = (c0 - Ca_raw(t-1)) / c1;
        pinf   = IP3 / (IP3 + d1);
        ninf   = Ca_raw(t-1) / (Ca_raw(t-1) + d5);
        Q2     = d2 * (IP3 + d1) / (IP3 + d3);
        tauh   = 1 / (a2 * (Q2 + Ca_raw(t-1)));
        hinf   = Q2 / (Q2 + Ca_raw(t-1));
        alphah = hinf / tauh;
        betah  = (1 - hinf) / tauh;

        dh  = alphah * (1 - h(t-1)) - betah * h(t-1);

        dCa = c1*v1*(pinf^3)*(ninf^3)*(h(t-1)^3)*(CaER - Ca_raw(t-1)) ...
            + c1*v2*(CaER - Ca_raw(t-1)) ...
            - v3*Ca_raw(t-1)^2 / (Ca_raw(t-1)^2 + k3^2);

        Ca_raw(t) = Ca_raw(t-1) + dCa * dt;
        h(t)      = h(t-1)      + dh  * dt;
    end

    % Normalise to [0, 1]
    lo = min(Ca_raw);
    hi = max(Ca_raw);
    if (hi - lo) < 1e-12
        Ca = Ca_raw;          % flat trace – return as-is
    else
        Ca = (Ca_raw - lo) / (hi - lo);
    end
end
