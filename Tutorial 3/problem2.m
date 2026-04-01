% Systems Biology: Calcium Oscillation (Problem 2)
clear; clc; close all;

% 1. Define Constant Parameters
p.v_M2  = 15.0;     % uM/s
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
p.m     = 2.2;

% 2. Initial Conditions [X; Y; Z]
% X = [Ca2+]cyt, Y = [Ca2+]ER, Z = [IP3]
y0 = [0.1; 1.5; 0.1]; 

% 3. Time Span for each simulation (200 seconds)
tspan = [0 200];

% 4. Parameter 'n' Range (2 to 4 with increment of 0.01)
n_values = 2:0.01:4;

% 5. Initialize Figure for Bifurcation Plot
figure('Name', 'Problem 2: Bifurcation Plot', 'Color', 'w');
xlabel('n');
ylabel('[Ca^{2+}]_{cyt}');
xlim([2 4]);
set(gca, 'FontSize', 10);

% 6. For loop to iterate through values of 'n'
for i = 1:length(n_values)
    % Update 'n' for the current iteration
    p.n = n_values(i);
    
    % Simulate the ODE system
    [~, y] = ode45(@(t,y) lh_model(t, y, p), tspan, y0);
    
    % Extract cytoplasmic calcium concentration (1st column)
    Ca_cyt = y(:,1);
    
    % Calculate the maximum and minimum values
    max_Ca = max(Ca_cyt);
    min_Ca = min(Ca_cyt);
    
    % Plotting in the loop 
    hold on; % Retain the plot to add new points
    plot(p.n, max_Ca, 'b.', 'MarkerSize', 5); % Max values in blue
    plot(p.n, min_Ca, 'r.', 'MarkerSize', 5); % Min values in red
    drawnow; % Force MATLAB to draw the plot progressively
end

hold off;

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
