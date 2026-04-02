
% Parametric Sensitivity Analysis (2D)
clear; clc; close all;

% 1. Define Constant Parameters
p.v_in  = 0.05;     
p.v_p   = 0.05;     
p.k_2   = 0.1;      
p.k_CaA = 0.15;     
p.k_CaI = 0.15;     
p.k_ip3 = 0.1;      
p.k_p   = 0.3;      
p.k_deg = 0.08;     
p.k_out = 0.5;      
p.k_f   = 0.5;      
p.n     = 2.02;
p.m     = 2.2;

% 2. Arrays based on the explicit hints in the document
v_M2_values = 1:3:28;       
v_M3_values = 5:2.5:27.5;   

% Initialize Matrix (rows = v_M2, cols = v_M3)
M = zeros(length(v_M2_values), length(v_M3_values));

% 3. Fixed time span to prevent adaptive solver artifacts
tspan = 0:1:1000; 

% 4. Nested loops 
for i = 1:length(v_M2_values)
    for j = 1:length(v_M3_values)
        y0 = [0.1; 1.5; 0.1]; 
        
        % Update parameters
        p.v_M2 = v_M2_values(i);
        p.v_M3 = v_M3_values(j);
        
        % Simulate the ODE system
        [~, y] = ode45(@(t,y) lh_model(t, y, p), tspan, y0);
        
        % Store max cytoplasmic calcium
        M(i,j) = max(y(:,1));
    end
end

% =========================================================================
% Plotting Figure 1: 3D Bar Chart
% =========================================================================
figure('Name', 'Problem 5: 3D Bar Chart');

% Plotting M directly maps the high peaks exactly as shown in the document
b = bar3(M);

% Apply gradient coloring
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

% Force labels to match the document's specific axis swap
set(gca, 'XTick', 1:length(v_M2_values), 'XTickLabel', v_M2_values); 
set(gca, 'YTick', 1:length(v_M3_values), 'YTickLabel', v_M3_values);

xlabel('v_{M2}', 'FontSize', 12);
ylabel('v_{M3}', 'FontSize', 12);
zlabel('Maximum Ca^{2+}', 'FontSize', 12);

% Document orientation view
view([-37.5 30]); 
grid on;

% =========================================================================
% Plotting Figure 2: Surface Plot
% =========================================================================
figure('Name', 'Problem 5: Surface Plot');

% To get the exact mismatched alignment seen in the tutorial's surface plot,
% X and Y must be mapped against M without transposition.
[X, Y] = meshgrid(v_M2_values, v_M3_values);
s = surf(X, Y, M);

colormap('parula'); 

% Set exact axis ticks from the document image
set(gca, 'XTick', [10 20 30]); 
set(gca, 'YTick', [10 15 20 25 30]);
xlim([1 30]); 
ylim([5 30]); 

xlabel('v_{M2}', 'FontSize', 12);
ylabel('v_{M3}', 'FontSize', 12);
zlabel('Maximum Ca^{2+}', 'FontSize', 12);

grid on;
view([-37.5 30]);

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
