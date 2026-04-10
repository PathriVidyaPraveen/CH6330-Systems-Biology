clear; clc; close all;

% =============================
% BASE PARAMETERS
% =============================
p.alpha1 = 0.1;
p.alpha2 = 3;
p.alpha3 = 3;

p.beta1  = 3;
p.beta2  = 1;
p.beta3  = 1;

p.K1 = 0.5;
p.K2 = 0.5;
p.K3 = 0.5;

p.n1 = 8;
p.n2 = 8;
p.n3 = 8;

% =============================
% TIME + INITIAL CONDITIONS
% =============================
tspan = [0 25];
y0 = [0.01 0.01 0.01];

% =============================
% SENSITIVITY ANALYSIS
% =============================
factors = [0.4, 1, 1.6];   % -60%, baseline, +60%

% ---------- n1 variation ----------
figure;
hold on;
for i = 1:length(factors)
    p_temp = p;
    p_temp.n1 = p.n1 * factors(i);

    [t, y] = ode45(@(t,y) model3(t,y,p_temp), tspan, y0);
    plot(t, y(:,1), 'LineWidth', 2);
end
xlabel('Time'); ylabel('CDK');
title('Sensitivity: n1 variation');
legend('-60%','baseline','+60%');
grid on;

% ---------- n2 variation ----------
figure;
hold on;
for i = 1:length(factors)
    p_temp = p;
    p_temp.n2 = p.n2 * factors(i);

    [t, y] = ode45(@(t,y) model3(t,y,p_temp), tspan, y0);
    plot(t, y(:,1), 'LineWidth', 2);
end
xlabel('Time'); ylabel('CDK');
title('Sensitivity: n2 variation');
legend('-60%','baseline','+60%');
grid on;

% ---------- n3 variation ----------
figure;
hold on;
for i = 1:length(factors)
    p_temp = p;
    p_temp.n3 = p.n3 * factors(i);

    [t, y] = ode45(@(t,y) model3(t,y,p_temp), tspan, y0);
    plot(t, y(:,1), 'LineWidth', 2);
end
xlabel('Time'); ylabel('CDK');
title('Sensitivity: n3 variation');
legend('-60%','baseline','+60%');
grid on;

% =============================
% MODEL FUNCTION
% =============================
function dydt = model3(~, y, p)

    C = y(1);
    P = y(2);
    A = y(3);

    % Hill functions
    x1 = A^p.n1 / (A^p.n1 + p.K1^p.n1);
    x2 = C^p.n2 / (C^p.n2 + p.K2^p.n2);
    x3 = P^p.n3 / (P^p.n3 + p.K3^p.n3);

    % ODEs
    dC = p.alpha1 - p.beta1*C*x1;
    dP = p.alpha2*(1-P)*x2 - p.beta2*P;
    dA = p.alpha3*(1-A)*x3 - p.beta3*A;

    dydt = [dC; dP; dA];
end
