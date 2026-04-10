clear; clc; close all;

% =============================
% PARAMETERS
% =============================
alpha1 = 0.1;
beta1  = 3;

K1 = 0.5;
n1 = 8;

% =============================
% TIME
% =============================
tspan = [0 20];

% =============================
% INITIAL CONDITION
% =============================
C0 = 0.01;

% =============================
% SOLVE ODE
% =============================
[t, C] = ode45(@(t,y) model1(t,y,alpha1,beta1,K1,n1), tspan, C0);

% =============================
% PLOT
% =============================
figure;
plot(t, C, 'LineWidth', 2);
xlabel('Time');
ylabel('CDK activity');
title('1-ODE Model with Hill Nonlinearity');
grid on;

% =============================
% ODE FUNCTION
% =============================
function dCdt = model1(~, C, alpha1, beta1, K1, n1)

    % Hill function
    hill = (C^n1) / (K1^n1 + C^n1);

    % ODE
    dCdt = alpha1 - beta1 * C * hill;
end
