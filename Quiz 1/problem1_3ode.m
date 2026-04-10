clear; clc; close all;

% =============================
% PARAMETERS 
% =============================
alpha1 = 0.1;   % CDK synthesis
alpha2 = 3;     % APC effect on CDK
alpha3 = 3;

beta1  = 3;     % CDK degradation
beta2  = 1;     % APC production
beta3 = 1;

K1 = 0.5;       % Hill constant
K2 = 0.5;
K3 = 0.5;
n1 = 8;         % Hill coefficient
n2 = 8;
n3 = 8;
% =============================
% TIME
% =============================
tspan = [0 25];

% =============================
% INITIAL CONDITIONS
% =============================
y0 = [0.01 0.01 0.01];   % [CDK, APC]

% =============================
% SOLVE ODE
% =============================
[t, y] = ode45(@(t,y) model3(t,y,alpha1,alpha2,alpha3,beta1,beta2,beta3,K1,K2,K3,n1,n2,n3), tspan, y0);

C = y(:,1);
P = y(:,2);
A = y(:,3);

% =============================
% PLOT (TEMPORAL PATTERN)
% =============================
figure;
plot(t, C, 'LineWidth', 2); hold on;
plot(t, A, 'LineWidth', 2);
plot(t, P, 'LineWidth', 2);

xlabel('Time');
ylabel('APC/CDK/Plk activity');
title('3-ODE Model');
legend('CDK','APC','Plk');
grid on;

% =============================
% ODE FUNCTION
% =============================
function dydt = model3(~, y, alpha1, alpha2,alpha3, beta1, beta2,beta3, K1, K2,K3, n1 , n2, n3)

    C = y(1);
    P = y(2);
    A = y(3);
    x1 = A^n1 / (A^n1 + K1^n1);
    dC = alpha1 - beta1*C*x1;
    x2 = C^n2 / (C^n2 + K2^n2);
    dP = alpha2*(1-P)*x2 - beta2*P;
    x3 = P^n3 / (P^n3 + K3^n3);
    dA = alpha3*(1-A)*x3 - beta3*A;
    dydt = [dC;dP;dA];
end
