clear; clc; close all;

% =============================
% PARAMETERS 
% =============================
alpha1 = 0.1;   % CDK synthesis
alpha2 = 3;     % APC effect on CDK

beta1  = 3;     % CDK degradation
beta2  = 1;     % APC production

K1 = 0.5;       % Hill constant
K2 = 0.5;
n1 = 8;         % Hill coefficient
n2 = 8;

% =============================
% TIME
% =============================
tspan = [0 25];

% =============================
% INITIAL CONDITIONS
% =============================
y0 = [0.01 0.01];   % [CDK, APC]

% =============================
% SOLVE ODE
% =============================
[t, y] = ode45(@(t,y) model2(t,y,alpha1,alpha2,beta1,beta2,K1,K2,n1,n2), tspan, y0);

C = y(:,1);
A = y(:,2);

% =============================
% PLOT (TEMPORAL PATTERN)
% =============================
figure;
plot(t, C, 'LineWidth', 2); hold on;
plot(t, A, 'LineWidth', 2);

xlabel('Time');
ylabel('APC/CDK activity');
title('2-ODE Model: Temporal Pattern of CDK and APC');
legend('CDK','APC');
grid on;

% =============================
% ODE FUNCTION
% =============================
function dydt = model2(~, y, alpha1, alpha2, beta1, beta2, K1, K2, n1 , n2)

    C = y(1);
    A = y(2);
    hill1 = A^n1/(A^n1 + K1^n1);

    dC = alpha1 - beta1*C*hill1;
    hill2 = C^n2/(C^n2 + K2^n2);
    dA = alpha2*(1-A)*hill2 - beta2*A;

    dydt = [dC; dA];
end
