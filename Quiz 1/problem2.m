clear; clc; close all;

% =============================
% PARAMETERS (same as yours)
% =============================
alpha1 = 0.1;
alpha2 = 3;
alpha3 = 3;

beta1  = 3;
beta2  = 1;
beta3  = 1;

K1 = 0.5;
K2 = 0.5;
K3 = 0.5;

n1 = 8;
n2 = 8;
n3 = 8;

% =============================
% TIME
% =============================
tspan = [0 50];   % longer to clearly see limit cycle

% =============================
% INITIAL CONDITION
% =============================
y0 = [0.01 0.01 0.01];

% =============================
% SOLVE ODE
% =============================
[t, y] = ode45(@(t,y) model3(t,y,alpha1,alpha2,alpha3,...
                            beta1,beta2,beta3,...
                            K1,K2,K3,n1,n2,n3), tspan, y0);

C = y(:,1);
P = y(:,2);
A = y(:,3);

% =============================
% PHASE PLANE PLOTS
% =============================
figure;

% ---- CDK vs APC ----
subplot(1,3,1);
plot(C, A, 'LineWidth', 2);
xlabel('CDK');
ylabel('APC');
title('Phase: CDK vs APC');
grid on;

% ---- CDK vs PIK ----
subplot(1,3,2);
plot(C, P, 'LineWidth', 2);
xlabel('CDK');
ylabel('PIK');
title('Phase: CDK vs PIK');
grid on;

% ---- PIK vs APC ----
subplot(1,3,3);
plot(P, A, 'LineWidth', 2);
xlabel('PIK');
ylabel('APC');
title('Phase: PIK vs APC');
grid on;

% =============================
% ODE FUNCTION
% =============================
function dydt = model3(~, y, alpha1, alpha2, alpha3, ...
                             beta1, beta2, beta3, ...
                             K1, K2, K3, n1, n2, n3)

    C = y(1);
    P = y(2);
    A = y(3);

    x1 = A^n1 / (A^n1 + K1^n1);
    x2 = C^n2 / (C^n2 + K2^n2);
    x3 = P^n3 / (P^n3 + K3^n3);

    dC = alpha1 - beta1*C*x1;
    dP = alpha2*(1-P)*x2 - beta2*P;
    dA = alpha3*(1-A)*x3 - beta3*A;

    dydt = [dC; dP; dA];
end
