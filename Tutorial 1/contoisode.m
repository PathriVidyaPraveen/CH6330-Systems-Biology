function dydt = contoisode(~,y,P)
dydt = zeros(3,1);
mumax = P(1,1);
kd = P(1,2);
YXS = P(1,3);
ms = P(1,4);
alpha = P(1,5);
kp = P(1,6);
kx = P(1,7);
k = P(1,8);
mu = mumax * y(2)/(kx * y(1)+y(2));
dydt(1) = mu*y(1)-kd*y(1);      %% This is the quation of dX/dt
dydt(2) = -(((1/YXS)*mu*y(1)) + ((ms*y(1)*y(2))/(k+y(2)))); %% dS/dt
%% Write the third equation of the system (~ 1 line of code)
dydt(3) = alpha*mu*y(1) - kp*y(3);
%% 
end
