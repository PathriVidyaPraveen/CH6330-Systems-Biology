tspan = [0:2:32];
y0 = [0.01 40 0];
P = [0.8400 0.0146 0.5826 0.0134 0.7772 0.0372 7.4455 1.1086];
[t,y] = ode45(@(t,y) contoisode(t,y,P), tspan, y0);
%% Write the syntax to plot the simulation results
plot(t,y(:,1),'r','LineWidth',2); hold on;
plot(t,y(:,2),'g','LineWidth',2);
plot(t,y(:,3),'b','LineWidth',2);

xlabel('Time');
ylabel('Concentration');
legend('Biomass (X)','Substrate (S)','Product (P)');
title('Contois Model Simulation');
%% end
