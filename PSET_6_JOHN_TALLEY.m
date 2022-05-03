% John Talley
% BE 306 PSET #6

%% 2.b. & c.

timespan = 0:0.1:20; %seconds

[t,Y] = ode45(@Gprotein, timespan, [0;0;0]);
complexes = Y(:,1);
G_active = Y(:,2);
Q_product = Y(:,3);

figure;
plot(G_active);
title('Active G-protein Generated in 20 Seconds');
xlabel('Time (0.1 sec)');
ylabel('Active G-protein (molecules)');
