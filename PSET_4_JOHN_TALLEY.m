% John Talley
% BE 306 PSET #4

%% 1.a.

Rt = 1;
k1 = 2; % min^-1
k2 = 1; %min^-1
Km1 = 0.05; %fixed
% Km2 = [0.01, 0.05, 0.1, 0.5];
Km2 = 0.01;

syms S Rp
eqn = k1*S*(Rt-Rp)/(Km1+Rt-Rp) - k2*Rp/(Km2+Rp) == 0;
soln_Rp = solve(eqn, Rp);

figure;
fplot(soln_Rp(2,1));
hold on;

Km2 = 0.05;

syms S Rp
eqn = k1*S*(Rt-Rp)/(Km1+Rt-Rp) - k2*Rp/(Km2+Rp) == 0;
soln_Rp = solve(eqn, Rp);

fplot(soln_Rp(2,1));
hold on;

Km2 = 0.1;

syms S Rp
eqn = k1*S*(Rt-Rp)/(Km1+Rt-Rp) - k2*Rp/(Km2+Rp) == 0;
soln_Rp = solve(eqn, Rp);

fplot(soln_Rp(2,1));
hold on;

Km2 = 0.5;

syms S Rp
eqn = k1*S*(Rt-Rp)/(Km1+Rt-Rp) - k2*Rp/(Km2+Rp) == 0;
soln_Rp = solve(eqn, Rp);

fplot(soln_Rp(2,1));
hold on;

xlim(0:1);
ylim(0:1);
xlabel('S');
ylabel('Rp');
legend('Km2 = 0.01', 'Km2 = 0.05', 'Km2 = 0.1', 'Km2 = 0.5', 'Location','Southeast');
title('Rp as Function of S @ Steady State');

%% 1.b.

Rt = 1;
k1 = 5; % min^-1
k2 = 1; %min^-1
S = 0.5; %fixed
% Km2 = [0.01, 0.05, 0.1, 0.5];
Km2 = 0.02;

syms Km1 Rp
eqn = k1*S*(Rt-Rp)/(Km1+Rt-Rp) - k2*Rp/(Km2+Rp) == 0;
soln_Rp = solve(eqn, Rp);

eqn1 = soln_Rp(2,1) == 0.5;
soln_Km1 = solve(eqn1, Km1);
% printing -> soln_Km1 = 4/5

%% 2.d.

tspan = 0:(1/1440):10;
[t,Y] = ode45(@virus, tspan, [5*10^3; 2.25*10^5; 7.885*10^5]);
Tcells = Y(:,1);
Tinfected = Y(:,2);
Virus = Y(:,3);

time = find(Virus < 7.885*10^4);
days = time(1)/1440;
%days = 3.56

%% 3.a.& b.
timespan = 0:0.01:50;
[t1,Y1] = ode45(@protein, timespan, [1;0;0]);
protein1 = Y1(:,1);
protein3 = Y1(:,2);
protein2 = Y1(:,3);

figure;
plot(protein1);
hold on;
plot(protein3);
hold on;
plot(protein2);
xlabel('Time');
ylabel('Protein');
title('Synthesis of Three Proteins Over Time');

%% 4.a.

N2t = 500; % molecules/um^2
Kd = 50; % molecules/um^2
N1t = linspace(100, 1000, 100);

Nb = (N1t+N2t+Kd)./2.*(1-sqrt(1-(4.*N1t.*N2t)./((N1t+N2t+Kd).^2)));

figure;
plot(N1t, Nb);
xlabel('Density of Receptors on Cell Surface (molecules/um^2)');
ylabel('Bond Density (molecules/um^2)');
title('Bond Density Between Two Cells at Steady State');

%% 4.b.

N2t = 500; % molecules/um^2
N1t = 500; % molecules/um^2

syms Kd
eqn2 = 100 == (N1t+N2t+Kd)./2.*(1-sqrt(1-(4.*N1t.*N2t)./((N1t+N2t+Kd).^2)));
soln_Kd = solve(eqn2, Kd);
% Kd = 1600 molecules/um^2