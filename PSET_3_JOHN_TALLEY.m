% John Talley
% BE 306 PSET #3

%% 1.a.

L = [10^-10, 3*10^-10, 10^-9]; % mol/L
n = [2*10^5, 10^6, 5*10^6]; % cells/mL
n = n.*10^3; % mL/L
Rt = 3*10^5; % receptors/cell
Nav = 6.022*10^23;
eta = zeros(3,3);

for i=1:3
    for j =1:3
        eta(i,j) = n(j)*Rt/(Nav*L(i));
    end
end

% Depletion occurs in all except eta(3,1) < 0.1 case

% NO DEPLETION:

% U = C/Rt, <- this is "binding level", as in complexes / total receptors
% C = Rt*Lt/(Lt+Kd)

U = L(3)/(L(3) + 10^-9);
% U = 0.500

% DEPLETION:

U = zeros(3,3); % ignore (3,1) where there is no depletion
Kd = 10^-9; % mol/L

for i=1:3
    for j =1:3
        U(i,j) = ((1 + eta(i,j) + Kd/L(i))/(2*eta(i,j)))*(1 - (1 - (4*eta(i,j)/((1 + eta(i,j) + Kd/L(i))^2)))^0.5);
    end
end

% SEE RESULTS IN WRITTEN TABLE

%% 1.b.
%Uss -> Co
% convert to moles/L
% U = C/Rt

Uss = U(2,2);
Co = Uss*Rt; % complexes / cell initially

tspan = 0:100;

[t,Y] = ode45(@binding, tspan, [Co; 0]);

complexes = Y(:,1); % complexes / cell
free_ligand = Y(:,2)*Nav/(n(2)*0.001); % free ligands / cell

figure;
plot(tspan, complexes);
hold on;
plot(tspan, free_ligand);
legend({'EGF-Receptor Complex', 'Free EGF Ligand'},'Location', 'east');
xlabel('Time (s)');
ylabel('Molecules per Cell');
title('Dissocation of EGF from Receptor');

% 26.1921 complexes / cell at S.S. (by inspection of column vector)

%% 2.a.

time = 0:800; % 15 min = 800 seconds

[t1, Y1] = ode45(@endocytosis, time, [3*10^5; 0; 0; 0]);

surf_receptors = Y1(:,1);
surf_complexes = Y1(:,2);
intra_complexes = Y1(:,3);
intra_ligand = Y1(:,4);

figure;
plot(time, surf_receptors);
hold on;
plot(time, surf_complexes);
plot(time, intra_complexes);
plot(time, intra_ligand);

%% 4.

timespan = 1:5;
surface = [4000, 5000, 5500, 5700, 6000];
internal = [1250, 2600, 4700, 5900, 7000];

mdl = fitlm(surface, internal);

figure;
plot(surface, internal);
xlabel('Surface Complexes');
ylabel('Internal Complexes');


fun = @(x) internal(timespan)./2.8883 + 10/2.8883;
int_Cs = integral(fun, 1, 5, 'ArrayValued', true);

mdl2 = fitlm(int_Cs, internal);