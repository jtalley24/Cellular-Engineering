function dydt = Gprotein(t,Y)
L = 10^-8; % M
P = 10^-7; % M
R = 10^5; % receptors
Kf = 10^7; % M^-1*S^-1
Kr = 0.1; % sec^-1
Gt = 10^4; % proteins
Ka = 10^-4; % cell/molecules-sec
Kd = 10^7; % M^-1*S^-1
K1 = 1.0; % sec^-1
Km = 2*10^-6; % M
S = 10^-6; % M
% Y(1) = C(t), Y(2) = G_active(t), Y(3) = Q(t)
dydt = zeros(3,1);
dydt(1) = Kf*L*R - Kr*Y(1); %dC/dT
dydt(2) = Ka*(Gt - Y(2))*Y(1) - Kd*P*Y(2); %dG_active/dt
dydt(3) = K1*S*Y(2)/(Km + S); %dQ/dt
end

