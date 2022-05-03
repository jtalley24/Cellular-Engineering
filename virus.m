function dydt = virus(t,Y)
S = 2*10^2; %T-cells/day
c = 4; %day^-1
p = 50; %day^-1
dt = 6; %day^-1
Tmax = 2*10^4; %day^-1
delta = 0.7; %day^-1
N = 20; %virions
% Y(1) = T(t), Y(2) = Tstar(t), Y(3) = V(t)
dydt = zeros(3,1);
dydt(1) = S + p*Y(1)*(1-(Y(1)/Tmax)) - dt*Y(1); %dt/dT
dydt(2) = -delta*Y(2); %dTstar/dt
dydt(3) = N*delta*Y(2) - c*Y(3);
end

