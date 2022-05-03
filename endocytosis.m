function dYdt = endocytosis(t1,Y1)
L = 10^-9; %mol/L
n = 10^9; %cells/L
ki = 0.3; %1/min
ki = ki/60; %1/sec
krec = 0.4; %1/min
krec = krec/60; %1/sec
kdeg = 0;
Rt = 3*10^5; %receptors/cell
kf = 2*10^8; %1/(M*s)
kr = 0.2; %1/s

% Y1(1) = Rs, Y1(2) = C, Y1(3) = Ci, Y1(4) = Li
dYdt = zeros(4,1);
dYdt(1) = -kf*Y1(1)*L + kr*Y1(2) + krec*Y1(3);
dYdt(2) = kf*Y1(1)*L - kr*Y1(2) - ki*Y1(2);
dYdt(3) = ki*Y1(2) - krec*Y1(3);
dYdt(4) = krec*Y1(3) - kdeg*Y1(4); % kdeg = 0

end

