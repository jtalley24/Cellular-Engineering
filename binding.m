function dydt = binding(t,Y)
kf = 2*10^8; 
kr = 0.2;
Rt = 3*10^5;
L = 3*10^-10;
n = 10^6;
% Y(1) = C(t), Y(2) = L(t)
dydt = zeros(2,1);
dydt(1) = kf*(Rt-Y(1))*Y(2) - kr*Y(1);
dydt(2) = (-n/(6.022*10^23))*(kf*(Rt-Y(1))*Y(2) - kr*Y(1));
end

