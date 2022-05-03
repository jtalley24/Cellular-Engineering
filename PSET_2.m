% plot X/Rt vs 2L/Kd

kxrt = 10^-2; % 1, 10^2
TwoL_Kd = 10^-3:10^3;
X_Rt = 0.5*kxrt*TwoL_Kd*((((1+TwoL_Kd)/(2*kxrt*TwoL_Kd))*(-1+(((1+4*kxrt*TwoL_Kd)/((1+TwoL_Kd).^2)).^0.5))).^2);
semilogx(TwoL_Kd, X_Rt);
hold on;
kxrt = 1;
X_Rt = 0.5*kxrt*TwoL_Kd*((((1+TwoL_Kd)/(2*kxrt*TwoL_Kd))*(-1+(((1+4*kxrt*TwoL_Kd)/((1+TwoL_Kd).^2)).^0.5))).^2);
semilogx(TwoL_Kd, X_Rt);
hold on;
kxrt = 10^2;
X_Rt = 0.5*kxrt*TwoL_Kd*((((1+TwoL_Kd)/(2*kxrt*TwoL_Kd))*(-1+(((1+4*kxrt*TwoL_Kd)/((1+TwoL_Kd).^2)).^0.5))).^2);
semilogx(TwoL_Kd, X_Rt);