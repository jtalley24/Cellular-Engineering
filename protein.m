function dYdt = protein(t1,Y1)
n = 3;
alpha = 2.0;
alpha0 = 0.025;
dYdt = zeros(3,1);
dYdt(1) = -Y1(1) + alpha/(1+(Y1(2)^n)) + alpha0; %protein1
dYdt(2) = -Y1(2) + alpha/(1+(Y1(3)^n)) + alpha0; %protein3
dYdt(3) = -Y1(3) + alpha/(1+(Y1(1)^n)) + alpha0; %protein2
end

