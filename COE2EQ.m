function EQ = COE2EQ(COE, I)
%Converts classical orbital elements to Equinoctial

e = COE(2);
i = COE(3);
omega = COE(4);
Omega = COE(5);
nu = COE(6);

h = e*sin(omega + I*Omega);
k = e*cos(omega + I*Omega);
p = (tan(i/2))^(I)*sin(Omega);
q = (tan(i/2))^(I)*cos(Omega);

E = nu2E(e,nu);
M = E - e*sin(E);
lambda = M + omega + I*Omega;

EQ(1) = COE(1);
EQ(2) = h;
EQ(3) = k;
EQ(4) = p;
EQ(5) = q;
EQ(6) = lambda;

end