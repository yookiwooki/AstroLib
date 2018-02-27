%Function to convert classical orbital elements to orbital elements with M
%Classical Orbital Elements - with True Anomaly (nu)

%Sean McArdle
%10/16/2016

function M_OE = COE2M_OE(COE)

%Calculate E

E = nu2E(COE(2),COE(6));

%Calculate M

M = E-COE(2)*sin(E);

%Build Orbital Element vector with M

M_OE = COE;
M_OE(6) = M;

end