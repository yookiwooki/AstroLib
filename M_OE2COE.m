%Function to convert orbital elements with M to classical orbital elements
%Classical Orbital Elements - with True Anomaly (nu)

%Sean McArdle
%10/16/2016



function COE = M_OE2COE(M_OE)

%Calculate eccentric anomaly
E = KepEqtnE(M_OE(6),M_OE(2));

%Calculate true anomaly
nu = E2nu(M_OE(2),E);

%Express input OE's as COE
COE = M_OE;
COE(6) = nu;

end