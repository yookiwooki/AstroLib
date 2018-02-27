%Sean McArdle
%10/16/2016

%True Anomaly to Eccentric Anomaly
%From Vallado 4th ed. pg 77
%Considering elliptical case only

function E = nu2E(e,nu)

SINE_E = (sin(nu)*sqrt(1-e^2))/(1+e*cos(nu));
COSINE_E = (e+cos(nu))/(1+e*cos(nu));

E = 2*atan2(SINE_E,COSINE_E+1);

end