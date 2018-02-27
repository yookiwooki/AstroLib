%Sean McArdle
%10/16/2016

%Eccentric Anomaly to True Anomaly
%From Vallado 4th ed. pg. 77
%Considering elliptical case only

function nu = E2nu(e,E)

SINE_NU = sin(E)*sqrt(1-e^2)/(1-e*cos(E));
COSINE_NU = (cos(E)-e)/(1-e*cos(E));

nu = 2*atan2(SINE_NU,COSINE_NU+1);

end