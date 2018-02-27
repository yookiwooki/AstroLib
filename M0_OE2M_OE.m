%Function to convert orbital elements with M0 to orbital elements with M

%Sean McArdle
%10/16/2016


function M_OE = M0_OE2M_OE(M0_OE, deltat, GM)

%Calculate current mean anomaly
n = sqrt(GM/M0_OE(1)^3);
M = M0_OE(6)+n*deltat;

%Express OE's with M instead of M0
M_OE = M0_OE;
M_OE(6) = M;

end