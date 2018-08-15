%Sean McArdle
%3/20/2017

%Converts orbital elements to position and velocity
%Orbital elements input must use true anomaly
%From Vallado 4th Ed. pg. 118

function RV = COE2RV(COE,GM)
%Orbital Element (OE) inputs are
%1 Semimajor Axis (LU)
%2 Eccentricity
%3 Inclination (radians)
%4 Argument of Periapsis (radians)
%5 Right Ascension of Ascending Node (radians)
%6 True Anomaly (radians)

%Check if Circular Equatorial
if (COE(2) == 0 && COE(3) == 0)
    COE(4) = 0;
    COE(5) = 0;
end

%Check if Circular Inclined
if (COE(2) == 0)
    COE(4) = 0;
end

%Check if Elliptical Equatorial
if (COE(3) == 0)
    COE(5) = 0;
end

e = COE(2); %extract eccentricity
nu = COE(6); %extract true anomaly
p = COE(1)*(1-COE(2)^2); %calculate semiparameter

%Calculate r and v in perifocal coordinates
rPQW = [p*cos(nu)/(1+e*cos(nu)) p*sin(nu)/(1+e*cos(nu)) 0]';
vPQW = [-sqrt(GM/p)*sin(nu) sqrt(GM/p)*(e+cos(nu)) 0]';

%Calculate Rotation Matrix
if COE(3) == 0 %Equatorial case (i = 0)
    R_TOTAL = eye(3);
else %General case
    R_TOTAL = ROT3(-COE(5))*ROT1(-COE(3))*ROT3(-COE(4));
end

%Rotate r and v
r = R_TOTAL*rPQW;
v = R_TOTAL*vPQW;

%Concatenate into RV vector

RV = [r;v];

end