%Sean McArdle
%10/16/2016

%Converts position and velocity to orbital elements
%Orbital elements output uses true anomaly
%From Vallado 4th Ed. pg 113

function COE = RV2COE(RV,mu)

%Orbital Element (OE) outputs are
%1 Semimajor Axis (LU)
%2 Eccentricity
%3 Inclination (radians)
%4 Argument of Periapsis (radians)
%5 Right Ascension of Ascending Node (radians)
%6 True Anomaly (radians)

%Split RV vector into position and velocity vectors and find magnitudes
r_vector = RV(1:3);
v_vector = RV(4:6);
v_mag = norm(v_vector);
r_mag = norm(r_vector);

%Calculate angular momentum vector and magnitude
h_vector = cross(r_vector,v_vector);
h_mag = norm(h_vector);

%Calculate node vector
K_vector = [0 0 1]; %unit vector K
n_vector = cross(K_vector,h_vector);
n_mag = norm(n_vector);

%Calculate eccentricity vector and magnitude
e_vector = ((v_mag^2-mu/r_mag)*r_vector-dot(r_vector,v_vector)*v_vector)/mu;
e_mag = norm(e_vector);

%Calculate specific mechanical energy (ksi)
ksi = v_mag^2/2-mu/r_mag;

%Calculate semimajor axis and semiparameter
if (e_mag ~= 1.0) %Non-parabolic orbits
    a = -mu/(2*ksi);
    p = a*(1-e_mag^2);
else %Parabolic orbits
    p = h_mag^2/mu;
    a = inf;
end

%Calculate inclination
i = acos(h_vector(3)/norm(h_vector));

%Calculate right ascension of the asecnding node
OMEGA = acos(n_vector(1)/norm(n_vector));
%Resolve quadrant: If J component of n vector is negative, pi<OMEGA<2*pi
if (n_vector(2)<0)
    OMEGA = 2*pi-OMEGA;
end

%Calculate argument of periapsis
omega = acos(dot(n_vector, e_vector)/(norm(n_vector)*norm(e_vector)));
%Resolve quadrant: If K component of e vector is negative, pi<omega<2*pi
if (e_vector(3)<0)
    omega = 2*pi-omega;
end

%Calculate true anomaly
nu = acos(dot(e_vector,r_vector)/(norm(e_vector)*norm(r_vector)));
%Resolve quadrant: if flight path angle is negative, pi<nu<2*pi
if (dot(r_vector,v_vector) < 0)
    nu = 2*pi-nu;
end

COE = [a e_mag i omega OMEGA nu]';

% %Special Cases
% 
% %Elliptical equatorial
% omega_true = acos(e_vector(1)/e_mag);
% if (e_vector(2)<0)
%     omega_true = 2*pi-omega_true;
% end
% 
% %Circular Inclined
% u = acos(dot(n_vector,r_vector)/(n_mag*r_mag));
% if (r_vector(3) < 0)
%     u = 2*pi-u;
% end
% 
% %Circular Equatorial
% lambda_true = r_vector(1)/r_mag;
% if (r_vector(2)<0)
%     lambda_true = 2*pi-lambda_true;
% end

end