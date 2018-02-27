function EQ = RV2EQ( RV, I, GM)
R = RV(1:3);
R_dot = RV(4:6);

%Compute semimajor axis by inverting energy integral
a = 1/(2/norm(R)-norm(R_dot)^2/GM);

%Compute equionoctial basis vector w
w = cross(R, R_dot)/norm(cross(R, R_dot));

%Compute equinoctial elements p and q
p = w(1)/(1+I*w(3));
q = -w(2)/(1+I*w(3));

%Compute equinoctial basis vectors f and g
f = 1/(1+p^2+q^2)*[1-p^2+q^2; 2*p*q; -2*I*p];
g = 1/(1+p^2+q^2)*[2*I*p*q; (1+p^2-q^2)*I; 2*q];

%Find eccentricity vector
e = -R/norm(R) + cross(R_dot, cross(R, R_dot))/GM;

%Compute equinoctial elements h and k
h = dot(e, g);
k = dot(e, f);

%Compute mean longitude by finding equinoctial pos and solving Kepler's eq
X = dot(R, f);
Y = dot(R, g);

b = 1/(1+sqrt(1-h^2-k^2));

SIN_F = h + ((1-h^2*b)*Y-h*k*b*X)/(a*sqrt(1-h^2-k^2));
COS_F = k + ((1-k^2*b)*X-h*k*b*Y)/(a*sqrt(1-h^2-k^2));

F = atan2(SIN_F, COS_F);

lambda = F + h*cos(F) - k*sin(F);

%Collect components into equionoctial element set
EQ = [a; p; q; h; k; lambda];
end