function RV = EQ2RV(EQ, I, GM)
%Convert Equinoctial Oribtal Elements to Cartesian position and velocity
%Reference: Sagovac, Neta, Early - Semianalytic Satellite Theory
%Last Updated: 4/7/2017
%
%INPUT:
%EQ ----- Equinoctial Orbtial Elements (a, h, k, p, q, lambda)
%I  ----- Retrograde Factor (+1 for direct, -1 for retrograde)
%GM ----- Gravitational Parameter
%
%OUTPUT:
%RV ----- Cartesian position and velocity (x, y, z, x_dot, y_dot, z_dot)


%Extract equinoctal elements from input vector
a = EQ(1);
h = EQ(2);
k = EQ(3);
p = EQ(4);
q = EQ(5);
lambda = EQ(6);

%I Retrograde Factor (+1 for direct, -1 for retrograde)

%Determine equinoctial reference frame basis
f = 1/(1+p^2+q^2)*[1-p^2+q^2; 2*p*q; -2*I*p];
g = 1/(1+p^2+q^2)*[2*I*p*q; (1+p^2-q^2)*I; 2*q];

%Solve equinoctial form of Kepler's equation using Newton's Method
F = lambda;
F_iplus1 = Inf;
while (abs(F_iplus1-F) > 1e-10)
    F_iplus1 = F - ...
        (F+h*cos(F)-k*sin(F)-lambda)/(1-h*sin(F)-k*cos(F));
    F = F_iplus1;
end

%Define auxiliary quantities
n = sqrt(GM/a^3);
b = 1/(1+sqrt(1-h^2-k^2));

%Calculate true longitude
SIN_L = ((1-h^2*b)*sin(F)+h*k*b*cos(F)-h)/(1-h*sin(F)-k*cos(F));
COS_L = ((1-h^2*b)*cos(F)+h*k*b*sin(F)-k)/(1-h*sin(F)-k*cos(F));

L = atan2(SIN_L, COS_L);

%Calculate position and velocity in equionctial reference frame
r = a*(1-h*sin(F)-k*cos(F));

X = r*cos(L);
Y = r*sin(L);

X_dot = -n*a*(h+sin(L))/sqrt(1-h^2-k^2);
Y_dot = n*a*(k+cos(L))/sqrt(1-h^2-k^2);

R = X*f + Y*g;
V = X_dot*f+Y_dot*g;

RV = [R;V];

end