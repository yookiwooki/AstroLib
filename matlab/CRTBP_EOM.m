%Sean McArdle
%11/29/2016

%Circular Restricted Three Body Problem Equations of Motion


function dY = CRTBP_EOM(t,Y)

n = 6;

global mu_star

%Initalize dY
dY = zeros(6+36,1);

x = Y(1);
y = Y(2);
z = Y(3);
x_dot = Y(4);
y_dot = Y(5);
z_dot = Y(6);

dY(1) = Y(4);
dY(2) = Y(5);
dY(3) = Y(6);

%Maple generated solution for Fx+2y_dot
t2 = x + mu_star;
t3 = t2 ^ 2;
t4 = (y ^ 2);
t5 = (z ^ 2);
t6 = t3 + t4 + t5;
t7 = sqrt(t6);
t13 = x - 1 + mu_star;
t14 = t13 ^ 2;
t15 = t14 + t4 + t5;
t16 = sqrt(t15);
dY(4) = x - (1 - mu_star) / t7 / t6 * t2 - mu_star / t16 / t15 * t13 + (2 * y_dot);



%Maple generated solution for Fy-2x_dot
t3 = (x + mu_star) ^ 2;
t4 = (y ^ 2);
t5 = (z ^ 2);
t6 = t3 + t4 + t5;
t7 = sqrt(t6);
t13 = (x - 1 + mu_star) ^ 2;
t14 = t13 + t4 + t5;
t15 = sqrt(t14);
dY(5) = y - (1 - mu_star) / t7 / t6 * y - mu_star / t15 / t14 * y - (2 * x_dot);


%Maple generated solution for Fz
t3 = (x + mu_star) ^ 2;
t4 = (y ^ 2);
t5 = (z ^ 2);
t6 = t3 + t4 + t5;
t7 = sqrt(t6);
t13 = (x - 1 + mu_star) ^ 2;
t14 = t13 + t4 + t5;
t15 = sqrt(t14);
dY(6) = -(1 - mu_star) / t7 / t6 * z - mu_star / t15 / t14 * z;

%Extract STM
Phi = reshape(Y(7:end),n,n);

%Calculate new Phi
A = getA(x,y,z);
dPhi = A*Phi;

%Rearrange Phi as column and insert into output
dY(7:end) = reshape(dPhi,36,1);

end