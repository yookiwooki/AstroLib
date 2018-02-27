%Generate Poincare Map for planar Circular Restricted Three Body Problem

%Sean McArdle
%12/01/2016

%Ouptut: [time, state vector, event time, state vector at event time]
%Input: (initial x pos., initial x vel., Jacobi Constant, # of Crossings)
function [t,Y,te,Ye] = poincare(x0,u0,C,N,want_traj)

global mu_star

%Take y0 on the x axis
y0 = 0; %LU

%Solve for v0 using x0, y0, u0, C, and mu
r1 = sqrt((x0+mu_star)^2+y0^2);
r2 = sqrt((x0+mu_star-1)^2+y0^2);
v0 = sqrt(x0^2+y0^2+2*(1-mu_star)/r1+2*mu_star/r2-C-u0^2);

%Planar problem so z position and velocity are zero
z0 = 0;
w0 = 0;

Y0 = [x0 y0 z0 u0 v0 w0]';

tspan = [0 5000];

if ~isreal(sum(Y0))
    t = 0;
    Y = 0;
    te = 0;
    Ye = 0;
    disp(['ERROR: Imaginary component in initial conditions at (' num2str(x0) ',' num2str(u0) ').'])
else    
    %Propagate Equations of Motion for 5000 ascending x axis crossings
    options = odeset('AbsTol', 1e-10, 'Events', @poincareEventsFunction, 'RelTol', 1e-10);
    [t,Y,te,Ye,~] = ode45(@CRTBP_EOM,tspan,Y0,options);

    %Consider trajectory on first N ascending x axis crossings only
    if length(te) >= N
    t_end = te(N);
    index_t_end = find(t >= t_end, 1) - 1;
    t = t(1:index_t_end);
    Y = Y(1:index_t_end,:);
    te = te(1:N);
    Ye = Ye(1:N,:);
    
    else
    disp('ERROR: Distance to Moon is < 1e-4 LU.  Stopped Propagation.')
    end
end

end