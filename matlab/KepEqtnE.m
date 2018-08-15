%Sean McArdle
%10/16/2016

%Kepler's Equation Solver (Elliptical Solutions)
%Uses Newton-Raphson Method to find eccentric anomaly E
%Inputs are Mean Anomaly and Eccentricity
%Algorithm from Vallado 4th Ed. pg. 65

function E = KepEqtnE(M,e)

tolerance = 1E-8; %set convergence tolerance

%Choose appropriate first guess for E
if ((-pi<M && M<0) || M>pi)
    E = M-e;
else
    E = M+e;
end

i = 0; %Initialize iteration counter
error = tolerance+1; %Initialize error value

%Counts iterations 
while(error>tolerance && i<100)
    Etest = E+(M-E+e*sin(E))/(1-e*cos(E)); %Calculate new E value
    error = abs(Etest-E); %Calculate error value on current iteration
    i=i+1; %Add 1 to iteration counter
    E = Etest; %Set E to new value
    if i == 100
        disp('Error: Maximum iterations reached for KepEqtnE.')
    end
end

end