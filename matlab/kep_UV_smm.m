% Sean McArdle
% Universal Variables Kepler solver
% 2/6/2018

function [RV, histStore, finalOut] = kep_UV_smm(RV0, TOF, mu, outFlag)

%outFlag - (1 = RV), (2 = RV, histStore), (3 = RV, histStore, finalOut)

% Step 1: get position magnitude and SMA
r0Mag = sqrt(RV0(1)^2 + RV0(2)^2 + RV0(3)^2);
v0Mag = sqrt(RV0(4)^2 + RV0(5)^2 + RV0(6)^2);

a = 0.5*(-mu/(v0Mag^2/2-mu/r0Mag));
alpha = -v0Mag^2/mu + 2/r0Mag; % 1/a

% Step 2: Solve UV TOF Equation with Newton step root solve

% Initial guess - see Vallado KEPLER algorithm

if (abs(TOF) < 1.0e-6)  % Guess zero if TOF close to zero
    X = 0;
else
    if (alpha > 1.0e-6) % Elliptical
        X = sqrt(mu)*TOF*alpha;
    elseif (abs(alpha) < 1.0e-6) % Parabolic
        p = (norm(cross(RV0(1:3),RV0(4:6))))^2/mu;
        s = acot(3*sqrt(mu/(p^3))*TOF)/2;
        w = atan(tan(s)^(1/3));
        X = sqrt(p)*2*cot(2*w);
    else % Hyperbolic
        X = sign(TOF)*sqrt(-a)*log(-2*mu*alpha*TOF/...
            (dot(RV0(1:3),RV0(4:6))+sign(TOF)*sqrt(-mu*a)*(1-r0Mag*alpha)));
    end
end

DeltaX = 1; % Something greater than tol

tol = 1.0e-12;

iter = 0;
itermax = 100;

% Initalize storage for Newton step history output
if (outFlag > 1)
    histStore = zeros(itermax,3);
end

while(abs(DeltaX) > tol && iter < itermax)
    z = X^2*alpha;
    
    [rMag,t] = UV_t_r(X, RV0, z, mu);
    
    DeltaX = sqrt(mu)*(TOF-t)/rMag;
    X = X + DeltaX;
    
    iter = iter + 1;
    
    if (outFlag > 1)
        histStore(iter,1) = iter;
        histStore(iter,2) = X;
        histStore(iter,3) = TOF-t;
    end
end

% Prune history output storage
if(outFlag > 1)
    histStore = histStore(any(histStore,2),:);
end

if (iter == itermax)
    disp('ERROR: Maximum Newton step iterations reached. (kep_UV_smm)')
end

% Step 3: Calculate f, g, fDot, gDot
[C,S] = get_Stumpff(z);
f = 1 - X^2*C/r0Mag;
g = t - X^3*S/sqrt(mu);
fDot = sqrt(mu)*X*(z*S - 1)/(r0Mag*rMag);
gDot = 1 - X^2*C/rMag;

% Check for accuracy condition
checkTol = 1.0e-5;
if (abs(1-(f*gDot-fDot*g)) > checkTol)
    disp('WARNING:  1 = f*gDot - fDot*g not being met. (kep_UV_smm)')
end

% Step 4: Compute final position and velocity
RV(1:3) = f*RV0(1:3) + g*RV0(4:6);
RV(4:6) = fDot*RV0(1:3) + gDot*RV0(4:6);

% Store final values for all parameters
if (outFlag > 2)
    finalOut = cell(10,1);
    finalOut{1} = C;
    finalOut{2} = S;
    finalOut{3} = X;
    finalOut{4} = z;
    finalOut{5} = f;
    finalOut{6} = g;
    finalOut{7} = fDot;
    finalOut{8} = gDot;
    finalOut{9} = RV(1:3); %r2vec
    finalOut{10} = RV(4:6); %v2vec
end

end


% Calculates r and t for Universal Variable approach
function [rMag,t] = UV_t_r(X, RV0, z, mu)

[C,S] = get_Stumpff(z);

r0Mag = (RV0(1)^2 + RV0(2)^2 + RV0(3)^2)^(1/2);

rMag = X^2*C + (dot(RV0(1:3),RV0(4:6))/sqrt(mu))*X*(1-z*S) + r0Mag*(1-z*C);
t = (1/sqrt(mu))*((X^3)*S + (dot(RV0(1:3),RV0(4:6))/sqrt(mu))*X^2*C + r0Mag*X*(1-z*S));
end


% Calculates Stumpff function values
function [C,S] = get_Stumpff(z)

eps = 1.0e-6;

% Ellipse
if ((z > 0) && (abs(z) > eps))
    C = (1-cos(sqrt(z)))/z;
    S = (sqrt(z)-sin(sqrt(z)))/sqrt(z^3);
end

% Parabola (Series solution)
if (abs(z) <= eps)
    C = 1/2 - z/factorial(4) + z^2/factorial(6) - z^3/factorial(8);
    S = 1/factorial(3) - z/factorial(5) + z^2/factorial(7) - z^3/factorial(9);
end

% Hyperbola
if ((z < 0) && (abs(z) > eps))
    C = (1-cosh(sqrt(-z)))/z;
    S = (sinh(sqrt(-z))-sqrt(-z))/sqrt((-z)^3);
end
end