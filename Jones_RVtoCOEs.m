function [COEs] = RVtoCOEs(R,V,MU)
%
% function [COEs] = RVtoCOEs(R,V,MU)
% ---------------------------------------------------------------------
% 
% Description:
%
%   Function to convert the provided radius and velocity vector to 
%   classical orbital elements (COEs).
% 
% Inputs:
%
%   R  - 3x1 vector representing position
%   V  - 3x1 vector representing velocity
%   MU - scalar gravitation parameter (G*M)
% 
% Outputs:
% 
%   A - 6x1 vector including:
%
%                   Semimajor Axis
%                   Eccentricity
%                   Inclination
%                   Right Ascention of the Ascending Node
%                   Argument of Periapse
%                   True Anomaly
%
% Assumptions/References:
%
%   NONE
%
% Dependencies:
%
%   NONE
%
% Modification History:
% 
%   17mar08     Brandon A. Jones      original version (header added)

%  Just make sure we have column vectors for future processing.
R = R(:);
V = V(:);
dotRV = R'*V;

%  Get the energy and the SMA
energy=norm(V)^2/2 - MU/norm(R);
a = -MU/2/energy;

%  Get the angular momentum vector
H = cross(R,V);

%  Get the eccentricity
%E = cross(V,H)/MU - R/norm(R)
E = ((norm(V)*norm(V)-MU/norm(R))*R - (dotRV)*V)/MU';
ecc = norm(E);

%  Get the inclination
inclination = myacos(H(3)/norm(H));

%  Now we get the right ascension of the ascending node
k = [0;0;1];
N = cross(k,H);
normN = norm(N);
node = myacos(N(1)/normN);
if N(2)<0
    node=2*pi-node;
end


arg = myacos(  dot(N,E) / ( normN*ecc )  );
if E(3)<0
    arg = 2*pi - arg;
end;

argument = dot(E,R)/(ecc*norm(R));
if abs(argument-1) < 1e-12
    argument = 1;
elseif abs(argument+1) < 1e-12
    argument = -1;
end
true = myacos(argument);
if dotRV < 0
    true = 2*pi - true;
end;

if ecc < 1e-10 && inclination < 1e-10
    
    %  Set these values to zero.
    node = 0.0;
    arg  = 0.0;
    
    %  Compute the true longitude of periapsis
    true = myacos( R(1)/norm(R) );
    if R(2) < 0
        true = 2*pi - true;
    end
    
elseif ecc < 1e-10
    
    %  Set the singular elements to zero
    arg = 0.0;
    
    %  Compute the argument of latitude
    true = myacos( dot(N,R)/(normN*norm(R)) );
    if R(3) < 0
        true = 2*pi - true;
    end
    
elseif inclination < 1e-10
    
    %  Set the singular elements to zero
    node = 0.0;
    
    %  Compute the Longitude of periapsis
    arg = myacos( E(1)/ecc );
    if E(2) < 0
        arg = 2*pi - arg;
    end
    
end

% eccentric = myacos( (ecc+cos(true))/(1+ecc*cos(true)) );
% if true>pi
%     eccentric = 2*pi - eccentric;
% end;


% mean = eccentric - ecc*sin(eccentric);


COEs = [a; ...
        ecc; ...
        inclination; ...
        node; ...
        arg; ...
        true ];
    
    function angle = myacos( arg )
        
        if abs(arg-1.0) < 1e-13 && arg > 1.0
            arg = 1.0;
        elseif abs(arg+1.0) < 1e-13 && arg < -1.0
            arg = 1.0;
        end
        angle = acos( arg );
        
    end

end
