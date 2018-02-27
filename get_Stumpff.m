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
