%Events function for numerical integration in poincare function

%Sean McArdle
%12/01/2016

function [event,isterminal,direction] = poincareEventsFunction(t,y)

global mu_star

% Event 1 is when y position is zero
% Event 2 is when distance from the moon is less than 
event1 = y(2);

rM = norm([y(1) y(2)]-[1-mu_star 0]);
event2 = (rM < 1e-4);

event = [event1 event2];

% Do not halt integration for Event 1
% Halt integration for Event 2
isterminal = [0 1];

% Event 1 is only when the fucntion is increasing
% Event 2 is all zeros
direction = [1 0]; 
end