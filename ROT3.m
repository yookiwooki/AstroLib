%Sean McArdle
%10/16/2016

%Rotation 3 Function
%Takes an angle in radians and outputs the corresponding 3 Rotation Matrix

function ROT3_MATRIX = ROT3(angle)

ROT3_MATRIX = [[cos(angle) sin(angle) 0];[-sin(angle) cos(angle) 0];[0 0 1]];

end
