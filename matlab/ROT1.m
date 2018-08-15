%Sean McArdle
%10/16/2016

%Rotation 1 Function
%Takes an angle in radians and outputs the corresponding 1 Rotation Matrix

function ROT1_MATRIX = ROT1(angle)

ROT1_MATRIX = [[1 0 0];[0 cos(angle) sin(angle)];[0 -sin(angle) cos(angle)]];

end