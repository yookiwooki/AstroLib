%Sean McArdle
%10/16/2016

%Rotation 2 Function
%Takes an angle in radians and outputs the corresponding 2 Rotation Matrix

function ROT2_MATRIX = ROT2(angle)

ROT2_MATRIX = [[cos(angle) 0 -sin(angle)];[0 1 0];[sin(angle) 0 cos(angle)]];

end
