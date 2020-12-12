%ECE 4560 - Homework 4, Problem 4
%Caitlyn Caggia

%part a
parta = RotX(pi/4)

%part b
partb = RotY(pi/6)*RotZ(3*pi/4)

%part c
partc = RotX(pi/3)*RotZ(-pi/4)

%part d
partd1 = RotY(pi/2)*RotY(pi/2)
partd2 = RotY(pi)
%partd1 == partd2



function Rx = RotX(theta)

Rx = [ 1 0 0;
       0 cos(theta) -sin(theta);
       0 sin(theta) cos(theta)];
end

function Ry = RotY(theta)

Ry = [ cos(theta) 0 sin(theta);
       0 1 0;
       -sin(theta) 0 cos(theta)];
end

function Rz = RotZ(theta)

Rz = [ cos(theta) -sin(theta) 0;
       sin(theta) cos(theta) 0;
       0 0 1];
end