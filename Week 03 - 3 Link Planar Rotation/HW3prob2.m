%ECE 4560 - Homework 4, Problem 2
%Caitlyn Caggia

gOA = [R(-3*pi/4) [7;2]; 0 0 1];
gOB = [R(pi/2) [0;8]; 0 0 1];

gAB = inv(gOA) * gOB;

theta = 5*pi/4;
d = gAB(:,3);
d(3) = [];

q = inv(R(theta) - [1 0; 0 1]) * (-d)