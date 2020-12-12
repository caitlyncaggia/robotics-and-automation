% ECE 4560 - Homework 5, Problem 2
% Caitlyn Caggia

gOA = SE2([5 2], pi/3);
gOB = SE2([2 7], pi/2);
gBC = SE2([0 3], -pi/6);

%part a
zA1 = SE2([2 5], -pi/12);
gBA = inv(gOB) * gOA;
zB1 = adjoint(zA1,inv(gBA))
ZO1 = adjoint(zA1,inv(gOA))

%part b
zC2 = SE2([1 -4], pi/15);
zB2 = adjoint(zC2,inv(gBC))
gOC = gOB * gBC;
zO2 = adjoint(zC2,inv(gOC))