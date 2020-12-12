% ECE 4560 - Homework 5, Problem 3
% Caitlyn Caggia

gOA = SE2([5 2], pi/3);
gOB = SE2([2 7], pi/2);
gBC = SE2([0 3], -pi/6);
zbB = SE2([2 -3], pi/9);


zbC = adjoint(zbB, inv(gBC))