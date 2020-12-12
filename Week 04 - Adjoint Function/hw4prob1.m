% ECE 4560 - Homework 4, Problem 1
% Caitlyn Caggia

g1 = SE2([1;2], pi/3);
g2 = SE2([-2,1], pi/6);

z = adjoint(g2, g1)
