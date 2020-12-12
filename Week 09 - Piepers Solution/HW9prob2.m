% ECE 4560 - Homework 9.2
% Caitlyn Caggia

g1 = SE3([-5;6;-1] , [1 0 0;0 0 -1;0 1 0]);
g2 = SE3([ 0;0; 0] , [ 1/sqrt(2) 0 1/sqrt(2);0 1 0;-1/sqrt(2) 0 1/sqrt(2)]);
adjoint(g1, g2)

g = SE3([1;0;-1],[sqrt(2)/2 0 -sqrt(2)/2 ; 0 1 0; sqrt(2)/2 0 sqrt(2)/2]);
xi = [0; 1; 0; pi/10; 0; pi/12];
adjoint(g, xi)

xi = [1; 0; 1; 0; pi/6; -pi/12];
adjoint(g, xi)
