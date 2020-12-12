%ECE 4560 - Homework 4, Problem 3
%Caitlyn Caggia

a1 = pi; a2 = pi/8; a3 = -pi/4; a4 = pi/8;
l1 = 1; l2 = 0.75; l3 = 0.75;

%syms a1 a2 a3 a4 l1 l2 l3;

d1 = [l1;0]; d2 = [l2;0]; d3 = [l3;0];
a = [a1; a2; a3; a4];

gOE = [R(a1+a2+a3+a4) R(a1)*d1 + R(a1+a2)*d2 + R(a1+a2+a3)*d3;
    0 0 1]


