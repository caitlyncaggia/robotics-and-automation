%ECE 4560 - Homework 4, Problem 1
%Caitlyn Caggia

syms a1 a2 a3 l1 l2 l3;

d1 = [l1;0]; d2 = [l2;0]; d3 = [l3;0];
a = [a1; a2; a3];

fprintf("part a: ")
gOEa = [R(a1+a2+a3) R(a1)*d1 + R(a1+a2)*d2 + R(a1+a2+a3)*d3;
    0 0 1]

a1 = -pi/12; a2 = pi/6; a3 = -2*pi/3;
l1 = 1; l2 = 1/2; l3 = 1/4;
adot = [1/6; -1/2; 1/4];
d1 = [l1;0]; d2 = [l2;0]; d3 = [l3;0];

fprintf("part b: ")
gOEb = [R(a1+a2+a3) R(a1)*d1 + R(a1+a2)*d2 + R(a1+a2+a3)*d3;
    0 0 1]

fprintf("part c:")
J = jacobian(gOEa(:,3), a)

fprintf("part d: ")
vsymbolic = J * adot
v1 =  (l2*sin(a1 + a2))/3 - (l1*sin(a1))/6 + (l3*sin(a1 + a2 + a3))/12;
v2 = (l1*cos(a1))/6 - (l2*cos(a1 + a2))/3 - (l3*cos(a1 + a2 + a3))/12;
v3 = 0;

v = [v1; v2; v3]
