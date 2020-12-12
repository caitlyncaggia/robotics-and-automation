%ECE 4560 - Homework 4, Problem 4
%Caitlyn Caggia

theta1 = pi/3; d1 = [1;2];
theta2 = pi/6; d2 = [-2;1];

g1 = SE2(d1, theta1);
g2 = SE2(d2, theta2);

g1inv = inv(g1);
fprintf("inverse: \n")
display(g1inv)

product = g1 * g2;
fprintf("product: \n")
display(product)