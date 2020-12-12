% ECE 4560 - Homework 7, Problem 2
% Caitlyn Caggia

disp('exponent function in SE(3)')
g = SE3.exp([2; 1; -3; pi/10; -pi/4; pi/9], 2)

disp('logarithm function in SE(3)')
ln = log(g, 2)