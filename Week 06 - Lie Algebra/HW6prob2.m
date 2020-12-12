% ECE 4560 - Homework 6, Problem 2
% Caitlyn Caggia

gi = SE2([-3; 2], pi/4);
gf = SE2([3; -3], -pi/3);
tau = 1;

disp('lie group displacement (g):');
g = inv(gi) * gf;
display(g);

disp('lie algebra element (xi):');
xi = SE2.ln(g, tau);
disp(xi);


