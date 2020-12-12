% ECE 4560 - Homework 6, Problem 1
% Caitlyn Caggia

za = [1; -1; 1/4];
zb = [1; -1; 0];
gi = SE2([0; 0.5], pi/6);
tau = pi;
J = [0 1; -1 0];

%part a
disp('part a:')
parta = SE2.exp(za,tau);
gfa = gi .* parta

%part b
disp('part b:')
partb = SE2.exp(zb, tau);
gfb = gi .* partb


  