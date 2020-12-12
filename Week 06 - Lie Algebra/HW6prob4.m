% ECE 4560 - Homework 6, Problem 4
% Caitlyn Caggia

gi = SE3([3; 5; 2], [0 1 0; 0 0 1; 1 0 0]);
gf = SE3([6; -3; 9], [sqrt(3)/2 1/4 sqrt(3)/4;
                      1/2 -sqrt(3)/4 -3/4;
                      0 sqrt(3)/2 -1/2]);
tau = 1;
                  
g = inv(gi) * gf;
xi = log(g, tau)

expXi = SE3.exp(xi, tau);
disp('gi*expXi')
exp = gi * expXi



