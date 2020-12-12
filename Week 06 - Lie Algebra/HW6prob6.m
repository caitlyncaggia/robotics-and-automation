% ECE 4560 - Homework 6, Problem 6
% Caitlyn Caggia

g1 = SE3([1;2;3],[1 0 0;0 0 -1; 0 1 0]);
g2 = SE3([2;-1;-1], [sqrt(2)/2 0 -sqrt(2)/2;0 1 0;sqrt(2)/2 0 sqrt(2)/2]);
p1 = [4;5;6];
p2 = [4;5;6;1];

disp('g3')
g3 = g1*g2

disp('g3.*p1')
ans1 = g3.*p1;
disp(ans1)

disp('g3.*p2')
ans2 = g3.*p2;
disp(ans2)

disp('inv(g3)')
ans3 = inv(g3)
