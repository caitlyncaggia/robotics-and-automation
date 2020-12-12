% ECE 4560 - Homework 13.3
% Caitlyn Caggia

g = [-3; 1; -pi/4];
h = [0; 1/8; -pi/4];

%part a
Fea = [1; 0; 0];
Rot = R(h(3));
adjTh = [Rot' zeros(2,1); h(2) h(1) 1];
Fca = inv(adjTh)*Fea

%part b
Fcb = [-3/(4*sqrt(2)); 1/(4*sqrt(2)); 3/(32*sqrt(2))];
Feb = adjTh*Fcb
