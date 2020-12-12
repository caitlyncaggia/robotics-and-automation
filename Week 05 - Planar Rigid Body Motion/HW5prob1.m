% ECE 4560 - Homework 5, Problem 1
% Caitlyn Caggia

gOA = SE2([9 4], pi/2);
gOB = SE2([3 6], pi/3);
gBC = SE2([4 0], pi/12);

%part a
vA1 = [1; 1; 1];
vO1 = gOA * vA1; 
sprintf('v1 in frame O: [%f %f]', vO1(1), vO1(2))
vB1 = inv(gOB) * vO1;
sprintf('v1 in frame B: [%f %f]', vB1(1), vB1(2))


%part b
vC2 = [sqrt(3)/2; -1/2; 1];
vB2 = gBC * vC2;
sprintf('v2 in frame B: [%f %f]', vB2(1), vB2(2))
vO2 = gOB * vB2;
sprintf('v2 in frame O: [%f %f]', vO2(1), vO2(2))