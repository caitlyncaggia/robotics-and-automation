% ECE 4560 - Homework 4, Problem 2
% Caitlyn Caggia

%from HW2prob8:
Arot = R(pi/3);
Atrans = [5;12];
A = [Arot Atrans];
A = [A; 0 0 1];
Brot = R(pi);
Btrans = [2;1];
B = [Brot Btrans];
B = [B; 0 0 1];

%part a
Aprime = R(-pi/2)*Atrans + [1; 0]
Bprime = R(-pi/2)*Btrans + [1; 0]

%part b
gOAprime = [Arot Aprime];
gOAprime = [gOAprime; 0 0 1];
gOBprime = [Brot Bprime];
gOBprime = [gOBprime; 0 0 1];

gAprimeBprime = gOAprime * inv(gOBprime)

%part c
gAB = A * inv(B);
h = [R(-pi/2) [1;0]];
h = [h; 0 0 1];

gApBp = h * gAB * inv(h)
% gApBp == gAprimeBprime 


