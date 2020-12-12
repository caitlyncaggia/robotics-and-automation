% ECE 4560 - Homework 6, Problem 3
% Caitlyn Caggia

what = [0     -1   -0.25; 
        1      0    0.5; 
        0.25  -0.5  0];
tau = 2;

w1 = what(3,2); 
w2 = what(1,3);
w3 = what(2,1);
w = [w1; w2; w3];
wmag = sqrt(w1^2 + w2^2 + w3^2);

%Rodrigues' formula
expW = eye(3) + (what./wmag)*sin(wmag*tau) + (what^2./wmag^2)*(1-cos(wmag*tau))