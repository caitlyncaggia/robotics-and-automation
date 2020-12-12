% ECE 4560 - Homework 7, Problem 1
% Caitlyn Caggia

%part a ==================================================================
disp('part a:')
syms a1 a2 l1 l2
de = [l1*cos(a1) + l2*cos(a1+a2);
      l1*sin(a1) + l2*sin(a1+a2)]

  
%part b ==================================================================
disp('part b:')
l1 = 1; l2 = 0.5;
desol = [1.3595; 0.2113];

%using law of cosines...
r = sqrt(desol(1)^2 + desol(2)^2);
a1 = acos( (l1^2 + r^2 - l2^2) / (2*l1*r) ) + atan2(desol(2), desol(1))
a2 = acos( (l1^2 + l2^2 - r^2) / (2*l1*l2) ) - pi

%check solutions for alphas by plugging a1 and a2 back into equation from
%part a:
decheck = [l1*cos(a1) + l2*cos(a1+a2);
      l1*sin(a1) + l2*sin(a1+a2)]
sprintf('results from part a and b match, so alphas are correct')

%part c ==================================================================
sprintf('part c:', ...
    'acos() requires an input between 0 and 1 (inclusive)', ...
    'this condition is not met when:', ...
    'l1^2 + r^2 - l2^2 < 2*l1^2*r', ...
    '2*l1*l2 < l1^2 + l2^2 -r^2', ...
    'l2^2 < l1^2 +r^2', ...
    '0 < l1^2 - r^2 + l2^2')



  
  