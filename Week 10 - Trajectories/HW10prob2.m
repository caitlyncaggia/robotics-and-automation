% ECE 4560 - Homework 10.2
% Caitlyn Caggia

l1 = 1; l2 = 0.5; l3 = 0.25;
l = [l1;l2;l3];
ai = [-pi/6; pi/4; -pi/3];
af = [0; -pi/12; pi/4];

syms t alphas;
tvec = [1; t; t^2; t^3];
% p = a0 + a1*t + a1*t^2 + a1*t^3
% pdot = a1 + 2*a2*t + 3*a3*t^2 
coeffs = [-pi/6 0  pi/50 -pi/375;
           pi/4 0 -pi/25 2*pi/375;
          -pi/3 0 7*pi/100 -7*pi/750];
alphas = coeffs*tvec

%reusing code to plot from Homework 9.3...

%plot a: joint angles
tlong = linspace(0,5,100);
alphalong = zeros(length(tlong), 3);
for i = 1:length(tlong)
    t = tlong(i);
    alphalong(i,1) = [-pi/6 0  pi/50 -pi/375] * [1; t; t^2; t^3];
    alphalong(i,2) = [pi/4 0 -pi/25 2*pi/375] * [1; t; t^2; t^3];
    alphalong(i,3) = [-pi/3 0 7*pi/100 -7*pi/750] * [1; t; t^2; t^3];
end
figure
plot(tlong, alphalong);
title('Joint Angles');
legend('alpha1', 'alpha2', 'alpha3');

%plot b: end effector configuration
e = zeros(length(tlong), 3);
for i = 1:length(tlong)
    a1 = alphalong(i, 1); a2 = alphalong(i,2); a3 = alphalong(i,3);
    e(i,:) = [l1*cos(a1) + l2*cos(a1+a2) + l3*cos(a1+a2+a3);
      l1*sin(a1) + l2*sin(a1+a2) + l3*sin(a1+a2+a3);
      a1+a1+a3];
end
figure
plot(tlong, e);
title('End Effector Configuration');
legend('x', 'y', 'theta');

%plot c: parametric plot
 %use display functions from class wiki to plot manipulator:
figure
hold on
planarR3_display(alphalong(1,:), l); %initial position
planarR3_display(alphalong(20,:), l);
planarR3_display(alphalong(40,:), l);
planarR3_display(alphalong(60,:), l); 
planarR3_display(alphalong(80,:), l); 
planarR3_display(alphalong(100,:), l); %final position
title('Parametric Plot');
