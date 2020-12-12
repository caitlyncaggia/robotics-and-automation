% ECE 4560 - Homework 9.3
% Caitlyn Caggia

tspan = [0 2];
alpha0 = [-pi/4; 2*pi/3; -pi/5];
l1 = 1; l2 = 0.5; l3 = 0.35;
l = [l1; l2; l3];

%part a: joint angles
[t, alphas] = ode45(@HW9prob3func, tspan, alpha0);
figure
plot(t, alphas);
title('Problem 3a: Joint Angles');
legend('alpha1', 'alpha2', 'alpha3');

%part b: end effector configuration
e = zeros(length(alphas(:,1)), 3);
for i = 1:length(alphas(:,1))
    a1 = alphas(i, 1); a2 = alphas(i,2); a3 = alphas(i,3);
    e(i,:) = [l1*cos(a1) + l2*cos(a1+a2) + l3*cos(a1+a2+a3);
      l1*sin(a1) + l2*sin(a1+a2) + l3*sin(a1+a2+a3);
      a1+a1+a3];
end
figure
plot(t, e);
title('Problem 3b: End Effector Configuration');
legend('x', 'y', 'theta');

%part c: parametric plot
figure
 %plot path of end effector:
plot(e(:,1),e(:,2), 'Linewidth', 2, 'Color', [0.8 0 0]);
title('Problem 3c: Parametric Plot');
 %use display functions from class wiki to plot manipulator:
hold on
planarR3_display(alphas(1,:), l); %initial position
planarR3_display(alphas(14,:), l);
planarR3_display(alphas(28,:), l);
planarR3_display(alphas(41,:), l); %final position


