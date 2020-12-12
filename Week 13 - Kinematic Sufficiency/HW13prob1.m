% ECE 4560 - Homework 13.1
% Caitlyn Caggia

clear all; close all;

syms a1 a2 a3 l1 l2 l3 t
T = [l1*cos(a1) + l2*cos(a1+a2) + l3*cos(a1+a2+a3);
    l1*sin(a1) + l2*sin(a1+a2) + l3*sin(a1+a2+a3)];
Jac = jacobian(T, [a1 a2 a3]);
x = 1 - 0.5*cos(t);
xdot = 0.5*sin(t);
y = -0.75 + 0.49*sin(t);
ydot = 0.49*cos(t);

%part a ------------------------------------------------------------------
ti = 0; tf = 2*pi; time = linspace(ti, tf, 50000);
l1 = 1; l2 = 0.5; l3 = 0.25;

ai = [-0.2987; -2.4189; 1.1468];
[alphaA, posA] = resolvedRate(ai, time, 2, 0);

%plot alphas vs time
figure
hold on
plot(time, alphaA)
plot(time, zeros(1, length(time))); %to find singularity graphically
title('Part A Joint Angles')
legend('alpha1', 'alpha2', 'alpha3', 'zero')


%plot position vs time
figure
plot(time, posA)
title('Part A End Effector Position')
legend('x','y')

%parametric plot
tend = length(time);
figure
hold on
planarR3_display([alphaA(1,1); alphaA(2,1); alphaA(3,1)], [l1; l2; l3])
planarR3_display([alphaA(1,0.2*tend); alphaA(2,0.2*tend); alphaA(3,0.2*tend)], [l1; l2; l3])
planarR3_display([alphaA(1,0.4*tend); alphaA(2,0.4*tend); alphaA(3,0.4*tend)], [l1; l2; l3])
planarR3_display([alphaA(1,0.6*tend); alphaA(2,0.6*tend); alphaA(3,0.6*tend)], [l1; l2; l3])
planarR3_display([alphaA(1,0.8*tend); alphaA(2,0.8*tend); alphaA(3,0.8*tend)], [l1; l2; l3])
planarR3_display([alphaA(1,tend); alphaA(2,tend); alphaA(3,tend)], [l1; l2; l3])
plot(posA(1,:), posA(2,:)) %compare to desired path
title('Part A Parametric Plot')

%part b ------------------------------------------------------------------
[alphaB1, posB1] = resolvedRate(ai, time, 'd', 0.05);
%plot alphas vs time
figure
plot(time, alphaB1)
title('Part B Joint Angles, p^2 = 0.05')
legend('alpha1', 'alpha2', 'alpha3')
%plot position vs time
figure
plot(time, posB1)
title('Part B End Effector Position, p^2 = 0.05')
legend('x','y')
%parametric plot
figure
hold on
planarR3_display([alphaB1(1,1); alphaB1(2,1); alphaB1(3,1)], [l1; l2; l3])
planarR3_display([alphaB1(1,0.2*tend); alphaB1(2,0.2*tend); alphaB1(3,0.2*tend)], [l1; l2; l3])
planarR3_display([alphaB1(1,0.4*tend); alphaB1(2,0.4*tend); alphaB1(3,0.4*tend)], [l1; l2; l3])
planarR3_display([alphaB1(1,0.6*tend); alphaB1(2,0.6*tend); alphaB1(3,0.6*tend)], [l1; l2; l3])
planarR3_display([alphaB1(1,0.8*tend); alphaB1(2,0.8*tend); alphaB1(3,0.8*tend)], [l1; l2; l3])
planarR3_display([alphaB1(1,tend); alphaB1(2,tend); alphaB1(3,tend)], [l1; l2; l3])
plot(posB1(1,:), posB1(2,:)) %actual path
plot(posA(1,:), posA(2,:)) %compare to desired path
title('Part B Parametric Plot, p^2 = 0.05')

[alphaB2, posB2] = resolvedRate(ai, time, 'd', 0.005);
%plot alphas vs time
figure
plot(time, alphaB2)
title('Part B Joint Angles, p^2 = 0.005')
legend('alpha1', 'alpha2', 'alpha3')
%plot position vs time
figure
plot(time, posB2)
title('Part B End Effector Position, p^2 = 0.005')
legend('x','y')
%parametric plot
figure
hold on
planarR3_display([alphaB2(1,1); alphaB2(2,1); alphaB2(3,1)], [l1; l2; l3])
planarR3_display([alphaB2(1,0.2*tend); alphaB2(2,0.2*tend); alphaB2(3,0.2*tend)], [l1; l2; l3])
planarR3_display([alphaB2(1,0.4*tend); alphaB2(2,0.4*tend); alphaB2(3,0.4*tend)], [l1; l2; l3])
planarR3_display([alphaB2(1,0.6*tend); alphaB2(2,0.6*tend); alphaB2(3,0.6*tend)], [l1; l2; l3])
planarR3_display([alphaB2(1,0.8*tend); alphaB2(2,0.8*tend); alphaB2(3,0.8*tend)], [l1; l2; l3])
planarR3_display([alphaB2(1,tend); alphaB2(2,tend); alphaB2(3,tend)], [l1; l2; l3])
plot(posB2(1,:), posB2(2,:)) %actual path
plot(posA(1,:), posA(2,:)) %compare to desired path
title('Part B Parametric Plot, p^2 = 0.005')

%part c ------------------------------------------------------------------
finalposA = posA(:, end)
finalposB1 = posB1(:, end)
finalposB2 = posB2(:, end)

%extra credit ------------------------------------------------------------
disp('Singularity occurs when a2 and a3 are zero, as rank is lost.')
disp('Graphically, we can examine joint angles in Figure 1.')
disp('This shows we are closest to singularity just before 4 seconds.')
disp('More specifically, we get closest to singularity when a2 is at a maximum.')
disp('Thus, singularity occurs at...')
[value, index] = max(alphaA(2,:));
singularity_time = index*(time(2) - time(1))


% function to integrate over resolved rate -------------------------------
function [a, pos] = resolvedRate(ai, time, invtype, p)

a1 = ai(1); a2 = ai(2); a3 = ai(3);
a = zeros(3,length(time));
a(1,1) = a1; a(2,1) = a2; a(3,1) = a3;

pos = zeros(2,length(time));
pos(:,1) = forwardkin(a(:,1));

deltat = time(2) - time(1);

for i = 2:length(time)
    
    %     x = 1 - 0.5*cos(t); xdot = 0.5*sin(t);
    %     y = -0.75 + 0.49*sin(t); ydot = 0.49*cos(t);
    
    t = time(i);
    
    %store old values
    aold = [a1; a2; a3];
    a(:,i) = aold;
    pos(:,i) = forwardkin(aold);
    
    %calculate velocity
    xi = [0.5*sin(time(i)); 0.49*cos(time(i))];
    
    %calculate new alphas
    Jp = pseudoInv(aold, invtype, p);
    alphanew = aold + deltat*Jp*xi;
    a1 = alphanew(1); a2 = alphanew(2); a3 = alphanew(3);
    
end

end

% function to compute pseudoinverse --------------------------------------
function Jp = pseudoInv(alphas, type, p2)

l1 = 1; l2 = 0.5; l3 = 0.25;
a1 = alphas(1); a2 = alphas(2); a3 = alphas(3);

%calculate Jacobian
J11 = -l2*sin(a1+a2) - l1*sin(a1) - l3*sin(a1+a2+a3);
J12 = -l2*sin(a1+a2) - l3*sin(a1+a2+a3);
J13 = -l3*sin(a1+a2+a3);
J21 = l2*cos(a1+a2) + l1*cos(a1) + l3*cos(a1+a2+a3);
J22 = l2*cos(a1+a2) + l3*cos(a1+a2+a3);
J23 = l3*cos(a1+a2+a3);
Jac = [J11 J12 J13; J21 J22 J23];

% m = dimension of task space
% n = number of independent variables (joints)

%find pseudoinverse
if(type == 1)
    %true inverse m = n, kinematially sufficient
    Jp = inv(Jac);
elseif (type == 2)
    %pseudoinverse m < n, kinematically redundant
    Jp = Jac' * inv(Jac*Jac');
elseif (type == 3)
    %pseudoinverse m > n, kinematically insufficient
    Jp = inv(Jac' * Jac) * Jac';
elseif(type == 'd')
    %damped pseudo inverse
    Jp = Jac' * inv(Jac*Jac' + p2*eye(2));
end

end

%function to compute forward kinematics ----------------------------------
function T = forwardkin(a)

l1 = 1; l2 = 0.5; l3 = 0.25;
g1 = SE2([0,0],a(1));
g2 = SE2([l1,0],a(2));
g3 = SE2([l2,0],a(3));
g4 = SE2([l3,0],0);
ge = g1 * g2 * g3 * g4;
T = getTranslation(ge);

end

