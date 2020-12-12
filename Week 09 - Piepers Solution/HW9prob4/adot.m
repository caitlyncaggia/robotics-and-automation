function alphadot = adot(t, alpha)

% Make sure to convert angular rates from radians per second
% to degrees per second.  Use diag function to help you.

%since x and y are controlled by alphas 2-4 we can reuse code from HW9.3
a1 = alpha(2); a2 = alpha(3); a3 = alpha(4); %alphas in deg
l1 = 1; l2 = 0.5; l3 = 0.25;
vdes = [-0.15; 0.7; -0.25];

J11 = -l1*sind(a1) - l2*sind(a1+a2) - l3*sind(a1+a2+a3);
J12 = -l2*sind(a1+a2) - l3*sind(a1+a2+a3);
J13 = -l3*sind(a1+a2+a3);
J21 = l1*cosd(a1) + l2*cosd(a1+a2) + l3*cosd(a1+a2+a3);
J22 = l2*cosd(a1+a2) + l3*cosd(a1+a2+a3);
J23 = l3*cosd(a1+a2+a3);

J = [J11 J12 J13; J21 J22 J23];
        
alpharad = pinv(J)*vdes(1:2);
%convert from rad back to deg
alphadeg = 180/pi * alpharad;

%z is controlled only by alpha 1
a1init = 1.5;
a1final = a1init + a1init*vdes(3)*3; %use constant velocity to find final pos
adot1 = (a1final - a1init) / 3; %rate of change is (final-initial)/time

%alpha 5 is gripper width and not considered here

alphadot = [adot1; alpharad; alpha(5)];

end