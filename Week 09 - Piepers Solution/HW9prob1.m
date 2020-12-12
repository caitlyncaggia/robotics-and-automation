% ECE 4560 - Homework 9.1
% Caitlyn Caggia

%forward kinematics:
syms l1 l2 l3 a1 a2 a3;
ge = [l1*cos(a1) + l2*cos(a1+a2) + l3*cos(a1+a2+a3);
      l1*sin(a1) + l2*sin(a1+a2) + l3*sin(a1+a2+a3);
      a1+a1+a3];
gw = [l1*cos(a1) + l2*cos(a1+a2);
      l1*sin(a1) + l2*sin(a1+a2);
      a1+a2];
      
%PART A ================================================================
xw = gw(1);
yw = gw(2);
thetae = ge(3);

gamma = atan2(yw, xw);
r = sqrt(xw^2 + yw^2);
delta = acos((l1^2 + l2^2 - r^2)/(2*l1*l2)) - pi;
alpha = acos((l1^2 + l2^2 - r^2)/(2*l1*l2));
beta = acos((l1^2 + r^2 - l2^2)/(2*l1*r));

a1 = simplify(gamma + beta);
a2 = simplify(alpha - pi);
a3 = simplify(thetae - a1 - a2);
alphaparta = [a1; a2; a3]

%PART B =================================================================
l1 = 1; l2 = 0.5; l3 = 0.25;
ge = [1.5560; 0.7288; 0.7854];
geh = [R(ge(3)) [ge(1); ge(2)]; 0 0 1];
g4 = [l3*cos(ge(3)); l3*sin(ge(3)); ge(3)];
g4h = [R(ge(3)) [l3*cos(ge(3)); l3*sin(ge(3))]; 0 0 1];
gw = geh *inv(g4h);

xw = gw(1,3);
yw = gw(2,3);
thetae = ge(3);

gamma = atan2(yw, xw);
r = sqrt(xw^2 + yw^2);
delta = acos((l1^2 + l2^2 - r^2)/(2*l1*l2)) - pi;
alpha = acos((l1^2 + l2^2 - r^2)/(2*l1*l2));
beta = acos((l1^2 + r^2 - l2^2)/(2*l1*r));

a1 = gamma + beta;
a2 = alpha - pi;
a3 = thetae - a1 - a2;

% all angles should fall in the range [-pi, pi]
if (a3 > pi)
    a3 = a3 - 2*pi;
elseif (a3 < -pi)
    a3 = a3 + 2*pi;
end

alphapartb = [a1; a2; a3]
