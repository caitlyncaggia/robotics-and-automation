% ECE 4560 - Homework 11.1
% Caitlyn Caggia

l1 = 1; l2 = 0.5; l3 = 0.25;
l = [l1;l2;l3];
a = [pi/3; pi/4; pi/12];
adot = [pi/12; pi/8; 0];

%part a
syms l1 l2 l3 a1 a2 a3;
a = [a1; a2; a3];
ge = [l1*cos(a(1)) + l2*cos(a(1)+a(2)) + l3*cos(a(1)+a(2)+a(3));
      l1*sin(a(1)) + l2*sin(a(1)+a(2)) + l3*sin(a(1)+a(2)+a(3));
      a(1) + a(2) + a(3)];
Jac = jacobian(ge, a)

%part b
adot = [pi/12; pi/8; 0];
l1 = 1; l2 = 0.5; l3 = 0.25;
a1 = pi/3; a2 = pi/4; a3 = pi/12;
J11 = - l2*sin(a1+a2) - l1*sin(a1) - l3*sin(a1+a2+a3);
J12 = - l2*sin(a1+a2) - l3*sin(a1+a2+a3);
J13 = -l3*sin(a1+a2+a3);
J21 = l2*cos(a1 + a2) + l1*cos(a1) + l3*cos(a1 + a2 + a3);
J22 = l2*cos(a1 + a2) + l3*cos(a1 + a2 + a3);
J23 = l3*cos(a1 + a2 + a3);
J31 = 1;
J32 = 1; 
J33 = 1;
Jac = [J11 J12 J13; J21 J22 J23; J31 J32 J33];

v = Jac * adot
disp('theta (v(3)) is positive so CCW rotation')

%part c
ge = [l1*cos(a1) + l2*cos(a1+a2) + l3*cos(a1+a2+a3);
      l1*sin(a1) + l2*sin(a1+a2) + l3*sin(a1+a2+a3);
      a1 + a2 + a3];
gc = SE2([ge(1); ge(2)], ge(3));
planarR3_display([a1; a2; a3], [l1; l2; l3])
velocityPlot(gc, v)

