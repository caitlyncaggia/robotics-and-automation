%ECE 4560 - Homework 11.3
%Caitlyn Caggia

%link lengths from lynx6.m from class wiki:
mm2in   = 1/25.4;
linklen = [110.6 120 120 130 20]*mm2in;
l0 = linklen(1);
l1 = linklen(2);
l2 = linklen(3);
l3 = linklen(4);
l4 = linklen(5);

%choose alphas from HW10.3
a1 = -pi/6;
a2 = -pi/3;
a3 = pi/6;
a4 = -pi/2;
a5 = pi/3;
a6 = 0;

%syms a1 a2 a3 a4 a5;

%forward kinematics as a product of Lie groups from HW8.3:
disp('g using product of lie groups')
g1 = SE3([0;0;l0], SE3.RotZ(a1));
g2 = SE3([0;0;0],SE3.RotX(a2));
g3 = SE3([0;0;l1], SE3.RotX(a3));
g4 = SE3([0;0;l2],SE3.RotX(a4));
g5 = SE3([0;0;l3],SE3.RotZ(a5));
g6 = SE3([0;0;l4],eye(3));
glie = g1*g2*g3*g4*g5*g6

%forward kinematics as a product of exponentials:
disp('g using product of exponentials')
g0 = SE3([0;0;l0+l1+l2+l3+l4], eye(3));
w1 = [0;0;1]; q1 = [0;0;l0];
z1 = [cross(q1,w1); w1]; exp1 = SE3.exp(z1*a1);
w2 = [1;0;0]; q2 = [0;0;l0]; 
z2 = [cross(q2,w2); w2]; exp2 = SE3.exp(z2*a2);
w3 = [1;0;0]; q3 = [0;0;l0+l1]; 
z3 = [cross(q3,w3); w3]; exp3 = SE3.exp(z3*a3);
w4 = [1;0;0]; q4 = [0;0;l0+l1+l2]; 
z4 = [cross(q4,w4); w4]; exp4 = SE3.exp(z4*a4);
w5 = [0;0;1]; q5 = [0;0;l0+l1+l2+l3]; 
z5 = [cross(q5,w5); w5]; exp5 = SE3.exp(z5*a5);
w6 = [0;0;0]; q6 = [0;0;l0+l1+l2+l3+l4]; 
z6 = [cross(q6,w6); w6]; exp6 = SE3.exp(z6*a6);
gexp = exp1 * exp2 * exp3 * exp4 * exp5 * exp6 * g0
