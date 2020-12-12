% ECE 4560 - Homework 12.1
% Caitlyn Caggia


%part a 
syms l1 l2 l3 a1 a2 a3

g1 = [R(a1) [0;0]; 0 0 1];
g2 = [R(a2) [l1;0]; 0 0 1];
g3  =[R(a3) [l2;0]; 0 0 1];
g4 = [eye(2) [l3;0]; 0 0 1];

J1 = [0; 0; 1];
J2 = [0; 0; 1];
J3 = [0; 0; 1];

ad1 = [R(-a2-a3) [l2*sin(a3) + l1*sin(a2+a3); l3+l2*cos(a3)+l1*cos(a2+a3)]; 0 0 1]
ad2 = [R(-a3) [l2*sin(a3); l3+l2*cos(a3)]; 0 0 1];
ad3 = [eye(2) [0;l3]; 0 0 1];   

Jb1 = ad1*J1;
Jb2 = ad2*J2;
Jb3 = ad3*J3;

Jbody = [Jb1 Jb2 Jb3]

%part b
ad1 = [R(a1) [0; 0]; 0 0 1];
ad2 = [R(a1+a2) [l1*sin(a1); l1*cos(a1)]; 0 0 1];
ad3 = [R(a1+a2+a3) [l1*sin(a1)+l2*sin(a1+a2); -l1*cos(a1)-l2*cos(a1+a2)]; 0 0 1];

Js1 = ad1*J1;
Js2 = ad2*J2;
Js3 = ad3*J3;

Jspatial = [Js1 Js2 Js3]

%part c
l1 = 1; l2 = 2/3; l3 = 1/2; 
a1 = pi/3; a2 = -pi/4; a3 = pi/6;
adot = [-pi/12; pi/8; -pi/12];

ad1 = [R(-a2-a3) [l2*sin(a3) + l1*sin(a2+a3); l3+l2*cos(a3)+l1*cos(a2+a3)]; 0 0 1];
ad2 = [R(-a3) [l2*sin(a3); l3+l2*cos(a3)]; 0 0 1];
ad3 = [eye(2) [0;l3]; 0 0 1];   

Jb1 = ad1*J1;
Jb2 = ad2*J2;
Jb3 = ad3*J3;

Jbody = [Jb1 Jb2 Jb3];

vbody = Jbody * adot

%part d
ad1 = [R(a1) [0; 0]; 0 0 1];
ad2 = [R(a1+a2) [l1*sin(a1); l1*cos(a1)]; 0 0 1];
ad3 = [R(a1+a2+a3) [l1*sin(a1)+l2*sin(a1+a2); -l1*cos(a1)-l2*cos(a1+a2)]; 0 0 1];

Js1 = ad1*J1;
Js2 = ad2*J2;
Js3 = ad3*J3;

Jspatial = [Js1 Js2 Js3];

vspatial = Jspatial*adot

%part e
h = [0.25 0 0];
adj = [eye(2) [0; -0.25]; 0 0 1];
vtool = adj * vbody





