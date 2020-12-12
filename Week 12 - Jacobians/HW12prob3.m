% ECE 4560 - Homework 12.3
% Caitlyn Caggia

l0 = 4.35; l1 = 4.725; l2 = l1; l3 = 5.12;
syms a1 a2 a3 a4 a5

%forward kinematics
g1 = SE3([0;0;l0], SE3.RotZ(a1));
g2 = SE3([0;0;0],SE3.RotX(a2));
g3 = SE3([0;0;l1], SE3.RotX(a3));
g4 = SE3([0;0;l2],SE3.RotX(a4));
g5 = SE3([0;0;0],SE3.RotZ(a5));
g6 = SE3([0;0;l3],eye(3));
ge = g1*g2*g3*g4*g5*g6;

%individual Jacobians
J1 = [0;0;0;0;0;1]
J2 = [0;0;0;1;0;0]
J3 = [0;0;0;1;0;0]
J4 = [0;0;0;1;0;0]
J5 = [0;0;0;0;0;1]

%inverted adjoints
g23456 = g2*g3*g4*g5*g6;
R1 = getRotationMatrix(g23456);
T1 = getTranslation(g23456);
ad1 = [transpose(R1) -transpose(R1)*SE3.hat(T1); zeros(3) transpose(R1)];
g3456 = g3*g4*g5*g6;
R2 = getRotationMatrix(g3456);
T2 = getTranslation(g3456);
ad2 = [transpose(R2) -transpose(R2)*SE3.hat(T2); zeros(3) transpose(R2)];
g456 = g4*g5*g6;
R3 = getRotationMatrix(g456);
T3 = getTranslation(g456);
ad3 = [transpose(R3) -transpose(R3)*SE3.hat(T3); zeros(3) transpose(R3)];
g56 = g5*g6;
R4 = getRotationMatrix(g56);
T4 = getTranslation(g56);
ad4 = [transpose(R4) -transpose(R4)*SE3.hat(T4); zeros(3) transpose(R4)];
R5 = getRotationMatrix(g6);
T5 = getTranslation(g6);
ad5 = [transpose(R5) -transpose(R5)*SE3.hat(T5); zeros(3) transpose(R5)];


%part a: body Jacobian ---------------------------------------------------
Jb1 = simplify(ad1 * J1);
Jb2 = simplify(ad2 * J2);
Jb3 = simplify(ad3 * J3);
Jb4 = simplify(ad4 * J4);
Jb5 = ad5 * J5;
Jbody = [Jb1 Jb2 Jb3 Jb4 Jb5]

%part b ------------------------------------------------------------------
a = [-pi/3; -pi/6; pi/4; -pi/2; 0];
adot = [pi/10; -pi/20; pi/20; 0; -pi/16];

%forward kinematics
g1 = SE3([0;0;l0], SE3.RotZ(a(1)));
g2 = SE3([0;0;0],SE3.RotX(a(2)));
g3 = SE3([0;0;l1], SE3.RotX(a(3)));
g4 = SE3([0;0;l2],SE3.RotX(a(4)));
g5 = SE3([0;0;0],SE3.RotZ(a(5)));
g6 = SE3([0;0;l3],eye(3));
ge = g1*g2*g3*g4*g5*g6;

%inverted adjoints
g23456 = g2*g3*g4*g5*g6;
R1 = getRotationMatrix(g23456);
T1 = getTranslation(g23456);
ad1 = [transpose(R1) -transpose(R1)*SE3.hat(T1); zeros(3) transpose(R1)];
g3456 = g3*g4*g5*g6;
R2 = getRotationMatrix(g3456);
T2 = getTranslation(g3456);
ad2 = [transpose(R2) -transpose(R2)*SE3.hat(T2); zeros(3) transpose(R2)];
g456 = g4*g5*g6;
R3 = getRotationMatrix(g456);
T3 = getTranslation(g456);
ad3 = [transpose(R3) -transpose(R3)*SE3.hat(T3); zeros(3) transpose(R3)];
g56 = g5*g6;
R4 = getRotationMatrix(g56);
T4 = getTranslation(g56);
ad4 = [transpose(R4) -transpose(R4)*SE3.hat(T4); zeros(3) transpose(R4)];
R5 = getRotationMatrix(g6);
T5 = getTranslation(g6);
ad5 = [transpose(R5) -transpose(R5)*SE3.hat(T5); zeros(3) transpose(R5)];

%body jacobian
Jb1 = ad1 * J1;
Jb2 = ad2 * J2;
Jb3 = ad3 * J3;
Jb4 = ad4 * J4;
Jb5 = ad5 * J5;
Jbody = [Jb1 Jb2 Jb3 Jb4 Jb5];

%body velocity
vbody = Jbody*adot

%part c ------------------------------------------------------------------
Rc = getRotationMatrix(ge);
Tc = getTranslation(ge);
adj = [Rc SE3.hat(Tc)*Rc; zeros(3) Rc];

vspatial = adj*vbody
