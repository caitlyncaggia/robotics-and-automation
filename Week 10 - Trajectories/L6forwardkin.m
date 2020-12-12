% ECE 4560 - Homework 8.3
% Caitlyn Caggia

%syms a1 a2 a3 a4 a5 l0 l1 l2 l3 l4;
  %link lengths from lynx6.m from class wiki:
  mm2in   = 1/25.4;
  linklen = [110.6 120 120 130 20]*mm2in; 
  l0 = linklen(1);
  l1 = linklen(2);
  l2 = linklen(3);
  l3 = linklen(4);
  l4 = linklen(5);
  
  alphas = deg2rad([0 -90 0 0 0 0]);
  
  a1 = alphas(1);
  a2 = alphas(2);
  a3 = alphas(3);
  a4 = alphas(4);
  a5 = alphas(5);
  a6 = alphas(6);

%translations
T1 = [0; 0; l0];
T2 = [0; 0; 0];
T3 = [0; 0; l1];
T4 = [0; 0; l2];
T5 = [0; 0; l3];
T6 = [0; 0; l4];

%rotations
R1 = SE3.RotZ(a1);
R2 = SE3.RotX(a2);
R3 = SE3.RotX(a3);
R4 = SE3.RotX(a4);
R5 = SE3.RotZ(a5);
R6 = eye(3);

g1 = SE3(T1, R1);
g2 = SE3(T2, R2);
g3 = SE3(T3, R3);
g4 = SE3(T4, R4);
g5 = SE3(T5, R5);
g6 = SE3(T6, R6);
ge = g1*g2*g3*g4*g5*g6;

de = getTranslation(ge)