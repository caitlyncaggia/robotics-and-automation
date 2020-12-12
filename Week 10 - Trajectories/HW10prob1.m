% ECE 4560 - Homework 10.1
% Caitlyn Caggia

% given link lengths
l0 = 0.5; l1 = 1; l2 = 1; l3 = 0.5;
 
%part a: inverse kinematics ===============================================
%forward kinematics...
syms a1 a2 a3 a4 a5 a6 w xw yw zw;
g1 = SE3([0;0;l0], SE3.RotZ(a1));
g2 = SE3([0;0;0],SE3.RotX(a2));
g3 = SE3([0;l1;0], SE3.RotX(a3));
g4 = SE3([0;l2;0],SE3.RotX(a4));
g5 = SE3([0;0;0],SE3.RotZ(a5));
g6 = SE3([0;0;0],SE3.RotY(a6));
g7 = SE3([0;l3;0],eye(3));
ge = g1*g2*g3*g4*g5*g6*g7;

%Pieper's Approach...
%xw = simplify(-sin(a1)*(l1*cos(a2) + l2*cos(a2+a3)));
%yw = simplify(cos(a1)*(l1*cos(a2) + l2*cos(a2+a3)));
%zw = simplify(l0 + l1*sin(a2)+l2*sin(a2+a3));

a1 = simplify(atan2(-xw, yw));
acosarg = (xw^2 + yw^2 + (zw - l0)^2 - l1^2 -l2^2)/(2*l1*l2);
a3 = simplify(acos(acosarg));

r = simplify(sqrt(xw^2 + yw^2));
A = r + l1 + l2*cos(a3);
B = 2*l2*sin(a3);
C = r - l1 - l2*cos(a3);
w = simplify(solve(A*w^2 + B*w + C, w));
a2 = simplify(2*atan(w));

Rh = simplify(getRotationMatrix(g4)*getRotationMatrix(g5)*getRotationMatrix(g6));
a5 = simplify(asin(-Rh(1,2)));

if (cos(a5) ~= 0)
    a6 = simplify(atan2(Rh(1,3), Rh(1,1)));
    a4 = simplify(atan2(Rh(3,2)/sqrt(Rh(1,1)^2 + Rh(1,3)^2),Rh(2,2)/sqrt(Rh(1,1)^2 + Rh(1,3)^2)));
else
    if (sin(a5) == 1)
        a4 = 0; a5 = pi/2; a6 = simplify(atan2(-R(3,1),R(3,3)));
    elseif (sin(a5) == -1)
        a4 = 0; a5 = -pi/2; a6 = simplify(atan2(-R(3,1),R(3,3)));
    end
end

alphas = [a1;a2;a3;a4;a5;a6]


%part b: calculations =====================================================
%initial configuration ----------------------------------------------------
giR = [0.9280 0.3536 0.1174;
      -0.3245 0.6124 0.7209;
      0.1830 -0.7071 0.6830];
giT = [1.01; 1.7551; 0.5947];
gi = SE3(giT, giR);

xwi = giT(1);
ywi = giT(2);
zwi = giT(3);
a1 = atan2(-xwi, ywi);
acosarg = (xwi^2 + ywi^2 + (zwi - l0)^2 - l1^2 -l2^2)/(2*l1*l2);
a3 = acos(acosarg);

r = sqrt(xwi^2 + ywi^2);
A = r + l1 + l2*cos(a3);
B = 2*l2*sin(a3);
C = r - l1 - l2*cos(a3);
w = (-B - sqrt(B^2 - 4*A*C)) / (2*A);
a2 = 2*atan(w);

a5 = asin(-giR(1,2));

if (cos(a5) ~= 0)
    a6 = atan2(giR(1,3), giR(1,1));
    a4 = atan2(giR(3,2)/sqrt(giR(1,1)^2 + giR(1,3)^2),giR(2,2)/sqrt(giR(1,1)^2 + giR(1,3)^2));
else
    if (sin(a5) == 1)
        a4 = 0; a5 = pi/2; a6 = atan2(-giR(3,1),giR(3,3));
    elseif (sin(a5) == -1)
        a4 = 0; a5 = -pi/2; a6 = atan2(-giR(3,1),giR(3,3));
    end
end

alphai = [a1;a2;a3;a4;a5;a6]

%final configuration ------------------------------------------------------
gfR = [0.2500 -0.6424 -0.7244;
      -0.4330  0.5950 -0.6771;
       0.8660  0.4830 -0.1294];
gfT = [-1.071; 1.5966; 1.6075];
gf = SE3(gfT, gfR);

xwf = gfT(1);
ywf = gfT(2);
zwf = gfT(3);
a1 = atan2(-xwf, ywf);
acosarg = (xwf^2 + ywf^2 + (zwf - l0)^2 - l1^2 -l2^2)/(2*l1*l2);
a3 = acos(acosarg);

r = sqrt(xwf^2 + ywf^2);
A = r + l1 + l2*cos(a3);
B = 2*l2*sin(a3);
C = r - l1 - l2*cos(a3);
w = (-B - sqrt(B^2 - 4*A*C)) / (2*A);
a2 = 2*atan(w);

a5 = asin(-gfR(1,2));

if (cos(a5) ~= 0)
    a6 = atan2(gfR(1,3), gfR(1,1));
    a4 = atan2(gfR(3,2)/sqrt(gfR(1,1)^2 + gfR(1,3)^2),gfR(2,2)/sqrt(gfR(1,1)^2 + gfR(1,3)^2));
else
    if (sin(a5) == 1)
        a4 = 0; a5 = pi/2; a6 = atan2(-gfR(3,1),gfR(3,3));
    elseif (sin(a5) == -1)
        a4 = 0; a5 = -pi/2; a6 = atan2(-gfR(3,1),gfR(3,3));
    end
end

alphaf = [a1;a2;a3;a4;a5;a6]
