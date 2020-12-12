% ECE 4560 - Homework 13.4
% Caitlyn Caggia

%given
linklen = [4.35433 4.72441 4.72441 5.11811 0.787402];
ai = rad2deg([-pi/6; -pi/3; pi/6; -pi/2; 0]);

%original end effector position
af = rad2deg([pi/6; 0; -pi/12; pi/3; -pi/3]);
figure
lynx6.displayArm(af)
title('Original End Effector Position')

%resolved rate end effector position
gi = forwardkin([ai; 0] * 180/pi);
gf = forwardkin([af; 0] * 180/pi);
T1 = getTranslation(gi);
T2 = getTranslation(gf);
xi = (T2 - T1);

a = ai;
tspan = [0 8];
pos = struct('pvec', a, 'pdot', xi, 'tspan', tspan);
aRR = genPositionTrajectory(pos, linklen);
afinalRR = aRR.alpha(:,end);
figure
lynx6.displayArm(afinalRR)
title('Resolved Rate End Effector Position')

%error
georiginal = forwardkin(af);
geresolvedrate = forwardkin(afinalRR);
disp('gerr = ')
gerr = inv(geresolvedrate) * georiginal
xierr = SE3.log(gerr)
xinorm = norm(xi)

% generate trajectory ----------------------------------------------------
function jtraj = genPositionTrajectory(ptraj, linklen)

d = .1/3;
v = ptraj.pdot;
a(1:1:5,1) = ptraj.pvec;
a(6,1) = 0;

d1 = [0;0;linklen(1)];
d2 = [0;0;linklen(2)];
d3 = [0;0;linklen(3)];
d4 = [0;0;linklen(4)];

for j = 0:.1:ptraj.tspan(2);
    i = uint64(10*j+1);
    
    R1 = SE3.RotZ(a(1,i)); R2 = SE3.RotX(a(2,i));
    R3 = SE3.RotX(a(3,i));
    R4 = SE3.RotX(a(4,i)); R5 = SE3.RotZ(a(5,i));
    
    g1 = SE3(d1,R1);
    g2 = SE3([0;0;0],R2);
    g3 = SE3(d2,R3);
    g4 = SE3(d3,R4);
    g5 = SE3([0;0;0],R5);
    g6 = SE3(d4,eye(3));
    
    ge = g1 * g2 * g3 * g4 * g5 * g6;
    Re = getRotationMatrix(ge);
    
    Jb = Jacobian(a(:,i), linklen);
    vb = inv(Re) * v;
    Jb1 = Jb(1:1:3,:);
    
    adot(:,i) = pinv(Jb1) * vb;
    a(1:1:5,i+1) = a(1:1:5,i) + adot(1:1:5,i) * d;
    a(6,i+1) = 0;
end

jtraj.alpha = a*180/pi;
jtraj.adot = adot*180/pi;
jtraj.time  = ptraj.tspan;
end

% compute body jacobian -------------------------------------------------
function JBody = Jacobian(alpha, linklen)

% Construction the Jacobian matrix.
a = alpha;
Jz = [0;0;0;0;0;1];
Jx = [0;0;0;1;0;0];

d1 = [0;0;linklen(1)];
d2 = [0;0;linklen(2)];
d3 = [0;0;linklen(3)];
d4 = [0;0;linklen(4)];

R1 = SE3.RotZ(a(1)); 
R2 = SE3.RotX(a(2)); 
R3 = SE3.RotX(a(3));
R4 = SE3.RotX(a(4)); 
R5 = SE3.RotZ(a(5));
g1 = SE3(d1,R1);
g2 = SE3([0;0;0],R2);
g3 = SE3(d2,R3);
g4 = SE3(d3,R4);
g5 = SE3([0;0;0],R5);
g6 = SE3(d4,eye(3));

JBody = [adjoint(inv(g2*g3*g4*g5*g6),Jz),adjoint(inv(g3*g4*g5*g6),Jx), adjoint(inv(g4*g5*g6),Jx),adjoint(inv(g5*g6),Jx),adjoint(inv(g6),Jz)];

end


%forward kinematics ------------------------------------------------------
function ge = forwardkin(alphas)

linklen = [4.35433 4.72441 4.72441 5.11811 0.787402];

d1 = [0;0;linklen(1)];
d2 = [0;0;linklen(2)];
d3 = [0;0;linklen(3)];
d4 = [0;0;linklen(4)];
d5 = [0;0;linklen(5)];

R1 = SE3.RotZ(alphas(1));
R2 = SE3.RotX(alphas(2));
R3 = SE3.RotX(alphas(3));
R4 = SE3.RotX(alphas(4));
R5 = SE3.RotZ(alphas(5));
g1 = SE3(d1,R1*R2);
g2 = SE3(d2,R3);
g3 = SE3(d3,R4);
g4 = SE3(d4,R5);
ge = g1*g2*g3*g4;

end
