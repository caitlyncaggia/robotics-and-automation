% ECE 4560 - Homework 6, Problem 5
% Caitlyn Caggia

%PART A ===================================================================
%extend SE2 leftact to work for point (2x1) or vector (3x1) velocity
ga = SE2([1 1], pi/4);
T = [-2; 3];
v = [9; -3; pi/8];
fprintf('Part A \n\n')

%verify with SE2 * homogeneous matrix translation
disp('Verify leftact works with SE2 * homogeneous matrix translation:')
parta1 = ga.leftact(T)

%verify with SE2 * velocity vector
disp('Verify leftact works with SE2 * velocity vector:')
parta2 = ga.leftact(v)


%PART B ===================================================================
%amend adjoint to operate on vectors
%show SE2 works with adjoint * vector-velocity problem : HW5prob3

gCB = SE2([0,3], -pi/6); %given in HW5.2
zBB = [2; -3; pi/9]; %given in HW5.2
gBC = inv(gCB);
disp('Part B')

zCC = adjoint(zBB, gBC)
zCCsoln = [2.3252; -2.127; 0.3491] %from HW5.3

%PART C ===================================================================
%show that (Adgxi)^ = Adg(xi)^ (in vector form) and vice versa with unhat
%(in homogeneous form)
g = SE2([4;5], pi/3); 
disp('Part C')

%verify hat
disp('verify (Adgxi)^ = Adg(xi)^:')
xiVec1 = [4; 6; pi/12];
xiHat1 = SE2.hat(xiVec1);
adj1 = adjoint(xiVec1, g);
left1 = SE2.hat(adj1)
right1 = adjoint(xiHat1, g)


%verify unhat
disp('verify (Adgxihat)v = Adg(xihat)v:')
xiHat2 = [0      -pi/8   7;
         pi/8   0       -3;
         0      0       0];
xiVec2 = SE2.unhat(xiHat2);
adj2 = adjoint(xiHat2, g);
left2 = SE2.unhat(adj2)
right2 = adjoint(xiVec2, g)


