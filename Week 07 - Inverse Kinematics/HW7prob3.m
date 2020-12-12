% ECE 4560 - Homework 7, Problem 3
% Caitlyn Caggia

l1 = 4.5; l2 = 4.0;

%Only alpha2 (shoulder) and alpha3 (elbow) contribute to translation
%Total rotation = theta from R matrix ...
% = alpha2 (shoulder) + alpha3 (elbow) + alpha4 (wrist)

%want alphas to be in degrees, not radians

%first configuration
config1 = SE2([3.165; -7.811], atan2(-0.707, 0.707));
T1 = getTranslation(config1);
r1 = sqrt(T1(1)^2 + T1(2)^2);
c1a2 = acos( (l1^2 + r1^2 - l2^2) / (2*l1*r1) ) + atan2(T1(2), T1(1));
c1a3 = acos( (l1^2 + l2^2 - r1^2) / (2*l1*l2) ) - pi;
c1a4 = getTheta(config1) - c1a2 - c1a3;
disp('calculate angles for configuration 1:')
config1alphas = [0.5 rad2deg(c1a2) rad2deg(c1a3) rad2deg(c1a4) 0.5]

%second configuration
config2 = SE2([7.328; 2.828], atan2(0.966, 0.259));
T2 = getTranslation(config2);
r2 = sqrt(T2(1)^2 + T2(2)^2);
c2a2 = acos( (l1^2 + r2^2 - l2^2) / (2*l1*r2) ) + atan2(T2(2), T2(1));
c2a3 = acos( (l1^2 + l2^2 - r2^2) / (2*l1*l2) ) - pi;
c2a4 = getTheta(config2) - c2a2 - c2a3;
disp('calculate angles for configuration 2:')
config2alphas = [0.25 rad2deg(c2a2) rad2deg(c2a3) rad2deg(c2a4) 0.5]


%check calculated alphas with forward kinematics
disp('verify angles for configuration 1:')
T1check = [(l1*cos(c1a2) + l2*cos(c1a2 + c1a3));
    (l1*sin(c1a2) + l2*sin(c1a2 + c1a3))];
g1 = SE2(T1check, c1a2 + c1a3 + c1a4)

disp('verify angles for configuration 2:')
T2check = [(l1*cos(c2a2) + l2*cos(c2a2 + c2a3));
    (l1*sin(c2a2) + l2*sin(c2a2 + c2a3))];
g2 = SE2(T2check, c2a2 + c2a3 + c2a4)

