% myparms = load('lynxparms.mat');
% l6 = lynx6(myparms);

ai = rad2deg([-pi/6; -pi/3; pi/6; -pi/2; 0]);
af = rad2deg([pi/6; 0; -pi/12; pi/3; -pi/3]);
linklen = [4.35433 4.72441 4.72441 5.11811 0.787402];

figure
lynx6.displayArm(ai)
figure
lynx6.displayArm(af)