lynx6.closeSerialPorts();
myparms = load('lynxparms.mat');

tspan = [0,8];
linklen = [4.35433 4.72441 4.72441 5.11811 0.787402];
ai = [-pi/6;-pi/3;pi/6;-pi/2;0];
af = [pi/6;0;-pi/12;pi/3;-pi/3];
a = ai;
l6 = lynx6(myparms);

gi = l6.forwardkin([ai; 0] * 180/pi, 1);
gf = l6.forwardkin([af; 0] * 180/pi, 1);
p1 = getTranslation(gi);
p2 = getTranslation(gf);
v = (p2 - p1);

traj = struct('pvec', a, 'pdot', v, 'tspan', tspan);
jtraj = l6.genPositionTrajectory(traj);
l6.followJointTrajectory(jtraj, tspan);


pTraj = [];
for i = 1:1:30
    g = l6.forwardkin(jtraj.alpha(:,i), 1);
    pTraj = [pTraj, getTranslation(g)];
end

figure
l6.displayArm(jtraj.alpha(:,1))
figure
l6.displayArm(jtraj.alpha(:,10))
figure
l6.displayArm(jtraj.alpha(:,20))
figure
l6.displayArm(jtraj.alpha(:,30))

figure
t = linspace(0,8,82);
plot(t, jtraj.alpha)
