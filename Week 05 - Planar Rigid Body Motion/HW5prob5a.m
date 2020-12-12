function zdot1 = f(t, z)

theta = z(3);
xdot1 = 1*cos(theta) - (-1)*sin(theta);
ydot1 = 1*sin(theta) + (-1)*cos(theta);
thetadot1 = 1*1/4;

zdot1 = [xdot1; ydot1; thetadot1];


end