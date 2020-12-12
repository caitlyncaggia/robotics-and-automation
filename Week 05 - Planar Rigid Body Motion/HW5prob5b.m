function zdot2 = f(t, z)

theta = z(3);
xdot2 = 1*cos(theta) - (-1)*sin(theta);
ydot2 = 1*sin(theta) + (-1)*cos(theta);
thetadot2 = 1*0;

zdot2 = [xdot2; ydot2; thetadot2];


end