% ECE 4560 - Homework 5, Problem 5
% Caitlyn Caggia

% gdot = g*z
% A = [cos(theta) -sin(theta) 0;
%      sin(theta)  cos(theta) 0;
%      0           0          1];

%Case 1
z = [1 -1 1/4];
g0 = [0 0.5 pi/6];
T0 = 0; Tfinal = pi;
[t1, z1] = ode45(@HW5prob5a, [T0 Tfinal], g0);
x1 = z1(:,1); y1 = z1(:,2); theta1 = z1(:,3);

%Case 2
z = [1 -1 0];
[t2, z2] = ode45(@HW5prob5b, [T0 Tfinal], g0);
x2 = z2(:,1); y2 = z2(:,2); theta2 = z2(:,3);


%Case 1 Plot
figure
hold on
title('Case 1 Velocity')
xlabel('x position')
ylabel('y position')
plot(x1, y1)

a1 = SE2([x1(10) y1(10)], theta1(10));
plot(a1,'','r')
b1 = SE2([x1(20) y1(20)], theta1(20));
plot(b1,'','y')
c1 = SE2([x1(30) y1(30)], theta1(30));
plot(c1,'','g')
d1 = SE2([x1(40) y1(40)], theta1(40));
plot(d1,'','b')

%Case 2 Plot
figure
hold on
title('Case 2 Velocity')
xlabel('x position')
ylabel('y position')
plot(x2, y2)

a2 = SE2([x2(10) y2(10)], theta2(10));
plot(a2,'','r')
b2 = SE2([x2(20) y2(20)], theta2(20));
plot(b2,'','y')
c2 = SE2([x2(30) y2(30)], theta2(30));
plot(c2,'','g')
d2 = SE2([x2(40) y2(40)], theta2(40));
plot(d2,'','b')


%DIFFERENCE BETWEEN CASE 1 AND 2
% Case 1 experiences changes in x, y, and theta while Case 2 only
% experiences changes in x and y. Thus in Case 1 the object rotates and
% translates, but in Case 2 the object only translates. This makes sense as
% change in theta (represented by the last element of z) is nonzero in Case
% 1, and 0 in Case 2.


