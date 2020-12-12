% ECE 4560 - Homework 12.2
% Caitlyn Caggia

l1 = 1; l2 = 0.5; l3 = 0.25;

%use values from HW 11.2 so new values can be compared to previous spline
ai = [-pi/6; pi/4; -pi/3];
af = [0; -pi/12; pi/4];
T = 5;

%part a ------------------------------------------------------------------
gi = forwardkin(ai);
gf = forwardkin(af);

Ti = getTranslation(gi);
Tf = getTranslation(gf);
xi = Ti(1);
yi = Ti(2);
xf = Tf(1);
yf = Tf(2);
thetai = getTheta(gi);
thetaf = getTheta(gf);

%create new plots
tlong = linspace(0,5,100);
pos = zeros(length(tlong), 3);
alphas = zeros(length(tlong), 3);
for i = 1:length(tlong)
    t = tlong(i);
    pos(i,1) = xi + t/T*(xf - xi);
    pos(i,2) = yi + t/T*(yf - yi);
    pos(i,3) = thetai + t/T*(thetaf - thetai);
    alphas(i,:) = inversekin([pos(i,1); pos(i,2); pos(i,3)]);
end
figure
plot(tlong, alphas)
title('Part A Joint Angles')
legend('alpha1', 'alpha2', 'alpha3')

figure
plot(tlong, pos)
title('Part A End Effector Configuration')
legend('x', 'y', 'theta')

%old spline plots
alphaold = zeros(length(tlong),3);
posold = zeros(length(tlong), 3);
for i = 1:length(tlong)
    t = tlong(i);
    alphaold(i,1) = [-pi/6 0  pi/50 -pi/375] * [1; t; t^2; t^3];
    alphaold(i,2) = [pi/4 0 -pi/25 2*pi/375] * [1; t; t^2; t^3];
    alphaold(i,3) = [-pi/3 0 7*pi/100 -7*pi/750] * [1; t; t^2; t^3];
    ge = forwardkin([alphaold(i,1); alphaold(i,2); alphaold(i,3)]);
    posold(i,1:2) = getTranslation(ge);
    posold(i,3) = getTheta(ge);
end
figure
plot(tlong, alphaold);
title('Joint Angles from Spline');
legend('alpha1', 'alpha2', 'alpha3');

figure 
plot(tlong, posold)
title('End Effector from Spline');
legend('x', 'y', 'theta')

%part b ------------------------------------------------------------------
z = SE2.log(inv(gi)*gf);
expi = gi * SE2.exp(z*0);
alphab = zeros(length(tlong), 3);
posb = zeros(length(tlong), 3);
for i = 1:length(tlong)
    t = tlong(i);
    exp = expi * SE2.exp(z*t);
    gexp = [getTranslation(exp); getTheta(exp)];
    posb(i,:) = gexp;
    alphab(i,:) = inversekin([gexp(1); gexp(2); gexp(3)]);
end
figure
plot(tlong, alphab)
title('Part B Joint Angles')
legend('alpha1', 'alpha2', 'alpha3')

figure
plot(tlong, posb)
title('Part B End Effector Configuration')
legend('x','y','theta')


%function to compute inverse kinematics ----------------------------------
function alphas = inversekin(ge)

l1 = 1; l2 = 0.5; l3 = 0.25;

geh = [R(ge(3)) [ge(1); ge(2)]; 0 0 1];
g4 = [l3*cos(ge(3)); l3*sin(ge(3)); ge(3)];
g4h = [R(ge(3)) [l3*cos(ge(3)); l3*sin(ge(3))]; 0 0 1];
gw = geh *inv(g4h);

xw = gw(1,3);
yw = gw(2,3);
thetae = ge(3);

gamma = atan2(yw, xw);
r = sqrt(xw^2 + yw^2);
delta = acos((l1^2 + l2^2 - r^2)/(2*l1*l2)) - pi;
alpha = acos((l1^2 + l2^2 - r^2)/(2*l1*l2));
beta = acos((l1^2 + r^2 - l2^2)/(2*l1*r));

a1 = gamma - beta;
a2 = pi - alpha;
a3 = thetae - a1 - a2;

% all angles should fall in the range [-pi, pi]
if (a3 > pi)
    a3 = a3 - 2*pi;
elseif (a3 < -pi)
    a3 = a3 + 2*pi;
end

alphas = [a1; a2; a3];

end

%function to compute forward kinematics ----------------------------------
function ge = forwardkin(a)

l1 = 1; l2 = 0.5; l3 = 0.25;
g1 = SE2([0,0],a(1));
g2 = SE2([l1,0],a(2));
g3 = SE2([l2,0],a(3));
g4 = SE2([l3,0],0);
ge = g1 * g2 * g3 * g4;

end
