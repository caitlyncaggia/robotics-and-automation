% ECE 4560 - Homework 8.2
% Caitlyn Caggia

%part a: forward kinematics
syms a1 a2 a3 a4 l0 l1 l2;
g1 = [SE3.RotZ(a1) [0; 0; l0]; 0 0 0 1];
g2 = [SE3.RotX(a2) [0; 0; 0]; 0 0 0 1];
g3 = [SE3.RotZ(a3) [0; l1; 0]; 0 0 0 1];
g4 = [eye(3) [0; l2; 0]; 0 0 0 1];
geparta = g1*g2*g3*g4

%part b: visualization
l0 = 1; l1 = 0.75; l2 = 0.5;
alpha1 = (-90:5:90) * pi()/180;
alpha2 = (-45:5:60) * pi()/180;
alpha3 = (-90:5:90) * pi()/180;
alpha1Cnt = numel(alpha1);
alpha2Cnt = numel(alpha2);
alpha3Cnt = numel(alpha3);
total = alpha1Cnt*alpha2Cnt*alpha3Cnt;
x = zeros(1,total);
y = zeros(1,total);
z = zeros(1,total);
ll = 1;

for ii=1:alpha1Cnt
    for jj=1:alpha2Cnt
        for kk=1:alpha3Cnt
            x(ll) = (l1+l2*cos(alpha3(kk)))*sin(alpha1(ii))*cos(alpha2(jj)) ...
                -l2*cos(alpha1(ii))*sin(alpha3(kk));
            y(ll) = (l1+l2*cos(alpha3(ii)))*cos(alpha1(ii))*cos(alpha2(jj)) ...
                -l2*sin(alpha1(ii))*sin(alpha3(kk));
            z(ll) = l0 + (l1+l2*cos(alpha3(kk)))*sin(alpha2(jj));
            ll = ll+1;
        end
    end
end

x = reshape(x, [alpha2Cnt*alpha3Cnt, alpha1Cnt]);
y = reshape(y, [alpha2Cnt*alpha3Cnt, alpha1Cnt]);
z = reshape(z, [alpha2Cnt*alpha3Cnt, alpha1Cnt]);
surf(x,y,z);

%part c: end-effector configuration
a = [pi/3; pi/3; -pi/4];
g1 = [SE3.RotZ(a(1)) [0; 0; l0]; 0 0 0 1];
g2 = [SE3.RotX(a(2)) [0; 0; 0]; 0 0 0 1];
g3 = [SE3.RotZ(a(3)) [0; l1; 0]; 0 0 0 1];
g4 = [eye(3) [0; l2; 0]; 0 0 0 1];
gepartc = g1*g2*g3*g4



