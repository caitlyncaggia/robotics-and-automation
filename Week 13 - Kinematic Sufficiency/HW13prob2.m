%ECE 4560 - Homework 13.2
%Caitlyn Caggia

%part a ------------------------------------------------------------------
syms a1 a2 a3 a4 l1 l2
ge = forwardkin([a1 a2 a3 a4], [l1 l2]);
bodyJac = bodyJacobian([a1 a2 a3 a4], [l1 l2])

%part b ------------------------------------------------------------------
l1 = 1; l2 = 2/3;
alphai = [pi/4; 5*pi/6; -pi/4; 1.5];
a1 = alphai(1); a2 = alphai(2); a3 = alphai(3); a4 = alphai(4);
xi = [-0.415; -2.263; -1];
T = 2.618; time = linspace(0,T,10000); deltat = time(2) - time(1);

ab = zeros(4,length(time));
posb = zeros(3, length(time));
for i = 1:length(time)
    
    t = time(i);
    aold = [a1; a2; a3; a4];
    ab(:,i) = aold;
    posb(:,i) = forwardkin(aold, [l1 l2]);
    
    Jb = bodyJacobian(aold, [l1 l2]);
    Jp = Jb' * inv(Jb * Jb'); %m = 2, n = 4
    alphanew = aold + deltat*Jp*xi;
    a1 = alphanew(1); a2 = alphanew(2); a3 = alphanew(3); a4 = alphanew(4);
    
end

abfinal = alphanew
gbfinal = forwardkin(abfinal, [l1 l2])

%part c ------------------------------------------------------------------
a1 = alphai(1); a2 = alphai(2); a3 = alphai(3); a4 = alphai(4);
W = [100 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 20];
winv = inv(W);

ac = zeros(4,length(time));
posc = zeros(3,length(time));
for i = 1:length(time)
    
    t = time(i);
    aold = [a1; a2; a3; a4];
    ac(:,i) = aold;
    posc(:,i) = forwardkin(aold, [l1 l2]);
    
    Jb = bodyJacobian(aold, [l1 l2]);
    Jp = winv * Jb' * inv(Jb * winv * Jb'); 
    alphanew = aold + deltat*Jp*xi;
    a1 = alphanew(1); a2 = alphanew(2); a3 = alphanew(3); a4 = alphanew(4);
    
end

acfinal = alphanew
gcfinal = forwardkin(acfinal, [l1 l2])

%part d -----------------------------------------------------------------
figure
plot(time, ab)
legend('a1', 'a2', 'a3', 'a4')
title('Part B Alphas')

figure
plot(time, ac)
legend('a1', 'a2', 'a3', 'a4')
title('Part C Alphas')

figure
plot(time, posb)
legend('x', 'y', 'theta')
title('Part B End Effector Config')

figure
plot(time, posc)
legend('x','y','theta')
title('Part C End Effector Config')

figure
plot(posb(1,:), posb(2,:))
title('Part B Parametric Plot')

figure
plot(posc(1,:), posc(2,:))
title('Part C Parametric Plot')


%functions ---------------------------------------------------------------
function Jbody = bodyJacobian(alphas, len)

a1 = alphas(1); a2 = alphas(2); a3 = alphas(3); a4 = alphas(4);
l1 = len(1); l2 = len(2);

J1 = [0; 0; 1];
J2 = [0; 0; 1];
J3 = [0; 0; 1];
J4 = [1; 0; 0];

ad1 = [R(-a2-a3) [l2*sin(a3) + l1*sin(a2+a3); a4+l2*cos(a3)+l1*cos(a2+a3)]; 0 0 1];
ad2 = [R(-a3) [l2*sin(a3); a4+l2*cos(a3)]; 0 0 1];
ad3 = [eye(2) [0; a4]; 0 0 1];
ad4 = [eye(2) [0; 0]; 0 0 1];

Jb1 = ad1*J1;
Jb2 = ad2*J2;
Jb3 = ad3*J3;
Jb4 = ad4*J4;

Jbody = [Jb1 Jb2 Jb3 Jb4];

end

function ge = forwardkin(a, l)

l1 = l(1); l2 = l(2);
a1 = a(1); a2 = a(2); a3 = a(3); a4 = a(4);

ge = [l1*cos(a1) + l2*cos(a1+a2) + a4*cos(a1+a2+a3);
    l1*sin(a1) + l2*sin(a1+a2) + a4*sin(a1+a2+a3);
    a1+a2+a3];

end
