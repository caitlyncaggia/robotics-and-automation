function alphadot = HWprob3func(t, alpha)

a1 = alpha(1); a2 = alpha(2); a3 = alpha(3);
l1 = 1; l2 = 0.5; l3 = 0.25;

J11 = -l1*sin(a1) - l2*sin(a1+a2) - l3*sin(a1+a2+a3);
J12 = -l2*sin(a1+a2) - l3*sin(a1+a2+a3);
J13 = -l3*sin(a1+a2+a3);
J21 = l1*cos(a1) + l2*cos(a1+a2) + l3*cos(a1+a2+a3);
J22 = l2*cos(a1+a2) + l3*cos(a1+a2+a3);
J23 = l3*cos(a1+a2+a3);

J = [J11 J12 J13; J21 J22 J23];
alphadot = pinv(J)*[-0.5;0.5];

end