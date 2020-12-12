%ECE 4560 - Homework 4, Problem 3
%Caitlyn Caggia

gH21 = SE2([-3.151; -3.906], -7*pi/12);
gH2O = SE2([3.243; 2.512], acos(-0.131));
gKH = SE2([3;-1], -pi/4);

%part a
gH1O = adjoint(gH2O, gH21)

%part b
gK2O = adjoint(gKH, gH2O)

%part c
gK1O = adjoint(gKH, gH1O);
gK2K1 = gK2O * inv(gK1O)
