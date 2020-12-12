% ECE 4560 - Homework 5, Problem 4
% Caitlyn Caggia

%part a
za = [4 6 pi/12];
J = [0 -1; 1 0];
zahat = [J*pi/12 [4;6]; 0 0 0]

%part b
zbhat = [0      -pi/8   7;
         pi/8   0       -3;
         0      0       0];
x = zbhat(1,3);
y = zbhat(2,3);
w = zbhat(2,1);
zb = [x y w]

%part c
zc = [1 -9 3 1 -0.8 -0.5];
x = zc(1); y = zc(2); z = zc(3);
w1 = zc(4); w2 = zc(5); w3 = zc(6);
what = [0 -w3 w2;
        w3 0 -w1;
        -w2 w1 0];
zchat = [what [x;y;z]; 0 0 0 0]

%part d
zdhat = [0   -0.9  -0.1    3;
         0.9  0    -0.2   -2;
         0.1  0.2   0      1;
         0    0     0      0];
x = zdhat(1,4); y = zdhat(2,4); z = zdhat(3,4);
w1 = zdhat(3,2); w2 = zdhat(1,3); w3 = zdhat(2,1);
zd = [x y z w1 w2 w3]
