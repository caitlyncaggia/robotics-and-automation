figure;
hold on;

O = SE2();
O.plot('O', '-k');

A = SE2([5,12], pi/3);
A.plot('A', '-b');

B = SE2([2,-1], pi);
B.plot('B', '-g');

%from HW1prob6
C = SE2([5.7321, 9.2679], pi/2);
C.plot('C', '-m');

%from HW1prob7, after nudge
A2 = SE2([cosd(50) + 5, sind(50) + 12], deg2rad(50));
A2.plot('A2', '-r');

%from HW1prob4
qax = 5.5; qay = 12.866;

%from HW1prob5
qbx = 1; qby = -1;

%from HW1prob6
qcx = 5.7321; qcy = 10.2679; 

plot(qax, qay, 'ob', qbx, qby, 'og', qcx, qcy, 'om');