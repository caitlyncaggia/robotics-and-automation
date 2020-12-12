tspan = [0, 40];
a01 = -3*pi/2;
a02 = pi/6;

[t1, a1] = ode45(@HW2prob7a1, tspan, a01);
[t2, a2] = ode45(@HW2prob7a2, tspan, a02);

%alpha 1 plot
figure;
plot(t1, a1);
xlabel('t1');
ylabel('a1');
grid on;

%alpha 2 plot
figure;
plot(t2, a2);
xlabel('t2');
ylabel('a2');
grid on;