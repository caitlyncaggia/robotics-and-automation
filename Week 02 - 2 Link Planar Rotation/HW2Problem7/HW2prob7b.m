tspan = linspace(0,0.5,40);
a01 = -3*pi/2;
a02 = pi/6;

%adot1 = (1/3)*cos(t); - value calculated by @HW2prob7a1
%adot2 = (-1/4)*sin(t); - value calculated by @HW2prob7a2

[t1, a1] = ode45(@HW2prob7a1, tspan, a01);
[t2, a2] = ode45(@HW2prob7a2, tspan, a02);

x = 1*cos(a1) + 0.5*cos(a2);
y = 1*sin(a1) + 0.5*sin(a2);

%x plot
figure;
plot(t1, x);
xlabel('t');
ylabel('x');
grid on;

%y plot
figure;
plot(t2, y);
xlabel('t');
ylabel('y');
grid on;