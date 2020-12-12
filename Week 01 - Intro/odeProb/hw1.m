%============================== Homework #1 ==============================
%
%
%  script hw1.m
%
%
%  Script to perform the work associated with Homework #1.  There is 
%  enough to get the code running for the sample problem.  It is up to
%  you to generate the remaining code.
%
%============================== Homework #1 ==============================


doPart = 'b';			% Set to 'a' or 'b' as desired.

switch doPart

  case 'a'
    tspan = [0, 40];
    x0 = 1;
    [t, x] = ode45( @hw1a , tspan, x0);

    figure(1);
    plot(t, x);
    xlabel('t');
    ylabel('x');
    grid on;

  case 'b'
    tspan = [0, 40];
    x0 = zeros(3,1);			% Should be a column vector.
    [t, x] = ode45( @hw1b, tspan, x0 );

    figure(1);
    plot(x(:,1), x(:,2));
    xlabel('x');
    ylabel('y');
    grid on;

    figure(2);
    plot(t, x(:,3) * 180/pi);	% Factor is to convert to degrees from radians.
    xlabel('t');
    ylabel('\theta');
end

