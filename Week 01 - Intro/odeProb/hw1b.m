%================================== hw1b =================================
%
%  function xdot = f(t, x)
%
%
%  Performs integration of the dynamics associated to Homework 1, 
%  Problem b.  It needs your help to get coded.
%
%  If Matlab complains about column versus row vector then transpose the
%  result, xdot, before returning.
%
%================================== hw1b =================================
function xdot = f(t, x)

u = 0.5*cos(t);
v = 0.3*sin(t);

xdot = zeros(3,1);			% Force to be a 3x1 (column) vector.
xdot(1) = cos(t) * u;				% Add in differential equations.
xdot(2) = sin(t) * u;
xdot(3) = v;

end
