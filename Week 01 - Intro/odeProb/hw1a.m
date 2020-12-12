%================================== hw1a =================================
%
%  function xdot = f(t, x)
%
%
%  Performs integration of the dynamics associated to Homework 1.
%  It needs your help to incorporate the input signal.
%
%================================== hw1a =================================
function xdot = f(t, x)

u = 8*exp(-t);

xdot = -2.5*x + 0.75*u;

end
