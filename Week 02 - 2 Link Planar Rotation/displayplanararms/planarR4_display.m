%================================ planarR4 ===============================
%
%  function planarR4_display(alpha, ll)
%
%  Plot the planar R4 manipulator given the joint configuration alpha,
%  the link lengths ll, and the gripper size, gl.  Also plots the base
%  frame and the circular boundary of the manipulator's reachable points
%  on the plane.  An optional 5th joint angle will plot the gripper open
%  to the amount specified (no checks on reasonableness of gripper
%  open width).
%
%  The last link length can be zero, in which case, the gripper is
%  colocated with the last joint.  The display is a bit awkward since it
%  is tough to visualize such a thing in 2D.  The link lengths and
%  gripper length are optional.  They get set to default values.
%
%  INPUT:
%    alpha      - The four joint angles plus optional gripper open width.
%                   Open width can't be bigger than 3*gripper length.
%                   Default gripper open width is (2/3)*gripper length.
%    ll         - The four link lengths (last one may be zero).
%                   Optional, default is [1; 1; 1/2; 1/4].
%    gl         - The gripper length.
%                   Optional, default is mean(ll)/5;
%
%  OUTPUT:
%     N/A. Displays the manipulator in current figure.
%
%
%  EXAMPLE:
%
%  >> planarR4_display([pi/10; pi/12; -pi/3; pi/12],  [1; 1; 1/2; 1/2], 1/4);
%
%
%================================ planarR4 ===============================
