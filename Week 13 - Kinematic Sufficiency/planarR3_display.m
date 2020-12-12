%================================ planarR4 ===============================
%
%  function planarR3_display(alpha, ll)
%
%  Plot the planar R3 manipulator given the joint configuration alpha
%  and the link lengths ll.  Also plots the base frame and the circular
%  boundary of the manipulator's reachable points on the plane.
%
%  The last link length can be zero, in which case, the gripper is
%  colocated with the last joint.  The display is a bit awkward since it
%  is tough to visualize such a thing in 2D.  The link lengths and
%  gripper length are optional.  They get set to default values.
%
%  INPUT:
%    alpha      - The three joint angles plus optional gripper open width.
%                   Open width can't be bigger than 3*gripper length.
%                   Default gripper open width is (2/3)*gripper length.
%    ll         - The three link lengths (last one may be zero).
%                   Optional, default is [1; 1/2; 1/4].
%    gl         - The gripper length.
%                   Optional, default is mean(ll)/6;
%
%  OUTPUT:
%     N/A. Displays the manipulator in current figure.
%
%
%  Example:
%
%  >> planarR3_display([pi/10; pi/12; -pi/3],  [1; 1; 1/2], 1/5);
%
%
%================================ planarR3 ===============================
