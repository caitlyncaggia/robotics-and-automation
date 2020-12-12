

giR = [0.8660 0.5000 0.0000;
      0.5000 -0.866 0.0000;
      0.0000  0.000 -1.000];
giT = [-4.4079; 7.6348; 0.8110];
gi = SE3(giT, giR);
ai = inversekin(gi);

gfR = [0.7071 0.0000 -0.7071;
      0.7071 0.0000  0.7071;
      0.0000 -1.000  0.0000];
gfT = [-8.7393; 8.7393; 6.0836];
gf = SE3(gfT, gfR);
af = inversekin(gf);


  %----------------------------- inversekin ----------------------------
  %
  %  Compute the inverse kinematics of the manipulator.  Either to get
  %  the finger tips to the desired configuration or the middle of the
  %  grip region (second is default).  Also, solPick is an optional
  %  argument specifying which of the various solutions to choose from.
  %  At minimum, it should specify the upper or the lower solution.
  %  A second option would be to specify if it will be a turn around
  %  and reach backwards solution or a reach forwards.
  %
  %  Gripper width is ignored in these inverse kinematics, so only
  %  a 5x1 column vector is returned.
  %
  function alpha = inversekin(gdes, totip, solPick)
  
  if ( (nargin <= 2) || isempty(totip) )
    totip = false;
  end
  
  if (nargin < 4)
    solPick = [];       % Up to you to decide how to use.
  end
  
  %link lengths from lynx6.m from class wiki:
  mm2in   = 1/25.4;
  linklen = [110.6 120 120 130 20]*mm2in; 
  l0 = linklen(1);
  l1 = linklen(2);
  l2 = linklen(3);
  l3 = linklen(4);
  l4 = linklen(5);
  
  g6 = SE3([0;0;l3+l4],eye(3));
  gwdes = gdes*inv(g6);
  
  T = getTranslation(gwdes);
  xwdes = T(1);
  ywdes = T(2); 
  zwdes = T(3);
  
  %there are 2 solutions for a3...could be +/- acos
  a3 = acos((xwdes^2 + ywdes^2 + (zwdes - l0)^2 - l1^2 - l2^2)/(2*l1*l2));
  
  %there are 4 solutions for a2
  r = sqrt(xwdes^2 + ywdes^2);
  A1 = r + l2*sin(a3);
  B1 = -2*(l1 + l2*cos(a3));
  C1 = r - l2*sin(a3);
  w11 = (-B1 + sqrt(B1^2 - 4*A1*C1))/(2*A1);
  w12 = (-B1 - sqrt(B1^2 - 4*A1*C1))/(2*A1);
  A2 = r - l2*sin(a3);
  B2 = 2*(l1 + l2*cos(a3));
  C2 = r + l2*sin(a3);
  w21 = (-B2 + sqrt(B2^2 - 4*A2*C2))/(2*A2);
  w22 = (-B2 - sqrt(B2^2 - 4*A2*C2))/(2*A2);
  
  %check values...
  a211 = 2*atan(w11); a212 = 2*atan(w12);
  a221 = 2*atan(w21); a222 = 2*atan(w22);
  %a212 and a221 are valid for a3 and -a3
  
  a2 = a221;
  a1 = atan2(xwdes/(l1*sin(a2)+l2*sin(a2+a3)), -ywdes/(l1*sin(a2)+l2*sin(a2+a3)));
  
  %a4 and a5 are chosen using the Rh matrix
  g1 = SE3([0;0;l0], SE3.RotZ(a1));
  g2 = SE3([0;0;0], SE3.RotX(a2));
  g3 = SE3([0;0;l1],SE3.RotX(a3));
  g4bar = SE3([0;0;l2],eye(3));
  gw = g1*g2*g3*g4bar;
  ghdes = inv(gw)*gdes*inv(g6);
  Rhdes = getRotationMatrix(ghdes);
  a4 = atan2(-Rhdes(2,3), Rhdes(3,3));
  a5 = atan2(-Rhdes(1,2), Rhdes(1,1));
  
  alpha = rad2deg([a1; a2; a3; a4; a5]);
  
  end