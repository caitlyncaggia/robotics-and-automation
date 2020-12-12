%================================= lynx6 =================================
%
%
%
%
%================================= lynx6 =================================

%
%  Name:	    lynx6.m
%
%  Author:	    Patricio A. Vela, pvela@gatech.edu
%
%  Created:	    2007/08/09
%  Modified:	2014/03/26
%
%================================= lynx6 =================================
classdef lynx6 < handle


properties
  alpha;

  alphaIds    = [   0    1    2    3    4    5];		
  				% Servo ID numbers on the SSC-32 control board.
  alpha2musec = diag(1./[0.09 0.09 0.09 0.09 0.09 0.09/75]);	
  				% Converts degrees to microsecs for manipulator.
        				%  This may be only kinda accurate since it 
        				%  assumes the same swing for all servos, which
        				%  is not necessarily true.  
  				%  I don't really use this since the musec
  				%  commands are interpolated from the xxLims
  				%  variables.
        				%  Need to investigate.
  
  alphaHome   = [   0;   0;   0;   0;   0; 1.0];		
  				% Home position as a joint configuration.
  %alphaSleep  = [   0,  75, -120, -75,   0, 1.0];	% Folded up.
  alphaSleep  = [   0;  35;   -90;  -50;   0;  1.0];	% Leaning on block.
  				% Sleep position as a joint configuration.
  				%   This is what you put it to before powering 
        				%   down.
  
  alphaOrient = [ 1,   1,   -1,  -1,   -1,  -1];		
  				% To switch orientation in case servo rotates 
        				%   in the wrong direction 
        				%   (e.g. need right-hand rule).
  
  alphaLims = [ -90, -90, -45, -90, -90, 0.00;		
                  0,   0,   0,   0,   0, 3/4;
        	       90,  90,  90,  90,  80, 1.125];
  				% Limits of the joints, either angular values 
      			%   or linear values as in the gripper.
  				%   The middle value is some measured
  				%   intermediate location.
  musecLims = [ 525  640  500  580  790 1475;			
               1460 1520 868 1455 1580 1760;			
               2320 2400 1630 2340 2500 2500];
  				% How these limits translate to microseconds 
        				%   for the servo commands.
  
  linklen;
  serialid;
end

%)
%
%============================ Member Functions ===========================
%
%(

methods

  %------------------------------- lynx6 -------------------------------
  %
  %  Constructor
  %
  function this = lynx6(setup)

  mm2in   = 1/25.4;
  this.linklen = [4.35433 4.72441 4.72441 5.11811 0.787402];			
  				% Measured link lengths of the manipulator in
  				%   millimeters but converted to inches.

  % Serial port setup defaults:
  serport.com  = 'COM1';		% COM port is COM1.
  serport.baud = 115200;		% Baud rate is 115200 bps.
  serport.terminator = 'CR';	% Message terminator is Carriage Return.

  % Parse the setup argument now.
  if (nargin > 0)
    fnames = fieldnames(setup);
    for ii = 1:length(fnames)
      fnames{ii}
      switch fnames{ii}
        case {'linklen','alphaHome','alphaOrient'},
          eval(['this.' fnames{ii} ' = getfield(setup, fnames{ii});']);

          case 'musecLims',
          disp('Here 2');
          if (size(setup.musecLims,2) == 6)
  	        if (size(setup.musecLims, 1) == 2)
  	          this.musecLims(:,1) = setup.musecLims(:,1);
  	          this.musecLims(:,3) = setup.musecLims(:,2);
  	          this.musecLims(:,2) = mean(setup.musecLims,1)
  	        elseif (size(setup.musecLims, 1) == 3)
  	          this.musecLims = setup.musecLims;
  	        end
          end

        case 'alphaLims'
          disp('Here 3');
          if (size(setup.alphaLims,2) == 6)
            if (size(setup.alphaLims,1) == 2)
  	          this.alphaLims(:,1) = setup.alphaLims(:,1);
  	          this.alphaLims(:,3) = setup.alphaLims(:,2);
  	          this.alphaLims(:,2) = mean(setup.alphaLims,1)
  	        elseif (size(setup.alphaLims, 1) == 3)
  	          this.alphaLims = setup.alphaLims;
  	        end
          end

        case 'COM','baud','terminator',
          eval(['serport.' fnames{ii} ' = getfield(setup, fnames{ii});']);
      end
    end
  end
  this.musecLims
  this.alphaLims
  this.alphaOrient
  % Open up the serial port.
  if (ismac)
    this.serialid = serial(serport.com,'BaudRate',serport.baud, ...
                                       'Terminator',serport.terminator);
    fopen(this.serialid);

  elseif (ispc)
    this.serialid = serial(serport.com,'BaudRate',serport.baud, ...
                                       'Terminator',serport.terminator)
    fopen(this.serialid)
    % Note that in Windows, if the program crashes then the com port is
    % still open.  You will not be able to reopen the port and you won't
    % be able to access the serialid you had before.  The only option I
    % know of is to restart Matlab.  Or you could try to do an 'fclose
    % all' and see if that works.  I haven't tested it out.
    %
    %   Read documentation: 'help fclose'
  
  elseif isunix
    device = '/dev/ttyS0'
    this.serialid = fopen(device,'w+')
  end

  end

  %----------------------------- gotoSleep -----------------------------
  %
  %  The first and second to last member function to call.  This puts
  %  the manipulator into a rest position for shutdown.  At startup, it
  %  returns the manipulator to the rest position.  This is to prevent
  %  huge jerky movements as the motor goes to the first commanded
  %  position from wherever it is located at.
  %
  function gotoSleep(this, time)
  
  if (nargin == 1)
    time = 4;
  end
  
  this.setArm(this.alphaSleep, time);
  
  end

  %------------------------------ gotoHome -----------------------------
  %
  %  This is the default position after achieving the sleep position.
  %  The home position is typically where it will go when powered up,
  %  but not being used.
  %
  function gotoHome(this, time)
  
  if (nargin == 1)
    time = 4;
  end
  
  this.setArm(this.alphaHome, time);
  
  end
  
  
  %----------------------------- setServos -----------------------------
  %
  %  Send raw motor commands to the servomotors.
  %
  function setServos(this, pwmsigs, time, defaultLims)
  
  if ((nargin < 3) || isempty(time))
    time = 2;
  end
  
  if (nargin < 4)
    defaultLims = true;
  end
    
  if ( (size(pwmsigs,1) ~= 6) || (size(pwmsigs,2) ~= 1) )
    disp('Incorrect dimensions, ignoring command [expect: 6x1 column vector].');
    return;
  end
    
  if (defaultLims)
    ticks = min(this.musecLims(3,:),max(this.musecLims(1,:), pwmsigs'));
  else
    ticks = min(2500, max(500, pwmsigs'));
  end
  
  cmdstr = sprintf('#%d P%d ',[this.alphaIds ; ticks]);		
  cmdstr = [cmdstr sprintf(' T%d \n', round(time*1000))];
  
  %fprintf(cmdstr); 			% Print actual command to STDOUT. To debug.
  fprintf(this.serialid, cmdstr);
 
  end
  
  %------------------------------- setArm ------------------------------
  %
  %  Send the arm to some position, with an optional time duration of
  %  movement.  If the gripper width is not passed, then defaults to the
  %  home position gripper width.
  %
  function setArm(this, alpha, time)

  if (nargin == 2)
    time = 4;
  end

  if ( (size(alpha,1) == 5) && (size(alpha,2) == 1) )
    alpha(6) = alphaHome(6);
  elseif ( (size(alpha,1) ~= 6) || (size(alpha,2) ~= 1) )
    disp('Incorrect dimensions, ignoring command');
    disp('    [expect: 5x1 or 6x1 column vector].');
    return;
  end

  alpha = min(this.alphaLims(3,:),max(this.alphaLims(1,:), alpha'));
  for ii = 1:length(this.alphaIds)
    if (this.alphaOrient(ii) == 1)
      ticks(ii) = interp1(this.alphaLims(:,ii), this.musecLims(:,ii), ...
                                                                   alpha(ii));
    else
      ticks(ii) = interp1(this.alphaLims(end:-1:1,ii), this.musecLims(:,ii), ...
                                                                   alpha(ii));
    end
  end
  ticks = round(ticks);
  cmdstr = sprintf('#%d P%d ',[this.alphaIds ; ticks]);
  cmdstr = [cmdstr sprintf(' T%d \n', round(time*1000))];

  %fprintf(cmdstr); 			% Print actual command to STDOUT. To debug.
  fprintf(this.serialid, cmdstr);
  this.alpha = alpha;           % Keep track of last command.
  
  end

  %----------------------------- displayArm ----------------------------
  %
  %  Display the manipulator.  Joint angles should be in degrees since
  %  the are converted.  (Or remove the line giving the conversion and
  %  let it accept radians.  Up to you.)
  %
  %(
  function displayArm(this, alphadisp)

  if (nargin == 1)
    alphadisp = this.alpha;
  end
  linklenght = [4.35433 4.72441 4.72441 5.11811 0.787402];
  vparms.home = 'straight-up';
  alphadisp = diag([(pi/180)*ones([1 5]), 1])*alphadisp;	% Convert to rads.
  lynx6_display(alphadisp, linklenght, vparms);

  end


  %)
  %----------------------------- forwardkin ----------------------------
  %
  %  Compute the forward kinematics of the manipulators.  If the
  %  manipulator was calibrated and the measurements taken so that it is
  %  possible to select either the middle of the grip region, or just
  %  the tip of the fingers, then totip would toggle between forward
  %  kinematics for one and the other.
  %
  %  Make sure that the SE3 class is in your path.
  %
  function [g, grip] = forwardkin(this, a, totip)
  a = a*pi/180;
  if (nargin == 2)
    totip = false;
  end
    d1 = [0;0;this.linklen(1)]; 
    d2 = [0;0;this.linklen(2)]; 
    d3 = [0;0;this.linklen(3)]; 
    d4 = [0;0;this.linklen(4)];
    d5 = [0;0;this.linklen(5)];
  if totip
    R1 = SE3.RotZ(a(1)); R2 = SE3.RotX(a(2)); R3 = SE3.RotX(a(3)); 
    R4 = SE3.RotX(a(4)); R5 = SE3.RotZ(a(5));
    g1 = SE3(d1,R1*R2);
    g2 = SE3(d2,R3);
    g3 = SE3(d3,R4);
    g4 = SE3(d4,R5);
    g5 = SE3(d5,eye(3));
    g = g1*g2*g3*g4*g5;
  else
    R1 = SE3.RotZ(a(1)); R2 = SE3.RotX(a(2)); R3 = SE3.RotX(a(3)); 
    R4 = SE3.RotX(a(4)); R5 = SE3.RotZ(a(5)); 
    g1 = SE3(d1,R1*R2);
    g2 = SE3(d2,R3);
    g3 = SE3(d3,R4);
    g4 = SE3(d4,R5);
    g = g1*g2*g3*g4;
  end
  grip = a(6);
  end
  
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
  function alpha = inversekin(this, gdes, totip, solPick)
  
  if ( isempty(totip) )
    totip = false;
  end
  
  if (nargin < 4)
    solPick = [];       % Up to you to decide how to use.
  end
  
  alpha = zeros(5,1);
  d1 = [0;0;this.linklen(1)];
  d2 = [0;0;this.linklen(2)];
  d3 = [0;0;this.linklen(3)];
  if totip
      d4 = [0;0;this.linklen(4)+this.linklen(5)];
  else
      d4 = [0;0;this.linklen(4)];
  end
  
  gh = SE3(d4,eye(3));
  gw = gdes*inv(gh); display(gw);
  dw = getTranslation(gw); x = dw(1); y = dw(2); z = dw(3);
  
  l1 = this.linklen(2); l2 = this.linklen(3); 
  a1 = atan2(-x,y); xp = sqrt(x^2 + y^2); yp = z-this.linklen(1); disp(xp); disp(yp);
  a3 = acos((xp^2 + yp^2 - l1^2 - l2^2)/(2*l1*l2));
  a2 = atan2(xp,yp)-atan2(l2*sin(a3),l1+l2*cos(a3));
  a2 = -a2;
  a3 = -a3;
  g1 = SE3(d1,SE3.RotZ(a1));
  g2 = SE3([0;0;0],SE3.RotX(a2));
  g3 = SE3(d2,SE3.RotX(a3));
  g4 = SE3(d3,eye(3));
  
  gRot = (inv(g4)*inv(g3)*inv(g2)*inv(g1)*gdes*inv(gh));
  Rmat = getRotationMatrix(gRot); display(Rmat);
  
  a4 = -atan2(Rmat(2,3),Rmat(3,3));
  a5 = atan2(-Rmat(1,2),Rmat(1,1));
  
  alpha = [a1;a2;a3;a4;a5]*(180/pi);
  
  end
  
  
  
  %---------------------------- posJacobian ----------------------------
  %
  %  Computes the manipulator Jacobian for the position part only.
  %  The frame of reference is the world coordinate system (it is
  %  not done in body coordinates).
  %(
  function pJac = posJacobian(this, alpha)

  % Construction the Jacobian matrix.
  pJac = [];

  end

  %)
  %----------------------- genPositionTrajectory -----------------------
  %
  %  Given a desired position trajectory, generate the joint angles
  %  that will follow it.  Also assume that the initial joint angles
  %  are availalbe as an argument and that they map to the initial
  %  value of the position trajectory.
  %
  %  ptraj tspan should be valid over the duration of the actual
  %  tspan.  There are no checks here for that case, just assumed
  %  to be the case.
  %
function jtraj = genPositionTrajectory(this, ptraj)

    d = .1/3;
    v = ptraj.pdot;
    a(1:1:5,1) = ptraj.pvec;
    a(6,1) = 0;
    
    d1 = [0;0;this.linklen(1)]; 
    d2 = [0;0;this.linklen(2)]; 
    d3 = [0;0;this.linklen(3)]; 
    d4 = [0;0;this.linklen(4)];
    
    for j = 0:.1:ptraj.tspan(2);
        i = uint64(10*j+1);
        
        R1 = SE3.RotZ(a(1,i)); R2 = SE3.RotX(a(2,i)); 
        R3 = SE3.RotX(a(3,i));
        R4 = SE3.RotX(a(4,i)); R5 = SE3.RotZ(a(5,i));
        
        g1 = SE3(d1,R1);     
        g2 = SE3([0;0;0],R2); 
        g3 = SE3(d2,R3);     
        g4 = SE3(d3,R4);     
        g5 = SE3([0;0;0],R5);   
        g6 = SE3(d4,eye(3));
        
        ge = g1 * g2 * g3 * g4 * g5 * g6;  
        Re = getRotationMatrix(ge); 
        
        Jb = Jacobian(this, a(:,i));  
        vb = inv(Re) * v;      
        Jb1 = Jb(1:1:3,:);
        
        adot(:,i) = pinv(Jb1) * vb;  
        a(1:1:5,i+1) = a(1:1:5,i) + adot(1:1:5,i) * d; 
        a(6,i+1) = 0;
    end
    
    jtraj.alpha = a*180/pi;  
    jtraj.adot = adot*180/pi;  
    jtraj.time  = ptraj.tspan;
end


%     %------------------------- posODE ------------------------
%     %
%     %  this function can access the variables in 
%     %  genPositionTrajectory, so it is all good.
%     %  Access ptraj properly here.
%     %
%     function alphaDot = posODE(alpha)
% 
%     alphaDot = zeros(5,1);
%     Jb = lynx6.Jacobian(alpha);
%     adot = pinv(Jb)*vb;
%     alphaDot = conversionHere;
% 
%     end
% 
%   end

  %------------------------------ Jacobian -----------------------------
  %
  %  Computes the manipulator Jacobian for the lynx6 manipulator.
  %  In this case, it should be the body Jacobian.
  %
  %(
 %------------------------------ Jacobian -----------------------------
  %
  %  Computes the manipulator Jacobian for the lynx6 manipulator.
  %  In this case, it should be the body Jacobian.
  %
  %(
  function JBody = Jacobian(this, alpha)

  % Construction the Jacobian matrix.
    a = alpha;
    Jz = [0;0;0;0;0;1];
    Jx = [0;0;0;1;0;0];

    d1 = [0;0;this.linklen(1)]; 
    d2 = [0;0;this.linklen(2)]; 
    d3 = [0;0;this.linklen(3)]; 
    d4 = [0;0;this.linklen(4)];

    R1 = SE3.RotZ(a(1)); R2 = SE3.RotX(a(2)); R3 = SE3.RotX(a(3)); 
    R4 = SE3.RotX(a(4)); R5 = SE3.RotZ(a(5));
    g1 = SE3(d1,R1);
    g2 = SE3([0;0;0],R2);
    g3 = SE3(d2,R3);
    g4 = SE3(d3,R4);
    g5 = SE3([0;0;0],R5);
    g6 = SE3(d4,eye(3));

    Jb = [adjoint(inv(g2*g3*g4*g5*g6),Jz),adjoint(inv(g3*g4*g5*g6),Jx), adjoint(inv(g4*g5*g6),Jx),adjoint(inv(g5*g6),Jx),adjoint(inv(g6),Jz)];
    
    JBody = [Jb];

  end
  
  
  function Js = SJacobian(this, alpha)
    a = alpha;
    Jz = [0;0;0;0;0;1];
    Jx = [0;0;0;1;0;0];

    d1 = [0;0;this.linklen(1)]; 
    d2 = [0;0;this.linklen(2)]; 
    d3 = [0;0;this.linklen(3)]; 
    d4 = [0;0;this.linklen(4)];

    R1 = SE3.RotZ(a(1)); 
    R2 = SE3.RotX(a(2)); 
    R3 = SE3.RotX(a(3)); 
    R4 = SE3.RotX(a(4)); 
    R5 = SE3.RotZ(a(5));
    
    g1 = SE3(d1,R1);
    g2 = SE3([0;0;0],R2);
    g3 = SE3(d2,R3);
    g4 = SE3(d3,R4);
    g5 = SE3([0;0;0],R5);
    g6 = SE3(d4,eye(3));

    J1 = adjoint(g1, Jz); 
    J2 = adjoint(g2, Jx);
    J3 = adjoint(g3, Jx); 
    J4 = adjoint(g4, Jx); 
    J5 = adjoint(g5, Jz);
    
    Js = [J1, adjoint(g1,J2), adjoint(g1*g2,J3),adjoint(g1*g2*g3,J4),adjoint(g1*g2*g3*g4,J5)];
  
  end

  %)
  %--------------------------- genTrajectory ---------------------------
  %
  %  Given a desired position trajectory, generate the joint angles
  %  that will follow it.  Also assume that the initial joint angles
  %  are availalbe as an argument and that they map to the initial
  %  value of the position trajectory.
  %
  %  ptraj tspan should be valid over the duration of the actual
  %  tspan.  There are no checks here for that case, just assumed
  %  to be the case.
  %
  function alphaTraj = genTrajectory(this, ptraj, tspan)

  % Resolved rate code here.  Use groupODE as the ode function.

  [alphaTraj.time, alphaTraj.alpha] = odeCallHere;

    %------------------------ groupODE -----------------------
    %
    %  this function can access the variables in 
    %  genPositionTrajectory, so it is all good.
    %  Access ptraj properly here.
    %
    function alphaDot = groupODE(t, alpha)

    alphaDot = zeros(5,1);

    % Fill in the ODE here.

    % Make sure to convert angular rates from radians per second
    % to degrees per second.  Use diag function to help you.
    alphaDot = conversionHere;

    end

  end

  %--------------------- followJointTrajectory --------------------
  %
  %  Given a desired joint trajectory, command the joint angles
  %  that will follow it for the desired time span passed as a
  %  second argument.
  %
  %  jtraj tspan should be valid over the duration of the desired
  %  tspan.  There are no checks here for that case, just assumed
  %  to be the case.
  %
  function followJointTrajectory(this, jtraj, time)

  % Code to command trajectory to follow here.
  %
  % 1. Generate linear spacing in time across tspan.
  % 2. Interpolate the alpha values for these times (use interp1).
  % 3. Compute the time differences between times.
  % 4. Using for loop from 2nd element to last element (1st is initial
  %     condition which we are at, so can ignore it), ...
  %    a] send manipulator to new joint witth time duration slightly
  %       more than delta t (say by 5-10%).
  %    b] wait for actual delta t.
  %    c] go to the next joint angle.
  
    a = jtraj.alpha;
    for ii = 1:5            	% Send several commands until it listens.
      this.gotoSleep();       	%  Usually responds by third command.
      pause(0.05);
    end
    pause(1); 
    this.gotoHome(); % Send it to the home position.
    disp('At home position, press return to start ...');
    pause();
    this.setArm(a(:,1), 1.5);
    pause();
    disp('At start position, press return to move ...');
    for i = 1:1:30;
        this.setArm(a(:,i),1.1*.1);
        pause(.1);
    end
    v = this.Jacobian(a)*jtraj.adot; 
    %display(v);

  end

  %------------------------------ shutdown -----------------------------
  %
  %  Send to sleep then close the serial port.
  %
  function shutdown(this)
  
  this.gotoSleep();
  fclose(this.serialid);
  
  end

end

%
%
%============================ Helper Functions ===========================
%
%

methods(Static)

  %-------------------------- closeSerialPort --------------------------
  %
  %  Find an open port using instrfind and then close it.  If the serial
  %  port is not closed for some funny reason, then call this helper
  %  method, as per:
  %
  %  lynx6.closeSerialPorts();
  %
  %  Should work.  Not tested.
  %
  function closeSerialPorts

  ports = instrfind();
 
  for ii=1:length(ports)
    if ( strcmp(ports(ii).Status,'open') )
      fclose(ports(ii));
    end
  end

  end

end


end
%
%================================= lynx6 =================================
