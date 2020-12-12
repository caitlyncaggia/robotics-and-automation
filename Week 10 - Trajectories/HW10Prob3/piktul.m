%================================= piktul ================================
%
%
%
%
%================================= piktul ================================

%
%  Name:		lynx6.m
%
%  Author:		Patricio A. Vela, pvela@gatech.edu
%
%  Created:		2013/01/28
%  Modified:	2013/02/01
%
%
%  NOTES:
%    indent is equal to 2, tabstop is equal to 4.
%
%================================= piktul ================================
classdef piktul < handle


properties

  alphaIds    = [   0    1    2    3    4   ];		
						% Servo ID numbers on the SSC-32 control board.
  alphaHome   = [   1.0;  0;  0;  0;  0.5];		
						% Home position as a joint configuration.
  alphaSleep  = [   1.0;  0;  0;  0;  0.5];	% Straight out slightly up.
				% Sleep position as a joint configuration.
				%   This is where to send the manipulator before powering 
   				%   down.
  alphaLims = [  1, -90, -90, -90, 1/8;
                 1.75 ,   0,   0,   0, 5.5/8;
      	         2.5,  90,  90,  90, 5/4];
						% Limits of the joints, either angular values 
      	    			%   or linear values as in the gripper.
						%   The middle value is some measured
						%   intermediate location.
  musecLims = [920  660  589 621 1064;			
               1580 1581 1560 1520 1527;			
               2240 2370 2324 2350  1990];
						% How these limits translate to microseconds 
      					%   for the servo commands.
  alphaOrient = [ -1,   -1,   1,   1,   -1];		
						% To switch orientation in case servo rotates 
      					%   in the wrong direction 
      					%   (e.g. need right-hand rule).
  mm2in   = 1/25.4;
  linklen;				% Measured link lengths of the manipulator in
						%   millimeters but converted to inches.

  % Serial port setup defaults:
  serport = [];
  serialid = [];

end

%)
%
%============================ Member Functions ===========================
%
%(

methods

  %------------------------------- piktul ------------------------------
  %
  %  Constructor
  %
  %(
  function this = piktul(setup)

  this.linklen = [4.5 4];	

  % Parse the setup argument now.
  if (nargin > 0)
    fnames = fieldnames(setup)
    for ii = 1:length(fnames)
      switch fnames{ii}
        case {'linklen','alphaHome','alphaOrient'},
          eval(['this.' fnames{ii} ' = getfield(setup, fnames{ii});']);
        case 'musecLims',
          if (size(setup.musecLims,2) == 5)
  	        if (size(setup.musecLims, 1) == 2)
  	          this.musecLims(:,1) = setup.musecLims(:,1);
  	          this.musecLims(:,3) = setup.musecLims(:,2);
  	          this.musecLims(:,2) = mean(setup.musecLims,1)
  	        elseif (size(setup.musecLims, 1) == 3)
  	          this.musecLims = setup.musecLims;
  	        end
          end
        case 'alphaLims',
          if (size(setup.alphaLims,2) == 5)
  	        if (size(setup.alphaLims, 1) == 2)
  	          this.alphaLims(:,1) = setup.alphaLims(:,1);
  	          this.alphaLims(:,3) = setup.alphaLims(:,2);
  	          this.alphaLims(:,2) = mean(setup.alphaLims,1)
  	        elseif (size(setup.alphaLims, 1) == 3)
  	          this.alphaLims = setup.alphaLims;
  	        end
          end
        case 'COM','baud','terminator',
          eval(['this.serport.' fnames{ii} ' = getfield(setup, fnames{ii});']);
      end
    end
  end

  this.serport.com  = 'COM1';		% COM port is COM1.
  this.serport.baud = 115200;		% Baud rate is 115200 bps.
  this.serport.terminator = 'CR';	% Message terminator is Carriage Return.

  % Open up the serial port.
  ispc
  if isunix
    device = '/dev/ttyS0';
    this.serialid = fopen(device,'w+');
  elseif ispc
    this.serialid = serial(this.serport.com,'BaudRate',this.serport.baud, ...
                                       'Terminator',this.serport.terminator)
    fopen(this.serialid)
    % Note that in Windows, if the program crashes then the com port is still 
    %   open.  You will not be able to reopen the port and you won't be able to 
    %   access the serialid you had before.  If you can still run the
	%   manipulator class variable, then run the free routine.  If you no 
	%   longer have access to the class variable, the only option I know of
	%   is to restart Matlab.  Or you could try to do an 'fclose all' and
	%   see if that works.  I haven't tested it out.
    %
    %   Read documentation: 'help fclose'
  end

  end

  %)
  %------------------------------ gotoHome -----------------------------
  %
  %  Send manipulator to home position.
  %
  %(
  function gotoHome(this, time)
  
  if (nargin == 1)
    time = 3;
  end
  
  this.setArm(this.alphaHome, time);
  
  end
  
  %)
  %----------------------------- gotoSleep -----------------------------
  %
  %  Send manipulator to sleep position.
  %
  %(
  function gotoSleep(this, time)
  
  if (nargin == 1)
    time = 3;
  end
  
  this.setArm(this.alphaSleep, time);
  
  end
  
  %)
  %----------------------------- setServos -----------------------------
  %
  %  Send low level servo command to control manipulator.
  %
  %(
  function setServos(this, pwmsigs, time, defaultLims)
  
    if (nargin < 3)
      time = 2;
    elseif (isempty(time))
      time = 2;
    end
  
    if (nargin < 4)
      defaultLims = true;
    end
    
    if ( (size(pwmsigs,1) ~= length(this.alphaIds)) || (size(pwmsigs,2) ~= 1) )
      disp('Incorrect dimensions, ignoring command.');
      return;
    end
    
    if (defaultLims)
      ticks = min(this.musecLims(3,:),max(this.musecLims(1,:), pwmsigs'));
    else
      ticks = min(2500, max(500, pwmsigs'));
    end
  
    cmdstr = sprintf('#%d P%d ',[this.alphaIds ; ticks]);		
    cmdstr = [cmdstr sprintf(' T%d \n', round(time*1000))];
  
    %fprintf(cmdstr); 			% Print actual command to STDOUT.
    fprintf(this.serialid, cmdstr);
  
  end
  
  %)
  %------------------------------- setArm ------------------------------
  %
  %  Send joint angle command to control manipulator.
  %
  %(
  function setArm(this, alpha, time)
  
    if ( (nargin < 3) || isempty(time) )
      time = 2;
    end
  
    if ( (size(alpha,1) == length(this.alphaIds)-1) && (size(alpha,2) == 1) )
      alpha(5) = alphaHome(5);
    elseif ( (size(alpha,1) ~= length(this.alphaIds)) || (size(alpha,2) ~= 1) )
      disp('Incorrect dimensions, ignoring command.');
      return;
    end
  
    alpha = min(this.alphaLims(3,:),max(this.alphaLims(1,:), alpha'));
    for ii = 1:length(this.alphaIds)
      if (this.alphaOrient(ii) == 1)
        ticks(ii) = interp1(this.alphaLims(:,ii), this.musecLims(:,ii), ...
		                                                           alpha(ii));
      else
        ticks(ii) = interp1(this.alphaLims(end:-1:1,ii), ...
		                                      this.musecLims(:,ii), alpha(ii));
      end
    end
    ticks = round(ticks);
    cmdstr = sprintf('#%d P%d ',[this.alphaIds ; ticks]);
    cmdstr = [cmdstr sprintf(' T%d \n', round(time*1000))];
  
    %fprintf(cmdstr); 			% Print actual command to STDOUT.
    fprintf(this.serialid, cmdstr);
  
    end
  
    %)
    %---------------------------- displayArm ---------------------------
    %
    %  Display the manipulator.  Joint angles should be in degrees since
    %  the are converted.  (Or remove the line giving the conversion and
    %  let it accept radians.  Up to you.)
    %
    %(
    function displayArm(this, alphadisp)
    if (nargin == 0)
      alphadisp = alpha;
    end

    vparms.home = 'straight-up';
    alphadisp = diag([(pi/180)*ones([1 5]), 1])*alphadisp;	% Convert to rads.
    piktul_display(alphadisp, this.linklen, vparms);

    end

    %)
    %---------------------------- forwardkin ---------------------------
    %
    %  Compute the forward kinematics: get end-effector configurations
	%  given joint angle configuration.
	%
	%(
    function [g, z] = forwardkin(this, alpha1, alpha2, alpha3)
  
%     if (nargin == 1)
%       totip = false;
%     end
  
%     g = forwardkin_straightup((pi/180)*alpha(1:5), this.linklen, totip);
% 	  NEED TO REDO THIS PART!
%     Rz = @(alpha)[cos(alpha), -sin(alpha), 0; ...
%                   sin(alpha), cos(alpha), 0; ...
% 	              0, 0, 1];
    
% note we changed linklen to [4.5 4]
    T = [(4.5*cosd(alpha1) + 4*cosd(alpha1 + alpha2));
         (4.5*sind(alpha1) + 4*sind(alpha1 + alpha2))];
    g = SE2(T, degtorad(alpha1 + alpha2 + alpha3));
    z = 1;

%     %-- This is the straight out home configuration: 
%     ge = SE3([0;0;alpha(1)] , eye(3)) ...
%          * SE3([this.linklen(2); 0; 0] , Rx(alpha(2))) ...
%          * SE3([this.linklen(3); 0; 0] , Rx(alpha(3))) ...
%          * SE3([this.linklen(4); 0; 0] , Rx(alpha(4))) ...
%          * SE3([0; 0; this.linklen(5)] , Rz(alpha(5)));
% 
%     if (totip)
%       ge = ge*SE3([0;0;this.linklen(6)], eye(3));
%     end
  
    end
  
    %)
    %---------------------------- inversekin ---------------------------
    %
	%  Compute the inverse kinematics: get joint configuration given
	%  the desired end-effector configuration.
    %
	%(
    function alpha = lynx6_inversekin(this, gdes)
  
    %if (nargin < 3)
    %  solfact = [];
    %end
    l1 = 4.5;
    l2 = 4.0;
    
    T1 = getTranslation(gdes);
    r1 = sqrt(T1(1)^2 + T1(2)^2);
    c1a2 = acos( (l1^2 + r1^2 - l2^2) / (2*l1*r1) ) + atan2(T1(2), T1(1));
    c1a3 = acos( (l1^2 + l2^2 - r1^2) / (2*l1*l2) ) - pi;
    c1a4 = getTheta(gdes) - c1a2 - c1a3;
    alpha = [0.5 rad2deg(c1a2) rad2deg(c1a3) rad2deg(c1a4) 0.5];

    %alpha = (180/pi)*inversekin_straightup(gdes, linklen, totip, solfact);
  
    end
  
	%)
    %------------------------------- free ------------------------------
    %
	%  Shutdown the system.
    %
	%(
    function shutdown(this)
  
    fclose(this.serialid);
  
    end

	%)
    
    %---------------------- genPositionTrajectory ----------------------
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
    function [time, alphaTraj] = genPositionTrajectory(this, ptraj, tspan, alpha0)
        [time, alphaTraj] = ode45(@adot, tspan, alpha0);
        alphas = alphaTraj(:, 2:4);
        l1 = 1; l2 = 0.5; l3 = 0.35;
        %from HW9.3:
        for i = 1:length(alphas(:,1))
            a1 = alphas(i, 1); a2 = alphas(i,2); a3 = alphas(i,3);
            ptraj(i,:) = [alphaTraj(i,1); 
                         l1*cos(a1) + l2*cos(a1+a2) + l3*cos(a1+a2+a3);
                         l1*sin(a1) + l2*sin(a1+a2) + l3*sin(a1+a2+a3);
                         a1+a1+a3;
                         alphaTraj(i,5)];
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
    function followJointTrajectory(this, time, jtraj, tspan, nsteps)
        % Code to command trajectory to follow here.
        
        % 1. Generate linear spacing in time across tspan.
        t = linspace(tspan(1), tspan(end), nsteps);
        
        % 2. Interpolate the alpha values for these times (use interp1).
        alphas = interp1(time, jtraj, t)
        
        % 3. Compute the time differences between times.
        deltat = t(2) - t(1);
        
        % 4. Using for loop from 2nd element to last element (1st is initial
        %     condition which we are at, so can ignore it), ...
        for i = 2:nsteps
            %  a] send manipulator to new joint with time duration slightly
            %  more than delta t (say by 5-10%).
            duration = deltat * 1.10;
            this.setArm(transpose(alphas(i,:)), duration);
            %displayArm(this, alphas(i,:))
            %  b] wait for actual delta t.
            pause(deltat);  
            %  c] go to the next joint angle.
            % continue loop
        end
        
    end
  
    %--------------------------- posJacobian ---------------------------
    %
    %  Computes the manipulator Jacobian for the position part only.
    %  The frame of reference is the world coordinate system (it is
    %  not done in body coordinates).
    
    function mJ = posJacobian(this, alpha)
        %from HW 9.3: Construction the Jacobian matrix.
        a1 = alpha(1); a2 = alpha(2); a3 = alpha(3);
        l1 = 1; l2 = 0.5; l3 = 0.25;

        J11 = -l1*sin(a1) - l2*sin(a1+a2) - l3*sin(a1+a2+a3);
        J12 = -l2*sin(a1+a2) - l3*sin(a1+a2+a3);
        J13 = -l3*sin(a1+a2+a3);
        J21 = l1*cos(a1) + l2*cos(a1+a2) + l3*cos(a1+a2+a3);
        J22 = l2*cos(a1+a2) + l3*cos(a1+a2+a3);
        J23 = l3*cos(a1+a2+a3);

        Jac = [J11 J12 J13; J21 J22 J23];
        
        % Multiply by alphaDot to get the estimated linear velocity.
         alphaTraj = this.genPositionTrajectory(ptraj, tspan, alpha0);
         mJ = Jac * alphaTraj;
    end

end

%)


end

%
%================================= piktul ================================
