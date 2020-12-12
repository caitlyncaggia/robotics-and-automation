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
    function alpha = lynx6_inversekin(gdes, totip, solfact)
  
    if (nargin < 3)
      solfact = [];
    end
  
    alpha = (180/pi)*inversekin_straightup(gdes, linklen, totip, solfact);
  
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

end

%)


end

%
%================================= piktul ================================
