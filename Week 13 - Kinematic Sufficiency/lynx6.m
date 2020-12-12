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
        alphaSleep  = [   0;  35;   70;  70;   0;  1.0];	% Leaning on block.
        % Sleep position as a joint configuration.
        %   This is what you put it to before powering
        %   down.
        
        alphaOrient = [ -1,   1,   1,  -1,   1,  -1];
        % To switch orientation in case servo rotates
        %   in the wrong direction
        %   (e.g. need right-hand rule).
        
        alphaLims = [ -90, -90, -90, -90, -90, 0.00;
            0,   0,   0,   0,   0, 3/4;
            90,  90,  90,  90,  80, 1.125];
        % Limits of the joints, either angular values
        %   or linear values as in the gripper.
        %   The middle value is some measured
        %   intermediate location.
        musecLims = [ 525  640  860  580  790 1475;		%using L6-08
            1530 1480 881 1455 1480 1760;
            2320 2400 2500 2340 2500 2500];
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
            
            %mm2in   = 1/25.4;
            %linklen = [110.6 120 120 130 20]*mm2in;
            % Measured link lengths of the manipulator in
            %   millimeters but converted to inches.
            linklen = [2.75 4.75 5.0 3.0 1.25];
            
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
        function g = forwardkin(this, alpha, totip)
            
            if (nargin == 2)
                totip = false;
            end
            
            alpha = deg2rad(alpha);
            
            a1 = alpha(1);
            a2 = alpha(2);
            a3 = alpha(3);
            a4 = alpha(4);
            a5 = alpha(5);
            
            %link lengths from lynx6.m from class wiki:
            mm2in   = 1/25.4;
            linklen = [110.6 120 120 130 20]*mm2in;
            l0 = linklen(1);
            l1 = linklen(2);
            l2 = linklen(3);
            l3 = linklen(4);
            l4 = linklen(5);
            
            g1 = SE3([0;0;l0], SE3.RotZ(a1));
            g2 = SE3([0;0;0],SE3.RotX(a2));
            g3 = SE3([0;0;l1], SE3.RotX(a3));
            g4 = SE3([0;0;l2],SE3.RotX(a4));
            g5 = SE3([0;0;l3],SE3.RotZ(a5));
            g6 = SE3([0;0;l4],eye(3));
            g = g1*g2*g3*g4*g5*g6;
            
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
            a3 = -acos((xwdes^2 + ywdes^2 + (zwdes - l0)^2 - l1^2 - l2^2)/(2*l1*l2));
            
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
        
        %---------------------------- posJacobian ----------------------------
        %
        %  Computes the manipulator Jacobian for the position part only.
        %  The frame of reference is the world coordinate system (it is
        %  not done in body coordinates).
        %(
        function pJac = posJacobian(this, alpha)
            
            % Construction the Jacobian matrix.
            s1 = sind(alpha(2));
            c1 = cosd(alpha(2));
            12
            s12 = sind(alpha(2)+alpha(3));
            c12 = cosd(alpha(2)+alpha(3));
            mJ = [ 0, -this.linklen(1)*s1 - this.linklen(2)*s12 , ...
                -this.linklen(2)*s12 , 0 , 0 ; ...
                0, this.linklen(1)*c1 + this.linklen(2)*c12 , ...
                this.linklen(2)*c12 , 0 , 0 ; ...
                1, 0 , 0 , 0 , 0 ];
            
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
        function alphaTraj = genPositionTrajectory(this, ptraj, tspan)
            
            % Resolved rate code here.  Use posODE as the ode function.
            
            [alphaTraj.time, alphaTraj.alpha] = ode45(@posODE, tspan, alpha0);
            
            %------------------------- posODE ------------------------
            %
            %  this function can access the variables in
            %  genPositionTrajectory, so it is all good.
            %  Access ptraj properly here.
            %
            function alphaDot = posODE(t, alpha)
                
                alphaDot = zeros(5,1);
                
                % Fill in the ODE here.
                
                % Make sure to convert angular rates from radians per second
                % to degrees per second.  Use diag function to help you.
                J = posJacobian(alpha);
                alphaDot = pinv(J)*ptraj.velocity(t);
                % Make sure to convert angular rates from radians per second
                % to degrees per second. Use diag function to help you.
                alphaDot = diag([1, 180/pi, 180/pi, 180/pi, 1])*alphaDot;
                
            end
            
        end
        
        %------------------------------ Jacobian -----------------------------
        %
        %  Computes the manipulator Jacobian for the lynx6 manipulator.
        %  In this case, it should be the body Jacobian.
        %
        %
        function mJacBody = Jacobian(this, alpha)
            
            alpha = ((pi/180)*[1;1;1;1;1]).*alpha;
            
            Jbjoint = [ [0;0;0;0;0;1], [0;0;0;1;0;0], [0;0;0;1;0;0], ...
                [0;0;0;1;0;0], [0;0;0;0;0;1] ];
            
            g(1) = SE3([0;0;l0], RotZ(alpha(1)));
            g(2) = SE3([0;0;0] , RotX(alpha(2)));
            g(3) = SE3([0;0;l1], RotX(alpha(3)));
            g(4) = SE3([0;0;l2], RotX(alpha(4)));
            g(5) = SE3([0;0;l3], RotZ(alpha(5)));
            g(6) = SE3([0;0;l4], eye(3));
            
            adInv = eye(6);
            Jb = zeros(6, 5);
            for ii=size(Jb,2):-1:1
                adInv = adInv*inv(g(ii+1).adjoint());
                Jb(:,ii) = adInv * Jbjoint(:,ii);
            end
            
            
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
            
            [alphaTraj.time, alphaTraj.alpha] = ode45(@groupODE, tspan, alpha0);
            
            %------------------------ groupODE -----------------------
            %
            %  this function can access the variables in
            %  genPositionTrajectory, so it is all good.
            %  Access ptraj properly here.
            %
            function alphaDot = groupODE(t, alpha)
                
                % Get end-effector frame and map velocity into body coords.
                g = this.forwardkin(alpha);
                vBody = inv(g).*[ptraj.velocity(t);0;0;0];
                vBody = vBody(1:3);
                % Get the body Jacobian, strip the angular part, and then
                % compute the pseudo-inverse.
                J = this.Jacobian(alpha);
                J = J(1:3,:); % Take only the position part.
                pinvJ = transpose(J)*inv(J*transpose(J) + 0.001^2*eye(size(J,1)));
                % Apply pseudo-inverse to desired end-effector linear body velocity.
                alphaDot = pinvJ*diag([1,1,1,pi/180,pi/180,pi/180])*vBody;
                % Make sure to convert angular rates from radians per second
                % to degrees per second. Use diag function to help you.
                alphaDot = diag((180/pi)*[1, 1, 1, 1, 1])*alphaDot;
                
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
        function followJointTrajectory(this, jtraj, tspan)
            
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
            %
            
            % 1. Generate linear spacing in time across tspan.
            tsteps = linspace(tspan(1),tspan(2),nsteps);
            
            % 2. Compute the time differences between times.
            tdiffs = diff(tsteps);
            
            % 3. Interpolate the alpha values for these times (use interp1).
            alphaSpace = transpose(interp1(jtraj.time, jtraj.alpha, tsteps));
            
            % 4. Using for loop from 2nd element to last element (1st is initial
            % condition which we are at, so can ignore it), ...
            % a] send manipulator to new joint witth time duration slightly
            % more than delta t (say by 5-10%).
            % b] wait for actual delta t.
            % c] go to the next joint angle.
            %
            figure(4);
            clf; hold off;
            for ii=2:length(tsteps) % Assuming already at first step.
                this.setArm(alphaSpace(:,ii));
                pause(tdiffs(ii-1));
            end
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
    
    
            %----------------------------- displayArm ----------------------------
        %
        %  Display the manipulator.  Joint angles should be in degrees since
        %  the are converted.  (Or remove the line giving the conversion and
        %  let it accept radians.  Up to you.)
        %
        %(
        function displayArm(alphadisp)
            
%             if (nargin == 1)
%                 alphadisp = this.alpha;
%             end
            
            linklen = [4.35433 4.72441 4.72441 5.11811 0.787402];
            
            vparms.home = 'straight-up';
            %alphadisp = diag([(pi/180)*ones([1 5])])*alphadisp;	% Convert to rads.
            lynx6_display(alphadisp, linklen, vparms);
            
        end
    
end

end
%
%================================= lynx6 =================================
