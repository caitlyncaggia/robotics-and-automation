%============================== verifyCalib ==============================
%
%  script verifyCalib
%
%  This script moves to some pre-programmed positions and outputs the
%  end-effector configuration associated to the joint positions.
%  The manipulator will be given a series of joint angle commands 
%  that it will be commanded to go to.  If the calibration is off, then
%  the manipulator will not achieve the proper manipulator joint
%  configurations, nor will the forward kinematics accurately predict
%  the end-effector configuration.
%
%============================== verifyCalib ==============================

clear all
close all

calib_param = load('lynxparms.mat');

l6 = lynx6(calib_param);
for ii = 1:5			% Send several commands until it listens.
  l6.gotoSleep();		%  Usually responds by third command.
  pause(0.05);
end
pause(1);
l6.gotoHome();			% Send it to the home position.
disp('At home position, press return to start ...');
pause();

				% Initial and final joint configurations.
alpha_traj    = transpose([ 0 ,  0 ,  0 ,  0 ,  0 , 0.5  ;		% part (a)
                            0 , -90,  0 ,  0 ,  0 , 1.0  ;		% part (b)
                           -45,  0 ,  30,  60,  30, 1.0  ;		% part (c)
                           -90,  0 , -90,  0 , -45, 1.0  ;
                           -90, -60,  30,  30, -90, 0.1  ;
                            0 ,  30, -30,  45,  0 , 1.5  ;
                            30, -45, -45,  60,  45, 0.25 ]);

pauses = [2, 2, 2, 2, 2, 2, 2];
				% Time between waypoints (in seconds).
				% You can adjust a bit if you'd like,
				% so long as it's clear what is going on.

for ii = 1:size(alpha_traj,2)				% Run the routine.
  gend = l6.forwardkin(alpha_traj(:,ii));	% Get forward kinematics.
  disp('=========');
  disp([ 'alpha = [ ' num2str(transpose(alpha_traj(:,ii)), ' %5.2f ') ' ]' ]);
  disp('  -----  ');
  gend										% Display g_e.
  disp('=========');
  l6.setArm(alpha_traj(:,ii), pauses(ii));	% Go to position.
  %pause(pauses(ii)*1.1);					% Wait.
  pause();
end
disp('Done, shutting down.');

l6.gotoSleep();			% Goto sleep.
l6.shutdown();			% Shutdown communication with manipulator.
clear  ii l6;

%
%============================== verifyCalib ==============================