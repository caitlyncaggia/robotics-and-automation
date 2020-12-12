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
%============================== splineL6 ==============================


l6 = lynx6;
for ii = 1:5			% Send several commands until it listens.
  l6.gotoSleep();		%  Usually responds by third command.
  pause(0.05);
end
pause(1);
l6.gotoHome();			% Send it to the home position.
disp('At home position, press return to start ...');
pause();

alpha_i = (180/pi)*[-pi/6 ; -pi/3 ; pi/6 ; -pi/2 ; 0];
alpha_f = (180/pi)*[ pi/6 ; 0 ; -pi/12 ; pi/3 ; -pi/3];
% Generate trajectory.
T = 8;
pptraj = pathspline(alpha_i, alpha_f, T);
alpha_traj = zeros(8,5);
for i = 1:8
    alpha_traj(i,:) = pptraj * [1; i; i^2; i^3];
end

pauses = [2, 2, 2, 2, 2, 2, 2];
				% Time between waypoints (in seconds).
				% You can adjust a bit if you'd like,
				% so long as it's clear what is going on.

for ii = 1:length(alpha_traj(1,:))				% Run the routine.
  gend = l6.forwardkin(transpose(alpha_traj(ii,:)));	% Get forward kinematics.
  disp('=========');
  disp([ 'alpha = [ ' num2str(alpha_traj(ii,:), ' %5.2f ') ' ]' ]);
  disp('  -----  ');
  gend										% Display g_e.
  disp('=========');
  l6.setArm(transpose(alpha_traj(ii,:)), pauses(ii));	% Go to position.
  %pause(pauses(ii)*1.1);					% Wait.
  pause();
end
disp('Done, shutting down.');

l6.gotoSleep();			% Goto sleep.
l6.shutdown();			% Shutdown communication with manipulator.
clear  ii l6;

%
%============================== verifyCalib ==============================