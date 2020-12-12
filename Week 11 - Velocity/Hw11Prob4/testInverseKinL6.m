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
%============================== testInverseKin ==============================

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

giR = [0.8660 0.5000 0.0000;
      0.5000 -0.866 0.0000;
      0.0000  0.000 -1.000];
giT = [-4.4079; 7.6348; 0.8110];
gi = SE3(giT, giR);
ai = l6.inversekin(gi);
ai = [ai; 0.8];
aifml = l6.inversekin(SE3([-4.4079; 7.6348; 1.8110], giR));
aifml = [aifml; 0.8];

gfR = [0.7071 0.0000 -0.7071;
      0.7071 0.0000  0.7071;
      0.0000 -1.000  0.0000];
gfT = [-8.7393; 8.7393; 6.0836];
gf = SE3(gfT, gfR);
af = l6.inversekin(gf);
af = [af; 0.8];

alpha_traj = [transpose(aifml); transpose(af)];
alpha = [transpose(ai); transpose(af)];

pauses = [2, 2, 2, 2, 2, 2, 2];
				% Time between waypoints (in seconds).
				% You can adjust a bit if you'd like,
				% so long as it's clear what is going on.

for ii = 1:2				% Run the routine.
  gend = l6.forwardkin(transpose(alpha(ii,:)));	% Get forward kinematics.
  disp('=========');
  disp([ 'alpha = [ ' num2str(alpha(ii,:), ' %5.2f ') ' ]' ]);
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

