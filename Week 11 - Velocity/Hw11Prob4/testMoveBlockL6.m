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
%============================== testMoveBlockL6 ==============================

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

giR = [0.7071 -0.7071 0.0000;
      -0.7071 -0.7071 0.0000;
      0.0000  0.000 -1.000];
giT = [5; 5; 0.35];
gi = SE3(giT, giR);
ai = l6.inversekin(gi);
aimove = [ai(1:5); 1.25];
aipick = [ai(1:5); 0.62];

gfR = [1  0  0;
       0 -1  0;
       0  0 -1];
gfT = [-4.4; 5.5; 0.25];
gf = SE3(gfT, gfR);
af = l6.inversekin(gf);
afmove = [af(1:5); 0.62];
afplace = [af(1:5); 1.25];

alpha_traj = [transpose(aimove); 
              transpose(aipick);
              transpose(afmove);
              transpose(afplace)];

pauses = [2, 2, 2, 2, 2, 2, 2];
				% Time between waypoints (in seconds).
				% You can adjust a bit if you'd like,
				% so long as it's clear what is going on.

for ii = 1:length(alpha_traj(:,1))				% Run the routine.
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