%============================== testCalib01 ==============================
%
%  script testCalib
%
%  This script moves to some pre-programmed positions and outputs the
%  end-effector configuration associated to the joint positions.
%  The manipulator will be given a series of joint angle commands 
%  that it will be commanded to go to.  If the calibration is off, then
%  the manipulator will not achieve the proper manipulator joint
%  configurations, nor will the forward kinematics accurately predict
%  the end-effector configuration.
%
%============================== testCalib01 ==============================

%
%  Name:		testCalib01
%
%  Author:		Patricio A. Vela,		pvela@gatech.edu
%  Created:		2015/03/03
%  Modified:	2016/10/04
%
%============================= testJac =============================
function testJac

myparms = load('pikparms.mat');
pikh = piktul(myparms);
fingerWidth = 0.69;

for ii = 1:5            	% Send several commands until it listens.
  pikh.gotoSleep();       	%  Usually responds by third command.
  pause(0.05);
end
pause(1);
pikh.gotoHome();          	% Send it to the home position.
disp('At home position, press return to start ...');
pause();


% Joint configurations.
alpha0 = [1.75; -30; 18; -25; 0.45];
ptraj = [];
tspan = [0 3];
[time, alpha_traj] = pikh.genPositionTrajectory(ptraj, tspan, alpha0);
pikh.followJointTrajectory(time, alpha_traj, tspan, 4)

                    
                            
% for ii = 1:size(alpha_traj,2)   				% Run the routine.
%   pikh.setArm(alpha_traj(:,ii), 3);
%   disp('Verify manipulator joint configuration, then press return.');
%   %alphas = alpha_traj(:,ii);
%   %[g, ~] = pikh.forwardkin(alphas(2), alphas(3), alphas(4))
%   pause();      
% end             
disp('Done, shutting down.');

%alpha_fk = [1.2; 30; -60; -10; 0.75];    
%pikh.setArm(alpha_fk);
disp('Done, hit return to end.');    
pause();
  
pikh.gotoSleep();         	% Goto sleep,
pikh.shutdown();			%    then shutdown.


end


%============================= testJac ==============================
