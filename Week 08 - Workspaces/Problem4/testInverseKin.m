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
%============================= testInverseKin =============================
function testInverseKin

myparms = load('pikparms.mat');
pikh = piktul(myparms);
fingerWidth = 1;

for ii = 1:5            	% Send several commands until it listens.
  pikh.gotoSleep();       	%  Usually responds by third command.
  pause(0.05);
end
pause(1);
pikh.gotoHome();          	% Send it to the home position.
disp('At home position, press return to start ...');
pause();


gpick = SE2([6.893; -3.447], pi/4);
gplace = SE2([6.468; 5.315], atan2(-0.866, 0.5));

%inverse kinematics
alphapick = pikh.lynx6_inversekin(gpick);
alphaplace = pikh.lynx6_inversekin(gplace);


% Initial and final joint configurations.
alpha_traj = transpose([1.75 ,  alphapick(2),  alphapick(3), alphapick(4),  1.0; 
                        0.25 ,  alphapick(2),  alphapick(3), alphapick(4),  1.0;
                        0.25 ,  alphapick(2),  alphapick(3), alphapick(4),  0.5;
                        1.75 ,  alphapick(2),  alphapick(3), alphapick(4),  0.5;
                        1.75 ,  alphaplace(2),  alphaplace(3), alphaplace(4), 0.5;
                        0.20 ,  alphaplace(2),  alphaplace(3), alphaplace(4), 0.5;
                        0.20 ,  alphaplace(2),  alphaplace(3), alphaplace(4), 1.0;
                        1.75 ,  alphaplace(2),  alphaplace(3), alphaplace(4), 1.0]);
                    
                            
for ii = 1:size(alpha_traj,2)   				% Run the routine.
  pikh.setArm(alpha_traj(:,ii), 3);
  disp('Verify manipulator joint configuration, then press return.');
  alphas = alpha_traj(:,ii);
  [g, ~] = pikh.forwardkin(alphas(2), alphas(3), alphas(4))
  pause();      
end             
disp('Done, shutting down.');

alpha_fk = [1.2; 30; -60; -10; 0.75];    
pikh.setArm(alpha_fk);
disp('Done, hit return to end.');    
pause();
  
pikh.gotoSleep();         	% Goto sleep,
pikh.shutdown();			%    then shutdown.


end


% Other options for alpha_traj
%
%   1, -30, -30, -30, 0.75 ;
%   Check out the arm. Gripper should be x-axis aligned.
%
%   2, -30, -60, 45,  0.25 ;
%   Checkout the arm.  Arm link should be aligned with x-axis.
%   Gripper should be diagonal.
%
%   0.25, 45, -45, -45, 0.25 ;
%   Arm should align with y-axis.
%   Gripper should be diagonal, other way now.
%
%============================= testForwardKin==============================
