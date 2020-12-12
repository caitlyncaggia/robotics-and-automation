% ECE 4560 - Homework 6, Problem 7
% Caitlyn Caggia

%modified from testCalib01.m

function HW6prob7

myparms = load('pikparms.mat');
pikh = piktul(myparms);
fingerWidth = 1;
z = [1 -1 1/4];
g0 = [0 0.5 pi/6];
T0 = 0; Tfinal = 4;

[t, z] = ode45(@HW6prob7func, [T0 Tfinal], g0);
a2 = z(:,1); a3 = z(:,2); a4 = z(:,3);
dtNow = t(2) - t(1);

for ii = 1:5            	% Send several commands until it listens.
  pikh.gotoSleep();       	%  Usually responds by third command.
  pause(0.05);
end
pause(1);
pikh.gotoHome();          	% Send it to the home position.

disp('At home position, press return to start ...');
pause();
                            
for ii = 1:t  				% Run the routine.
  pikh.setArm(alphaNow, dtNow*1.1);
  pause(dtNow); 
  alphaNow = [1; a2(ii); a3(ii); a4(ii); 3/4];
end             
disp('Done, shutting down.');

alpha_fk = [1.2; 30; -60; -10; 0.75];    
pikh.setArm(alpha_fk);
disp('Done, hit return to end.');    
pause();
  
pikh.gotoSleep();         	% Go to sleep,
pikh.shutdown();			%    then shutdown.


end

