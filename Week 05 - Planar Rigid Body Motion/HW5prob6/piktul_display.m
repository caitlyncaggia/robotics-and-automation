%============================= piktul_display ============================
%
%  function piktul_display(alpha, mgeom, viewparms)
%
%
%  Displays a 3D version of the piktul SCARA-type manipulator.
%
%
%  INPUT:
%    alpha      -  the joint angles, (z, a1, a2, a3, gw)
%               z   - height of manipulator.
%               ai  - one of the three revolute joints (a1, a2, a3).
%               gw  - the gripper width.
%    mgeom      - the manipulator geometry. optional.
%               default values will depict the piktul.
%    viewparms  - optional structure that specifies what view to take.
%               fields:
%                   view    - specify lat/long view angle.
%
%  OUTPUT:
%    The current figure should have a 3D rendering of the piktul.
%
%
%  REQUIRES:
%    SE3 class should have minimal functionality, which is:
%       times, mtimes, leftact, and inverse.
%
%============================= piktul_display ============================

%
%  Name:		piktul_display.m
%
%  Author:		Patricio A. Vela,			pvela@gatech.edu
%
%  Created:		2013/02/01
%  Modified:	2014/02/25
%
%============================= piktul_display ============================
function piktul_display(alpha, mgeom, viewparms)

if ( (nargin < 2) || isempty(mgeom) )
  mgeom.base  = [4.5, 3, 1.5, 1.5];
  mgeom.tower = [1.5, 2, 4.5];			% (x, y, z) dims of tower.
  mgeom.stage = [2.5, 1.5, 0.25, 1, 0.5];
  mgeom.link0 = [1.5, 1.5,4];
  mgeom.link1 = [4, 1, 0.5];
  mgeom.link2 = [3, 1, 0.25];
  mgeom.link3 = [0, 0, -1.5];
  mgeom.wrist = [0.75, 1.75, 0.5];
  mgeom.gripper = [0.75, 0.2, 1.00];
end

try
  g0 = SE3();
catch
  error('Requires SE3 class to function. Should have basic operations.');
end

maxlen = mgeom.link0(1) + mgeom.link1(1) + mgeom.link2(1) + mgeom.wrist(2)/2;

if ( (nargin < 3) )
  viewparms = [];
end
if (~isfield(viewparms,'view'))
  viewparms.view = {-20, 50};
  %viewparms.view = {90, 0}; 	%Front
  %viewparms.view = {0, 0};		%Side.
  %viewparms.view = {0, 90};	%Top.
end

wasHeld = ishold;

deg2rad = pi/180;
alpha = diag([1, deg2rad, deg2rad, deg2rad, 1])*alpha;

Rx = @(alpha)[ 1, 0, 0; 0, cos(alpha), -sin(alpha); 0, sin(alpha), cos(alpha)];
Ry = @(alpha)[cos(alpha), 0, sin(alpha); 0, 1, 0; -sin(alpha), 0, cos(alpha)];
Rz = @(alpha)[cos(alpha), -sin(alpha), 0; sin(alpha), cos(alpha), 0; 0, 0, 1];

gTower = SE3([0;0;mgeom.tower(3)/2], eye(3));
g1 = SE3([0; 0; alpha(1)], eye(3));
gStage = SE3([mgeom.link0(1)/2; 0; mgeom.link0(2)], eye(3) );

g2 = SE3([mgeom.link0(1); 0; mgeom.link0(2)], Rz(alpha(2)) );
g3Half = SE3([mgeom.link1(1)/2; 0; 0], eye(3) );

g3 = SE3([mgeom.link1(1); 0; 0], Rz(alpha(3)) );
g4Half = SE3([mgeom.link2(1)/2; 0; 0], eye(3) ); 

g4 = SE3([mgeom.link2(1)-mgeom.wrist(1)/2; 0; 0], Rz(alpha(4)) ); 
g5Half = SE3([0; 0; -mgeom.wrist(3)/2], eye(3));

g5 = SE3([0; 0; mgeom.link3(3)], eye(3));

axfact = 10;

plot3(axfact*[-1 1],[0 0],[0 0],'b-.');
hold on;
plot3([0 0], axfact*[-1 1],[0 0],'b-.');
plot3([0 0], [0 0], axfact*[0 1],'b-.');

gStage = g1*gStage;
gLink1 = g1*g2*g3Half;
gLink2 = g1*g2*g3*g4Half;
gWrist = g1*g2*g3*g4*g5Half;
gGripper = g1*g2*g3*g4*g5Half;
gGripL = gGripper*SE3([0; ...
                    alpha(5)/2+mgeom.gripper(2)/2;-mgeom.gripper(3)/2], eye(3));
gGripR = gGripper*SE3([0; ...
                   -alpha(5)/2-mgeom.gripper(2)/2;-mgeom.gripper(3)/2], eye(3));
hold on;
shBlock(gTower, mgeom.tower, [1, 1, 1; 0.5, 0.5, 0.5]);
useTrans = [0.5 1];
shBlock(gStage, mgeom.stage(1:3),  [1, 1, 1; 0.5, 0.5, 0.5], useTrans,1);
shBlock(gLink1, mgeom.link1(1:3),  [1, 1, 1; 0.5, 0.5, 0.5], useTrans,1);
shBlock(gLink2, mgeom.link2(1:3),  [1, 1, 1; 0.5, 0.5, 0.5], useTrans,1);
shBlock(gWrist, mgeom.wrist(1:3),  [0.5, 0.5, 1.0; 0.5, 0.5, 0.5], useTrans,1);
shBlock(gGripL, mgeom.gripper(1:3), [0.5, 0.5, 1.0; 0.5, 0.5, 0.5], useTrans,1);
shBlock(gGripR, mgeom.gripper(1:3), [0.5, 0.5, 1.0; 0.5, 0.5, 0.5], useTrans,1);

hold on;
outerRing = [];
innerRing = [];
leftRing = [];
rightRing = [];
manRad = mgeom.link1(1) + mgeom.link2(1) - mgeom.gripper(1)/2;
minRad = sqrt(mgeom.link1(1)^2 + (mgeom.link2(1) - mgeom.gripper(1)/2)^2);
edgRad = mgeom.link2(1) - mgeom.gripper(1)/2;
manOrig = [mgeom.link0(1);0;0];
leftOrig = manOrig + [0;mgeom.link1(1);0];
rightOrig = manOrig + [0;-mgeom.link1(1);0];
theta = linspace(-pi/2,pi/2,100);
for th = theta
  outerRing = [outerRing, Rz(th) * [manRad; 0; 0] + manOrig];
end
angOffset = atan( (mgeom.link2(1) - mgeom.gripper(1)/2) / mgeom.link1(1));
theta = linspace( -angOffset-pi/2, angOffset+pi/2, 100 );
for th = theta
  innerRing = [innerRing, Rz(th) * [minRad; 0; 0] + manOrig];
end
theta = linspace(3*pi/2,2*pi,100);
for th = theta
  leftRing  = [leftRing, Rz(th)*[-edgRad;0;0] + leftOrig];
  rightRing = [rightRing, Rz(-th)*[-edgRad;0;0] + rightOrig];
end
plot3(outerRing(1,:), outerRing(2,:), outerRing(3,:),'r-.');
plot3(innerRing(1,:), innerRing(2,:), innerRing(3,:),'r-.');
plot3(leftRing(1,:), leftRing(2,:), leftRing(3,:),'r-.');
plot3(rightRing(1,:), rightRing(2,:), rightRing(3,:),'r-.');

axis equal;
set(gca,'Color',[0, 0, 0]);
set(gca,'XColor',[1, 1, 1]);
set(gca,'YColor',[1, 1, 1]);
set(gca,'ZColor',[1, 1, 1]);
xlabel('x');
ylabel('y');
axis([-maxlen*2/3, maxlen, -maxlen, maxlen, 0, 5]);
view(viewparms.view{:});
wasHeld
if ~wasHeld
  hold off;
end

end

%-------------------------------- shBlock --------------------------------
%
%  function blockh = shBlock(g, blens, colors, transp, drawLines)
%
%
%
%  Inputs:
%    g			- SE(3) description of block location.
%    blens		- Lengths of the block sides (along body x, y, and z axes)
%    colors		- Vector describing the color of the faces [1x3 vector]
%					If desired, the edge color can be specified.  In 
%					this case, then input is a [2x3] vector.
%					Default is same as face color.
%    transp		- The transparency level of the block.
%    drawLines	- Draw the edge lines?  [default is yes]
%					This overrides the colors specification.
%
%-------------------------------- shBlock --------------------------------

%
%  Name:		shBlock.m
%
%  Author:		Patricio A. Vela, 				pvela@gatech.edu
%
%  Created:		2011/05/11
%  Modified:	2011/05/11
%
%-------------------------------- shBlock --------------------------------
function blockh = shBlock(g, blens, color, transp, drawLines)

if (size(color,1) == 1)
  color = repmat(color,[2,1]);
end

if (nargin <  5)
  drawLines = true;
end

if (nargin < 4)
  transp = [1 1];
elseif (isscalar(transp))
  transp = [transp transp];
end

facepts = blens([1, 1, 1, 1; 2, 2, 2, 2; 3, 3, 3, 3])/2;

f{1} = facepts .* [ -1, -1,  1,  1 ; -1, -1, -1, -1 ;  1, -1, -1,  1 ];
f{2} = facepts .* [  1,  1,  1,  1 ; -1, -1,  1,  1 ; -1,  1,  1, -1 ];
f{3} = facepts .* [  1,  1, -1, -1 ;  1,  1,  1,  1 ;  1, -1, -1,  1 ];
f{4} = facepts .* [ -1, -1, -1, -1 ;  1,  1, -1, -1 ;  1, -1, -1,  1 ];
f{5} = facepts .* [ -1, -1,  1,  1 ;  1, -1, -1,  1 ;  1,  1,  1,  1 ];
f{6} = facepts .* [ -1, -1,  1,  1 ;  1, -1, -1,  1 ; -1, -1, -1, -1 ];

blockHold = ishold;
hold on;
for ii=1:6
  f{ii} = g .* [f{ii}; ones(1, 4)];
  f{ii} = transpose(f{ii}(1:3,:));
  fh(ii) = fill3(f{ii}(:,1), f{ii}(:,2), f{ii}(:,3), color(1,:));

  if (drawLines)
    set(fh(ii),'EdgeColor',color(2,:));
  else
    set(fh(ii),'LineStyle','none');
  end
  set(fh(ii),'EdgeAlpha', transp(1), 'FaceAlpha', transp(2));
end
if (~ishold)
  hold off;
end

if (nargout)
  blockh = fh;
end

end

%
%-------------------------------- shBlock --------------------------------
