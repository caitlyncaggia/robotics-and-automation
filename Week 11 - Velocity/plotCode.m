
  %============================ velocityPlot ===========================
  %
  %  Plots a vector velocity of SE(2) as a vector and a rotation.
  %  Assumes that this is not given in body coordinates, but in the
  %  world frame (so it is in mixed frames for velocities, so to speak).
  %
  function velocityPlot(g, vect, rad)

  if (nargin < 3)
    rad = 0.5;
  end

  basePt = g.getTranslation();      % Get base point of twist.

  rotAngle = g.getAngle();          % Get rotation angle for velocity.
  thetaVals = rotAngle + linspace(0, 3*pi/6, 20);
  if (vect(3) < 0)
    thetaVals = -thetaVals;
  end
  arcPts = rad*[cos(thetaVals) ; sin(thetaVals)];
  arcVec = diff(arcPts(:,end-1:end),1,2);

  wasHeld = ishold;

  hold on;
  if ((vect(1) ~= 0) || (vect(2) ~= 0))
    quiver(basePt(1), basePt(2), vect(1), vect(2), ...
                                        'LineWidth',2, 'Color', [0, 0, 1]);
  end

  if (vect(3) ~= 0)
    plot(arcPts(1,1:end-1), arcPts(2,1:end-1), 'LineWidth', 2);
    qh = quiver(arcPts(1,end-1), arcPts(2,end-1), arcVec(1), arcVec(2), ...
           0, 'LineWidth', 2, 'MaxHeadSize', 20, 'Color', [0, 0, 1]);
    get(qh)
  end

  if (~wasHeld)
    hold off;
  end

  end

  %============================= twistPlot =============================
  %
  %  Plots a vector velocity of SE(2) as a vector and a rotation at
  %  the object frame defined by g, presuming that the group element
  %  is given in terms of the frame of reference to plot in.  Typically
  %  will be the body frame in the world frame if it is the body
  %  velocity.  Otherwise, should be the identity element to plot as the
  %  spatial velocity.
  %
  %  When rotation is positive, arc will be in first quadrant.
  %  When negative, arc will be in second quadrant.
  %  When zero, no arc.
  %
  %  Likewise, if linear velocity is zero, then no vector.
  %
  %  TODO: Modify so that linespec can be adjusted.
  %
  function twistPlot(g, xi, rad)

  if (nargin < 3)
    rad = 0.5;
  end

  basePt = g.getTranslation();      % Get base point of twist.
  rotMat = g.getRotation();         % Get rotation matrix to transform ...
  vect = rotMat * xi(1:2);          %   twist vector to obs frame.

  rotAngle = g.getAngle();          % Get rotation angle for velocity.
  thetaVals = rotAngle + linspace(0, 3*pi/6, 20);
  if (xi(3) < 0)
    thetaVals = -thetaVals;
  end
  arcPts = rad*[cos(thetaVals) ; sin(thetaVals)];
  arcVec = diff(arcPts(:,end-1:end),1,2);

  wasHeld = ishold;

  hold on;
  if ((xi(1) ~= 0) || (xi(2) ~= 0))
    quiver(basePt(1), basePt(2), vect(1), vect(2), ...
                                        'LineWidth',2, 'Color', [0, 0, 1]);
  end

  if (xi(3) ~= 0)
    plot(arcPts(1,1:end-1), arcPts(2,1:end-1), 'LineWidth', 2);
    qh = quiver(arcPts(1,end-1), arcPts(2,end-1), arcVec(1), arcVec(2), ...
           0, 'LineWidth', 2, 'MaxHeadSize', 20, 'Color', [0, 0, 1]);
    get(qh)
  end

  if (~wasHeld)
    hold off;
  end

  end

