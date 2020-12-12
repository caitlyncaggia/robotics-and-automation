    %--------------------------- posJacobian ---------------------------
    %
    %  Computes the manipulator Jacobian for the position part only.
    %  The frame of reference is the world coordinate system (it is
    %  not done in body coordinates).
    %(
    function pDot = posJacobian(this, alpha, alphaDot)

    % Construction the Jacobian matrix.
    pJac = [];

    % Multiply by alphaDot to get the estimated linear velocity.
    pDot = [];

    end

    %)
    %---------------------- genPositionTrajectory ----------------------
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

    [alphaTraj.time, alphaTraj.alpha] = odeCallHere;

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
      alphaDot = conversionHere;

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
    function followPositionTrajectory(this, jtraj, tspan)

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

    end


