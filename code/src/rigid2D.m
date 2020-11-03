function [Xfinal, trajectory] = rigid2D(options,prams,xc,tau)
ttotal = tic; % start a timer

om = monitor(options,prams);
% build object for doing I/O

tt = tstep(options,prams);

% build initial condition of rigid body particle
geom = capsules(prams,xc,tau);
% plot geometry if usePlot == true
om.plotData(geom);

om.initializeFiles(geom);

time = 0;
step = 0;

xc0  = xc;
tau0 = tau;
om.writeData(time,geom.center,geom.tau,geom.X);
% write the velocity to a file,initial velocities are 0
om.writeVelData(time,0*geom.center,0*geom.tau);  


trajectory = [xc(1,:) xc(2,:) tau];

% begin time loop
%while time < prams.T
while step < prams.m
  % start a timer for this particular time step
  tSingleStep = tic;   

%  msg = ['Yukwa reguired  ' num2str(iterYukawa) ' iterations'];
%  om.writeMessage(msg);
%  msg = ['Stokes reguired ' num2str(iterStokes) ' iterations'];
%  om.writeMessage(msg);

  % comupute density function, translational velocity, angular velocity,
  % and the GMRES output

  % Regular Forward Euler
  if prams.order == 1
    geom = capsules(prams,xc,tau);
    [Up, wp] = tt.timeStep(geom);
    xc  = xc  + tt.dt*Up;
    tau = tau + tt.dt*wp;
    geom2 = geom;
  
  %Adams-Bashforth 
  elseif prams.order == 2
    if step == 0
      % update geometry
      geom0 = capsules(prams,xc0,tau0);
      [Up0, wp0,~,~,etaY0,etaS0] = tt.timeStep(geom0,geom0.X,geom0.X);
      
      % causes xc2, tau2 to be forward Euler, at step 0
      xc1  = xc0;
      tau1 = tau0;

      % update time
      time = time + tt.dt;

      % write current time to console
      message = ['Completed time ', num2str(time,'%4.2e'), ...
          ' of time horizon ' , num2str(prams.T,'%4.2e')];
      om.writeMessage(message);
      
      % write the shape to a file
      om.writeData(time,geom0.center,geom0.tau,geom0.X);

      % write the velocity to a file
      om.writeVelData(time,Up0,wp0);
      
      % update step counter
      step = step + 1;
    end
      
    geom1 = capsules(prams,xc1,tau1);
    [Up1, wp1,~,~,etaY,etaS] = tt.timeStep(geom1,etaY0,etaS0);
    etaY0 = etaY; etaS0 = etaS;
      
    % Applying two-step Adams-Bashforth
    xc2  = xc1  + tt.dt*(1.5*Up1 - 0.5*Up0 );
    tau2 = tau1 + tt.dt*(1.5*wp1 - 0.5*wp0 );
      
    xc0 = xc1; tau0 = tau1; Up0 = Up1; wp0 = wp1;
    xc1 = xc2; tau1 = tau2;
      
    % update geometry
    geom2 = capsules(prams,xc2,tau2);
  end

  om.plotData(geom2);
    
  % update time
  time = time + tt.dt;

  % write current time to console
  message = ['Completed time ', num2str(time,'%4.2e'), ...
      ' of time horizon ' , num2str(prams.T,'%4.2e')];
  om.writeMessage(message);
  
  % write the shape to a file
  om.writeData(time,geom2.center,geom2.tau,geom2.X);

  % write the velocity to a file
  om.writeVelData(time,Up0,wp0);  
  
  % update step counter
  step = step + 1;
end

% save final time step
Xfinal = geom2.X;

om.writeStars();
message = ['Finished entire simulation in ', ...
      num2str(toc(ttotal),'%4.2e'), ' seconds'];
om.writeMessage(message)


end %rigid2D
