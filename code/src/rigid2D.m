function [Xfinal, trajectory] = rigid2D(options,prams,xc,tau)
ttotal = tic; % start a timer

om = monitor(options,prams);
% build object for doing I/O

tt = tstep(options,prams);

% build initial condition of rigid body particle
geom = capsules(prams,xc,tau);
% plot geometry if usePlot == true
om.plotData(geom);
drawnow

time = 0;
step = 0;

xc0  = xc;
tau0 = tau;

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
  if 1 
      geom = capsules(prams,xc,tau);
      [Up, wp] = tt.timeStep(geom);
      xc  = xc  + tt.dt*Up;
      tau = tau + tt.dt*wp;
      geom2 = geom;
  end
  
  %Adams-Bashforth 
  if 0 
  
      if step == 0
          
          % update geometry
          geom0 = capsules(prams,xc0,tau0);
          [Up0, wp0] = tt.timeStep(geom0);
          
          % causes xc2, tau2 to be forward Euler, at step 0
          xc1  = xc0;
          tau1 = tau0;
          
      end
      
      geom1 = capsules(prams,xc1,tau1);
      [Up1, wp1] = tt.timeStep(geom1);
      
      % Applying two-step Adams-Bashforth
      xc2  = xc1  + tt.dt*( 1.5*Up1 - 0.5*Up0 );
      tau2 = tau1 + tt.dt*( 1.5*wp1 - 0.5*wp0 );
      
      xc0 = xc1; tau0 = tau1; Up0 = Up1; wp0 = wp1;
      xc1 = xc2; tau1 = tau2;
      
      % update geometry
      geom2 = capsules(prams,xc2,tau2);
      
  end
  
  
  
  
  % update time
  time = time + tt.dt;
  
  % om.plotField(geom2,Unum,Xtest,Ytest);
  hold off
  om.plotData(geom2);
  rhs = reshape(yukawaRHS(geom2), geom2.N, geom2.nb);
  hold on
  x1   = geom2.X(1:geom2.N,:);         % grid points on curves 
  x2   = geom2.X(geom2.N+1:2*geom2.N,:);            

  for p = 1:geom2.nb
    z = rhs(:,p);      %# colors
    h = surface([x1(:,p), x1(:,p)], [x2(:,p), x2(:,p)], [z(:), z(:)], ...
    [z(:), z(:)], 'EdgeColor','flat', 'FaceColor','none');
    colormap( jet )
  end
  drawnow
  
  % write the shape to a file
    
  step = step + 1;
  %[step prams.m time prams.T]
  
end

% save final time step
Xfinal = geom2.X;

om.writeStars();
message = ['Finished entire simulation in ', ...
      num2str(toc(ttotal),'%4.2e'), ' seconds'];
om.writeMessage(message)


end %rigid2D
