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

time = tt.dt*tt.sstep;
step0 = 0;

step = step0;

xc0  = xc;
tau0 = tau;
om.writeData(time,geom.center,geom.tau,geom.X);

trajectory = [xc(1,:) xc(2,:) tau];

% begin time loop
%while time < prams.T
while step <= prams.m
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
    if step == step0
      % update geometry
      geom0 = capsules(prams,xc0,tau0);
      [Up0, wp0,~,~,etaY0,etaS0,force,torque] = tt.timeStep(geom0,geom0.X,geom0.X);

      % write the velocity to a file
      om.writeVelData(time,Up0,wp0);
      
      % causes xc2, tau2 to be forward Euler, at step 0
      xc1  = xc0;
      tau1 = tau0;

      % update time
      time = time + tt.dt;

      % write current time to console
      %message = ['Completed time ', num2str(time,'%4.2e'), ...
      %           ' of time horizon ' , num2str(prams.T,'%4.2e')];
      %om.writeMessage(message);
      
      % write the shape to a file
%       om.writeData(time,geom0.center,geom0.tau,geom0.X);

      % write the velocity to a file
%       om.writeVelData(time,Up0,wp0);
      
      % update step counter
      step = step + 1;
    end
      
    geom1 = capsules(prams,xc1,tau1);
    [Up1, wp1,~,~,etaY,etaS,force,torque] = tt.timeStep(geom1,etaY0,etaS0);
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

  % write the shape to a file
  om.writeData(time,geom2.center,geom2.tau,geom2.X);

  % write the velocity to a file
  om.writeVelData(time,Up0,wp0); 
  
  
  % update time
  time = time + tt.dt;

  % write current time to console
  message = ['Completed time ', num2str(time,'%4.2e'), ...
      ' of time horizon ' , num2str(prams.T,'%4.2e')];
  om.writeMessage(message);
  

  
  % update tracers
  if options.tracer
  XX = load("../examples/tracers.dat");
  UU = load("../examples/tracer_vel.dat");
   
  XX = XX + tt.dt*UU;
  
  % pacman correction for points wandering out-of-bounds
  il = find(XX(:,1) < options.plotAxis(1)); 
  ir = find(XX(:,1) > options.plotAxis(2)); 
  ib = find(XX(:,2) < options.plotAxis(3)); 
  it = find(XX(:,2) > options.plotAxis(4)); 
  
  XX(il,1) = options.plotAxis(2);
  XX(ir,1) = options.plotAxis(1);
  XX(ib,2) = options.plotAxis(4);
  XX(it,2) = options.plotAxis(3);
  
  save("-ascii", "../examples/tracers.dat", "XX");
  end
  
  % output data
  
  DATA = [xc1(1,:)' xc1(2,:)' tau1' force(1,:)' force(2,:)' torque];
  fileName = sprintf("../output/data/frames/N%d_%f_%d.dat", geom2.nb, options.shearRate, step+tt.sstep);
  save("-ascii", fileName, "DATA");
 
  if options.tracer
  fileName = sprintf("../output/data/frames/N%d_%f_%d.tracer", geom2.nb, options.shearRate, step+tt.sstep);
  save("-ascii", fileName, "XX");  
  end
  
  
  % stress calculations
  if geom2.nb==58
       xo = xc1(1,27:end);
       yo = xc1(2,27:end);
       dist = sqrt(diff(xo).^2+diff(yo).^2);
       tmdist = cumsum([0 dist]);
       tInt = linspace(0,tmdist(end));
       xInt = interp1(tmdist,xo,tInt,"spline");
       yInt = interp1(tmdist,yo,tInt,"spline");
       Xtar = [xInt;yInt];
       Ntar = length(xInt);
[stress, pressure, velocity] = fluidstress(geom1,etaS0,force,torque,Xtar);
%   
SPV = [stress(1:Ntar) stress(Ntar+1:2*Ntar) stress(2*Ntar+1:end) ...
       pressure velocity(1:Ntar) velocity(Ntar+1:end)];
% size(stress)   % 3*Ntar
% size(pressure) % 1*Ntar
% size(velocity) % 2*Ntar  

fileName = sprintf("../output/data/frames/N%d_%f_%d.stress", geom2.nb, options.shearRate, step+tt.sstep);
  save("-ascii", fileName, "SPV");  
  end
  
  
  
  
  
  
  
  
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
