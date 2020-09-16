function Xfinal = rigid2D(options,prams,xc,tau)
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
% begin time loop
while time < prams.T
  % start a timer for the this particular time step
  tSingleStep = tic;   

  % comupute density function, translational velocity, angular velocity,
  % and the GMRES output
%  [eta,eta2,Up,wp,iter,iter2,iflag,iflag2,res,...
%      res2,Unum,Xtest,Ytest,force,torque] = tt.timeStep(geom);
  [Up,wp,iterYukawa,iterStokes] = tt.timeStep(geom);
  msg = ['Yukwa reguired  ' num2str(iterYukawa) ...
    ' iterations'];
  om.writeMessage(msg);
  msg = ['Stokes reguired ' num2str(iterStokes) ...
    ' iterations'];
  om.writeMessage(msg);
  
%  om.plotField(geom,Unum,Xtest,Ytest);  

  % update centers and angles with forward Euler
  xc = xc + tt.dt*Up;
  tau = tau + tt.dt*wp;
  geom = capsules(prams,xc,tau);
  % update geometry

%   om.plotData(geom);
  % plot geometry

  % update time
  time = time + tt.dt;    
  % update the shape
   
  % write the shape to a file
  om.writeData(time,xc,tau,geom.X);
end

% save final time step
Xfinal = geom.X;

om.writeStars();
message = ['Finished entire simulation in ', ...
      num2str(toc(ttotal),'%4.2e'), ' seconds'];
om.writeMessage(message)


end %rigid2D
