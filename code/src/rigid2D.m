function Xfinal = rigid2D(options,prams,xc,tau)
ttotal = tic; % start a timer

om = monitor(options,prams);
% build object for doing I/O

tt = tstep(options,prams);

% build initial condition of rigid body particle
geom = capsules(prams,xc,tau);
% plot geometry if usePlot == true
om.plotData(geom);

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
    
  % write the CPU time and the output of GMRES in Yukawa
%  om.writeMessage(....
%    ['Yukawa Finished t=', num2str(time, '%4.2e'), ' in ' ...
%        num2str(iter2) ' iterations after ', ...
%        num2str(toc(tSingleStep), '%4.2e'), ' seconds (residual ', ...
%        num2str(res2,'%4.2e'), ')']);  
  
  % write the CPU time and the output of GMRES in mobility problem
%  om.writeMessage(....
%    ['Mobility Finished t=', num2str(time, '%4.2e'), ' in ' ...
%        num2str(iter) ' iterations after ', ...
%        num2str(toc(tSingleStep), '%4.2e'), ' seconds (residual ', ...
%        num2str(res,'%4.2e'), ')']);
    
    % plot the shape
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
