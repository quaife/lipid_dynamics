classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Routines that we may need to add are different
% integrators such as semi-implicit, SDC, Runge-Kutta, adaptivity

properties

dt               % time step size
sstep            % starting step
gmresTol         % GMRES tolerance
farField         % background flow
janusbc          % particle boundary condition
tracer           % flag for tracer
confined         % flag to determine if geometry is confined or not
precoYukawa      % block-diagonal preconditioner for Yukawa on bodies
precoStokes      % block-diagonal preconditioner for Stokes on bodies
precoYukawaWalls % block-diagonal preconditioner for Yukawa on walls
precoStokesWalls % block-diagonal preconditioner for Stokes on walls

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options,prams)
% o.tstep(options,prams): constructor.  Initialize class.  Take all
% elements of options and prams needed by the time stepper

o.dt       = prams.T/prams.m;
o.sstep    = prams.sstep;
o.gmresTol = options.gmresTol;
o.farField = @(X) o.bgFlow(X,options);
o.tracer   = options.tracer;
o.confined = options.confined;

% for screen laplace BVP
o.janusbc  = @(X,tau,center) o.bcfunc(X,tau,center,options);

end % constructor: tstep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Up,wp,iterYukawa,iterStokes,etaYukawaJanus,etaStokesJanus,...
      fforce,torque,Energy] = timeStep(o,geom,etaY0,etaS0,walls)
% Main time stepping routine
oc     = curve;
N      = geom.N;
Nwall  = walls.N;
nb     = geom.nb;
nbwall = walls.nb;
X      = geom.X;
tau    = geom.tau;
center = geom.center;
rho    = geom.rho;
gam    = geom.gam;
radii  = geom.radii;

% CREATE NEAR SINGULAR INTEGRATION STRUCTURES
if ~o.confined
  geom.nearStructB2B = geom.getZone([],1);
else
  % contribution beween Janus particles and themselves, and between
  % particles and the solid walls
  [geom.nearStructB2B,geom.nearStructB2T] = ...
        geom.getZone(walls,3);
  % contribution between solid walls and Janus particles
  [~,walls.nearStructB2T] = walls.getZone(geom,3);
end

op = poten(geom.N,geom.rho);
% START OF REPULSION BETWEEN JANUS PARTICLES CALCULATION
[RE, R1, R2, RTq, pp] = op.Repul(geom);
% END OF REPULSION BETWEEN JANUS PARTICLES CALCULATION

% START OF SCREENED LAPLACE SOLVE TO FIND HYDROPHOBIC FORCE AND TORQUE
% right hand side for the screened Laplace solver
yukawaRHS = geom.yukawaRHS;
% append zeros for the boudnary condition for the Yukawa equation
yukawaRHS = [yukawaRHS;zeros(Nwall*nbwall,1)];

% build the DLP Yukawa matrix due to the Janus particles
geom.DLPYukawa = op.yukawaDLmatrix(geom); 
o.precoYukawa.L = zeros(N,N,nb);
o.precoYukawa.U = zeros(N,N,nb);
for k = 1:nb
  [o.precoYukawa.L(:,:,k),o.precoYukawa.U(:,:,k)] = ...
    lu(0.5*eye(N) + geom.DLPYukawa(:,:,k));
end
% build the DLP Yukawa matrix due to the solid walls if confined == true
if o.confined
  op2 = poten(Nwall,walls.rho);
  walls.DLPYukawa = op2.yukawaDLmatrix(walls); 
  o.precoYukawaWalls.L = zeros(Nwall,Nwall,nbwall);
  o.precoYukawaWalls.U = zeros(Nwall,Nwall,nbjwall);
  for k = 1:nbwall
    [o.precoYukawaWalls.L(:,:,k),o.precoYukawaWalls.U(:,:,k)] = ...
      lu(0.5*eye(Nwall) + walls.DLPYukawa(:,:,k));
  end
end

% build precomputed matrix with all necessary bessel functions that are
% needed at every GMRES iteration
geom.YukawaKernelMatrix;

% Solve for the density function using GMRES
[sigma,iflagYukawa,resYukawa,iterYukawa] = gmres(...
      @(etaY0) o.timeMatVecYukawa(etaY0,geom,walls),...
      yukawaRHS, [], o.gmresTol, N*nb + Nwall*nbwall,...
      @(etaY0) o.precoYukawaBD(etaY0)); 
iterYukawa = iterYukawa(2);
fprintf('Yukawa required %i iterations\n',iterYukawa);

% Unstack the Yukawa density function so that it is arranged as columns
% for each body
etaYukawaJanus = zeros(N,nb);
etaYukawaWalls = zeros(Nwall,nbwall);
for k = 1:nb
  etaYukawaJanus(:,k) = sigma((k-1)*N+1:k*N);
end
for k = 1:nbwall
  etaYukawaWalls(:,k) = sigma(nb*N + (k-1)*Nwall+1:nb*N + k*Nwall);
end

if 0
  [xtar,ytar] = meshgrid(-7:1.01:+7,-7:1.01:7);
  [nx,ny] = size(xtar);
  targets.N = numel(xtar);
  targets.nb = 1;
  targets.X = [xtar(:);ytar(:)];
  [~,NearJanusTargets] = geom.getZone(targets,2);

  kernel = @op.exactYukawaSL;
  kernelDirect = @op.exactYukawaSL;
  kernelSelf = @(z) op.exactYukawaSLdiag(geom,...
        S,z);

  potJanus = op.nearSingInt(geom,etaYukawaJanus,kernelSelf,...
    NearJanusTargets,kernel,kernelDirect,targets,false,false);
  potJanus = potJanus(1:targets.N,:);
end

% Finite difference way of computing energy
%Energy = geom.computeEnergy(etaYukawaJanus);

% generate Yukawa single-layer potential matrix for
% self-contribution of each individual Janus particle
geom.SLPYukawa = op.yukawaSLmatrix(geom); 
% Layer potential way of computing energy
Energy = geom.computeEnergyNew(etaYukawaJanus);

%[Energy Energy2]
%abs(Energy - Energy2)/Energy2
%pause

%Energy = geom.computeEnergy(etaYukawaJanus) + RE;

if 0
  [xtar,ytar] = meshgrid(-9:0.01:-7,-1:0.01:1);
  [nx,ny] = size(xtar);
  ind = [];
  for k = 1:geom.nb
    ind = [ind;find((xtar - geom.center(1,k)).^2 + ...
          (ytar - geom.center(2,k)).^2 < 1.0001*geom.radii(k).^2)];
  end
  targets.N = numel(xtar);
  targets.nb = 1;
  targets.X = [xtar(:);ytar(:)];
  [~,NearJanusTargets] = geom.getZone(targets,2);
  [~,NearWallsTargets] = walls.getZone(targets,2);

  kernel = @op.exactYukawaDL;
  kernelDirect = @op.exactYukawaDL;
  kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,...
        geom.DLPYukawa,z);

  potJanus = op.nearSingInt(geom,etaYukawaJanus,kernelSelf,...
    NearJanusTargets,kernel,kernelDirect,targets,false,false);
  potJanus = potJanus(1:targets.N,:);

  potJanus = reshape(potJanus,nx,ny);

  kernel = @op2.exactYukawaDL;
  kernelDirect = @op2.exactYukawaDL;
  kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(walls,walls.DLPYukawa,z);
  potWalls = op2.nearSingInt(walls,etaYukawaWalls,kernelSelf,NearWallsTargets,...
      kernel,kernelDirect,targets,false,false);
  potWalls = potWalls(1:targets.N,:);
  potWalls = reshape(potWalls,nx,ny);

  pot = potWalls + potJanus;
  pot(ind) = 0;
  clf;
  surf(xtar,ytar,pot)
  view(2); axis equal;
  shading interp;
  colorbar
  hold on;
%  plot3(geom.X(1:end/2),geom.X(end/2+1:end),geom.yukawaRHS,'k','linewidth',2)
  for k = 1:geom.nb
    fill3(geom.X(1:end/2,k),geom.X(end/2+1:end,k),10*ones(geom.N,1),'k')
  end
  plot(walls.X(1:end/2),walls.X(end/2+1:end),'k');
  axis([min(xtar(:)) max(xtar(:)) min(ytar(:)) max(ytar(:))]);
  caxis([0 max(yukawaRHS)]);
  pause
end


% Compute force and torques due to contribution of Janus particles and
% solid walls using QBX and Tpq + Tqp identity.
[F1,F2,Tq] = op.evalHydroForces(geom,etaYukawaJanus,walls,etaYukawaWalls);
% lift Tq so that it is for rotation about center (not origin)

% tst = fopen("HAP_force_distance.dat", "r");
% if tst == -1
%   fid = fopen("HAP_force_distance.dat", "w");
% else 
%   fclose(tst);
%   fid = fopen("HAP_force_distance.dat", "a");
% end
% fprintf(fid,"%d %d\n", [geom.center(1,2) - geom.center(1,1) - 2,  F1(1)]');
% fclose(fid);
%fprintf("%d %d\n", [geom.center(1,2) - geom.center(1,1) - 2,  F1(1)]');

if 0
  th = linspace(0, 2*pi)'; 
  hold off
  plot(center(1,:), center(2,:), 'ob');
  hold on
  plot(center(1,:) + radii.*cos(th), center(2,:) + radii.*sin(th), 'b');
  l0 = geom.RepulLength;
  plot(center(1,:) + (radii+l0).*cos(th), center(2,:) + (radii+l0).*sin(th), 'r:');
  quiver(center(1,:), center(2,:), R1', R2', 0, 'r');
  % quiver(center(1,:), center(2,:), F1', F2', 'm');
  if pp
    plot([pp(:,1) pp(:,3)], [pp(:,2) pp(:,4)], 'k*');
  end

  axis equal
  plot(walls.X(1:end/2),walls.X(end/2+1:end),'k','linewidth',2)
  axis([-16 -10 -3.1 3.1])
  pause
  pause(0.01)
end

%outputs: 
fforce  = [F1 + R1, F2 + R2].';
force  = fforce(:);
torque = Tq + RTq;
% END OF SCREENED LAPLACE SOLVE TO FIND HYDROPHOBIC FORCE AND TORQUE

% START OF MOBILITY PROBLEM
% far field condition plus the stokeslet and rotlet terms
if ~o.confined
  rhsJanus = o.farField(X) + op.StokesletRotlet(geom,force,torque);
  rhsWalls = [];
else
  % solve stokes in bounded domain
  rhsJanus = op.StokesletRotlet(geom,force,torque);
  [~,rhsWalls] = op.StokesletRotlet(geom,force,torque,walls.X,(1:nb));
  rhsWalls = rhsWalls - o.farField(walls.X);
end

% append body force and torque to the background flow to complete the
% right hand side for GMRES
rhs = [-rhsJanus(:); force; torque; -rhsWalls(:)];

% build double-layer potential matrix
geom.DLPStokes = op.stokesDLmatrix(geom);
if o.confined
  walls.DLPStokes = op2.stokesDLmatrix(walls);
  walls.N0Stokes = op2.stokesN0matrix(walls);
end

% Build LU factorization of the completed Stokes DLP system of equations
o.precoStokesMatrix(geom);

% max GMRES iterations
maxit = 2*N*nb + 2*Nwall*nbwall; 

% Build the upsampled matrix for doing Stokes Layer potentials
geom.StokesDLMatrix;

[sigma, iflagStokes, resStokes, iterStokes] = gmres(...
      @(etaS0) o.timeMatVecStokes(etaS0,geom,walls),...
      rhs,[],o.gmresTol,maxit,...
      @(etaS0) o.precoStokesBD(etaS0,walls));
iterStokes = iterStokes(2);
fprintf('Stokes required %i iterations\n',iterStokes);

% REORGANIZE COLUMN VECTOR INTO MATRIX EXTRACT DENSITY FUNCITONS ON
% FIBRES AND WALLS
% each column of etaJanus corresponds to the density function of a
% Janus particle
etaStokesJanus = zeros(2*N,nb);
for k = 1:nb
  etaStokesJanus(:,k) = sigma((k-1)*2*N+1:k*2*N);
end

% EXTRACT TRANSLATIONAL VELOCITIES
Up = zeros(2,nb);
for k = 1:nb
  Up(:,k) = sigma(2*N*nb+(k-1)*2+1:2*N*nb+2*k);
end

% EXTRACT ROTATIONAL VELOCITIES
wp = zeros(1,nb);
for k = 1:nb
  wp(k) = sigma(2*N*nb+2*nb+k);
end

if o.confined
  etaStokesWalls = sigma(2*N*nb + 3*nb + 1:end);
end


if 0
  [xtar,ytar] = meshgrid(-1.5:.1:4,2.999:.0001:2.9999);
  [nx,ny] = size(xtar);
  ind = [];
  for k = 1:geom.nb
    ind = [ind;find((xtar - geom.center(1,k)).^2 + ...
          (ytar - geom.center(2,k)).^2 < 1.0001*geom.radii(k)^2)];
  end
  targets.N = numel(xtar);
  targets.nb = 1;
  targets.X = [xtar(:);ytar(:)];
  [~,NearJanusTargets] = geom.getZone(targets,2);
  [~,NearWallsTargets] = walls.getZone(targets,2);

  kernel = @op.exactStokesDL;
  kernelDirect = @op.exactStokesDL;
  kernelSelf = @(z) +0.5*z + op.exactStokesDLdiag(geom,...
        geom.DLPStokes,z);

  StokesJanus = op.nearSingInt(geom,etaStokesJanus,kernelSelf,...
    NearJanusTargets,kernel,kernelDirect,targets,false,false);

  StokesJanusx = reshape(StokesJanus(1:end/2),nx,ny);
  StokesJanusy = reshape(StokesJanus(end/2+1:end),nx,ny);

  kernel = @op2.exactStokesDL;
  kernelDirect = @op2.exactStokesDL;
  kernelSelf = @(z) +0.5*z + op.exactStokesDLdiag(walls,walls.DLPStokes,z);
  StokesWalls = op2.nearSingInt(walls,etaStokesWalls,kernelSelf,...
      NearWallsTargets,kernel,kernelDirect,targets,false,false);
  StokesWallsx = reshape(StokesWalls(1:targets.N,:),nx,ny);
  StokesWallsy = reshape(StokesWalls(targets.N+1:end,:),nx,ny);


  [~,StokesLets] = op.StokesletRotlet(geom,force,torque,targets.X,(1:nb));
  StokesLetsx = reshape(StokesLets(1:end/2),nx,ny);
  StokesLetsy = reshape(StokesLets(end/2+1:end),nx,ny);

  Stokesx = StokesJanusx + StokesWallsx + StokesLetsx;
  Stokesy = StokesJanusy + StokesWallsy + StokesLetsy;
  Stokesx(ind) = 0;
  Stokesy(ind) = 0;
  clf;
%  quiver(xtar,ytar,Stokesx,Stokesy)
  surf(xtar,ytar,log10(sqrt(Stokesx.^2 + Stokesy.^2)))
  shading interp; view(2); colorbar
%  axis equal;
  hold on;
%  for k = 1:geom.nb
%    fill3(geom.X(1:end/2,k),geom.X(end/2+1:end,k),10*ones(geom.N,1),'k')
%  end
  plot(walls.X(1:end/2),walls.X(end/2+1:end),'k');
  axis([-1.5 4 2.999 3]);
  pause
end


%[mean(Up(1,:)); mean(Up(2,:))]
%Up = Up - [mean(Up(1,:)); mean(Up(2,:))];

% Routine to compute the velocity in the bulk and use the rigid body
% motion to define a velocity in rigid bodies. The goal is to see
% something continuous
if o.tracer
  op.bulkVelocity(geom,etaStokesJanus,Up,wp,force,torque,...
      @(X) o.farField(X));
end

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVecStokes(o,Xn,geom,walls)
% Tx = timeMatVecStokes(Xn,geom,walls) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% janus1: eta_x, eta_y janus2: eta_x, eta_y ... janusnv: eta_x, eta_y
% janus1: u, v janus2: u, v ... janusnv: u, v
% janus: omega janus: omega ... janusnv: omega
% walls: eta_x, eta_y

N = geom.N;   % points per body
nb = geom.nb; % number of bodies
rho = geom.rho; % screen length
if o.confined
  Nwall = walls.N; % number of points on outer wall
  nbwall = walls.nb; % number of components of outer wall
else
  Nwall = 0; nbwall = 0;
end

op = poten(N,rho);
op2 = poten(Nwall,rho);

% Output of Tx that corresponds to the velocity of the Janus particles 
valJanus = zeros(2*N,nb);
% output of Tx that corresponds to force on fibres
valForce = zeros(2,nb);
% output of Tx that corresponds to torque on fibres
valTorque = zeros(1,nb);
% output of Tx that corresponds to the velocity on the solid wall
valWalls = zeros(2*Nwall,nbwall);

% BEGIN FORMATTING UNKNOWN VECTOR
% EXTRACT DENSITY FUNCTION FOR FIBRES
etaJanus = zeros(2*N,nb);
for k = 1:nb
  etaJanus(:,k) = Xn((k-1)*2*N+1:k*2*N);
end

% EXTRACT TRANSLATIONAL AND ROTATIONAL VELOCITIES OF FIBRES
Up = zeros(2,nb);
for k = 1:nb
  Up(:,k) = Xn(2*N*nb+1+2*(k-1):2*N*nb+2*k);
end
wp = zeros(1,nb);
for k = 1:nb
  wp(k) = Xn(2*N*nb+2*nb+k);
end

% EXTRACT DENSITY FUNCTION FOR SOLID WALLS
etaWalls = zeros(2*Nwall,nbwall);
for k = 1:nbwall
  istart = 2*N*nb + 3*nb + 1;
  iend = istart + 2*Nwall - 1;
  etaWalls(:,k) = Xn(istart:iend);
end
% END FORMATTING UNKNOWN VECTOR

% ADD JUMP IN DLP
jump = +0.5;
valJanus = valJanus + jump*etaJanus;

% START OF CONTRIBUTIONS EVALUATED ON JANUS PARTICLES
% ADD SELF CONTRIBUTION
valJanus = valJanus + ...
    op.exactStokesDLdiag(geom,geom.DLPStokes,etaJanus);

% Matrix-free version
%% DEFINE STOKES DLP KERNELS
%kernel = @op.exactStokesDL;
%kernelSelf = @(z) +jump*z + op.exactStokesDLdiag(geom,geom.DLPStokes,z);
%
%% ADD CONTRIBUTION FROM OTHER JANUS PARTICLES
%stokesDLP = op.nearSingInt(geom,etaJanus,kernelSelf,...
%    geom.nearStructB2B,kernel,kernel,geom,true,false);

% Precomputing and using matrix (memory-intensive)
kernel = @op.exactStokesDLMatFree;
kernelDirect = @op.exactStokesDL;
kernelSelf = @(z) jump*z + op.exactStokesDLdiag(geom,geom.DLPStokes,z);

stokesDLP = op.nearSingInt(geom,etaJanus,kernelSelf,...
    geom.nearStructB2B,kernel,kernelDirect,geom,true,false);

valJanus = valJanus + stokesDLP;

% ADD TRANSLATIONAL VELOCITY CONTRIBUTION
for k = 1:nb
  valJanus(1:end/2,k) = valJanus(1:end/2,k) - Up(1,k);
  valJanus(end/2+1:end,k) = valJanus(end/2+1:end,k) - Up(2,k);
end

% ADD ROTATIONAL VELOCITY CONTRIBUTION
for k = 1:nb
  valJanus(1:end/2,k) = valJanus(1:end/2,k) - ...
                (-(geom.X(end/2+1:end,k) - geom.center(2,k)))*wp(k);
  valJanus(end/2+1:end,k) = valJanus(end/2+1:end,k) - ...
                (+(geom.X(1:end/2,k) - geom.center(1,k)))*wp(k);
end

if o.confined
  kernel = @op2.exactStokesDL;
  kernelSelf = @(z) +jump*z + op2.exactStokesDLdiag(walls,walls.DLPStokes,z);

  % ADD CONTRIBUTION FROM SOLID WALLS
  stokesDLPwalls = op2.nearSingInt(walls,etaWalls,kernelSelf,...
      walls.nearStructB2T,kernel,kernel,geom,false,false);

  valJanus = valJanus + stokesDLPwalls;
end


% EVALUTATE FORCES ON JANUS PARTICLES
for k = 1:nb
  valForce(1,k) = 1/2/pi*sum(etaJanus(1:N,k).*geom.sa(:,k))*2*pi/N;
  valForce(2,k) = 1/2/pi*sum(etaJanus(N+1:2*N,k).*geom.sa(:,k))*2*pi/N;
end

% EVALUATE TORQUES ON JANUS PARTICLES
for k = 1:nb
  valTorque(k) = 1/2/pi*sum((-...
      (geom.X(N+1:2*N,k)-geom.center(2,k)).*etaJanus(1:N,k) + ...
      (geom.X(1:N,k)-geom.center(1,k)).*etaJanus(N+1:2*N,k)).*...
                     geom.sa(:,k))*2*pi/N;
end

% START OF CONTRIBUTIONS EVALUATED ON SOLID WALLS
if o.confined
  % Jump condition due to solid walls
  valWalls = valWalls + jump*etaWalls;
  valWalls = valWalls + ...
    op2.exactStokesDLdiag(walls,walls.DLPStokes,etaWalls);
  % add in rank one modification so that null space due to
  % incompressibility condition is removed.
  valWalls = valWalls + ...
    op.exactStokesDLdiag(walls,walls.N0Stokes,etaWalls);

  kernel = @op.exactStokesDL;
  kernelSelf = @(z) +jump*z + op.exactStokesDLdiag(geom,geom.DLPStokes,z);

  % ADD CONTRIBUTION FROM JANUS PARTICLES
  stokesDLPJanus = op.nearSingInt(geom,etaJanus,kernelSelf,...
      geom.nearStructB2T,kernel,kernel,walls,false,false);

  valWalls = valWalls + stokesDLPJanus;
end


% CONSTRUCT OUTPUT VECTOR
Tx = [valJanus(:);+valForce(:);+valTorque(:);valWalls(:)];

end % timeMatVecStokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVecYukawa(o,Xn,geom,walls)
% Tx = timeMatVecYukawa(Xn,geom) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% fiber1: eta1 fiber2: eta2 ... fibernv: etan

N = geom.N;
nb = geom.nb;
rho = geom.rho;
op = poten(N,rho);

Nwall = walls.N;
nbwall = walls.nb;
op2 = poten(Nwall,rho);

% Store with each column being one part of the matvec. Will columnize
% the vector at the end of this routine
TxJanus = zeros(N,nb);
TxWalls = zeros(Nwall,nbwall);

etaJanus = zeros(N,nb);
etaWalls = zeros(Nwall,nbwall);
% unstack the density function
for k = 1:nb
  etaJanus(:,k) = Xn((k-1)*N+1:k*N);
end
for k = 1:nbwall
  etaWalls(:,k) = Xn(N*nb + (k-1)*Nwall+1:N*nb + k*Nwall);
end

% JUMP IN DLP ON JANUS PARTICLES
TxJanus = TxJanus + 0.5*etaJanus;

% JANUS PARTICLE SELF CONTRIBUTION
TxJanus = TxJanus + op.exactYukawaDLdiag(geom,geom.DLPYukawa,etaJanus);

% JANUS PARTICLE TO JANUS PARTICLE CONTRIBUTION (NON-SELF)
% Matrix-free version
%% DEFINE Yukawa DLP KERNELS
%kernel = @op.exactYukawaDL;
%kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,geom.DLPYukawa,z);
%
%yukawaDLP = op.nearSingInt(geom,etaJanus,kernelSelf,geom.nearStructB2B,...
%    kernel,kernel,geom,true,false);

% Precomputing and using matrix (memory-intensive)
kernel = @op.exactYukawaDLMatFree;
kernelDirect = @op.exactYukawaDL;
kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,geom.DLPYukawa,z);

yukawaDLP = op.nearSingInt(geom,etaJanus,kernelSelf,geom.nearStructB2B,...
    kernel,kernelDirect,geom,true,false);
TxJanus = TxJanus + yukawaDLP(1:geom.N,:);

% CONTRIBUTION DUE TO SOLID WALLS EVALUATED ON JANUS PARTICLES
% DEFINE Yukawa DLP KERNELS
if o.confined
  kernel = @op2.exactYukawaDL;
  kernelSelf = @(z) +0.5*z + op2.exactYukawaDLdiag(walls,walls.DLPYukawa,z);
  yukawaDLP = op2.nearSingInt(walls,etaWalls,kernelSelf,walls.nearStructB2T,...
      kernel,kernel,geom,false,false);

  TxJanus = TxJanus + yukawaDLP(1:geom.N,:);
end

if o.confined
  % JUMP IN DLP ON SOLID WALL
  TxWalls = TxWalls + 0.5*etaWalls;

  % SOLID WALL SELF CONTRIBUTION
  TxWalls = TxWalls + op2.exactYukawaDLdiag(walls,walls.DLPYukawa,etaWalls);

  % CONTRIBUTION DUE TO JANUS PARTICLES EVALUATED ON SOLID WALLS
  % DEFINE Yukawa DLP KERNELS
  kernel = @op.exactYukawaDL;
  kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,geom.DLPYukawa,z);
  yukawaDLP = op.nearSingInt(geom,etaJanus,kernelSelf,geom.nearStructB2T,...
      kernel,kernel,walls,false,false);
  TxWalls = TxWalls + yukawaDLP(1:walls.N,:);
else
  TxWalls = [];
end

Tx = [TxJanus(:);TxWalls(:)];

end % timeMatVecYukawa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vInf = bgFlow(~,X,options)
    
N = size(X,1)/2;
nb  = size(X,2);
oc = curve;

[x,y] = oc.getXY(X);
shearRate = options.shearRate; 
% shear rate

switch options.farField
  case 'shear'
    vInf = [shearRate*y;zeros(N,nb)];
  case 'extensional'
    vInf = shearRate*[-x;y];
  case 'taylorgreen'
    vInf = shearRate*[-cos(x/2).*sin(y/2);sin(x/2).*cos(y/2)];
  case 'parabolic'
    vInf = [shearRate*(1-(y/80).^2);zeros(N,nb)];
  case 'channel'
    vInf = zeros(2*N,nb);
    vx = (1 - (y/max(y)).^2).*(abs(x) > 8);
    vy = zeros(N,nb);
    vInf = shearRate*[vx;vy];
  otherwise
    vInf = zeros(2*N,nb);
end
    
end % bgFlow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bc = bcfunc(~,X,tau,center,options)
    
N = size(X,1)/2;
nb  = size(X,2);
oc = curve;

[x,y] = oc.getXY(X);
bc = zeros(N,nb);

pow = options.janusbc;
for i = 1:nb    
  xc = center(:,i);
  th = atan2(y(:,i)-xc(2),x(:,i)-xc(1));
  bc(:,i) = cos(0.5*(th-tau(i))).^pow;         
end


end % bcfunc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = precoYukawaBD(o,eta); 

N = size(o.precoYukawa.U,1);
nb = size(o.precoYukawa.U,3);
if o.confined
  Nwall = size(o.precoYukawaWalls.U,1);
  nbwall = size(o.precoYukawaWalls.U,3);
else
  Nwall = 0;
  nbwall = 0;
end
z = zeros(N*nb + Nwall*nbwall,1);
for k = 1:nb
  istart = (k-1)*N + 1;
  iend = istart + N - 1;
  z(istart:iend) = o.precoYukawa.U(:,:,k)\...
    (o.precoYukawa.L(:,:,k)\eta(istart:iend));
end

for k = 1:nbwall
  istart = N*nb + (k-1)*Nwall + 1;
  iend = istart + Nwall - 1;
  z(istart:iend) = o.precoYukawaWalls.U(:,:,k)\...
    (o.precoYukawaWalls.L(:,:,k)\eta(istart:iend));
end



end % precoYukawaBD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function precoStokesMatrix(o,geom)
% Build the LU decomposition of the Stokes linear system. It is actually
% applied in the function precoStokesBD

N = geom.N;
nb = geom.nb;
o.precoStokes.L = zeros(2*N + 3, 2*N + 3,nb);
o.precoStokes.U = zeros(2*N + 3, 2*N + 3,nb);

M11 = zeros(2*N,2*N,nb);
M12 = zeros(2*N,3,nb);
M21 = zeros(3,2*N,nb);
% allocate space for 3 of the 4 block matricies coming from the Power
% and Miranda completed double-layer potential
jump = +0.5;
oc = curve;
for k = 1:nb
  [x,y] = oc.getXY(geom.X(:,k));
  [cx,cy] = oc.getXY(geom.center(:,k));
  sa = geom.sa(:,k)*2*pi/N;
  M11(:,:,k) = jump*eye(2*N) + geom.DLPStokes(:,:,k);

  M12(:,:,k) = [[-ones(N,1);zeros(N,1)] [zeros(N,1);-ones(N,1)] ...
      [(y - cy);-(x-cx)]];
  M21(:,:,k) = 1/2/pi*[[sa' zeros(1,N)];...
                      [zeros(1,N) sa'];...
                      [(y - cy)'.*sa' -(x - cx)'.*sa']];
end

for k = 1:nb
  [o.precoStokes.L(:,:,k),o.precoStokes.U(:,:,k)] = lu(...
      [M11(:,:,k) M12(:,:,k); M21(:,:,k) zeros(3)]);
end

end % precoStokesMatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = precoStokesBD(o,eta,walls) 
% Apply the preconditioner. This routine is called within GMRES

N = (size(o.precoStokes.U,1) - 3)/2;
nb = size(o.precoYukawa.U,3);
if o.confined
  Nwall = walls.N;
end

istart = 1;
iend = istart + 2*N*nb - 1;
etaDensity = eta(istart:iend);
% extract density function part of eta

istart = 2*N*nb+1;
iend = istart + 2*nb - 1;
etaTranslational = eta(istart:iend);
% extract translational velocity part of eta

istart = 2*N*nb + 2*nb + 1;
iend = istart + 1*nb - 1;
etaRotational = eta(istart:iend);
% extract rotational velocity part of eta

if o.confined
  istart = 2*N*nb + 3*nb + 1;
  iend = istart + 2*Nwall - 1;
  etaWalls = eta(istart:iend);
end

zDensity = zeros(2*N*nb,1);
zTranslational = zeros(2*nb,1);
zRotational = zeros(nb,1);

for k = 1:nb
  istart1 = (k-1)*2*N + 1;
  iend1 = istart1 + 2*N - 1;
  istart2 = (k-1)*2 + 1;
  iend2 = istart2 + 1;
  eta = [etaDensity(istart1:iend1);...
         etaTranslational(istart2:iend2);...
         etaRotational(k)];

  eta = o.precoStokes.U(:,:,k)\...
    (o.precoStokes.L(:,:,k)\eta);

  zDensity(istart1:iend1) = eta(1:2*N);
  zTranslational(istart2:iend2) = eta(2*N+1:2*N+2);
  zRotational(k) = eta(2*N + 3);
end

if o.confined
%  zWalls = etaWalls;
  zWalls = (0.5*eye(2*walls.N) + walls.DLPStokes + ...
        walls.N0Stokes)\etaWalls;
else
  zWalls = [];
end

z = [zDensity;zTranslational;zRotational;zWalls];

end % precoStokesBD


end % methods

end % classdef


