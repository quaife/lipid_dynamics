classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Routines that we may need to add are different
% integrators such as semi-implicit, SDC, Runge-Kutta, adaptivity

properties

dt              % time step size
sstep           % starting step
gmresTol        % GMRES tolerance
farField        % background flow
janusbc         % particle boundary condition
tracer          % flag for tracer
precoYukawa     % block-diagonal preconditioner for Yukawa
precoStokes     % block-diagonal preconditioner for Stokes on bodies
precoWalls      % block-diagonal preconditioner for walls

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

% for screen laplace BVP
o.janusbc  = @(X,tau,center) o.bcfunc(X,tau,center,options);

end % constructor: tstep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Up,wp,iterYukawa,iterStokes,etaYukawa,etaStokes,...
      fforce,torque] = timeStep(o,geom,etaY0,etaS0)
% Main time stepping routine
oc     = curve;
N      = geom.N;
nb     = geom.nb;
X      = geom.X;
tau    = geom.tau;
center = geom.center;
rho    = geom.rho;
gam    = geom.gam;
radii  = geom.radii;

% CREATE NEAR SINGULAR INTEGRATION STRUCTURES
geom.nearStruct = geom.getZone([],1);

op = poten(geom.N,geom.rho);
% START OF REPULSION CALCULATION
[R1, R2, RTq, pp] = op.Repul(geom);
% END OF REPULSION CALCULATION


% START OF SCREENED LAPLACE SOLVE USING GMRES
% right hand side for the screened Laplace solver
yukawaRHS = geom.yukawaRHS;

% build the DLP Yukawa matrix
geom.DLPYukawa = op.yukawaDLmatrix(geom); 
o.precoYukawa.L = zeros(N,N,nb);
o.precoYukawa.U = zeros(N,N,nb);
for k = 1:nb
  [o.precoYukawa.L(:,:,k),o.precoYukawa.U(:,:,k)] = ...
    lu(0.5*eye(N) + geom.DLPYukawa(:,:,k));
end

% build precomputed matrix with all necessary bessel functions that are
% needed at every GMRES iteration
geom.YukawaKernelMatrix;

% Solve for the density function using GMRES
[sigma,iflagYukawa,resYukawa,iterYukawa] = gmres(...
      @(etaY0) o.timeMatVecYukawa(etaY0,geom) ,...
      yukawaRHS, [], o.gmresTol, N*nb, ...
      @(etaY0) o.precoYukawaBD(etaY0)); 
%[sigma,iflagYukawa,resYukawa,iterYukawa] = gmres(...
%      @(etaY0) o.timeMatVecYukawa(etaY0,geom) ,...
%      yukawaRHS, [], o.gmresTol, N*nb); 
iterYukawa = iterYukawa(2);
fprintf('Yukawa required %i iterations\n',iterYukawa);

% Unstack the density function so that it is arranged as columns for
% each body
etaYukawa = zeros(N,nb);
for k = 1:nb
  etaYukawa(:,k) = sigma((k-1)*N+1:k*N);
end

% Bryan's rewrite of Rolf's code to evaluate forces using QBX and Tpq +
% Tqp identity
[F1,F2,Tq] = op.evalForcesQBX(geom,etaYukawa);
% lift Tq so that it is for rotation about center (not origin)

% tst = fopen("HAP_force_distance.dat", "r");
% if tst == -1
%     fid = fopen("HAP_force_distance.dat", "w");
% else 
%     fclose(tst);
%     fid = fopen("HAP_force_distance.dat", "a");
% end
% fprintf(fid,"%d %d\n", [geom.center(1,2) - geom.center(1,1) - 2,  F1(1)]');
% fclose(fid);
%fprintf("%d %d\n", [geom.center(1,2) - geom.center(1,1) - 2,  F1(1)]');

%{
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
pause(0.01)
%}

%outputs: 
fforce  = [F1 + R1, F2 + R2].';
force  = fforce(:);
torque = Tq + RTq;

%force = 0*[-2;1;0;0;2;-1];
%torque = 1*[-10;0;+10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve mobility problem here

% far field condition plus the stokeslet and rotlet terms
rhs = o.farField(X) + op.StokesletRotlet(geom,force,torque);

% append body force and torque to the background flow to complete the
% right hand side for GMRES
rhs = [-rhs(:); force; torque];

op = poten(N,rho);
% build double-layer potential matrix
geom.DLPStokes = op.stokesDLmatrix(geom);

% Build LU factorization of the completed Stokes DLP system of equations
o.precoStokesMatrix(geom);

% DEBUGGING CODE TO BUILD LINEAR SYSTEM AND PRECONDITIONER AS MATRICIES
%A = zeros(numel(rhs));
%P = zeros(numel(rhs));
%for k = 1:numel(rhs)
%  e = zeros(numel(rhs),1);
%  e(k) = 1;
%  A(:,k) = o.timeMatVecStokes(e,geom);
%  P(:,k) = o.precoStokesBD(e);
%end
%clf
%surf(A*P - eye(size(P)))
%invP = inv(P);
%[invP(end-1:end,end-10:end)' ...
%A(end-1:end,end-10:end)']
%pause

% max GMRES iterations
maxit = 2*N*nb; 


% Build the upsampled matrix for doing stokes Layer potentials
geom.StokesDLMatrix;

% SOLVE SYSTEM USING GMRES
%tic

[sigma, iflagStokes, resStokes, iterStokes] = gmres(...
      @(etaS0) o.timeMatVecStokes(etaS0,geom),...
      rhs,[],o.gmresTol,maxit,...
      @(etaY0) o.precoStokesBD(etaY0));
      
%[sigma, iflagStokes, resStokes, iterStokes] = gmres(...
%      @(etaS0) o.timeMatVecStokes(etaS0,geom),...
%      rhs,[],o.gmresTol,maxit);
%toc

iterStokes = iterStokes(2);
fprintf('Stokes required %i iterations\n',iterStokes);

% REORGANIZE COLUMN VECTOR INTO MATRIX
% EXTRACT DENSITY FUNCITONS ON FIBRES AND WALLS
% each column of etaStokes corresponds to the density function of a rigid body
etaStokes = zeros(2*N,nb);
for k = 1:nb
  etaStokes(:,k) = sigma((k-1)*2*N+1:k*2*N);
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

%[mean(Up(1,:)); mean(Up(2,:))]
%Up = Up - [mean(Up(1,:)); mean(Up(2,:))];

% Routine to compute the velocity in the bulk and use the rigid body
% motion to define a velocity in rigid bodies. The goal is to see
% something cotinuous
if o.tracer
  op.bulkVelocity(geom,etaStokes,Up,wp,force,torque,...
      @(X) o.farField(X));
end

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVecStokes(o,Xn,geom)
% Tx = timeMatVecStokes(Xn,geom) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% fiber1: eta_x, eta_y fiber2: eta_x, eta_y ... fibernv: eta_x, eta_y
% fiber1: u, v fiber2: u, v ... fibernv: u, v
% fiber1: omega fiber2: omega ... fibernv: omega

N = geom.N;   % points per body
nb = geom.nb; % number of bodies
rho = geom.rho; % screen length

op = poten(N,rho);

% Output of Tx that corresponds to the velocity of the fibers
valFibers = zeros(2*N,nb);
% output of Tx that corresponds to force on fibres
valForce = zeros(2,nb);
% output of Tx that corresponds to torque on fibres
valTorque = zeros(1,nb);

% BEGIN FORMATTING UNKNOWN VECTOR
% EXTRACT DENSITY FUNCTION FOR FIBRES
eta = zeros(2*N,nb);
for k = 1:nb
  eta(:,k) = Xn((k-1)*2*N+1:k*2*N);
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
% END FORMATTING UNKNOWN VECTOR

% ADD JUMP IN DLP
jump = +0.5;
valFibers = valFibers + jump*eta;

% ADD SELF CONTRIBUTION
valFibers = valFibers + op.exactStokesDLdiag(geom,geom.DLPStokes,eta);

% Matrix-free version
if 0
% DEFINE STOKES DLP KERNELS
kernel = @op.exactStokesDL;
kernelSelf = @(z) +jump*z + op.exactStokesDLdiag(geom,geom.DLPStokes,z);

% ADD CONTRIBUTION FROM OTHER FIBERS
stokesDLP = op.nearSingInt(geom,eta,kernelSelf,geom.nearStruct,...
    kernel,kernel,geom,true,false);
end

% Precomputing and using matrix (memory-intensive)
if 1
kernel = @op.exactStokesDLMatFree;
kernelDirect = @op.exactStokesDL;
kernelSelf = @(z) jump*z + op.exactStokesDLdiag(geom,geom.DLPStokes,z);

stokesDLP = op.nearSingInt(geom,eta,kernelSelf,geom.nearStruct,...
    kernel,kernelDirect,geom,true,false);
end

valFibers = valFibers + stokesDLP;

% ADD TRANSLATIONAL VELOCITY CONTRIBUTION
for k = 1:nb
  valFibers(1:end/2,k) = valFibers(1:end/2,k) - Up(1,k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k) - Up(2,k);
end

% ADD ROTATIONAL VELOCITY CONTRIBUTION
for k = 1:nb
  valFibers(1:end/2,k) = valFibers(1:end/2,k) - ...
                (-(geom.X(end/2+1:end,k) - geom.center(2,k)))*wp(k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k) - ...
                (+(geom.X(1:end/2,k) - geom.center(1,k)))*wp(k);
end

% EVALUTATE FORCES ON FIBERS
for k = 1:nb
  valForce(1,k) = 1/(2*pi)*sum(eta(1:N,k).*geom.sa(:,k))*2*pi/N;
  valForce(2,k) = 1/(2*pi)*sum(eta(N+1:2*N,k).*geom.sa(:,k))*2*pi/N;
end

% EVALUATE TORQUES ON FIBERS
for k = 1:nb
  valTorque(k) = 1/(2*pi)*sum((-(geom.X(N+1:2*N,k)-geom.center(2,k)).*eta(1:N,k) + ...
                     (geom.X(1:N,k)-geom.center(1,k)).*eta(N+1:2*N,k)).*...
                     geom.sa(:,k))*2*pi/N;
end

% CONSTRUCT OUTPUT VECTOR
Tx = [valFibers(:);+valForce(:);+valTorque(:)];

end % timeMatVecStokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVecYukawa(o,Xn,geom)
% Tx = timeMatVecYukawa(Xn,geom) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% fiber1: eta1 fiber2: eta2 ... fibernv: etan

N = geom.N;
nb = geom.nb;
rho = geom.rho;
op = poten(N,rho);

% Store with each column being one part of the matvec. Will columnize
% the vector at the end of this routine
Tx = zeros(N,nb);

eta = zeros(N,nb);
% unstack the density function
for k = 1:nb
  eta(:,k) = Xn((k-1)*N+1:k*N);
end

% ADD JUMP IN DLP
Tx = Tx + 0.5*eta;

% ADD SELF CONTRIBUTION
Tx = Tx + op.exactYukawaDLdiag(geom,geom.DLPYukawa,eta);

if 0
% DEFINE Yukawa DLP KERNELS
kernel = @op.exactYukawaDL;
kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,geom.DLPYukawa,z);

yukawaDLP = op.nearSingInt(geom,eta,kernelSelf,geom.nearStruct,...
    kernel,kernel,geom,true,false);
end

if 1
% START OF AN EXPERIMENT TO SPEED UP THE CODE BY PRECOMPUTING BESSELK OF
% ALL DISTANCES
kernel = @op.exactYukawaDLMatFree;
kernelDirect = @op.exactYukawaDL;
kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,geom.DLPYukawa,z);

yukawaDLP = op.nearSingInt(geom,eta,kernelSelf,geom.nearStruct,...
    kernel,kernelDirect,geom,true,false);
end
% END OF AN EXPERIMENT TO SPEED UP THE CODE BY PRECOMPUTING BESSELK OF
% ALL DISTANCES

Tx = Tx + yukawaDLP(1:geom.N,:);
Tx = Tx(:);

end % timeMatVecYukawa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    vInf = shearRate*[-cos(x).*sin(y);sin(x).*cos(y)];
  case 'parabolic'
    vInf = [shearRate*(1-(y/80).^2);zeros(N,nb)];
  otherwise
    vInf = zeros(2*N,nb);
end
    
end % bgFlow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = precoYukawaBD(o,eta); 

N = size(o.precoYukawa.U,1);
nb = size(o.precoYukawa.U,3);
z = zeros(N*nb,1);
for k = 1:nb
  istart = (k-1)*N + 1;
  iend = istart + N - 1;
  z(istart:iend) = o.precoYukawa.U(:,:,k)\...
    (o.precoYukawa.L(:,:,k)\eta(istart:iend));
end

end % precoYukawaBD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function precoStokesMatrix(o,geom)
% Build the LU decomposition of the Stokes linear system.

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = precoStokesBD(o,eta); 

N = (size(o.precoStokes.U,1) - 3)/2;
nb = size(o.precoYukawa.U,3);

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

z = [zDensity;zTranslational;zRotational];

end % precoStokesBD


end % methods

end % classdef


