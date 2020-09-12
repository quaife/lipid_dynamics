classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Routines that we may need to add are different
% integrators such as semi-implicit, SDC, Runge-Kutta, adaptivity

properties

timeOrder       % time stepping order
dt              % time step size
Dp              % Stokes DLP for fiber-fiber interaction
rhs             % Right hand side of mobility problem
rhs2            % Right hand side of screen laplace problem
inear           % flag for using near-singular integration
gmresTol        % GMRES tolerance
plotAxis        % plot axes
farField        % background flow
janusbc         % particle boundary condition
precop          % block-diagonal preconditioner for fibres
precow          % block-diagonal preconditioner for walls
potp            % class for fiber layer potentials
potw            % class for wall layer potentials
om              % monitor class

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options, prams)
% o.tstep(options,prams): constructor.  Initialize class.  Take all
% elements of options and prams needed by the time stepper

o.timeOrder = options.timeOrder;
o.inear = options.inear;
o.dt = prams.T/prams.m;
o.gmresTol = options.gmresTol;
o.plotAxis = options.plotAxis;
o.farField = @(X) o.bgFlow(X,options); 

% for screen laplace BVP
o.janusbc = @(X,tau,center) o.bcfunc(X,tau,center,options);

end % constructor: tstep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Up,wp,iterYukawa,iterStokes] = timeStep(o,geom)
% Main time stepping routine
oc = curve;
N = geom.N;
nb = geom.nb;
X = geom.X;
tau = geom.tau;
center = geom.center;
rho = geom.rho;
radii = geom.radii;

% CREATE NEAR SINGULAR INTEGRATION STRUCTURES
if o.inear
  geom.nearStruct = geom.getZone([],1);
end

% START OF SCREENED LAPLACE SOLVE USING GMRES

% right hand side for the screened Laplace solver
yukawaRHS = geom.yukawaRHS;

op = poten(geom.N,geom.rho);
% build the DLP Yukawa matrix
geom.DLPYukawa = op.yukawaDLmatrix(geom);

% Solve for the density function using GMRES
[sigma,iflagYukawa,resYukawa,iterYukawa] = gmres(...
      @(X) o.timeMatVecYukawa(X,geom),...
      yukawaRHS,[],o.gmresTol,N*nb);
iterYukawa = iterYukawa(2);

% Unstack the density function so that it is arranged as columns for
% each body
etaYukawa = zeros(N,nb);
for k = 1:nb
  etaYukawa(:,k) = sigma((k-1)*N+1:k*N);
end

om = misc;

% dS = velocity*dt = |dx/dt| dt
% sa = |dx/dt|
% dt = 2*pi/N

N    = geom.N;                % number of points per componenet
Nb   = geom.nb;               % number of rigid bodies
x1   = geom.X(1:N,:);         % grid points on curves 
x2   = geom.X(N+1:2*N,:);            
pc   = geom.center;           % center of each rigid body
rho  = geom.rho;              % screen length of particles
xt   = geom.xt;               % tangent unit vector
tau1 = xt(1:N,:);      
tau2 = xt(N+1:2*N,:);
nu1  = +tau2;                 % outward normal : (nu, tau) is right-handed
nu2  = -tau1;

dS   = geom.sa*2*pi/N;        % Jacobian
cur  = geom.cur;              % curvature
h    = etaYukawa;

[F1, F2, Tq] = op.evalForcesQBX(Nb, N, x1, x2, nu1, nu2, dS, rho, etaYukawa);

[F1, F2, Tq; sum(F1) sum(F2) sum(Tq)]


if 0
% plot solution field of screen laplace problem
NX = 151; NY = 151;
% 50 x 50 grid

xmin = o.plotAxis(1);xmax = o.plotAxis(2);
ymin = o.plotAxis(3);ymax = o.plotAxis(4);

xx = linspace(xmin,xmax,NX);
yy = linspace(ymin,ymax,NY);
[xx,yy] = meshgrid(xx,yy);

Xtar = [xx(:);yy(:)];
geomTar.X = Xtar;
geomTar.N = numel(xx);
geomTar.nb = 1;
[~,NearOther] = geom.getZone(geomTar,2);


kernel = @op.exactYukawaDL;
kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,geom.DLPYukawa,z);
% calculate the numerical solution  
yukawaDLPtar = op.nearSingInt(geom,etaYukawa,kernelSelf,...
    NearOther,kernel,kernel,geomTar,false,false);
% Only first half is meaningful since we are solving for a scalar
yukawaDLPtar = yukawaDLPtar(1:end/2);
% reshape so its the same size as the target points
yukawaDLPtar = reshape(yukawaDLPtar,NX,NY);

[yukawaExact, dyuk1, dyuk2] = geom.yukawaExact(xx, yy);

% zero out points that are in the interior of the geometry
% using winding number test that works for arbitrary curves 

CUT_OFF = 0*xx;

for j = 1:geom.nb

  v1 = geom.X(1:end/2,j); 
  v2 = geom.X(end/2+1:end,j);
  Nv = geom.N;        
  v1(N+1) = v1(1);
  v2(N+1) = v2(1);
   
  CUT_OFF = CUT_OFF + oc.wn_PnPoly(xx, yy, v1, v2, Nv);

end

yukawaDLPtar(find(CUT_OFF == 1)) = 0;
yukawaExact(find(CUT_OFF == 1)) = 0;
dyuk1(find(CUT_OFF == 1)) = 0;
dyuk2(find(CUT_OFF == 1)) = 0;

%RJR THERE IS A (-) ERROR IN sigma 

% calculate absolute and relative error of the numerical solution
% in the L1-norm 
abserr = om.trapz2(xx, yy, abs(yukawaExact - yukawaDLPtar))
relerr = om.trapz2(xx, yy, abs(yukawaExact - yukawaDLPtar))/om.trapz2(xx, yy, abs(yukawaExact))
%RJR norm(yukawaExact - yukawaDLPtar,inf)
%RJR relerr = norm(yukawaExact - yukawaDLPtar,inf)/norm(yukawaExact,inf)

% Plot solution
figure(1); clf; hold on;
contour(xx,yy,yukawaExact)
hold on
quiver(xx,yy,dyuk1,dyuk2);
hold off
%surf(xx,yy,yukawaDLPtar)
%shading interp;
%view(2);
%fill3(geom.X(1:end/2,:),geom.X(end/2+1:end,:),10*ones(geom.N,geom.nb),'k')
axis equal
axis([xmin xmax ymin ymax])

% Plot error on a log-base-10 scale
figure(2); clf; hold on;
surf(xx,yy,log10(abs(yukawaDLPtar - yukawaExact)))
shading interp;
view(2);
fill3(geom.X(1:end/2,:),geom.X(end/2+1:end,:),10*ones(geom.N,geom.nb),'k')
axis equal
axis([xmin xmax ymin ymax])
colorbar

end %if

%outputs: 
%force  = [F1; F2];
%torque = Tq;

force = zeros(2*geom.nb,1);
torque = zeros(geom.nb,1);

%eta = 0; eta2 = 0; % variables not used; 
%!!!! BLOCK COMMENTED OUT THROUGH LINE 246 !!!!
%{

rhs2 = o.janusbc(X,tau,center);  
% specify the boundary condition for Janus particles
rhs2 = rhs2(:);

op2 = poten_yukawa(N,rho);
% build double-layer potential matrix
geom.DLP2 = op2.yukawaDLmatrix(geom);

% max GMRES iterations
maxit2 = N*nb; 

% SOLVE SYSTEM USING GMRES
[Xn2,iflag2,res2,I2] = gmres(@(X) o.timeMatVecHalf(X,geom),rhs2,[],...
      o.gmresTol,maxit2);
  
iter2 = I2(2);

% REORGANIZE COLUMN VECTOR INTO MATRIX
% EXTRACT DENSITY FUNCITONS ON FIBRES AND WALLS
% each column of eta corresponds to the density function of a rigid body
eta2 = zeros(N,nb);
for k = 1:nb
  eta2(:,k) = Xn2((k-1)*N+1:k*N);
end


%%%%%%%%%%%%%%%%%
% plot solution field of screen laplace problem
NX = 500; NY = 500;
% 50 x 50 grid

xmin = o.plotAxis(1);xmax = o.plotAxis(2);
ymin = o.plotAxis(3);ymax = o.plotAxis(4);

xx = linspace(xmin,xmax,NX);
yy = linspace(ymin,ymax,NY);


% test 07/31/2020 
% rr = linspace(1.1, 2, NX);
% tt = linspace(0,2*pi,NY);
% 
% [RR, TT] = meshgrid(rr,tt);
% Uxact = 1/2*(besselk(0,RR/rho)/besselk(0,1/rho)+...
%             besselk(1,RR/rho)/besselk(1,1/rho).*cos(TT));
% Xtest = RR.*cos(TT);
% Ytest = RR.*sin(TT);

[Xtest, Ytest] = meshgrid(xx,yy);
Ztest = Xtest+1i*Ytest;
ind_int = []; %indices of grid points inside particles
% we only want data outside particles
zc = center(1,:) + 1i*center(2,:);
for j = 1:nb
ind_int = [ind_int; find(abs(Ztest-zc(j))<radii(j))];
end

Unum = Xtest; 

  [xsou,ysou] = oc.getXY(geom.X);
  sa = geom.sa(:,:);
  sa = sa(:);
  % Jacobian
  zt = xsou(:)+1i*ysou(:);
  normalx = geom.xt(N+1:2*N,:);
  normaly = -geom.xt(1:N,:);


for k = 1:NX
    for j = 1:NY       
        indtmp = (k-1)*NY+j;
        xtest = Xtest(indtmp);
        ytest = Ytest(indtmp);
        ztest = xtest + 1i*ytest;
% calculate the term  (xsou-xtar)\cdot normal(y)
        const = -((xsou(:)-xtest).*normalx(:) + (ysou(:)-ytest).*normaly(:));
        rr = abs(zt-ztest);
        tmp = 1i*rho*rr;
        kernel = besselh(1,tmp);       
% boundary integration
        Unum(indtmp)= Xn2'*(-1i/4*tmp.*kernel(:).*const./rr.^2.*sa)*2*pi/N;
    end
end
Unum(ind_int)=0;

% test 07/31/2020 
% figure(4);
% surf(Xtest,Ytest,Uxact-Unum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bug is still there
% compute the body force and torque

force = zeros(2*nb,1);  
torque = zeros(nb,1);   

GradUmat = op2.yukawaGradDLmatrix(geom);
% spy(GradUmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           |           |           | %
%           |x component|     O     | %
%           |           |           | %
%  GradUmat=|-----------|-----------| %
%           |           |           | %
%           |     O     |y component| %
%           |           |           | %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

sa = geom.sa(:,:);
% sa2 = [sa; sa];
% jacobian

% density = [eta2; eta2];

grad = zeros(2*N,nb);
sol = zeros(N,nb);
stress = zeros(2,2,N,nb);
TdotNorm = zeros(2,N,nb);
r0 = zeros(2,N,nb);

% everything is on the boundary
% First, obtain the solution and gradient on the boundary
% by using double layer potential
for k = 1:nb
  sol(:,k) = eta2(:,k)'*(geom.DLP2(:,:,k).*sa(:,k))*2*pi/N;
  grad(1:N,k) = eta2(:,k)'*(GradUmat(1:N,1:N,k).*sa(:,k))*2*pi/N;
  grad(N+1:2*N,k) = eta2(:,k)'*(GradUmat(N+1:2*N,N+1:2*N,k).*sa(:,k))*2*pi/N;
end

figure(3)   % test numerical sol on boundary
% truesol = 0.5*(1+cos(0:2*pi/512:2*pi-2*pi/512));
% plot(rhs2-truesol')
% truegradx = (0.5*sin(1+cos(0:2*pi/512:2*pi-2*pi/512)))'...
%      ./(1+X(N+1:2*N).^2./X(1:N).^2).*(X(N+1:2*N)./X(1:N).^2);
%  size(truegradx)
% plot(truegradx)
plot(grad(1:N))

normal = [geom.xt(N+1:2*N,:); -geom.xt(1:N,:)]; 
  
% Second, calculate the hydrophobic stress, force and, torque.  
for k=1:nb   
    for j=1:N
        gradvec = [grad(j,k); grad(j+N,k)];

% calculuate point stress
        stress(:,:,j,k) = 1/rho*sol(j,k)^2*eye(2) + ...
            2*rho*(0.5*norm(gradvec)^2*eye(2)-gradvec*gradvec');
       
        normali = [normal(j,k); normal(j+N,k)];
        
%         TdotNorm(:,j,k) = (rho*norm(gradvec)^2 + 1/rho*sol(j,k)^2)*normali ...
%              - 2*rho*gradvec*(gradvec'*normali);

        TdotNorm(:,j,k) = stress(:,:,j,k)'*normali;           
        TnuVec = TdotNorm(:,j,k);

% we need the vector from center to the correpsonding boundary
        r0(:,j,k) = [xsou(j,k)-center(1,k); ysou(j,k)-center(2,k)];
        r0vec = r0(:,j,k);
        
        % cross product r0 x TdotN
        cpval = cross ([r0vec;0], [TnuVec;0]);

% force and torque calculated by using integrations
        force(k) = force(k) - TnuVec(1)*sa(j,k)*2*pi/N;
        force(k+nb) = force(k+nb) - TnuVec(2)*sa(j,k)*2*pi/N;
        torque(k) = torque(k) - cpval(3)*sa(j,k)*2*pi/N;
    end
    
end
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve mobility problem here

% far field condition
rhs = o.farField(X);

% append body force and torque to the background flow to complete the
% right hand side fro GMRES
rhs = [-rhs(:); force; torque];

op = poten(N,rho);
% build double-layer potential matrix
geom.DLPStokes = op.stokesDLmatrix(geom);

% max GMRES iterations
maxit = 2*N*nb; 

% SOLVE SYSTEM USING GMRES
[sigma,iflagStokes,resStokes,iterStokes] = gmres(@(X) o.timeMatVecStokes(X,geom),...
      rhs,[],o.gmresTol,maxit);
iterStokes = iterStokes(2);

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
valFibers = valFibers - 0.5*eta;

% ADD SELF CONTRIBUTION
valFibers = valFibers + op.exactStokesDLdiag(geom,geom.DLPStokes,eta);

% DEFINE STOKES DLP KERNELS
kernel = @op.exactStokesDL;
kernelSelf = @(z) +0.5*z + op.exactStokesDLdiag(geom,geom.DLPStokes,z);
% BQ: NOT SURE WHY THIS IS +0.5 and not -0.5 like it is two lines above

% ADD CONTRIBUTION FROM OTHER FIBERS
stokesDLP = op.nearSingInt(geom,eta,kernelSelf,geom.nearStruct,...
    kernel,kernel,geom,true,false);
valFibers = valFibers + stokesDLP;

% ADD TRANSLATIONAL VELOCITY CONTRIBUTION
for k = 1:nb
  valFibers(1:end/2,k) = valFibers(1:end/2,k) - Up(1,k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k) - Up(2,k);
end

% ADD ROTATIONAL VELOCITY CONTRIBUTION
for k = 1:nb
  valFibers(1:end/2,k) = valFibers(1:end/2,k) ...
                + (geom.X(end/2+1:end,k) - geom.center(2,k))*wp(k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k)...
                - (geom.X(1:end/2,k) - geom.center(1,k))*wp(k);
end

% EVALUTATE FORCES ON FIBERS
for k = 1:nb
  valForce(1,k) = sum(eta(1:N,k).*geom.sa(:,k))*2*pi/N;
  valForce(2,k) = sum(eta(N+1:2*N,k).*geom.sa(:,k))*2*pi/N;
end

% EVALUATE TORQUES ON FIBERS
for k = 1:nb
  valTorque(k) = sum(((geom.X(N+1:2*N,k)-geom.center(2,k)).*eta(1:N,k) - ...
                     ((geom.X(1:N,k)-geom.center(1,k)).*eta(N+1:2*N,k))).*...
                     geom.sa(:,k))*2*pi/N;
end

% CONSTRUCT OUTPUT VECTOR
Tx = [valFibers(:); -valForce(:);-valTorque(:)];

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

% DEFINE Yukawa DLP KERNELS
kernel = @op.exactYukawaDL;
kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,geom.DLPYukawa,z);

yukawaDLP = op.nearSingInt(geom,eta,kernelSelf,geom.nearStruct,...
    kernel,kernel,geom,true,false);

Tx = Tx + yukawaDLP(1:geom.N,:);

Tx = Tx(:);

end % timeMatVecYukawa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVecHalf(o,Xn,geom)
% Tx = timeMatVecHalf(Xn,geom) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% fiber1: eta1 fiber2: eta2 ... fibernv: etan

N = geom.N;   % points per body
nb = geom.nb; % number of bodies
rho = geom.rho; 
op = poten_yukawa(N,rho);

% Output of Tx that corresponds to the shape of eta
Tx = zeros(N,nb);

% BEGIN FORMATTING UNKNOWN VECTOR
% EXTRACT DENSITY FUNCTION FOR FIBRES
eta2 = zeros(N,nb);
for k = 1:nb
  eta2(:,k) = Xn((k-1)*N+1:k*N);
end

% ADD JUMP IN DLP
Tx = Tx - 1/2*eta2;

% ADD SELF CONTRIBUTION
Tx = Tx + op.exactYukawaDLdiag(geom,geom.DLP2,eta2);
Tx = Tx(:);
end % timeMatVecHalf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vInf = bgFlow(~,X,options)
    
N = size(X,1)/2;
nb  = size(X,2);
oc = curve;

[x,y] = oc.getXY(X);
shearRate = 0.1; % manually set the shear rate here for now

switch options.farField
  case 'shear'
    vInf = [shearRate*y;zeros(N,nb)];
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


end % methods

end % classdef


