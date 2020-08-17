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
function [eta,eta2,Up,wp,iter,iter2,iflag,res, ...
        iflag2,res2,Unum,Xtest,Ytest,force,torque] = timeStep(o,geom)
% Main time stepping routine
oc = curve;
N = geom.N;
nb = geom.nb;
X = geom.X;
tau = geom.tau;
center = geom.center;
rho = geom.rho;
radii = geom.radii;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVE SCREENED LAPLACE SYSTEM USING GMRES

% Line 65 solves the screened Laplace problem with hydrophobic Dirichlet boundary
% conditions and then evaluates the force and torque [F1 F2 Tq]

addpath("../tests")
yukawa_force

% outputs: 
force  = [F1; F2];
torque = Tq;

eta = 0; eta2 = 0; % variables not used; 
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
NX = 50; NY = 50;
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


% CREATE NEAR SINGULAR INTEGRATION STRUCTURES
if o.inear
  geom.nearStruct = geom.getZone([],1);
end




% append body force and torque to the background flow to complete the
% right hand side fro GMRES
rhs = [-rhs(:); force; torque];

op = poten(N);
% build double-layer potential matrix
geom.DLP = op.stokesDLmatrix(geom);

% max GMRES iterations
maxit = 2*N*nb; 

% SOLVE SYSTEM USING GMRES
[Xn,iflag,res,I] = gmres(@(X) o.timeMatVec(X,geom),rhs,[],...
      o.gmresTol,maxit);
iter = I(2);

% REORGANIZE COLUMN VECTOR INTO MATRIX
% EXTRACT DENSITY FUNCITONS ON FIBRES AND WALLS
% each column of eta corresponds to the density function of a rigid body
eta = zeros(2*N,nb);
for k = 1:nb
  eta(:,k) = Xn((k-1)*2*N+1:k*2*N);
end

% EXTRACT TRANSLATIONAL VELOCITIES
Up = zeros(2,nb);
for k = 1:nb
  Up(:,k) = Xn(2*N*nb+(k-1)*2+1:2*N*nb+2*k);
end

% EXTRACT ROTATIONAL VELOCITIES
wp = zeros(1,nb);
for k = 1:nb
  wp(k) = Xn(2*N*nb+2*nb+k);
end

end % timeStep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tx = timeMatVec(o,Xn,geom)
% Tx = timeMatVec(Xn,geom) does a matvec for GMRES 
% Xn is a state vector that contains the following in order:
% fiber1: eta_x, eta_y fiber2: eta_x, eta_y ... fibernv: eta_x, eta_y
% fiber1: u, v fiber2: u, v ... fibernv: u, v
% fiber1: omega fiber2: omega ... fibernv: omega

N = geom.N;   % points per body
nb = geom.nb; % number of bodies
op = poten(N);

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
valFibers = valFibers - 1/2*eta;

% ADD SELF CONTRIBUTION
valFibers = valFibers + op.exactStokesDLdiag(geom,geom.DLP,eta);

for k = 1:nb
  valFibers(1:end/2,k) = valFibers(1:end/2,k) - Up(1,k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k) - Up(2,k);
end

for k = 1:nb
  valFibers(1:end/2,k) = valFibers(1:end/2,k) ...
                + (geom.X(end/2+1:end,k) - geom.center(2,k))*wp(k);
  valFibers(end/2+1:end,k) = valFibers(end/2+1:end,k)...
                - (geom.X(1:end/2,k) - geom.center(1,k))*wp(k);
end

% EVALUTATE FORCES ON FIBRES
for k = 1:nb
  valForce(1,k) = sum(eta(1:N,k).*geom.sa(:,k))*2*pi/N;
  valForce(2,k) = sum(eta(N+1:2*N,k).*geom.sa(:,k))*2*pi/N;
end

% EVALUATE TORQUES ON FIBRES
for k = 1:nb
  valTorque(k) = sum(((geom.X(N+1:2*N,k)-geom.center(2,k)).*eta(1:N,k) - ...
                     ((geom.X(1:N,k)-geom.center(1,k)).*eta(N+1:2*N,k))).*...
                     geom.sa(:,k))*2*pi/N;
end

% CONSTRUCT OUTPUT VECTOR
Tx = [valFibers(:); -valForce(:);-valTorque(:)];

end % timeMatVec


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
%     figure(3);
%     plot(bc)


end % bcfunc


end % methods

end % classdef

