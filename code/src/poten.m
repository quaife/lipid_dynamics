classdef poten 
% this class defines single and double layers for various kernels
% (stokes, laplace) on 2D periodic curves.
% Defines the matricies that map a density function defined on the
% boundary of a curve to the layer potential evaluated on the curve,
% and defines the operators that takes a density function defined on
% the boundary of a curve and returns the layer 
% potential at arbitrary target points.
% This class also has the main routine that evaluates layer
% potentials using near-singular integration.
    
properties
    
  N; % points per curve
  interpMat;  
  profile;
  rho % screen length
  om;
  
end % properties

methods 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = poten(N,rho)
% o = poten(N): constructor; N is the number of points per curve

o.interpMat = o.lagrangeInterp;
% load in the interpolation matrix which is precomputed with 7
% interpolation points

o.N = N;
o.rho = rho;

end % poten: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingInt(o,geomSou,f,selfMat,...
    NearStruct,kernel,kernelDirect,geomTar,tEqualS,idebug)
% LP = nearSingInt(geom,f,selfMat,NearStruct,kernel,kernelDirect,
% geomTar,tEqualS,idebug) computes a layer potential due to f at all
% points in geomTar.X. If tEqualS==true, then the geomTar == geomSou and
% the self-geom interaction is skipped. selfMat is the diagonal of the
% potential needed to compute the layer potential of each geom
% indepenedent of all others. kernel and kernelDirect are two (possibly
% the same) routines that compute the layer potential.  kernelDirect
% always uses the direct method whereas kernel may use an
% FMM-accelerated method. NearStruct is a structure containing the
% variables zone,dist,nearest,icp,argnear which are required by
% near-singular integration (they keep everything sorted and
% precomputed). Everything is in the 2*N x n format. Can pass a final
% argument if desired so that plots of the near-singular integration
% algorithm are displayed

if (tEqualS && size(geomSou.X,2) == 1)
  LP = zeros(size(geomSou.X));
  return
end
% only a single geom, so layer potential on all other geoms will always
% be zero

if (nargin == 9)
  idebug = false;
end

dist = NearStruct.dist;
zone = NearStruct.zone;
nearest = NearStruct.nearest;
icp = NearStruct.icp;
argnear = NearStruct.argnear;
nearFibers = NearStruct.nearFibers;

Xsou = geomSou.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
nsou = size(Xsou,2); % number of source 'geoms'
Xtar = geomTar.X; % target positions
Ntar = size(Xtar,1)/2; % number of target points
ntar = size(Xtar,2); % number of target 'geoms'
h = geomSou.length/Nsou; % arclength term

Nup = Nsou*ceil(sqrt(Nsou));

vself = selfMat(f);

% pad vself with zeros if dealing with a scalar layer potential
% and upsample to N^(3/2).  
if size(vself,1) == geomSou.N
  vself = [vself;zeros(geomSou.N,geomSou.nb)];
  fup = interpft(f,Nup);
else
  fup = [interpft(f(1:Nsou,:),Nup); interpft(f(Nsou+1:2*Nsou,:),Nup)];
end

prams.N = Nup;
prams.nb = geomSou.nb;
prams.tau = geomSou.tau;
prams.radii = geomSou.radii;
prams.rho = geomSou.rho;
prams.gam = geomSou.gam;
prams.ar = geomSou.ar;
prams.RepulLength = geomSou.RepulLength;
prams.RepulStrength = geomSou.RepulStrength;
prams.bcShift = geomSou.bcShift;
prams.bcType = geomSou.bcType;

xc = geomSou.center;
tau = geomSou.tau;
% BQ: tau and xc are implemented funny

geomUp = capsules(prams,xc,tau);
% Build an object with the upsampled geom

interpOrder = size(o.interpMat,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right

if tEqualS % sources == targets
  if nsou > 1
    for k = 1:nsou
      K = [(1:k-1) (k+1:nsou)];
      [~,farField(:,k)] = kernelDirect(geomUp,fup,Xtar(:,k),K);
    end
    % This is a huge savings if we are using a direct method rather
    % than the fmm to evaluate the layer potential.  The speedup is
    % more than N^{1/2}, where N is the resolution of the geoms
    % that we are computing with
  else
    farField = zeros(2*Ntar,ntar);
  end

else % sources ~= targets
  [~,farField] = kernel(geomUp,fup,Xtar,1:nsou);
  % evaluate layer potential due to all 'geoms' at all points in Xtar
end
% Use upsampled trapezoid rule to compute layer potential

if size(farField,1) == geomTar.N
  farField = [farField;zeros(geomTar.N,geomTar.nb)];
end

nearField = zeros(2*Ntar,ntar);

beta = 1.1;
% small buffer to make sure Lagrange interpolation points are not in the
% near zone
vel = zeros(2*Ntar, nsou, nsou); %allocate array for vel

for k1 = 1:nsou
  if tEqualS % sources == targets
    % skip diagonal geom
    K = nearFibers{k1};
  else % sources ~= targets
    K = (1:ntar);
    % consider all geoms
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    % set of points on geom k2 close to geom k1
    if (numel(J) ~= 0)
      indcp = icp{k1}(J,k2);
      % closest point on geom k1 to each point on geom k2 that is close
      % to geom k1
      for j = 1:numel(J)
        pn = mod((indcp(j)-p+1:indcp(j)-p+interpOrder)' - 1,Nsou) + 1;
        % index of points to the left and right of the closest point
        v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
              o.interpMat*vself(pn,k1));
        vel(J(j),k2,k1) = v(end);
        % x-component of the layer potential at the closest point
        v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
              o.interpMat*vself(pn+Nsou,k1));
        vel(J(j)+Ntar,k2,k1) = v(end);
        % y-component of the layer potential at the closest point
      end
      % compute values of layer  potential at required intermediate
      % points using local interpolant
            
      if ((numel(J) + numel(fup)) >= 512 && numel(J) > 32)
        [~,potTar] = kernel(geomUp,fup,...
              [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
      else
        [~,potTar] = kernelDirect(geomUp,fup,...
              [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
      end
      % Need to subtract off contribution due to geom k1 since its layer
      % potential will be evaulted using Lagrange interpolant of nearby
      % points
      nearField(J,k2) =  nearField(J,k2) - potTar(1:numel(J));
      if numel(f) == geomSou.N*geomSou.nb
        nearField(J+Ntar,k2) = 0;
      else
        nearField(J+Ntar,k2) =  nearField(J+Ntar,k2) - ...
              potTar(numel(J)+1:end);
      end
            
      XLag = zeros(2*numel(J),interpOrder - 1);
      % initialize space for initial tracer locations
      for i = 1:numel(J)
        nx = (Xtar(J(i),k2) - nearest{k1}(J(i),k2))/...
              dist{k1}(J(i),k2);
        ny = (Xtar(J(i)+Ntar,k2) - nearest{k1}(J(i)+Ntar,k2))/...
              dist{k1}(J(i),k2);
        XLag(i,:) = nearest{k1}(J(i),k2) + ...
              beta*h(k1)*nx*(1:interpOrder-1);
        XLag(i+numel(J),:) = nearest{k1}(J(i)+Ntar,k2) + ...
              beta*h(k1)*ny*(1:interpOrder-1);
        % Lagrange interpolation points coming off of geom k1 All points
        % are behind Xtar(J(i),k2) and are sufficiently far from geom k1
        % so that the Nup-trapezoid rule gives sufficient accuracy
      end
            
      if (numel(XLag)/2 > 100)
        [~,lagrangePts] = kernel(geomUp,fup,XLag,k1);
      else
        [~,lagrangePts] = kernelDirect(geomUp,fup,XLag,k1);
      end
      % evaluate layer potential at the lagrange interpolation points

      if numel(f) == geomSou.N*geomSou.nb
        lagrangePts = [lagrangePts;lagrangePts];
      end
            
      for i = 1:numel(J)
        Px = o.interpMat*[vel(J(i),k2,k1) ...
              lagrangePts(i,:)]';
        Py = o.interpMat*[vel(J(i)+Ntar,k2,k1) ...
              lagrangePts(i+numel(J),:)]';
        % Build polynomial interpolant along the one-dimensional points
        % coming out of the geom
        dscaled = full(dist{k1}(J(i),k2)/(beta*h(k1)*(interpOrder-1)));
        % Point where interpolant needs to be evaluated
                
        v = filter(1,[1 -dscaled],Px);
        nearField(J(i),k2) = nearField(J(i),k2) + v(end);
                
        v = filter(1,[1 -dscaled],Py);
        nearField(J(i)+Ntar,k2) = nearField(J(i)+Ntar,k2) + v(end);
                
        if idebug
%          if numel(f) == geomSou.N*geomSou.nb
%            f = [f;zeros(geomSou.N,geomSou.nb)];
%          end

          figure(2); clf; hold on;
          plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.','markersize',10)
          plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.','markersize',10)
          plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.','markersize',10)
          plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx',...    
              'markersize',10)
          plot(XLag(i,:),XLag(numel(J)+i,:),'gx','markersize',10)
          axis equal
                    
          figure(1); clf; hold on
          plot((0:interpOrder-1)*beta*h(k1),...
                  real([vel(J(i),k2,k1) lagrangePts(i,:)]),'g-o')
          if size(f,1) == 2*geomSou.N
            plot((0:interpOrder-1)*beta*h(k1),...
                  real([vel(J(i)+Ntar,k2,k1) lagrangePts(...
                  i+numel(J),:)]),'r--o')
          end
                    
          figure(3)
          clf
          hold on
          plot(f(1:Nsou,k1));
          if size(f,1) == 2*geomSou.N
            plot(f(Nsou+1:2*Nsou,k1));
          end
                    
          drawnow;
          pause
                    
        end
        % DEBUG: PASS IN idebug=true INTO THIS ROUTINE AND THEN YOU CAN
        % SEE THE INTERPOLATION POINTS AND CHECK THE SMOOTHNESS OF THE
        % INTERPOLANT
      end % i
    end % numel(J) ~= 0
    % Evaluate layer potential at Lagrange interpolation points if there
    % are any
  end % k2
end % k1
% farField

LP = farField + nearField;
% Add kernel due to far points and near points.  Far points were
% upsampled if source==geom so need to truncate here.  We are 
% only using Ntar target points.  Note that it is only the sources 
% that were upsampled

end % nearSingInt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT BUILD LAYER-POTENTIAL MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLmatrix(~,geom)
% D = stokesDLmatrix(geom), generate double-layer potential for 
% Stokes. geom is a data structure defined as in the capsules class
% D is (2N,2N,n) array where Np is the number of points per curve and 
% n is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(geom.X);
% Vesicle positions
Ng = geom.N;
% number of points per geoms
D = zeros(2*Ng,2*Ng,geom.nb);
% initialize space for double-layer potential matrix

for k=1:geom.nb  % Loop over curves
  xx = x(:,k);
  yy = y(:,k);
  % locations
  [tx,ty] = oc.getXY(geom.xt);
  tx = tx(:,k); ty = ty(:,k);
  % Vesicle tangent
  sa = geom.sa(:,k)';
  % Jacobian
  cur = geom.cur(:,k)';
  % curvature

  xtar = xx(:,ones(Ng,1));
  ytar = yy(:,ones(Ng,1));
  % target points

  xsou = xx(:,ones(Ng,1))';
  ysou = yy(:,ones(Ng,1))';
  % source points

  txsou = tx';
  tysou = ty';
  % tangent at sources
  sa = sa(ones(Ng,1),:);
  % Jacobian

  diffx = xtar - xsou;
  diffy = ytar - ysou;
  rho4 = (diffx.^2 + diffy.^2).^(-2);
  rho4(1:Ng+1:Ng^2) = 0;
  % set diagonal terms to 0

  kernel = diffx.*(tysou(ones(Ng,1),:)) - ...
            diffy.*(txsou(ones(Ng,1),:));
  kernel = kernel.*rho4.*sa;

  D11 = kernel.*diffx.^2;
  % (1,1) component
  D11(1:Ng+1:Ng^2) = -0.5*cur.*sa(1,:).*txsou.^2;
  % diagonal limiting term

  D12 = kernel.*diffx.*diffy;
  % (1,2) component
  D12(1:Ng+1:Ng^2) = -0.5*cur.*sa(1,:).*txsou.*tysou;
  % diagonal limiting term

  D22 = kernel.*diffy.^2;
  % (2,2) component
  D22(1:Ng+1:Ng^2) = -0.5*cur.*sa(1,:).*tysou.^2;
  % diagonal limiting term

  D(:,:,k) = [D11 D12; D12 D22];
  % build matrix with four blocks
  D(:,:,k) = 1/pi*D(:,:,k)*2*pi/Ng;
  % scale with the arclength spacing and divide by pi
end % k

end % stokesDLmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N0 = stokesN0matrix(~,wall)
% N0 = stokesN0matrix(vesicle) generates the the integral operator with 
% kernel normal(x) \otimes normal(y) which removes the rank one defficiency 
% of the double-layer potential.  Need this operator for solid walls

% Normal vector
normal = [wall.xt(wall.N+1:2*wall.N,:);-wall.xt(1:wall.N,:)]; 
normal = normal(:,ones(2*wall.N,1));

sa = [wall.sa(:,1);wall.sa(:,1)];
sa = sa(:,ones(2*wall.N,1));
N0 = zeros(2*wall.N,2*wall.N,wall.n);
N0(:,:,1) = normal.*normal'.*sa'*2*pi/wall.N;
% Use N0 if solving (-1/2 + DLP)\eta = f where f has no flux through
% the boundary.  By solving (-1/2 + DLP + N0)\eta = f, we guarantee
% that \eta also has no flux through the boundary.  This is not
% required, but it means we're enforcing one addition condition on eta
% which removes the rank one kernel.  DLP is the double-layer potential
% for stokes equation

end % stokesN0matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = yukawaDLmatrix(~,geom)
% D = yukawaDLmatrix(geom) generates double-layer potential for yukawa.
% geom is a data structure defined as in the capsules class D is (N,N,n)
% array where Np is the number of points per curve and n is the number
% of curves in X 

oc = curve;
[x,y] = oc.getXY(geom.X);
% Vesicle positions
Ng = geom.N;
% number of points per geoms
D = zeros(Ng,Ng,geom.nb);
% initialize space for double-layer potential matrix

for k=1:geom.nb  % Loop over curves
  % x and y coordinates of the sources/targets
  xx = x(:,k);
  yy = y(:,k);

  % Normal vector(y)  
  nx = -geom.xt(Ng+1:2*Ng,k);
  ny = geom.xt(1:Ng,k);
  % normal that points into each body (ie. outward normal)
  nx = nx(:,ones(Ng,1))';
  ny = ny(:,ones(Ng,1))';
%  kernel = diffx.*(tysou(ones(Ng,1),:)) - ...
%            diffy.*(txsou(ones(Ng,1),:));

  sa = geom.sa(:,k)'; % Jacobian
  cur = geom.cur(:,k)'; % curvature

  % target points
  xtar = xx(:,ones(Ng,1));
  ytar = yy(:,ones(Ng,1));

  % source points
  xsou = xtar';
  ysou = ytar';

  % Jacobian
  sa = sa(ones(Ng,1),:);

  diffx = xtar - xsou;
  diffy = ytar - ysou;
  
  % r dot n term
  rdotn = diffx.*nx + diffy.*ny;

  r2 = diffx.^2 + diffy.^2;
  r12 = r2.^(0.5);
  r12(1:Ng+1:Ng^2) = 0;

  % set diagonal terms to 0

  D11 = 1/2/pi/geom.rho.*besselk(1,r12/geom.rho).*...
      rdotn./r12.*sa;
%  omega = 1i*geom.rho;
%  kernel = -1i/4*(omega*r12).*besselh(1,omega*r12);  

  D11(1:Ng+1:Ng^2) = +0.25/pi*cur.*sa(1,:);
  D(:,:,k) = D11;
  
  D(:,:,k) = -D(:,:,k)*2*pi/Ng;
  % scale with the arclength spacing and divide by pi
end % k

end % yukawaDLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GradDmat = yukawaGradDLmatrix(~,geom)
% D = yukawaDLmatrix(geom), generates the gradient of the double-layer
% potential for yukawa. geom is a data structure defined as in the
% capsules class D is (N,N,n) array where Np is the number of points per
% curve and n is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(geom.X);
% Vesicle positions
Ng = geom.N;
% number of points per geoms
GradDmat = zeros(2*Ng,2*Ng,geom.nb);
% initialize space for double-layer potential matrix

for k=1:geom.nb  % Loop over curves
  % locations
  xx = x(:,k);
  yy = y(:,k);

  % normal vector
  normal = [geom.xt(Ng+1:2*Ng,k) -geom.xt(1:Ng,k)]; % Normal vector(y)  
  
%   [tx,ty] = oc.getXY(geom.xt);
%   tx = tx(:,k); ty = ty(:,k);
%   % Vesicle tangent

  % Jacobian
  sa = geom.sa(:,k)';

  % target points
  xtar = xx(:,ones(Ng,1));
  ytar = yy(:,ones(Ng,1));

  % source points
  xsou = xtar';
  ysou = ytar';

%   txsou = tx';
%   tysou = ty';
%   % tangent at sources
  
  % Jacobian
  sa = sa(ones(Ng,1),:);

  % rx and ry terms
  diffx = xtar - xsou;
  diffy = ytar - ysou;
  
  r2 = diffx.^2 + diffy.^2;
  r12 = r2.^(0.5);
  r12(1:Ng+1:Ng^2) = 0;
  % set diagonal terms to 0

  rdotn = (diffx.*normal(:,1) + diffy.*normal(:,2))./r2;

  omega = 1i*geom.rho;
  tmp = omega*r12;
  kernel1 = 1i/4*(omega*r12).^2.*besselh(0,omega*r12);  
  kernel2 = 1i/4*omega*besselh(1,omega*r12)./r12;
  kernel3 = -1i/2*omega.*besselh(1,omega*r12)./r12;

  GradD11 = ((kernel1+kernel3).*diffx.*rdotn + ...
              kernel2.*normal(:,1)).*sa;
  % (1,1) component
  GradD11(1:Ng+1:Ng^2) = 0;
  % diagonal limiting term

%   figure(2);
%   surf(GradD11*2*pi/Ng)
%   shading interp
% plot(GradD11(:,10))
%   pause  
  
  GradD12 = zeros(Ng);
  % (1,2) component

  GradD22 = ((kernel1+kernel3).*diffy.*rdotn + ...
              kernel2.*normal(:,2)).*sa;
  % (2,2) component
  GradD22(1:Ng+1:Ng^2) = 0;
  % diagonal limiting term

  GradDmat(:,:,k) = [GradD11 GradD12; GradD12 GradD22];
  % build matrix with four blocks
  
  GradDmat(:,:,k) = GradDmat(:,:,k)*2*pi/Ng;
  % scale with the arclength spacing
  
end % k

end % yukawaGradDLmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DLP = exactStokesDLdiag(~,geom,D,f)
% DLP = exactStokesDLdiag(geom,f,K) computes the diagonal term of the
% double-layer potential due to f around all geoms. Source and target
% points are the same. This uses trapezoid rule with the curvature at
% the diagonal in order to guarantee spectral accuracy.  This routine
% can either compute the double-layer potential matrix-free, which may
% upsample the number of source points. Or, if the matrix D is passed in
% and anti-aliasing is not requested, it will simply do the
% matrix-vector product with the precomputed matrix D.

DLP = zeros(2*geom.N,geom.nb);
for k = 1:geom.nb
  A = D(:,:,k);
  DLP(:,k) = A * f(:,k);
end

end % exactStokesDLdiag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N0 = exactStokesN0diag(~,wall,N0,f)
% DLP = exactStokesN0diag(vesicle,f) computes the diagonal term of the
% modification of the double-layer potential due to f around outermost
% wall.  Source and target points are the same.  This uses trapezoid
% rule
if isempty(N0)
  Nf = size(f,1)/2;
  oc = curve;
  [fx,fy] = oc.getXY(f(:,1));
  fx = fx.*wall.sa(:,1);
  fy = fy.*wall.sa(:,1);
  [tx,ty] = oc.getXY(wall.xt(:,1));
  % tangent vector
  const = sum(ty.*fx - tx.*fy)*2*pi/Nf;
  % function to be integrated is dot product of normal with density
  % function
  N0 = zeros(2*Nf,1);
  N0 = const*[ty;-tx];
else
  N0 = N0(:,:,1)*f(:,1);
end

end % exactStokesN0diag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DLP = exactYukawaDLdiag(~,geom,D,f)
% DLP = exactYukawaDLdiag(geom,f,K) computes the diagonal term of the
% double-layer potential due to f around all geoms. Source and target
% points are the same. This uses trapezoid rule with the curvature at
% the diagonal in order to guarantee spectral accuracy. This routine
% can either compute the double-layer potential matrix-free, which may
% upsample the number of source points. Or, if the matrix D is passed in
% and anti-aliasing is not requested, it will simply do the
% matrix-vector product with the precomputed matrix D.

DLP = zeros(geom.N,geom.nb);
for k = 1:geom.nb
  A = D(:,:,k);
  DLP(:,k) = A * f(:,k);
end

end % exactYukawaDLdiag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES == TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES ~= TARGETS.  CAN COMPUTE LAYER POTENTIAL ON EACH
% VESICLE DUE TO ALL OTHER VESICLES (ex. stokesSLP) AND CAN
% COMPUTE LAYER POTENTIAL DUE TO VESICLES INDEXED IN K1 AT 
% TARGET POINTS Xtar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = exactStokesDL(~,geom,f,Xtar,K1)
% [stokesDLP,stokesDLPtar] = exactStokesDL(geom,f,Xtar,K1) computes the
% double-layer potential due to f around all parts of the geometry
% except itself.  Also can pass a set of target points Xtar and a
% collection of geom K1 and the double-layer potential due to components
% of the geometry in K1 will be evaluated at Xtar.  Everything but Xtar
% is in the 2*N x n format Xtar is in the 2*Ntar x ncol format

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stokesDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, the user does not need the velocity at arbitrary
  % points
end
den = f.*[geom.sa;geom.sa]*2*pi/geom.N;
% jacobian term and 2*pi/N accounted for here

oc = curve;
[xsou,ysou] = oc.getXY(geom.X(:,K1));
xsou = xsou(:); ysou = ysou(:);
xsou = xsou(:,ones(Ntar,1))';
ysou = ysou(:,ones(Ntar,1))';

[denx,deny] = oc.getXY(den(:,K1));
denx = denx(:); deny = deny(:);
denx = denx(:,ones(Ntar,1))';
deny = deny(:,ones(Ntar,1))';

[normalx,normaly] = oc.getXYperp(geom.xt(:,K1));
normalx = normalx(:); normaly = normaly(:);
normalx = normalx(:,ones(Ntar,1))';
normaly = normaly(:,ones(Ntar,1))';

for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  xtar = xtar(:,ones(geom.N*numel(K1),1));
  ytar = ytar(:,ones(geom.N*numel(K1),1));
  
  diffx = xtar - xsou; diffy = ytar - ysou;
  dis2 = (diffx).^2 + (diffy).^2;
  % difference of source and target location and distance squared
  
  rdotnTIMESrdotf = (diffx.*normalx + diffy.*normaly)./dis2.^2 .* ...
      (diffx.*denx + diffy.*deny);
  % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term
  
  stokesDLPtar(1:Ntar,k) = stokesDLPtar(1:Ntar,k) + ...
      sum(rdotnTIMESrdotf.*diffx,2);
  stokesDLPtar(Ntar+1:end,k) = stokesDLPtar(Ntar+1:end,k) + ...
      sum(rdotnTIMESrdotf.*diffy,2);
  % r \otimes r term of the double-layer potential
end
stokesDLPtar = stokesDLPtar/pi;
% double-layer potential due to geometry components indexed over K1
% evaluated at arbitrary points

stokesDLP = zeros(2*geom.N,geom.nb);
if (nargin == 4 && geom.nb > 1)
  for k = 1:geom.n
    K = [(1:k-1) (k+1:geom.nb)];
    [x,y] = oc.getXY(geom.X(:,K));
    [nx,ny] = oc.getXYperp(geom.xt(:,K));
    [denx,deny] = oc.getXY(den(:,K));
    for j=1:geom.N
      diffxy = [geom.X(j,k) - x ; geom.X(j+geom.N,k) - y];
      dis2 = diffxy(1:geom.N,:).^2 + ...
          diffxy(geom.N+1:2*geom.N,:).^2;
      % difference of source and target location and distance squared

      rdotfTIMESrdotn = ...
        (diffxy(1:geom.N,:).*nx + ...
        diffxy(geom.N+1:2*geom.N,:).*ny)./dis2.^2 .* ...
        (diffxy(1:geom.N,:).*denx + ...
        diffxy(geom.N+1:2*geom.N,:).*deny);
      % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term

      stokesDLP(j,k) = stokesDLP(j,k) + ...
          sum(sum(rdotfTIMESrdotn.*diffxy(1:geom.N,:)));
      stokesDLP(j+geom.N,k) = stokesDLP(j+geom.N,k) + ...
          sum(sum(rdotfTIMESrdotn.*diffxy(geom.N+1:2*geom.N,:)));
      % double-layer potential for Stokes
    end
  end

  stokesDLP = stokesDLP/pi;
  % 1/pi is the coefficient in front of the double-layer potential
end
% double-layer potential due to all components of the geometry except
% oneself

end % exactStokesDL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLPpressure,stokesDLPpressuretar] = ...
      exactStokesDLpressure(~,geom,f,Xtar,K1)
% [stokesDLPpressure,stokesDLPpressuretar] =
% exactStokesDLpressure(geom,f,Xtar,K1) computes the pressure tensor of
% the double-layer potential due to f around all parts of the geometry
% except itself. Also can pass a set of target points Xtar and a
% collection of geom K1 and the double-layer potential due to components
% of the geometry in K1 will be evaluated at Xtar. Everything but Xtar
% is in the 2*N x n format Xtar is in the 2*Ntar x ncol format

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesDLPpresstar = zeros(Ntar,ncol);
else
  K1 = [];
  stokesDLPpressuretar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, the user does not need the pressure at arbitrary
  % points
end
den = f.*[geom.sa;geom.sa]*2*pi/geom.N;
% jacobian term and 2*pi/N accounted for here

oc = curve;
[xsou,ysou] = oc.getXY(geom.X(:,K1));
xsou = xsou(:); ysou = ysou(:);
xsou = xsou(:,ones(Ntar,1))';
ysou = ysou(:,ones(Ntar,1))';

[denx,deny] = oc.getXY(den(:,K1));
denx = denx(:); deny = deny(:);
denx = denx(:,ones(Ntar,1))';
deny = deny(:,ones(Ntar,1))';

[nx,ny] = oc.getXYperp(geom.xt(:,K1));
nx = nx(:); ny = ny(:);
nx = nx(:,ones(Ntar,1))';
ny = ny(:,ones(Ntar,1))';

for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  xtar = xtar(:,ones(geom.N*numel(K1),1));
  ytar = ytar(:,ones(geom.N*numel(K1),1));
  
  rx = xtar - xsou; ry = ytar - ysou;
  rho2 = rx.^2 + ry.^2;
  % difference of source and target location and distance squared

  ndotf = nx.*denx + ny.*deny;
  rdotn = rx.*nx + ry.*ny;
  rdotf = rx.*denx + ry.*deny;
  % dot products in the equation for the pressure of the DLP

  stokesDLPpressuretar(1:Ntar,k) = sum(...
      ndotf./rho2 - 2*rdotn.*rdotf./rho2.^2,2);
  % pressure
end
stokesDLPpressuretar = -stokesDLPpressuretar/pi;
% pressure of the double-layer potential due to geometry components
% indexed over K1 evaluated at arbitrary points

stokesDLPpressure = [];
% April 1, 2021, BQ: THIS CODE HAS NOT BEEN UPDATED SINCE I DON'T
% BELIEVE WE'LL USE IT
%if (nargin == 4 && geom.nb > 1)
%  for k = 1:geom.n
%    K = [(1:k-1) (k+1:geom.nb)];
%    [x,y] = oc.getXY(geom.X(:,K));
%    [nx,ny] = oc.getXYperp(geom.xt(:,K));
%    [denx,deny] = oc.getXY(den(:,K));
%    for j=1:geom.N
%      diffxy = [geom.X(j,k) - x ; geom.X(j+geom.N,k) - y];
%      dis2 = diffxy(1:geom.N,:).^2 + ...
%          diffxy(geom.N+1:2*geom.N,:).^2;
%      % difference of source and target location and distance squared
%
%      rdotfTIMESrdotn = ...
%        (diffxy(1:geom.N,:).*nx + ...
%        diffxy(geom.N+1:2*geom.N,:).*ny)./dis2.^2 .* ...
%        (diffxy(1:geom.N,:).*denx + ...
%        diffxy(geom.N+1:2*geom.N,:).*deny);
%      % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term
%
%      stokesDLP(j,k) = stokesDLP(j,k) + ...
%          sum(sum(rdotfTIMESrdotn.*diffxy(1:geom.N,:)));
%      stokesDLP(j+geom.N,k) = stokesDLP(j+geom.N,k) + ...
%          sum(sum(rdotfTIMESrdotn.*diffxy(geom.N+1:2*geom.N,:)));
%      % double-layer potential for Stokes
%    end
%  end
%
%  stokesDLP = stokesDLP/pi;
%  % 1/pi is the coefficient in front of the double-layer potential
%end
%% double-layer potential due to all components of the geometry except
%% oneself

end % exactStokesDLstress


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLPstress,stokesDLPstresstar] = ...
      exactStokesDLstress(~,geom,f,Xtar,K1)
% [stokesDLPstress,stokesDLPstresstar] =
% exactStokesDLstress(geom,f,Xtar,K1) computes the stress tensor of the
% double-layer potential due to f around all parts of the geometry
% except itself. Also can pass a set of target points Xtar and a
% collection of geom K1 and the double-layer potential due to components
% of the geometry in K1 will be evaluated at Xtar. Everything but Xtar
% is in the 2*N x n format Xtar is in the 2*Ntar x ncol format

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesDLPstresstar = zeros(3*Ntar,ncol);
else
  K1 = [];
  stokesDLPstresstar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, the user does not need the stress at arbitrary
  % points
end
den = f.*[geom.sa;geom.sa]*2*pi/geom.N;
% jacobian term and 2*pi/N accounted for here

oc = curve;
[xsou,ysou] = oc.getXY(geom.X(:,K1));
xsou = xsou(:); ysou = ysou(:);
xsou = xsou(:,ones(Ntar,1))';
ysou = ysou(:,ones(Ntar,1))';

[denx,deny] = oc.getXY(den(:,K1));
denx = denx(:); deny = deny(:);
denx = denx(:,ones(Ntar,1))';
deny = deny(:,ones(Ntar,1))';

[nx,ny] = oc.getXYperp(geom.xt(:,K1));
nx = nx(:); ny = ny(:);
% need a negative so normal points into the body
nx = nx(:,ones(Ntar,1))';
ny = ny(:,ones(Ntar,1))';

for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  xtar = xtar(:,ones(geom.N*numel(K1),1));
  ytar = ytar(:,ones(geom.N*numel(K1),1));
  
  rx = xtar - xsou; ry = ytar - ysou;
  rho2 = rx.^2 + ry.^2;
  % difference of source and target location and distance squared

  ndotf = nx.*denx + ny.*deny;
  rdotn = rx.*nx + ry.*ny;
  rdotf = rx.*denx + ry.*deny;
  % dot products in the equation for the stress of the DLP

  stokesDLPstresstar(1:Ntar,k) = sum(...
      +1*ndotf./rho2 ...
      -8*rdotn.*rdotf.*rx.^2./rho2.^3 ...
      +2*rdotn.*rx.*denx./rho2.^2 ...
      +2*rdotf.*rx.*nx./rho2.^2,2);
  % (1,1) component of stress tensor

  stokesDLPstresstar(Ntar+1:2*Ntar,k) = sum(...
      -8*rdotn.*rdotf.*rx.*ry./rho2.^3 ...
      +1*rdotn.*(rx.*deny + ry.*denx)./rho2.^2 ...
      +1*rdotf.*(rx.*ny + ry.*nx)./rho2.^2,2);
  % (2,1) or (1,2) (it's symmetric) component of stress tensor

  stokesDLPstresstar(2*Ntar+1:3*Ntar,k) = sum(...
      +1*ndotf./rho2 ...
      -8*rdotn.*rdotf.*ry.^2./rho2.^3 ...
      +2*rdotn.*ry.*deny./rho2.^2 ...
      +2*rdotf.*ry.*ny./rho2.^2,2);
  % (2,2) component of stress tensor
end
stokesDLPstresstar = stokesDLPstresstar/pi;
% stress of the double-layer potential due to geometry components
% indexed over K1 evaluated at arbitrary points

%stokesDLPstress = zeros(3*geom.N,geom.nb);
stokesDLPstress = [];
% April 1, 2021, BQ: THIS CODE HAS NOT BEEN UPDATED SINCE I DON'T
% BELIEVE WE'LL USE IT
%if (nargin == 4 && geom.nb > 1)
%  for k = 1:geom.n
%    K = [(1:k-1) (k+1:geom.nb)];
%    [x,y] = oc.getXY(geom.X(:,K));
%    [nx,ny] = oc.getXYperp(geom.xt(:,K));
%    [denx,deny] = oc.getXY(den(:,K));
%    for j=1:geom.N
%      diffxy = [geom.X(j,k) - x ; geom.X(j+geom.N,k) - y];
%      dis2 = diffxy(1:geom.N,:).^2 + ...
%          diffxy(geom.N+1:2*geom.N,:).^2;
%      % difference of source and target location and distance squared
%
%      rdotfTIMESrdotn = ...
%        (diffxy(1:geom.N,:).*nx + ...
%        diffxy(geom.N+1:2*geom.N,:).*ny)./dis2.^2 .* ...
%        (diffxy(1:geom.N,:).*denx + ...
%        diffxy(geom.N+1:2*geom.N,:).*deny);
%      % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term
%
%      stokesDLP(j,k) = stokesDLP(j,k) + ...
%          sum(sum(rdotfTIMESrdotn.*diffxy(1:geom.N,:)));
%      stokesDLP(j+geom.N,k) = stokesDLP(j+geom.N,k) + ...
%          sum(sum(rdotfTIMESrdotn.*diffxy(geom.N+1:2*geom.N,:)));
%      % double-layer potential for Stokes
%    end
%  end
%
%  stokesDLP = stokesDLP/pi;
%  % 1/pi is the coefficient in front of the double-layer potential
%end
%% double-layer potential due to all components of the geometry except
%% oneself

end % exactStokesDLstress

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vLets,vLetsTar] = StokesletRotlet(o,geom,force,torque,Xtar,K1)
% [vLets,vLetsTar] = StokesletRotlet(geom,force,torque,Xtar,K1) computes
% the velocity due to stokeslets and rotlets. Also can pass a set of
% target points Xtar and a collection of geom K1 and the Stokeslet and
% Rotlet velocities due to components of the geometry in K1 will be
% evaluated at Xtar. Everything but Xtar is in the 2*nb x n format Xtar
% is in the 2*Ntar x ncol format

if nargin == 6
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  vLetsTar = zeros(2*Ntar,ncol);
else
  K1 = [];
  vLetsTar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 6, the user does not need the velocity at arbitrary
  % points
end

oc = curve;
X = geom.X;
N = size(X,1)/2;
nb = size(X,2);

% stokeslet and torques at the target points (if there are any)
for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  for j = K1
    [cx,cy] = oc.getXY(geom.center(:,j));
    fx = force(2*(j-1) + 1);
    fy = force(2*(j-1) + 2);
    tor = torque(j);
    % force and torque due to body j

    rx = xtar - cx;
    ry = ytar - cy;
    rho2 = rx.^2 + ry.^2;
    rdotf = rx*fx + ry*fy;
    vLetsTar(:,k) = vLetsTar(:,k) + (1/4/pi)*...
      [-0.5*log(rho2)*fx + rdotf./rho2.*rx; ...
      -0.5*log(rho2)*fy + rdotf./rho2.*ry];

    vLetsTar = vLetsTar + (1/4/pi)*tor*[-ry./rho2;+rx./rho2]; 
  end
end


% stokeslet and torques along the bodies
vLets = zeros(2*N,nb);

[x,y] = oc.getXY(X);
for j = 1:nb
  [cx,cy] = oc.getXY(geom.center(:,j));

  fx = force(2*(j-1) + 1);
  fy = force(2*(j-1) + 2);
  tor = torque(j);
  % force and torque due to body j

  rx = x - cx;
  ry = y - cy;
  rho2 = rx.^2 + ry.^2;
  rdotf = rx*fx + ry*fy;
  vLets = vLets + 1/4/pi*...
    [-0.5*log(rho2)*fx + rdotf./rho2.*rx; ...
     -0.5*log(rho2)*fy + rdotf./rho2.*ry];

  vLets = vLets + (1/4/pi)*tor*[-ry./rho2;+rx./rho2]; 
end
%clf; hold on
%plot(vLets(1:end/2,2));
%plot(vLets(1+end/2:end,2),'r');
%quiver(x,y,vLets(1:end/2,:),vLets(end/2+1:end,:))
%[x(:,2) y(:,2) vLets(1:end/2,2) vLets(end/2+1:end,2)]
%pause

end % StokesletRotlet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pressure = RSpressure(~,geom,force,torque,Xtar,K1)
% pressure = RSpressure(geom,force,torque,Xtar,K1) computes the pressure
% due to the Stokeslets and Rotlets at a collection of target points.
% Only looks at contribution due to bodies indexed inside K1.

Ntar = size(Xtar,1)/2;
ncol = size(Xtar,2);
pressure = zeros(Ntar,ncol);

oc = curve;

for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  for j = K1
    [cx,cy] = oc.getXY(geom.center(:,j));
    fx = force(2*(j-1) + 1);
    fy = force(2*(j-1) + 2);
%    tor = torque(j);
    % force and torque due to body j

    rx = xtar - cx;
    ry = ytar - cy;
    rho2 = rx.^2 + ry.^2;
    rdotf = rx*fx + ry*fy;
    % pressure of the stokeslet
    pressure(:,k) = pressure(:,k) + 1/2/pi*rdotf./rho2;

    % ROTLET HAS NO PRESSURE.
  end
end
% pressure of the Rotlet and Stokeslets due to geometry components
% indexed over K1 evaluated at arbitrary points

end % RSpressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stress = RSstress(~,geom,force,torque,Xtar,K1)
% stress = RSstress(geom,force,torque,Xtar,K1) computes the stress due
% to the Stokeslets and Rotlets at a collection of target points. Only
% looks at contribution due to bodies indexed inside K1.

Ntar = size(Xtar,1)/2;
ncol = size(Xtar,2);
stress = zeros(3*Ntar,ncol);

oc = curve;

for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  for j = K1
    [cx,cy] = oc.getXY(geom.center(:,j));
    fx = force(2*(j-1) + 1);
    fy = force(2*(j-1) + 2);
    tor = torque(j);
    % force and torque due to body j

    rx = xtar - cx;
    ry = ytar - cy;
    rho2 = rx.^2 + ry.^2;
    rdotf = rx*fx + ry*fy;

    stress(1:Ntar,k) = stress(1:Ntar,k) - ...
        1/pi*rdotf./rho2.^2.*rx.^2 + ...
        1/pi*rx.*ry./rho2.^2*tor;
    stress(Ntar+1:2*Ntar,k) = stress(Ntar+1:2*Ntar,k) - ...
        1/pi*(rx.*ry.*rdotf)./rho2.^2 + ...
        1/2/pi*tor*(ry.^2 - rx.^2)./rho2.^2;
    stress(2*Ntar+1:3*Ntar,k) = stress(2*Ntar+1:3*Ntar,k) - ...
        1/pi*rdotf./rho2.^2.*ry.^2 - ...
        1/pi*rx.*ry./rho2.^2*tor;

  end
end
% pressure of the Rotlet and Stokeslets due to geometry components
% indexed over K1 evaluated at arbitrary points

end % RSstress



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yukawaDLP,yukawaDLPtar] = exactYukawaDL(~,geom,f,Xtar,K1)
% [yukawaDLP,yukawaDLPtar] = exactYukawaDL(geom,f,Xtar,K1) computes the
% Yukawa double-layer potential due to f around all parts of the
% geometry except itself. Also can pass a set of target points Xtar and
% a collection of geom K1 and the double-layer potential due to
% components of the geometry in K1 will be evaluated at Xtar.
% Everything but Xtar is in the 2*N x n format Xtar is in the 2*Ntar x
% ncol format

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  yukawaDLPtar = zeros(Ntar,ncol);
else
  K1 = [];
  yukawaDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, the user does not need the velocity at arbitrary
  % points
end
den = f.*geom.sa*2*pi/geom.N;
% jacobian term and 2*pi/N accounted for here

oc = curve;
[xsou,ysou] = oc.getXY(geom.X(:,K1));
xsou = xsou(:); ysou = ysou(:);
xsou = xsou(:,ones(Ntar,1))';
ysou = ysou(:,ones(Ntar,1))';

den = den(:,K1);
den = den(:);
den = den(:,ones(Ntar,1))';

[tx,ty] = oc.getXY(geom.xt(:,K1));
nx = -ty(:); ny = tx(:);
normalx = nx(:,ones(Ntar,1))';
normaly = ny(:,ones(Ntar,1))';

for k = 1:ncol % loop over columns of target points
  [xtar,ytar] = oc.getXY(Xtar(:,k));
  xtar = xtar(:,ones(geom.N*numel(K1),1));
  ytar = ytar(:,ones(geom.N*numel(K1),1));
  
  diffx = xtar - xsou; diffy = ytar - ysou;
  % difference of source and target location
  dis2 = diffx.^2 + diffy.^2;
  % distance squared
  dis = sqrt(dis2);
  % distance
  rdotn = diffx.*normalx + diffy.*normaly;
  % difference dotted with normal

  BK1 = besselk(1,dis/geom.rho);
  kernel = -1/2/pi/geom.rho.*BK1.*rdotn./dis.*den;
%  kernel = -1/2/pi/geom.rho.*besselk(1,dis/geom.rho).*...
%        rdotn./dis.*den;
  % BQ: NOT SURE WHY WE NEED A MINUS SIGN HERE, BUT IT MAKES IT WORK
  
  yukawaDLPtar(:,k) = yukawaDLPtar(:,k) + sum(kernel,2);
  % Yukawa DLP
end
% double-layer potential due to geometry components indexed over K1
% evaluated at arbitrary points

yukawaDLP = zeros(geom.N,geom.nb);
if (nargin == 3 && geom.N > 1)
  for k = 1:geom.nb
    K = [(1:k-1) (k+1:geom.nb)];
    [x,y] = oc.getXY(geom.X(:,K));
    [tx,ty] = oc.getXY(geom.xt(:,K));
    nx = -ty; ny = tx;
    den = f(:,K).*geom.sa(:,K)*2*pi/geom.N;
    % density including the jacobian and arclength term
    for j = 1:geom.N
      diffx = geom.X(j,k) - x; diffy = geom.X(j+geom.N,k) - y;
      % difference of source and target location
      dis2 = diffx.^2 + diffy.^2;
      % distance squared
      dis = sqrt(dis2);
      % distance
      rdotn = diffx.*nx + diffy.*ny;
      % difference dotted with normal

      % double-layer potential for Yukawa
      yukawaDLP(j,k) = yukawaDLP(j,k) - ...
          sum(sum(besselk(1,dis/geom.rho).*rdotn./dis.*den));

      % BQ: NOT SURE WHY WE NEED A MINUS SIGN HERE, BUT IT MAKES IT WORK
    end
  end

  yukawaDLP = yukawaDLP/2/pi/geom.rho;
  % 1/2/pi is the coefficient in front of the double-layer potential
  % 1/geom.rho comes out of the chain rule
end
% double-layer potential due to all components of the geometry except
% oneself

end % exactYukawaDL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES ~= TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT CALCULATE FORCES USING JUMP RELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F1, F2, Tq] = evalForcesTaylor(o,geom,eta)
% NOTE: THIS CODE IS CURRENTLY NOT BEING CALLED, SO NO NEED TO CLEAN IT
% UP. IDEA IS TO USE TAYLOR EXPANSIONS RATHER THAN QBX

% Uses * Tpq + Tqp identity, 
%      * Jpq jump value to evaluate on curve itself 
%      * Taylor expansion to derive uq and uq_n on curve p, p ~= q
%      * FFT for tangential derivatives     

oc = curve;
N = geom.N;
[x1,x2] = oc.getXY(geom.X);
[tau1,tau2] = oc.getXY(geom.xt);
nu1 = +tau2;
nu2 = -tau1;
dS = geom.sa*2*pi/N;

F1 = zeros(geom.nb,1);
F2 = zeros(geom.nb,1);
Tq = zeros(geom.nb,1);

tol   = 1.0;

p_max = 8;
r_max = 0.2;
r_m   = linspace(0, r_max, p_max+1);
r_m   = r_m(2:end)';
r_0   = r_m(p_max);
pw    = (0:p_max-1);
AA    = (r_m - r_0).^pw;
BB    = (0 - r_0).^pw;

IK = oc.modes(geom.N,1); % Fourier modes

for p = 1:geom.nb
  for q = [1:p-1, p+1:geom.nb]
    x1p   = x1(:,p);
    x2p   = x2(:,p);        
    dSp   = dS(:,p);
    nu1p  = nu1(:,p);
    nu2p  = nu2(:,p);
    tau1p = -nu2p;
    tau2p =  nu1p;
    etap    = eta(:, p);

    % decide which strategy to use : expansion or standard evaluation 
    D = min(pdist2([x1(:,p) x2(:,p)],[x1(:,q), x2(:,q)]),[],"all");
               
    if D < tol
      [Uq,Uq_x1,Uq_x2] = o.evalDL(x1p - nu1p*r_m',x2p - nu2p*r_m', ...
         1,N,x1(:,q),x2(:,q),nu1(:,q),nu2(:,q),dS(:,q),geom.rho,eta(:,q));

      coef    = AA\Uq.';
      coef_x1 = AA\Uq_x1.';
      coef_x2 = AA\Uq_x2.';                

      uq      = (BB*coef)';
      uq_x1   = (BB*coef_x1)';
      uq_x2   = (BB*coef_x2)';

    else

      [uq,uq_x1,uq_x2] = o.evalDL(x1p,x2p,1,N,x1(:,q),x2(:,q),...
          nu1(:,q),nu2(:,q),dS(:,q),geom.rho,eta(:,q));

    end % if D < tol                                               

%    etap_t = o.tanDeriv(etap,geom.sa(:,p),IK);
%    uq_t = o.tanDeriv(uq,geom.sa(:,p),IK);

    % compute tangent derivatives of etap and uq.
    etap_t = oc.diffFT(etap,IK)./geom.sa(:,p);
    uq_t = oc.diffFT(uq,IK)./geom.sa(:,p);
    uq_n = uq_x1.*nu1p + uq_x2.*nu2p;

    Jpq1 = 2.0/geom.rho*etap.*uq.*nu1p + ...
        2.0*geom.rho*etap_t.*uq_t.*nu1p - ...
        2.0*geom.rho*etap_t.*uq_n.*tau1p;
    Jpq2 = 2.0/geom.rho*etap.*uq.*nu2p + ...
        2.0*geom.rho*etap_t.*uq_t.*nu2p - ...
        2.0*geom.rho*etap_t.*uq_n.*tau2p;        

    F1(p) = F1(p) + sum(Jpq1.*dSp);

    F2(p) = F2(p) + sum(Jpq2.*dSp);

    Tq(p) = Tq(p) + sum( ( x1p.*Jpq2 - x2p.*Jpq1 ) .* dSp );        
             
  end    
end     

end % evalForcesTaylor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F1, F2, Tq] = evalForcesQBX(o,geom,eta)

% Uses * Tpq + Tqp identity, 
%      * Jpq jump value to evaluate on curve itself 
%      * QBX expansion to derive uq and uq_n on curve p, p ~= q
%      * FFT for tangential derivatives     

oc  = curve;
N   = geom.N;
gam = geom.gam;

[x1,x2]     = oc.getXY(geom.X);
[tau1,tau2] = oc.getXY(geom.xt);
nu1 = +tau2;
nu2 = -tau1;
dS = geom.sa*2*pi/N;

F1 = zeros(geom.nb,1);
F2 = zeros(geom.nb,1);
Tq = zeros(geom.nb,1);

pc1 = geom.center(1,:).';
pc2 = geom.center(2,:).';

rad   = 0.3;
m_max = 6;

if N <= 32
  rad = 0.4; m_max = 2;
elseif N <= 64
  rad = 0.3; m_max = 6;
elseif N <= 128
  rad = 0.3; m_max = 12;
else 
  rad = 0.2; m_max = 12;
end

tol   = rad;
B_m   = zeros(N,2*m_max + 1);

IK = oc.modes(N,1); % Fourier modes;

for p = 1:geom.nb
  for q = [1:p-1, p+1:geom.nb]
    x1pc = pc1(p);
    x2pc = pc2(p);    

    x1p  = x1(:,p);
    x2p  = x2(:,p);        
    dSp  = dS(:,p);
    
    nu1p = nu1(:,p);
    nu2p = nu2(:,p);
    
    tau1p = tau1(:,p);
    tau2p = tau2(:,p);

    % setting argin geom.nb = 1 tricks evaluators to using only one geometry
    % column
    etap = eta(:, p);

    % decide which strategy to use : expansion or standard evaluation 
              
    D = min(pdist2([x1(:,p) x2(:,p)],[x1(:,q), x2(:,q)]),[],"all");

    if D < tol
      uq = zeros(N,1);
      uq_x1 = zeros(N,1);
      uq_x2 = zeros(N,1);

      c1 = x1p - rad*nu1p;
      c2 = x2p - rad*nu2p;

      % Form QBX coefficients
      for m = -m_max:m_max
        B_m(:,m+m_max+1) = o.QBX_coeff(c1,c2,m,x1(:,q),x2(:,q), ...
              nu1(:,q),nu2(:,q),dS(:,q),geom.rho,eta(:,q));
      end                

      % Evaluate QBX expansion
      [uq,uq_x1,uq_x2] = o.QBX_exp(x1p,x2p,c1,c2,B_m,m_max,geom.rho); 

    else
      % Evaluate DL potential directly
      [uq,uq_x1,uq_x2] = o.evalDL(x1p,x2p,1,N,x1(:,q),x2(:,q), ...
            nu1(:,q),nu2(:,q),dS(:,q),geom.rho,eta(:,q));
    end % if D < tol                                   

%    etap_t = o.tanDeriv(etap,geom.sa(:,p),IK);
%    uq_t = o.tanDeriv(uq,geom.sa(:,p),IK);

    % compute tangent derivatives of etap and uq.
    etap_t = oc.diffFT(etap,IK)./geom.sa(:,p);
    uq_t = oc.diffFT(uq,IK)./geom.sa(:,p);
    uq_n = uq_x1.*nu1p + uq_x2.*nu2p;

    Jpq1 = 2.0/geom.rho*etap.*uq.*nu1p + ...
        2.0*geom.rho*etap_t.*uq_t.*nu1p - ...
        2.0*geom.rho*etap_t.*uq_n.*tau1p;
    Jpq2 = 2.0/geom.rho*etap.*uq.*nu2p + ...
        2.0*geom.rho*etap_t.*uq_t.*nu2p - ...
        2.0*geom.rho*etap_t.*uq_n.*tau2p;        

    F1(p) = F1(p) + sum(Jpq1.*dSp);
    F2(p) = F2(p) + sum(Jpq2.*dSp);
    Tq(p) = Tq(p) + sum(( + (x1p - pc1(p)).* Jpq2 - (x2p - pc2(p)).* Jpq1 ) .* dSp);
    
  end    
end     

%F1, F2 and Tq have numerically mean zero

F1 = gam*F1;
F2 = gam*F2;
Tq = gam*Tq;

end % evalForcesQBX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F1, F2, Tq] = evalForcesTaylorOld(o,Nb,N,x1,x2,nu1,nu2,dS,rho,h)
% NOTE: THIS CODE IS CURRENTLY NOT BEING CALLED, SO NO NEED TO CLEAN IT
% UP. IDEA IS TO USE TAYLOR EXPANSIONS RATHER THAN QBX

% Uses * Tpq + Tqp identity, 
%      * Jpq jump value to evaluate on curve itself 
%      * Taylor expansion to derive uq and uq_n on curve p, p ~= q
%      * FFT for tangential derivatives     

F1 = zeros(Nb,1);
F2 = zeros(Nb,1);
Tq = zeros(Nb,1);

tol   = 1.0;

p_max = 8;
r_max = 0.2;
r_m   = linspace(0, r_max, p_max+1);
r_m   = r_m(2:end)';
r_0   = r_m(p_max);
pw    = (0:p_max-1);
AA    = (r_m - r_0).^pw;
BB    = (0 - r_0).^pw;

%errors = zeros(Nb, 3);

IK = 1i*[(0:N/2-1) 0 (-N/2+1:-1)]'; % Fourier modes    
dT = 2*pi/N;     % to convert dS to sa

for p = 1:Nb
  for q = [1:p-1, p+1:Nb]

    x1p   = x1(:,p);
    x2p   = x2(:,p);        
    dSp   = dS(:,p);
    nu1p  = nu1(:,p);
    nu2p  = nu2(:,p);
    tau1p = -nu2p;
    tau2p =  nu1p;
    hp    = h(:, p);

    % decide which strategy to use : expansion or standard evaluation 
              
    D = pdist2([x1(:,p) x2(:,p)],[x1(:,q), x2(:,q)]);
    D = min(D, [], "all");
               
    if D < tol
      [Uq,Uq_x1,Uq_x2] = o.evalDL(x1p - nu1p*r_m',x2p - nu2p*r_m', ...
         1,N,x1(:,q),x2(:,q),nu1(:,q),nu2(:,q),dS(:,q),rho,h(:,q));

      coef    = AA\Uq.';
      coef_x1 = AA\Uq_x1.';
      coef_x2 = AA\Uq_x2.';                

      uq      = (BB*coef)';
      uq_x1   = (BB*coef_x1)';
      uq_x2   = (BB*coef_x2)';

    else

      [uq,uq_x1,uq_x2] = o.evalDL(x1p,x2p,1,N,x1(:,q),x2(:,q),...
          nu1(:,q),nu2(:,q),dS(:,q),rho,h(:,q));

    end % if D < tol                                               

    hp_t = o.tanDeriv(hp,dSp/dT,IK);
    uq_t = o.tanDeriv(uq,dSp/dT,IK);
    uq_n = uq_x1.*nu1p + uq_x2.*nu2p;

    Jpq1 = 2.0/rho*hp.*uq.*nu1p + 2.0*rho*hp_t.*uq_t.*nu1p - ...
        2.0*rho*hp_t.*uq_n.*tau1p;
    Jpq2 = 2.0/rho*hp.*uq.*nu2p + 2.0*rho*hp_t.*uq_t.*nu2p - ...
        2.0*rho*hp_t.*uq_n.*tau2p;        

    F1(p) = F1(p) + sum(Jpq1.*dSp);

    F2(p) = F2(p) + sum(Jpq2.*dSp);

    Tq(p) = Tq(p) + sum(( + (x1p - x1pc).*Jpq2 - (x2p - x2pc).*Jpq1 ) .* dSp);

  end    
end     

end % evalForcesTaylorOld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F1,F2,Tq] = evalForcesQBXOld(o, Nb, N, x1, x2, nu1, nu2, dS, rho, h)

% Uses * Tpq + Tqp identity, 
%      * Jpq jump value to evaluate on curve itself 
%      * QBX expansion to derive uq and uq_n on curve p, p ~= q
%      * FFT for tangential derivatives     

F1 = zeros(Nb,1);
F2 = zeros(Nb,1);
Tq = zeros(Nb,1);

rad   = 0.3;
m_max = 6;
tol   = rad;
B_m   = zeros(N,2*m_max + 1);

%errors = zeros(Nb, 3);

IK = 1i*[(0:N/2-1) 0 (-N/2+1:-1)]'; % Fourier modes
dT = 2*pi/N;     % to convert dS to sa

for p = 1:Nb
  for q = [1:p-1, p+1:Nb]

    x1p = x1(:,p);
    x2p = x2(:,p);        
    dSp = dS(:,p);
    nu1p = nu1(:,p);
    nu2p = nu2(:,p);
    tau1p = -nu2p;
    tau2p = +nu1p;

    % setting argin Nb = 1 tricks evaluators to using only one geometry
    % column
    hp = h(:, p);

    % decide which strategy to use : expansion or standard evaluation 
              
    D = pdist2([x1(:,p) x2(:,p)],[x1(:,q), x2(:,q)]);
    D = min(D, [], "all");

    if D < tol
      uq = zeros(N,1);
      uq_x1 = zeros(N,1);
      uq_x2 = zeros(N,1);

      c1 = x1p - rad*nu1p;
      c2 = x2p - rad*nu2p;

      for m = -m_max:m_max
        B_m(:,m+m_max+1) = o.QBX_coeff(c1,c2,m,x1(:,q),x2(:,q), ...
              nu1(:,q),nu2(:,q),dS(:,q),rho,h(:,q));
      end                

      [uq,uq_x1,uq_x2] = o.QBX_exp(x1p,x2p,c1,c2,B_m,m_max,rho); 

    else
      [uq,uq_x1,uq_x2] = o.evalDL(x1p,x2p,1,N,x1(:,q),x2(:,q), ...
            nu1(:,q),nu2(:,q),dS(:,q),rho, h(:,q));
    end % if D < tol                                   

    hp_t = o.tanDeriv(hp,dSp/dT,IK);
    uq_t = o.tanDeriv(uq,dSp/dT,IK);
    uq_n = uq_x1.*nu1p + uq_x2.*nu2p;

    Jpq1 = 2.0/rho*hp.*uq.*nu1p + 2.0*rho*hp_t.*uq_t.*nu1p - ...
        2.0*rho*hp_t.*uq_n.*tau1p;
    Jpq2 = 2.0/rho*hp.*uq.*nu2p + 2.0*rho*hp_t.*uq_t.*nu2p - ...
        2.0*rho*hp_t.*uq_n.*tau2p;        

    F1(p) = F1(p) + sum(Jpq1.*dSp);
    F2(p) = F2(p) + sum(Jpq2.*dSp);
    Tq(p) = Tq(p) + sum(( + (x1p - x1pc).*Jpq2 - (x2p - x2pc).*Jpq1 ) .* dSp);

  end    

end     
    
end % evalForcesQBXOld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dh,Dh_x1,Dh_x2] = QBX_exp(o,X1,X2,c1,c2,A_m,m_max,rho)
%
% (X1, X2)  : target points, arbitrary size 
% (c1, c2)  : expansion center, scalars
% A_m       : 2*m_max + 1 many expanson coefficients calculated elsewhere 
%

Dh = 0*X1;
Dh_x1 = 0*X1;
Dh_x2 = 0*X1;
rcX = sqrt( (X1 - c1).^2 + (X2 - c2).^2);
thcX = atan2(c2 - X2,c1 - X1);          %the order c minus x is important here
for m = -m_max:m_max    
  I0 = besseli( m-1, rcX/rho);
  I1 = besseli( m,   rcX/rho);
  I2 = besseli( m+1, rcX/rho);
  mm = m + m_max + 1;
  Dh = Dh + A_m(:,mm).*I1.*exp(1i*m*thcX);
  Dh_x1 = Dh_x1  + A_m(:,mm).*(1/2*(I0 + I2).*(X1-c1)./(rcX*rho) ...
          + 1i*m*I1.*(-(X2 - c2))./rcX.^2).*exp(1i*m*thcX);
  Dh_x2 = Dh_x2  + A_m(:,mm).*(  1/2*(I0 + I2).*(X2-c2)./(rcX*rho) ...
          + 1i*m*I1.*(+(X1 - c1))./rcX.^2  ).*exp(1i*m*thcX);
end

Dh    = real(Dh);
Dh_x1 = real(Dh_x1);
Dh_x2 = real(Dh_x2);
    
end % QBX_exp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a_m = QBX_coeff(o, c1, c2, m, y1, y2, nu1, nu2, dS, rho, h)

% evaluates the coefficient : vectorized version 
% alpha(m) = int_gamma (-1/2pi) dnu 
%    [ K_{-m}(|y-c|/rho)exp( im(pi - arg(y,c))) ]h(y) dS(y)
% (c1, c2) is the expansion center

% the arc gamma can be one portion of a closed curve, for example, or a
% dispersed set of points with normals, for example. 

% (c1, c2) and (y1, y2) point sets of differing size
% output has size of (c1, c2)

%c1 = reshape(c1, numel(c1),1); % m elements
%c2 = reshape(c2, numel(c2),1); 

%y1 = reshape(y1, 1, numel(y1)); % n elements
%y2 = reshape(y2, 1, numel(y2));
y1 = y1';
y2 = y2';
% transpose is much fater than reshape

%nu1 = reshape(nu1, 1, numel(nu1));
%nu2 = reshape(nu2, 1, numel(nu2));
%dS  = reshape(dS, 1, numel(dS));
%h   = reshape(h, 1, numel(h));
nu1 = nu1';
nu2 = nu2';
dS = dS';
h = h';

r1 = y1 - c1; r2 = y2 - c2; % differences like this are size m by n
r = sqrt(r1.^2 + r2.^2);
rdotnu = +r1.*nu1 + r2.*nu2;
rpotnu = -r2.*nu1 + r1.*nu2; % r-perp

th = atan2(y2-c2,y1-c1);

K0 = besselk(-m-1, r/rho);
K1 = besselk(-m,   r/rho);
K2 = besselk(-m+1, r/rho);
dK = -1/2*(K0 + K2);
    
f = 1/(2*pi)*(dK .* rdotnu./(rho*r) - 1i*m * K1 .* rpotnu./r.^2) ...
      .* exp(1i*m*( pi - th));
a_m = sum(f.*h.*dS, 2); % want m by 1 output, second dimension 

end % QBX_coeff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dh,Dh_X1,Dh_X2] = evalDL(o,X1,X2,Nb,N,x1,x2,nu1,nu2,dS,rho,h)
% TODO: This routine must already be in the code
% evaluates double layer potential at (X1, X2)

Dh = 0*X1;
Dh_X1 = 0*X1;
Dh_X2 = 0*X1; 

for q = 1:Nb
  for j = 1:N
    r1 = X1 - x1(j,q); r2 = X2 - x2(j,q);
    r = sqrt(r1.^2 + r2.^2);                        
    rdotnu = r1.*nu1(j,q) + r2.*nu2(j,q);            

    K0 = besselk(0,r/rho);
    K1 = besselk(1,r/rho);
    K2 = besselk(2,r/rho);
    dK1 = -0.5*(K0 + K2); 
    % see identity
    % https://functions.wolfram.com/Bessel-TypeFunctions/BesselK/20/01/02/

    % using r = x - y; hence lack of minus below
    Dh = Dh + 1/(2*pi)*(r/rho).*K1.*rdotnu./r.^2.*dS(j,q).*h(j,q);

    Dh_X1 = Dh_X1 + 1/(2*pi)*( ...
        + r1./(r*rho).*K1.*rdotnu./r.^2 ...
        + (r/rho).*dK1.*r1./(r*rho).*rdotnu./r.^2 ...
        + (r/rho).*K1.*(1./r.^2).*(nu1(j,q) - 2*r1.*rdotnu./r.^2)) ...
        * dS(j,q)*h(j,q);

    Dh_X2 = Dh_X2 + 1/(2*pi)*( ...
        + r2./(r*rho).*K1.*rdotnu./r.^2 ...
        + (r/rho).*dK1.*r2./(r*rho).*rdotnu./r.^2 ...
        + (r/rho).*K1.*(1./r.^2).*(nu2(j,q) - 2*r2.*rdotnu./r.^2)) ...
        *dS(j,q)*h(j,q);
  end
end

end % evalDL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_t = tanDeriv(o, f,dS,IK)
% calculate tangential derivative of along curve with arclength differential dS
% TODO: THIS ROUTINE IS ALREADY SOMEWHERE IN THE CODE
    
f_t = real(ifft(IK.*fft(f)))./dS;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT CALCULATE FORCES USING JUMP RELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF ROUTINES THAT CALCULATE REPULSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R1, R2, RTq, pp] = Repul(o, geom)

%pp are point pairs (p1, p2, q1, q2) of interacting points

N    = geom.N;                % number of points per componenet
Nb   = geom.nb;               % number of rigid bodies
x1   = geom.X(1:N,:);         % grid points on curves
x2   = geom.X(N+1:2*N,:);
pc   = geom.center;           % center of each rigid body

pc1  = pc(1,:);
pc2  = pc(2,:);

l0 = geom.RepulLength;   % repulsion length
M  = geom.RepulStrength; % repulsion strength

% temporarily alter geom.length parameter  
% for present purposes
h_store = geom.length;          

geom.length = l0*N/2;

NearOther = geom.getZone([],1); % get near structure 

% identifies which pairs of particles have distance < l0
nf = NearOther.nearFibers; 

% structure identifying point in discrete curve q
% that is closest to point in curve p
icp = NearOther.icp; 

R1 = zeros(geom.nb,1); % repulsive force 
R2 = zeros(geom.nb,1); % repulsize force  
RTq = zeros(geom.nb,1); % repulsive torque
pp  = [];

for p = 1:Nb
  [iq, qq, jp] = find(icp{p}); %an N by nb sparse matrix         
  % point iq(l) in curve qq(l) is within l0 of curve p, and 
  % the jp(l) is the point in curve p closest to iq(l) 
   
  q_list = unique(qq); %set of q within l0 of curve p

  for k = 1:length(q_list)
      
    q  = q_list(k); %curve q within l0 of p
    in = find(qq == q); 
      
    i = jp(in); %vertex (i, p) closest to following point 
    j = iq(in); %vertex (j, q) within l0 of p curve

    %goal: find point on q curve closest to current, p curve      
    dist = (x1(i,p) - x1(j,q)).^2 + (x2(i,p) - x2(j,q)).^2;
      
    [dist, arg] = min(dist);
      
    jjp = i(arg);
    jjq = j(arg); 

    [dist,x1nrt,x2nrt,y1nrt,y2nrt] = geom.closestPntPair(geom.X,p,q,jjp,jjq);      

    %(x1nrt, x2nrt) and (y1nrt, y2nrt) are the nearest points in curves
    %p, q respectively
    
    r1 = x1nrt - y1nrt;
    r2 = x2nrt - y2nrt;

    r1 = r1/(dist+eps);
    r2 = r2/(dist+eps);            

    % repulsion profile
    [~, dR] = o.Repul_profile(dist/l0);

    % weight and chain rule
    dR = M*dR/l0; %<---repulsive strength multiplied here
    r1 = -dR*r1;
    r2 = -dR*r2; 

    R1(p)  = R1(p)  + r1;
    R2(p)  = R2(p)  + r2;            
    RTq(p) = RTq(p) + r1.*(x2nrt-pc2(p)) - r2.*(x1nrt-pc1(p));
%{    
     hold off
     plot(x1, x2, 'k');
     hold on
     plot(pc1(p), pc2(p), 'o', pc1(q_list), pc2(q_list), '*')
     plot(x1(:,p), x2(:,p), 'r');
     plot(x1(:,q), x2(:,q), 'm'); 
     plot(x1nrt, x2nrt, 'r+');
     plot(y1nrt, y2nrt, 'm+');
     plot([pc1(p) pc1(q)], [pc2(p) pc2(q)])
     quiver( pc1(p), pc2(p), r1, r2, 'k') 
     pause
%}
    pp = [pp; x1nrt, x2nrt, y1nrt, y2nrt];

  end

end

%system computes force free
%[mean(R1)/mean(abs(R1)) mean(R2)/mean(abs(R2))]
%pause

%  hold off
%  plot(x1, x2, 'k.-');
%  hold on
%  quiver(pc1', pc2', R1, R2);
%  pause

% restore h value   
geom.length = h_store;

end % Repul

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R, dR] = Repul_profile(o, z)

% a repulsion profile with cut-off. Cut-off must occur at z = 1 (inside function)
% ie dist = l0 on outside call. 

%R  =  exp(-z); %
R  = (1 - sin(z*pi/2)).*( z < 1 );
%dR = -exp(-z); %
dR = -pi/2*cos(z*pi/2).*( z < 1 );
    
end % Repul_profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF ROUTINES THAT CALCULATE REPULSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bulkVelocity(o,geom,eta,Up,wp,force,torque,hFarField);
% form the velocity in the fluid bulk and along rigid bodies and plot
% them.

oc = curve;
X = geom.X;
[x,y] = oc.getXY(X);

% create mesh grid with target points
% xmin = min(min(x)) - 1;
% xmax = max(max(x)) + 1;
% ymin = min(min(y)) - 1;
% ymax = max(max(y)) + 1;
% 
% xmin = -3;
% xmax = 3;
% ymin = -3;
% ymax = 3;
% dx = 0.1;
% dy = 0.1;
% [xx,yy] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);

% load tracers as target points
X_tar = load("../examples/tracers.dat");
xx = X_tar(:,1);
yy = X_tar(:,2);

% put in capsules-like structure so that it can be passed into getZone
bd2.N = numel(xx);
bd2.nb = 1;
bd2.X = [xx(:); yy(:)];

% Construct near-singular integration structure so that velocity can be
% evaulated in bulk.
[~,NearOther] = geom.getZone(bd2,2);

jump = 0.5;
kernel = @o.exactStokesDL;
kernelSelf = @(z) +jump*z + o.exactStokesDLdiag(geom,geom.DLPStokes,z);

vel = o.nearSingInt(geom,eta,kernelSelf,NearOther,...
    kernel,kernel,bd2,false,false);

% Add on Rotlets and Stokeslets
for k = 1:geom.nb
  [cx,cy] = oc.getXY(geom.center(:,k));

  fx = force(2*(k-1) + 1);
  fy = force(2*(k-1) + 2);
  tor = torque(k);

  rx    = xx(:) - cx;
  ry    = yy(:) - cy;
  rho2  = rx.^2 + ry.^2;
  rdotf = rx*fx + ry*fy;
  
  vel = vel + (1/4/pi)*...
      [-0.5*log(rho2)*fx + rdotf./rho2.*rx;...
       -0.5*log(rho2)*fy + rdotf./rho2.*ry];

  vel = vel + (1/4/pi)*tor*[-ry./rho2;+rx./rho2];
  
end

% Add on Background flow 

vel = vel + hFarField(bd2.X);

[velx,vely] = oc.getXY(vel);
velx = reshape(velx,size(xx,1),size(xx,2));
vely = reshape(vely,size(xx,1),size(xx,2));

for k = 1:geom.nb

  [cx,cy] = oc.getXY(geom.center(:,k));
  rx = xx - cx;
  ry = yy - cy;

  bx = X(1:geom.N,k);     bx = [bx; bx(1)]; %and boundary curve
  by = X(geom.N+1:end,k); by = [by; by(1)];    
  wn = o.wn_PnPoly(xx, yy, bx, by, geom.N); %winding number to calculate interior of shapes 

  s = find(wn > 0.5);
  
  velx(s) = 1*(Up(1,k) - wp(k)*ry(s));
  vely(s) = 1*(Up(2,k) + wp(k)*rx(s));
  
end

VEL = [velx vely];
save("-ascii", "../examples/tracer_vel.dat", "VEL");

%clf;
%hold on
%quiver(xx,yy,velx,vely)
%surf(xx,yy,vely)
%plot(xx(ceil(end/20),:),vely(ceil(end/2),:),'b-o')
%plot(xx(ceil(end/20),:),velx(ceil(end/2),:),'r-o')
%axis equal
%pause

end % bulkVelocity



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wn = wn_PnPoly(o, p1, p2, v1, v2, N)

% wn_PnPoly(): winding number test for a point in a polygon
%      Input:   P = a point (p1, p2)
%               V = (v1, v2) = vertex points of a polygon with V(1) = V(N+1)
%      Return:  wn = the winding number (=0 only when P is outside)

    wn = 0*p1;

    for i = 1:N
    
        wn = wn + (v2(i) <= p2).*(v2(i+1)  > p2).*(o.isLeft(v1(i), v2(i), v1(i+1), v2(i+1), p1, p2) > 0);
        wn = wn - (v2(i)  > p2).*(v2(i+1) <= p2).*(o.isLeft(v1(i), v2(i), v1(i+1), v2(i+1), p1, p2) < 0);

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = isLeft(o, p1, p2, q1, q2, r1, r2)
    %    Input:  three points p, q, and r
    %    Return: >0 for r left of the line through p and q
    %            =0 for r  on the line
    %            <0 for r  right of the line     
    %http://geomalgorithms.com/a03-_inclusion.html

    out = (q1 - p1) .* (r2 - p2) - (r1 -  p1) .* (q2 - p2) ;

end
    
end % methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods(Static)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = lagrangeInterp(~)
% interpMap = lagrangeInterp builds the Lagrange interpolation
% matrix that takes seven function values equally distributed
% in [0,1] and returns the seven polynomial coefficients

LP(1,1) = 6.48e1;
LP(1,2) = -3.888e2;
LP(1,3) = 9.72e2;
LP(1,4) = -1.296e3;
LP(1,5) = 9.72e2;
LP(1,6) = -3.888e2;
LP(1,7) = 6.48e1;

LP(2,1) = -2.268e2;
LP(2,2) = 1.296e3;
LP(2,3) = -3.078e3;
LP(2,4) = 3.888e3;
LP(2,5) = -2.754e3;
LP(2,6) = 1.0368e3;
LP(2,7) = -1.62e2;

LP(3,1) = 3.15e2;
LP(3,2) = -1.674e3;
LP(3,3) = 3.699e3;
LP(3,4) = -4.356e3;
LP(3,5) = 2.889e3;
LP(3,6) = -1.026e3;
LP(3,7) = 1.53e2;

LP(4,1) = -2.205e2;
LP(4,2) = 1.044e3;
LP(4,3) = -2.0745e3;
LP(4,4) = 2.232e3;
LP(4,5) = -1.3815e3;
LP(4,6) = 4.68e2;
LP(4,7) = -6.75e1;

LP(5,1) = 8.12e1;
LP(5,2) = -3.132e2;
LP(5,3) = 5.265e2;
LP(5,4) = -5.08e2;
LP(5,5) = 2.97e2;
LP(5,6) = -9.72e1;
LP(5,7) = 1.37e1;

LP(6,1) = -1.47e1;
LP(6,2) = 3.6e1;
LP(6,3) = -4.5e1;
LP(6,4) = 4.0e1;
LP(6,5) = -2.25e1;
LP(6,6) = 7.2e0;
LP(6,7) = -1e0;

LP(7,1) = 1e0;
% rest of the coefficients are zero

end % lagrangeInterp

end % methods(Static)

end % classdef
