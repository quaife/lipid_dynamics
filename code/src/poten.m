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
prams.ar = geomSou.ar;
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
% START OF ROUTINES THAT BUILD LAYER-POTENTIAL MATRICIES
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
  
  diffx = xtar-xsou; diffy = ytar-ysou;
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
xx = xsou; yy = ysou;
xsou = xsou(:,ones(Ntar,1))';
ysou = ysou(:,ones(Ntar,1))';

den = den(:,K1);
den = den(:);
den = den(:,ones(Ntar,1))';

[tx,ty] = oc.getXY(geom.xt(:,K1));
nx = -ty(:); ny = tx(:);
normalx = nx(:,ones(Ntar,1))';
normaly = ny(:,ones(Ntar,1))';

%den = f.*geom.sa*2*pi/geom.N;
%[xtar,ytar] = oc.getXY(Xtar);
%rx = xtar - xsou'; ry = ytar - ysou';
%dis2 = rx.^2 + ry.^2; dis = sqrt(dis2);
%rdotn = rx.*nx + ry.*ny;
%kernel = -1/2/pi*besselk(1,dis).*rdotn./dis.*den;
%yukawaDLPtar = sum(kernel);

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

  kernel = -1/2/pi/geom.rho.*besselk(1,dis/geom.rho).*...
        rdotn./dis.*den;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % methods

methods(Static)

    
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
