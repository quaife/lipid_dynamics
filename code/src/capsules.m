classdef capsules < handle
% This class implements standard calculations that need to be done to a
% componenet of the solid wall, or a collection of arbitrary target
% points (such as tracers).  The main tasks that can be performed
% include constructing structures required for near-singluar integration

properties
N;            % number of points per componenet
nb;           % number of rigid bodies
X;            % positions of component
center;       % center of each rigid body
tau;          % orientation of each rigid body
radii;        % radius of each rigid body
rho;          % screen length of particles
gam;          % HAP strength
ar;           % aspect ratio of particles
xt;           % tangent unit vector
sa;           % Jacobian
cur;          % curvature
length;       % total length of each component
RepulLength   % repulsion length
RepulStrength % repulusion strength
bcShift       % boundary condition shifting parameter
bcType        % type of the boundary condition
shape         % determine the type of particle shape
petal         % number of petal if particle shape is star

% structure for near-singular integration from body to body
nearStructB2B;   
% structure for near-singular integration from body to target
nearStructB2T;   
DLPStokes;    % double-layer potential matrix for Stokes
% rank-one correction for Stokes to remove one dimensional null space
N0Stokes;     
DLPYukawa;    % double-layer potential matrix for Yukawa
upRate;       % upsampling rate that is used for NSI quadrature
NPeaks;       % change the number of peaks in cosine bc
% besselk(1,~) kernel for each target point with upsampled source points
% on all other bodies
BK1;          
% Stokes double-layer kernel for each target point with upsampled source
% points on all other bodies
StokesDL;

end %properties


methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = capsules(prams,xc,tau)
% capsules(X) constructs an object of class capsules.  Mostly, it
% takes the values of prams and options that it requires.
% This is the constructor

o.N = prams.N;
o.nb = prams.nb;
o.radii = prams.radii;
o.rho = prams.rho;
o.gam = prams.gam;
o.ar = prams.ar;
o.center = xc; % center
o.tau = tau; % inclination angle
o.bcShift = prams.bcShift;
o.bcType = prams.bcType;
o.upRate = ceil(sqrt(o.N));
o.NPeaks = prams.NPeaks;
o.shape  = prams.shape;
o.petal  = prams.petal;

theta = (0:o.N-1)'*2*pi/o.N;
X = zeros(o.N*2,o.nb);  

switch o.shape
  case 'circle'
    for i = 1:o.nb
      % shape of particle with rotation
      refX = kron([cos(tau(i)) -sin(tau(i)); +sin(tau(i)) cos(tau(i))],...
          eye(o.N))*[o.ar(i)*cos(theta);sin(theta)]*o.radii(i);  
      % shift to the correct center
      X(1:o.N,i) = refX(1:o.N) + xc(1,i);
      X(o.N+1:2*o.N,i) = refX(o.N+1:end) + xc(2,i);
    end
  case 'star'
    for i = 1:o.nb
      % shape of star
      rr = o.radii(i) + o.ar(i)*cos(o.petal*theta);
      % rotated star
      refX = kron([cos(tau(i)) -sin(tau(i));...
             +sin(tau(i)) cos(tau(i))],...
          eye(o.N))*[rr.*cos(theta);rr.*sin(theta)];
      % shift to the correct center
      X(1:o.N,i) = refX(1:o.N) + xc(1,i);
      X(o.N+1:2*o.N,i) = refX(o.N+1:end) + xc(2,i);
    end

  case 'tube'
    % parameters for the boundary:
    % a and b control the length and height and order controls the
    % regularity
    a = 10; b = 3; order = 8;
    t = (0:o.N-1)'*2*pi/o.N;
    % radius of the "ellipse"
    r = (cos(-t).^order + sin(-t).^order).^(-1/order);
    % rounded off cylinder. 
    x = a*r.*cos(-t); y = b*r.*sin(-t);

    X = [x;y];

  case 'choke'
    % parameters for the boundary:
    % a and b control the length and height and order controls the
    % regularity
    a = 20; b = 3; c = 0.6; order = 8;
    t = (0:o.N-1)'*2*pi/o.N;
    % radius of the "ellipse"
    r = (cos(-t).^order + sin(-t).^order).^(-1/order);
    x = a*r.*cos(-t); y = b*r.*sin(-t);
    ind = abs(x) < 10;
    scalRat = 2*c/(1+c)*(0.5-0.5*cos(pi*x(ind(1:end/2))/10)).^10 + ...
      (1-c)/(1+c);
    y(ind) = y(ind).*[scalRat;scalRat];
%    y(ind) = y(ind).*(1 - c*cos(x(ind)))/(1+c);
    % rounded off cylinder and pinched off in the midle
    X = [x;y];
end
o.X = X;

% repulsion length and strength for interactions between Janus
% particles
o.RepulLength = prams.RepulLength;
o.RepulStrength = prams.RepulStrength;

oc = curve;
% compute arlength, tangent, and curvature
[o.sa,o.xt,o.cur] = oc.diffProp(o.X); 
% compute total length
[~,o.length] = oc.geomProp(o.X);


% store as empty arrays. These are only constructed when needed
o.BK1 = [];
o.StokesDL = [];

end % capsules: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function YukawaKernelMatrix(geom)
% YukawaKernelMatrix(o,geom) computes the bessel function between all
% pairs of points on different bodies. This is used to accelerate the
% number of times besselk has to be evaulated when doing matrix-vector
% multiplication in GMRES for Yukawa solve

oc = curve;
N = geom.N;
Nup = N*geom.upRate;
nb = geom.nb;
geom.BK1 = zeros(N,Nup*(nb-1),nb);
for k = 1:nb
  Ksou = [(1:k-1) (k+1:nb)];
  [xtar,ytar] = oc.getXY(geom.X(:,k));
  xtar = xtar(:,ones(Nup*(nb-1),1));
  ytar = ytar(:,ones(Nup*(nb-1),1));

  [xsou,ysou] = oc.getXY(geom.X(:,Ksou));
  xsou = interpft(xsou,Nup);
  ysou = interpft(ysou,Nup);
  [~,tangent] = oc.diffProp([xsou;ysou]);
  [tx,ty] = oc.getXY(tangent);
  xsou = xsou(:); ysou = ysou(:);
  xsou = xsou(:,ones(N,1))';
  ysou = ysou(:,ones(N,1))';

  nx = -ty; nx = nx(:);
  ny = tx; ny = ny(:);
  nx = nx(:,ones(N,1))';
  ny = ny(:,ones(N,1))';

  rdotn = (xtar - xsou).*nx + (ytar - ysou).*ny;
  dis = sqrt((xtar - xsou).^2 + (ytar - ysou).^2);

  geom.BK1(:,:,k) = -1/2/pi/geom.rho*besselk(1,dis/geom.rho).*rdotn./dis;
end

end % YukawaKernelMatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StokesDLMatrix(geom)
% StokesDLMatrix(o,geom) computes the Stokes double-layer potential
% between all pairs of points on different bodies. This is used to
% accelerate the number of times exactStokesDL has to be evaulated when
% doing matrix-vector multiplication in GMRES for Stokes solve

oc = curve;
N = geom.N;
Nup = N*geom.upRate;
nb = geom.nb;
geom.StokesDL = zeros(2*N,2*Nup*(nb-1),nb);
for k = 1:nb
  % skip self-contribution
  Ksou = [(1:k-1) (k+1:nb)];
  % compute target points
  [xtar,ytar] = oc.getXY(geom.X(:,k));
  xtar = xtar(:,ones(Nup*(nb-1),1));
  ytar = ytar(:,ones(Nup*(nb-1),1));

  % compute source points
  [xsou,ysou] = oc.getXY(geom.X(:,Ksou));
  xsou = interpft(xsou,Nup);
  ysou = interpft(ysou,Nup);
  [~,tangent] = oc.diffProp([xsou;ysou]);
  [tx,ty] = oc.getXY(tangent);
  xsou = xsou(:); ysou = ysou(:);
  xsou = xsou(:,ones(N,1))';
  ysou = ysou(:,ones(N,1))';

  % compute upsampled normal
  nx = -ty; nx = nx(:);
  ny = tx; ny = ny(:);
  nx = nx(:,ones(N,1))';
  ny = ny(:,ones(N,1))';

  rdotn = (xtar - xsou).*nx + (ytar - ysou).*ny;
  dis2 = (xtar - xsou).^2 + (ytar - ysou).^2;

  geom.StokesDL(1:N,1:Nup*(nb-1),k) = ...
        rdotn./dis2.*(xtar - xsou).^2./dis2;
  geom.StokesDL(1:N,Nup*(nb-1)+1:2*Nup*(nb-1),k) = ...
        rdotn./dis2.*(xtar - xsou).*(ytar - ysou)./dis2;
  geom.StokesDL(N+1:2*N,1:Nup*(nb-1),k) = ...
        geom.StokesDL(1:N,Nup*(nb-1)+1:2*Nup*(nb-1),k);
  geom.StokesDL(N+1:2*N,Nup*(nb-1)+1:2*Nup*(nb-1),k) = ...
        rdotn./dis2.*(ytar - ysou).^2./dis2;
end

geom.StokesDL = -geom.StokesDL/pi;

end % StokesDLMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NearSelf,NearOther] = getZone(o,bd2,relate)
% [NearSelf,NearOther] = getZone(bd1,bd2,relate) constructs
% each boundary, index of the closest point, nearest point on a local
% interapolant, and argument of that nearest point. bd1 contains
% the source points (which are also target points) and bd2
% contains additional target points.  The
% values of relate corresond to
% relate == 1  => only require NearSelf  (ie. bd1 to bd1)
% relate == 2  => only require NearOther (ie. bd1 to bd2)
% relate == 3  => require both NearSelf and NearOther

NearSelf = [];
NearOther = [];

N1 = o.N;  % number of source/target points per boundary component
nb = o.nb; % number of source/target boundaries
X1 = o.X;  % source and target points
oc = curve;
% separate targets into x and y coordinates
[xsou,ysou] = oc.getXY(X1);

% smallest arclength over all boundaries
h = max(o.length)/N1;

% Estimate for number of points per box. This simply sets the number of
% uniformly refined boxes we take. Estimate is not very accurate. What
% ptsperbox represents is the total number of points that could be put
% in each two-dimensional bin where no two are less than distance h from
% one another. However, our points live on curves and thus will not fill
% up an entire bin
ptsperbox = 10;

% Add a buffer around the points so that it is easier to
% work with bd2
H = sqrt(ptsperbox)*h;
xmin = min(min(xsou));
xmax = max(max(xsou));
xmin = xmin - H;
xmax = xmax + H;
ymin = min(min(ysou));
ymax = max(max(ysou));
ymin = ymin - H;
ymax = ymax + H;

% Find bounds for box that contains all points and add a buffer so that
% all points are guaranteed to be in the box
Nx = ceil((xmax - xmin)/H);
Ny = ceil((ymax - ymin)/H);

Nbins = Nx * Ny; % Total number of bins

% Index in x and y direction of the box containing each point
ii = ceil((xsou - xmin)/H);
jj = ceil((ysou - ymin)/H);
% Find bin of each point using lexiographic ordering (x then y)
bin = (jj-1)*Nx + ii;

% allocate space for storing first and last points
fpt = zeros(Nbins,nb);
lpt = zeros(Nbins,nb);
% build permute. Need binsort to find first and last points in each box
[binsort,permute] = sort(bin);

% Construct first and last point in each box corresponding to each
% boundary. The order is based on permute. For example,
% permute(fpt(ibox,k)),...,permute(lpt(ibox,k)) is the set of all points
% from boundary k contained in box ibox
for k = 1:nb % Loop over boundaries
  for j = 1:N1 % Loop over bins
    ibox = binsort(j,k);
    if (fpt(ibox,k) == 0)
      fpt(ibox,k) = j;
      lpt(ibox,k) = 1;
    else
      lpt(ibox,k) = lpt(ibox,k) + 1;
    end
  end
  lpt(:,k) = fpt(:,k) + lpt(:,k) - 1;
end

neigh = zeros(Nbins,9);

%Do corners first
% bottom left corner
neigh(1,1:4) = [1 2 Nx+1 Nx+2];
% bottom right corner
neigh(Nx,1:4) = [Nx Nx-1 2*Nx 2*Nx-1];
% top left corner
neigh(Nbins-Nx+1,1:4) = [Nbins-Nx+1 Nbins-Nx+2 ...
Nbins-2*Nx+1 Nbins-2*Nx+2];
% top right corner
neigh(Nbins,1:4) = [Nbins Nbins-1 Nbins-Nx Nbins-Nx-1];

% neighbors of bottom row
for j = 2:Nx-1
  neigh(j,1:6) = j + [-1 0 1 Nx-1 Nx Nx+1];
end

% neighbors of top row
for j = Nbins-Nx+2:Nbins-1
  neigh(j,1:6) = j + [-1 0 1 -Nx-1 -Nx -Nx+1];
end

% neighbors of left column
for j=Nx+1:Nx:Nbins-2*Nx+1
  neigh(j,1:6) = j + [-Nx -Nx+1 0 1 Nx Nx+1];
end

% neighbors of right column
for j=2*Nx:Nx:Nbins-Nx
  neigh(j,1:6) = j + [-Nx-1 -Nx -1 0 Nx-1 Nx];
end

% J is the index of boxes that are not on the boundary
J = (Nx + 1:Nbins - Nx);
J = J(mod(J-1,Nx)~=0);
J = J(mod(J,Nx)~=0);
% neighbors of interior points
for j=J
  neigh(j,:) = j + [-Nx-1 -Nx -Nx+1 -1 0 1 Nx-1 Nx Nx+1];
end
% TREE STRUCTURE IS COMPLETE


if (relate == 1 || relate == 3)
  for k = 1:nb
    % dist(n,k,j) is the distance of point n on boundary k to boundary j
    distSS{k} = spalloc(N1,nb,0);
    % near or far zone
    zoneSS{k} = spalloc(N1,nb,0);
    % nearest point using local interpolant
    nearestSS{k} = spalloc(2*N1,nb,0);
    % index of closest discretization point
    icpSS{k} = spalloc(N1,nb,0);
    % argument in [0,1] of local interpolant
    argnearSS{k} = spalloc(N1,nb,0);
    nearFibersSS{k} = {};
  end
  % Represent near-singular integration structure in such a way that we
  % can use sparse matricies.


  % begin classifying points where we are considering
  % boundary to boundary relationships
  for k = 1:nb
    % Find all boxes containing points of boundary k
    boxes = unique(bin(:,k));
    % Look at all neighbors of boxes containing boundary k
    boxes = neigh(boxes,:);
    % Remove repetition
    boxes = unique(boxes(:));
    % Delete non-existent boxes that came up because of neigh
    boxes = boxes(boxes~=0);

    K = [(1:k-1) (k+1:nb)];
    for k2 = K
      % Find index of all points in neighboring boxes of boundary k that
      % are in boundary k2
      istart = fpt(boxes,k2);
      iend = lpt(boxes,k2);
      istart = istart(istart ~= 0);
      iend = iend(iend ~= -1);

      % Allocate space to assign possible near points
      neighpts = zeros(sum(iend-istart+1),1);
      is = 1;
      % neighpts contains all points on boundary k2 that are in
      % neighboring boxes to boundary k
      for j=1:numel(istart)
        ie = is + iend(j) - istart(j);
        neighpts(is:ie) = permute(istart(j):iend(j),k2);
        is = ie + 1;
      end

      % sorting should help speedup as we won't be jumping around
      % through different boxes
      neighpts = sort(neighpts);

      n = 0;
      for i=1:numel(neighpts)
        ipt = neighpts(i);
        % box containing ipt on boundary k2
        ibox = bin(ipt,k2);
        if (ibox ~= n)
          % Check if we need to move to a new box
          n = ibox;
          % neighbors of this box
          neighbors = neigh(ibox,:);
          % Remove non-existent neighbors
          neighbors = neighbors(neighbors~=0);
          % Find points on boundary k in neighboring boxes
          istart = fpt(neighbors,k);
          iend = lpt(neighbors,k);
          istart = istart(istart ~= 0);
          iend = iend(iend ~= -1);
          neighpts2 = zeros(sum(iend-istart+1),1);
          is = 1;
          % neighpts2 contains all points on boundary k that
          % are in neighboring box of ibox
          for j=1:numel(istart)
            ie = is + iend(j) - istart(j);
            neighpts2(is:ie) = permute(istart(j):iend(j),k);
            is = ie + 1;
          end
        end % decide if we need to switch boxes

        % Find minimum distance between ipt on boundary k2 to possible
        % closest points on boundary k
        [d0,d0loc] = min((xsou(ipt,k2) - xsou(:,k)).^2 + ...
                         (ysou(ipt,k2) - ysou(:,k)).^2);
        % Save on not taking the square root on a vector but instead on
        % a single real number
        d0 = sqrt(d0);

        icpSS{k}(ipt,k2) = d0loc;
        if (d0 < 2*h);
          [distSS{k}(ipt,k2),nearestx,nearesty,...
                argnearSS{k}(ipt,k2)] = ...
            o.closestPnt([xsou;ysou],xsou(ipt,k2),ysou(ipt,k2),...
                k,icpSS{k}(ipt,k2));
            nearestSS{k}(ipt,k2) = nearestx;
            nearestSS{k}(ipt+N1,k2) = nearesty;
            
          % Find closest point along a local interpolant using Newton's
          % method.
          nearFibersSS{k}= [nearFibersSS{k},k2];
            
          % Point ipt of boundary k2 is in the near zone of
          % boundary k
          if (distSS{k}(ipt,k2) < h)
            zoneSS{k}(ipt,k2) = 1;
          end
        end
      end % ipt
    end % k2
    nftmp = nearFibersSS{k};
    nearFibersSS{k} = unique([nftmp{:}]);
  end % k

  % Store everything in the structure NearSelf. This way it is much
  % cleaner to pass everything around
  NearSelf.dist = distSS;
  NearSelf.zone = zoneSS;
  NearSelf.nearest = nearestSS;
  NearSelf.icp = icpSS;
  NearSelf.argnear = argnearSS;
  NearSelf.nearFibers = nearFibersSS;

end % relate == 1 || relate == 3


% Bin target points with respect to the source points
if (relate == 2 || relate == 3)
  Np2 = bd2.N; % number of additional targets
  nb2 = bd2.nb; % number of additional boundaries
  X2 = bd2.X; % additional target points
  [xtar,ytar] = oc.getXY(X2);

  % Represent near-singular integration structure so that we can use
  % sparse matricies.
  for k = 1:nb
    % dist(n,k,j) is the distance of point n on boundary k to
    distST{k} = spalloc(N1,nb2,0);
    % near or far zone
    zoneST{k} = spalloc(N1,nb2,0);
    % nearest point using local interpolant
    nearestST{k} = spalloc(2*N1,nb2,0);
    % index of closest discretization point
    icpST{k} = spalloc(N1,nb2,0);
    % argument in [0,1] of local interpolant
    argnearST{k} = spalloc(N1,nb2,0);
  end

  % Only have to consider xx(ind),yy(ind) since all other points
  % are not contained in the box [xmin xmax] x [ymin ymax]
  itar = ceil((xtar - xmin)/H);
  jtar = ceil((ytar - ymin)/H);
  [indx,indy] = find((itar >= 1) & (itar <= Nx) & ...
  (jtar >= 1) & (jtar <= Ny));

  for k = 1:nb % loop over sources
    for nind = 1:numel(indx)
      % loop over points that are not outside the box that surrounds
      % all target points with a sufficiently large buffer
      ii = indx(nind);
      jj = indy(nind);
      binTar = (jtar(ii,jj) - 1)*Nx + itar(ii,jj);
      boxesTar = neigh(binTar,:);
      boxesTar = boxesTar(boxesTar~=0);
      istart = fpt(boxesTar,k);
      iend = lpt(boxesTar,k);
      istart = istart(istart ~= 0);
      iend = iend(iend ~= -1);

      % Allocate space to assign possible near points
      neighpts = zeros(sum(iend-istart+1),1);
      if numel(neighpts) > 0
        % it is possible of the neighboring boxes to contain no points.
        is = 1;
        % Find set of potentially nearest points to (xtar(jj),ytar(jj))
        for j = 1:numel(istart)
          ie = is + iend(j) - istart(j);
          neighpts(is:ie) = permute(istart(j):iend(j),k);
          is = ie + 1;
        end

        % find closest point and distance between (xtar(jj),ytar(jj))
        % and boundary k. Only need to look at points in neighboring
        % boxes
        [d0,d0loc] = min((xtar(ii,jj) - xsou(neighpts,k)).^2 + ...
                         (ytar(ii,jj) - ysou(neighpts,k)).^2);
        d0 = sqrt(d0);
        icpST{k}(ii,jj) = neighpts(d0loc);

        if d0 < 2*h
          [distST{k}(ii,jj),nearestx,nearesty,argnearST{k}(ii,jj)] = ...
            o.closestPnt([xsou;ysou],xtar(ii,jj),ytar(ii,jj),...
            k,icpST{k}(ii,jj));
          nearestST{k}(ii,jj) = nearestx;
          nearestST{k}(ii+Np2,jj) = nearesty;
          if distST{k}(ii,jj) < h
            % (xtar(ii,jj),ytar(ii,jj)) is in the near zone of boundary
            % k
            zoneST{k}(ii,jj) = 1;
          end
        end % d0 < 2*h
      end % numel(neighpts) > 0

    end % ii and jj
  end % k

  % store near-singluar integration requirements in structure NearOther
  NearOther.dist = distST;
  NearOther.zone = zoneST;
  NearOther.nearest = nearestST;
  NearOther.icp = icpST;
  NearOther.argnear = argnearST;
  NearOther.nearFibers = [];

end % relate == 2 || relate == 3

end % getZone


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E = computeEnergy(geom,eta)

oc = curve;
[x,y] = oc.getXY(geom.X);
[tx,ty] = oc.getXY(geom.xt);
nx = ty;
ny = -tx;

ds = 1e-2;
xtar = x + ds*nx;
ytar = y + ds*ny;
%clf; hold on;
%plot(x,y,'r')
%plot(xtar,ytar,'b.')
%pause

targets.N = geom.N;
targets.nb = geom.nb;
targets.X = [xtar;ytar];

[~,NearJanusTargets] = geom.getZone(targets,2);

op = poten(geom.N,geom.rho);
kernel = @op.exactYukawaDL;
kernelDirect = @op.exactYukawaDL;
kernelSelf = @(z) +0.5*z + op.exactYukawaDLdiag(geom,...
      geom.DLPYukawa,z);

potJanus = op.nearSingInt(geom,eta,kernelSelf,...
    NearJanusTargets,kernel,kernelDirect,targets,false,false);
% get rid of second half which is only used when solving a vector PDE
% (Stokes)
potJanus = potJanus(1:end/2,:);

bdCond = reshape(geom.yukawaRHS,geom.N,geom.nb);
integrand = bdCond.*(bdCond - potJanus)/ds;
E = sum(integrand(:).*geom.sa(:))*2*pi/geom.N;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,nearestx,nearesty,theta] = closestPnt(~,X,...
        xTar,ytar,k,icp)
% [dist,nearestx,nearesty,theta] = closestPnt(X,xTar,ytar,k,icp)
% computes the closest point on boundary k to (xTar,ytar)
% using a Lagrange interpolant.  icp is the index of the closest
% point on the discrete mesh which is used as an initial guess

N = size(X,1)/2; % Number of points per boundary 
A = poten.lagrangeInterp;
% need interpolation matrix and its size
interpOrder = size(A,1);

% Accommodate for either an even or odd number of interpolation
% points
p = ceil((interpOrder+1)/2);
% band of points around icp.  The -1,+1 combination sets index
% 0 to N as required by the code
pn = mod((icp-p+1:icp-p+interpOrder)' - 1,N) + 1;

px = A*X(pn,k); % polynomial interpolant of x-coordinate
py = A*X(pn+N,k); % polynomial interpolant of y-coordinate
% To do Newton's method, need two derivatives
Dpx = px(1:end-1).*(interpOrder-1:-1:1)';
Dpy = py(1:end-1).*(interpOrder-1:-1:1)';
D2px = Dpx(1:end-1).*(interpOrder-2:-1:1)';
D2py = Dpy(1:end-1).*(interpOrder-2:-1:1)';

% midpoint is a good initial guess
theta = 1/2;
for newton = 1:2
  % Using filter is the same as polyval, but it is much
  % faster when only requiring a single polyval such as here.
  zx = filter(1,[1 -theta],px);
  zx = zx(end);
  zy = filter(1,[1 -theta],py);
  zy = zy(end);
  Dzx = filter(1,[1 -theta],Dpx);
  Dzx = Dzx(end);
  Dzy = filter(1,[1 -theta],Dpy);
  Dzy = Dzy(end);
  D2zx = filter(1,[1 -theta],D2px);
  D2zx = D2zx(end);
  D2zy = filter(1,[1 -theta],D2py);
  D2zy = D2zy(end);

  % numerator of Newton's method
  newtonNum = (zx-xTar)*Dzx + (zy-ytar)*Dzy;
  % denominator of Newton's method
  newtonDen = (zx-xTar)*D2zx + (zy-ytar)*D2zy + ...
      Dzx^2 + Dzy^2;
  % one step of Newton's method
  theta = theta - newtonNum/newtonDen;
end
% Do a few (no more than 3) Newton iterations

% Compute nearest point and its distance from the target point
nearestx = filter(1,[1,-theta],px);
nearestx = nearestx(end);
nearesty = filter(1,[1,-theta],py);
nearesty = nearesty(end);
dist = sqrt((nearestx - xTar)^2 + (nearesty - ytar)^2);

end % closestPnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,x1,x2,y1,y2] = closestPntPair(~,X,...
        k,l,jk,jl)
% [dist,x1,x2,y1,y2] = closestPntPair(X,k,l,icp,jcp) computes the
% closest pair of points (x1, x2) and (y1, y2) between boundaries k and
% l using a Lagrange interpolant. jk, jl are the indices of points on
% the discrete mesh used as an initial guess

N = size(X,1)/2; % Number of points per boundary 
A = poten.lagrangeInterp;
% need interpolation matrix and its size
interpOrder = size(A,1);

% Accommodate for either an even or odd number of interpolation
% points
p = ceil((interpOrder+1)/2);
% band of points around icp, jcp respectively.  The -1,+1 combination sets index
% 0 to N as required by the code
pn = mod((jk-p+1:jk-p+interpOrder)' - 1,N) + 1;
qn = mod((jl-p+1:jl-p+interpOrder)' - 1,N) + 1;

p1 = A*X(pn,  k); % polynomial interpolant of x1-coordinate
p2 = A*X(pn+N,k); % polynomial interpolant of x2-coordinate

q1 = A*X(qn,  l); % polynomial interpolant of x1-coordinate
q2 = A*X(qn+N,l); % polynomial interpolant of x2-coordinate

Dp1 = p1(1:end-1).*(interpOrder-1:-1:1)';
Dp2 = p2(1:end-1).*(interpOrder-1:-1:1)';

Dq1 = q1(1:end-1).*(interpOrder-1:-1:1)';
Dq2 = q2(1:end-1).*(interpOrder-1:-1:1)';

% To do Newton's method, need two derivatives
DDp1 = Dp1(1:end-1).*(interpOrder-2:-1:1)';
DDp2 = Dp2(1:end-1).*(interpOrder-2:-1:1)';

DDq1 = Dq1(1:end-1).*(interpOrder-2:-1:1)';
DDq2 = Dq2(1:end-1).*(interpOrder-2:-1:1)';

s = 1/2;
t = 1/2;
% midpoint is a good initial guess
for newton = 1:2
  u1 = filter(1,[1,-s],p1); u1 = u1(end);  
  u2 = filter(1,[1,-s],p2); u2 = u2(end);  
    
  v1 = filter(1,[1,-t],q1); v1 = v1(end);  
  v2 = filter(1,[1,-t],q2); v2 = v2(end);  

  Du1 = filter(1,[1,-s],Dp1); Du1 = Du1(end);  
  Du2 = filter(1,[1,-s],Dp2); Du2 = Du2(end);  
    
  Dv1 = filter(1,[1,-t],Dq1); Dv1 = Dv1(end);  
  Dv2 = filter(1,[1,-t],Dq2); Dv2 = Dv2(end);  

  DDu1 = filter(1,[1,-s],DDp1); DDu1 = DDu1(end);  
  DDu2 = filter(1,[1,-s],DDp2); DDu2 = DDu2(end);  
    
  DDv1 = filter(1,[1,-t],DDq1); DDv1 = DDv1(end);  
  DDv2 = filter(1,[1,-t],DDq2); DDv2 = DDv2(end);  
  % Using filter is the same as polyval, but it is much
  % faster when only requiring a single polyval such as here.

  f   = 0.5*(u1 - v1)^2 + 0.5*(u2 - v2)^2;
  f_s = (u1 - v1)*Du1 + (u2 - v2)*Du2;
  f_t = (v1 - u1)*Dv1 + (v2 - u2)*Dv2;
  
  f_ss = Du1*Du1 + (u1 - v1)*DDu1 + Du2*Du2 + (u2 - v2)*DDu2;
  f_st = -Dv1*Du1 - Dv2*Du2;
  f_ts = -Du1*Dv1 - Du2*Dv2;
  f_tt = Dv1*Dv1 + (v1 - u1)*DDv1 + Dv2*Dv2 + (v2 - u2)*DDv2;

  w = -([f_ss f_st; f_ts f_tt])\[f_s; f_t];
  
  s = s + w(1);
  t = t + w(2);
  % one step of Newton's method
end
% Do a few (no more than 3) Newton iterations

x1 = filter(1,[1,-s],p1);
x1 = x1(end);
x2 = filter(1,[1,-s],p2);
x2 = x2(end);

y1 = filter(1,[1,-t],q1);
y1 = y1(end);
y2 = filter(1,[1,-t],q2);
y2 = y2(end);

dist = sqrt((x1 - y1)^2 + (x2 - y2)^2);
% Compute nearest point and its distance from the target point

end % closestPntPair

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rhs = yukawaRHS(geom)
% Build right-hand side for the Yukawa equation solver

rhs = zeros(geom.N,geom.nb);
oc = curve;
[xc,yc] = oc.getXY(geom.center);
[x,  y] = oc.getXY(geom.X);
tau     = geom.tau;
rad     = geom.radii;
NPeaks = geom.NPeaks;

theta   = atan2(y - yc, x - xc);

%reshape nb-many boundary parameters into N by nb array
bcs = meshgrid(geom.bcShift, ones(geom.N,1)); 

%rhs = geom.yukawaExact(x,y);
switch geom.bcType
  case 'cosine'  
    rhs  = 0.5*(1 + cos(NPeaks*(theta - tau))) + bcs;

    % The boundary condition 'cosine' is normalize by a factor.
    for i = 1:geom.nb
        rhs(:,i) = rhs(:,i)/sqrt(sum(rhs(:,i).^2.*geom.sa(:,i))*2*pi/geom.N);
    end    
    
  case 'vonMises'
  rhs  = exp(bcs.*cos(theta - tau))/2/pi./besseli(0,bcs);
end

rhs = rhs(:);

end % yukawaRHS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, dz1, dz2] = yukawaExact(geom, Xtar, Ytar)
% The exact solution used at various places 

h = 1e-5;

z = geom.primal(Xtar, Ytar);

dz1 = (geom.primal(Xtar+h, Ytar) - geom.primal(Xtar-h, Ytar))/(2*h);
dz2 = (geom.primal(Xtar, Ytar+h) - geom.primal(Xtar, Ytar-h))/(2*h);
  
end % yukawaExact 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = primal(geom,xx,yy)

oc = curve;
[xc,yc] = oc.getXY(geom.center);

z = 0*xx;
for i = 1:geom.nb
  th  = atan2(yy - yc(i), xx - xc(i));
  rad = sqrt((xx-xc(i)).^2 + (yy-yc(i)).^2);
  j = i;
  r_ref = min(geom.radii(i)*geom.ar(i), geom.radii(i));
  z = z + besselk(j,rad/geom.rho).*cos(j*th)/besselk(j,r_ref/geom.rho);
end

end % primal
  
end % methods

end %capsules

