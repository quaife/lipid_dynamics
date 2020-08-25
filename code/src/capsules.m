classdef capsules < handle
% This class implements standard calculations that need to be done to a
% componenet of the solid wall, or a collection of arbitrary target
% points (such as tracers).  The main tasks that can be performed
% include constructing structures required for near-singluar integration

properties
N;          % number of points per componenet
nb;         % number of rigid bodies
X;          % positions of component
center;     % center of each rigid body
tau;        % orientation of each rigid body
radii;      % radius of each rigid body
rho;        % screen length of particles
ar;         % aspect ratio of particles
xt;         % tangent unit vector
sa;         % Jacobian
cur;        % curvature
length;     % total length of each component
nearStruct; % structure for near-singular integration
DLPStokes;  % double-layer potential matrix
DLPYukawa;  % double-layer potential matrix for yukawa
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
o.ar = prams.ar;
o.center = xc; % center
o.tau = tau; % inclination angle

oc = curve;

theta = (0:o.N-1)'*2*pi/o.N;
X = zeros(o.N*2,o.nb);  

for i = 1:o.nb
  refX = kron([cos(tau(i)) -sin(tau(i)); +sin(tau(i)) cos(tau(i))],eye(o.N)) * ... 
      [o.ar(i)*cos(theta);sin(theta)]*o.radii(i);  
  % shape of particle

  % rotated circle
  X(1:o.N,i) = refX(1:o.N) + xc(1,i);
  X(o.N+1:2*o.N,i) = refX(o.N+1:end) + xc(2,i);
  % shift to the correct center
end

o.X = X;

[o.sa,o.xt,o.cur] = oc.diffProp(o.X); 
% compute arlenght, tangent, and curvature
[~,o.length] = oc.geomProp(o.X);
% compute total length

end % capsules: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [force,torque] = bodyForceTorque(o)

% will direct the calculated force, torque here
    
force = zeros(2*o.nb,1);  % original [1;2]
torque = zeros(o.nb,1);   % original -10


end % bodyForce Torque



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NearSelf,NearOther] = getZone(o,bd2,relate)
% [NearSelf,NearOther] = getZone(bd1,bd2,relate) constructs
% each boundary, index of the closest point, nearest point on a local
% interapolant, and argument of that nearest point.  bd1 contains
% the source points (which are also target points) and bd2
% contains additional target points.  The
% values of relate corresond to
% relate == 1  => only require NearSelf  (ie. bd1 to bd1)
% relate == 2  => only requirpe NearOther (ie. bd1 to bd2)
% relate == 3  => require both NearSelf and NearOther

NearSelf = [];
NearOther = [];

N1 = o.N; % number of source/target points
nb = o.nb; % number of source/target boundaries
X1 = o.X; % source and target points
oc = curve;
[xsou,ysou] = oc.getXY(X1);
% separate targets into x and y coordinates

h = max(o.length)/N1;
% smallest arclength over all boundaries
ptsperbox = 10;
% Estimate for number of points per box.  This simply sets the
% number of uniformly refined boxes we take.  Estimate is not very
% accurate.  What ptsperbox represents is the total number of points
% that could be put in each two-dimensional bin where no two are
% less than distance h from one another.  However, our points live
% on curves and thus will not fill up an entire bin

H = sqrt(ptsperbox)*h;
xmin = min(min(xsou));
xmax = max(max(xsou));
xmin = xmin - H;
xmax = xmax + H;
ymin = min(min(ysou));
ymax = max(max(ysou));
ymin = ymin - H;
ymax = ymax + H;
% Add a buffer around the points so that it is easier to
% work with bd2

Nx = ceil((xmax - xmin)/H);
Ny = ceil((ymax - ymin)/H);
% Find bounds for box that contains all points and add a buffer
% so that all points are guaranteed to be in the box

Nbins = Nx * Ny; % Total number of bins

ii = ceil((xsou - xmin)/H);
jj = ceil((ysou - ymin)/H);
% Index in x and y direction of the box containing each point
bin = (jj-1)*Nx + ii;

% bin = reshape(bin,N1,nb);

% Find bin of each point using lexiographic ordering (x then y)

fpt = zeros(Nbins,nb);
lpt = zeros(Nbins,nb);
% allocate space for storing first and last points
[binsort,permute] = sort(bin);
% build permute.  Need binsort to find first and last points in each box


%%%
% binsort = reshape(binsort,N1,nb);
% permute = reshape(permute,N1,nb);
% bin = reshape(bin,N1,nb);
% xsou = reshape(xsou,N1,nb);
% ysou = reshape(ysou,N1,nb);
%%%



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
% Construct first and last point in each box corresponding to each
% boundary.  The order is based on permute.  For example,
% permute(fpt(ibox,k)),...,permute(lpt(ibox,k)) is the set of all points
% from boundary k contained in box ibox

neigh = zeros(Nbins,9);

%Do corners first
neigh(1,1:4) = [1 2 Nx+1 Nx+2];
% bottom left corner
neigh(Nx,1:4) = [Nx Nx-1 2*Nx 2*Nx-1];
% bottom right corner
neigh(Nbins-Nx+1,1:4) = [Nbins-Nx+1 Nbins-Nx+2 ...
Nbins-2*Nx+1 Nbins-2*Nx+2];
% top left corner
neigh(Nbins,1:4) = [Nbins Nbins-1 Nbins-Nx Nbins-Nx-1];
% top right corner

for j = 2:Nx-1
  neigh(j,1:6) = j + [-1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of bottom row

for j = Nbins-Nx+2:Nbins-1
  neigh(j,1:6) = j + [-1 0 1 -Nx-1 -Nx -Nx+1];
end
% neighbors of top row

for j=Nx+1:Nx:Nbins-2*Nx+1
  neigh(j,1:6) = j + [-Nx -Nx+1 0 1 Nx Nx+1];
end
% neighbors of left column

for j=2*Nx:Nx:Nbins-Nx
  neigh(j,1:6) = j + [-Nx-1 -Nx -1 0 Nx-1 Nx];
end
% neighbors of right column

J = (Nx + 1:Nbins - Nx);
J = J(mod(J-1,Nx)~=0);
J = J(mod(J,Nx)~=0);
% J is the index of boxes that are not on the boundary
for j=J
  neigh(j,:) = j + [-Nx-1 -Nx -Nx+1 -1 0 1 Nx-1 Nx Nx+1];
end
% neighbors of interior points
% TREE STRUCTURE IS COMPLETE


if (relate == 1 || relate == 3)
  for k = 1:nb
    distSS{k} = spalloc(N1,nb,0);
    % dist(n,k,j) is the distance of point n on boundary k to boundary j
    zoneSS{k} = spalloc(N1,nb,0);
    % near or far zone
    nearestSS{k} = spalloc(2*N1,nb,0);
    % nearest point using local interpolant
    icpSS{k} = spalloc(N1,nb,0);
    % index of closest discretization point
    argnearSS{k} = spalloc(N1,nb,0);
    % argument in [0,1] of local interpolant
    nearFibersSS{k} = {};
  end
  % New way of representing near-singular integration structure so that
  % we can use sparse matricies.


  % begin classifying points where we are considering
  % boundary to boundary relationships
  for k = 1:nb
    boxes = unique(bin(:,k));
    % Find all boxes containing points of boundary k
    boxes = neigh(boxes,:);
    % Look at all neighbors of boxes containing boundary k
    boxes = unique(boxes(:));
    % Remove repetition
    boxes = boxes(boxes~=0);
    % Delete non-existent boxes that came up because of neigh

    K = [(1:k-1) (k+1:nb)];
    for k2 = K
      istart = fpt(boxes,k2);
      iend = lpt(boxes,k2);
      istart = istart(istart ~= 0);
      iend = iend(iend ~= -1);
      % Find index of all points in neighboring boxes of boundary k that
      % are in boundary k2

      neighpts = zeros(sum(iend-istart+1),1);
      % Allocate space to assign possible near points
      is = 1;
      for j=1:numel(istart)
        ie = is + iend(j) - istart(j);
        neighpts(is:ie) = permute(istart(j):iend(j),k2);
        is = ie + 1;
      end
      % neighpts contains all points on boundary k2 that are in
      % neighboring boxes to boundary k

      neighpts = sort(neighpts);
      % sorting should help speedup as we won't be jumping around
      % through different boxes

      n = 0;
      for i=1:numel(neighpts)
        ipt = neighpts(i);
        ibox = bin(ipt,k2);
        % box containing ipt on boundary k2
        if (ibox ~= n)
          n = ibox;
          % Check if we need to move to a new box
          neighbors = neigh(ibox,:);
          % neighbors of this box
          neighbors = neighbors(neighbors~=0);
          % Remove non-existent neighbors
          istart = fpt(neighbors,k);
          iend = lpt(neighbors,k);
          istart = istart(istart ~= 0);
          iend = iend(iend ~= -1);
          % Find points on boundary k in neighboring boxes
          neighpts2 = zeros(sum(iend-istart+1),1);
          is = 1;
          for j=1:numel(istart)
            ie = is + iend(j) - istart(j);
            neighpts2(is:ie) = permute(istart(j):iend(j),k);
            is = ie + 1;
          end
          % neighpts2 contains all points on boundary k that
          % are in neighboring box of ibox
        end % decide if we need to switch boxes

        [d0,d0loc] = min((xsou(ipt,k2) - xsou(:,k)).^2 + ...
                         (ysou(ipt,k2) - ysou(:,k)).^2);
        % Find minimum distance between ipt on boundary k2 to possible
        % closest points on boundary k
        d0 = sqrt(d0);
        % Save on not taking the square root on a vector but instead on
        % a single real number

        icpSS{k}(ipt,k2) = d0loc;
        if (d0 < 2*h);
          [distSS{k}(ipt,k2),nearestx,nearesty,...
                argnearSS{k}(ipt,k2)] = ...
            o.closestPnt([xsou;ysou],xsou(ipt,k2),ysou(ipt,k2),...
                k,icpSS{k}(ipt,k2));
            nearestSS{k}(ipt,k2) = nearestx;
            nearestSS{k}(ipt+N1,k2) = nearesty;
            
          nearFibersSS{k}= [nearFibersSS{k},k2];
          % Find closest point along a local interpolant using Newton's
          % method.
            
          if (distSS{k}(ipt,k2) < h)
            zoneSS{k}(ipt,k2) = 1;
          end
          % Point ipt of boundary k2 is in the near zone of
          % boundary k
        end
      end % ipt
    end % k2
    nftmp = nearFibersSS{k};
    nearFibersSS{k} = unique([nftmp{:}]);
  end % k

  NearSelf.dist = distSS;
  NearSelf.zone = zoneSS;
  NearSelf.nearest = nearestSS;
  NearSelf.icp = icpSS;
  NearSelf.argnear = argnearSS;
  NearSelf.nearFibers = nearFibersSS;
  % Store everything in the structure NearSelf.  This way it is
  % much cleaner to pass everything around

end % relate == 1 || relate == 3


% Bin target points with respect to the source points
if (relate == 2 || relate == 3)
  Np2 = bd2.N; % number of additional targets
  nb2 = bd2.nb; % number of additional boundaries
  X2 = bd2.X; % additional target points
  [xtar,ytar] = oc.getXY(X2);

  for k = 1:nb
    distST{k} = spalloc(N1,nb2,0);
    % dist(n,k,j) is the distance of point n on boundary k to
    zoneST{k} = spalloc(N1,nb2,0);
    % near or far zone
    nearestST{k} = spalloc(2*N1,nb2,0);
    % nearest point using local interpolant
    icpST{k} = spalloc(N1,nb2,0);
    % index of closest discretization point
    argnearST{k} = spalloc(N1,nb2,0);
    % argument in [0,1] of local interpolant
  end
  % New way of representing near-singular integration structure so that
  % we can use sparse matricies.

  itar = ceil((xtar - xmin)/H);
  jtar = ceil((ytar - ymin)/H);
  [indx,indy] = find((itar >= 1) & (itar <= Nx) & ...
  (jtar >= 1) & (jtar <= Ny));
  % Only have to consider xx(ind),yy(ind) since all other points
  % are not contained in the box [xmin xmax] x [ymin ymax]

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

      neighpts = zeros(sum(iend-istart+1),1);
      % Allocate space to assign possible near points
      if numel(neighpts) > 0
        % it is possible of the neighboring boxes to contain no points.
        is = 1;
        for j = 1:numel(istart)
          ie = is + iend(j) - istart(j);
          neighpts(is:ie) = permute(istart(j):iend(j),k);
          is = ie + 1;
        end
        % Set of potentially nearest points to (xtar(jj),ytar(jj))

        [d0,d0loc] = min((xtar(ii,jj) - xsou(neighpts,k)).^2 + ...
                         (ytar(ii,jj) - ysou(neighpts,k)).^2);
        % find closest point and distance between (xtar(jj),ytar(jj))
        % and boundary k.  Only need to look at points in neighboring
        % boxes
        d0 = sqrt(d0);
        icpST{k}(ii,jj) = neighpts(d0loc);

        if d0 < 2*h
          [distST{k}(ii,jj),nearestx,nearesty,argnearST{k}(ii,jj)] = ...
            o.closestPnt([xsou;ysou],xtar(ii,jj),ytar(ii,jj),...
            k,icpST{k}(ii,jj));
          nearestST{k}(ii,jj) = nearestx;
          nearestST{k}(ii+Np2,jj) = nearesty;
          if distST{k}(ii,jj) < h
            zoneST{k}(ii,jj) = 1;
            % (xtar(ii,jj),ytar(ii,jj)) is in the near zone of boundary
            % k
          end
        end % d0 < 2*h
      end % numel(neighpts) > 0

    end % ii and jj
  end % k

  NearOther.dist = distST;
  NearOther.zone = zoneST;
  NearOther.nearest = nearestST;
  NearOther.icp = icpST;
  NearOther.argnear = argnearST;
  NearOther.nearFibers = [];
  % store near-singluar integration requirements in structure NearOther

end % relate == 2 || relate == 3

end % getZone


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist,nearestx,nearesty,theta] = closestPnt(~,X,...
        xTar,ytar,k,icp)
% [dist,nearestx,nearesty,theta] = closestPnt(X,xTar,ytar,k,icp)
% computes the closest point on boundary k to (xTar,ytar)
% using a Lagrange interpolant.  icp is the index of the closest
% point on the discrete mesh which is used as an initial guess

N = size(X,1)/2; % Number of points per boundary 
A = poten.lagrangeInterp;
interpOrder = size(A,1);
% need interpolation matrix and its size

p = ceil((interpOrder+1)/2);
% Accommodate for either an even or odd number of interpolation
% points
pn = mod((icp-p+1:icp-p+interpOrder)' - 1,N) + 1;
% band of points around icp.  The -1,+1 combination sets index
% 0 to N as required by the code

px = A*X(pn,k); % polynomial interpolant of x-coordinate
py = A*X(pn+N,k); % polynomial interpolant of y-coordinate
Dpx = px(1:end-1).*(interpOrder-1:-1:1)';
Dpy = py(1:end-1).*(interpOrder-1:-1:1)';
D2px = Dpx(1:end-1).*(interpOrder-2:-1:1)';
D2py = Dpy(1:end-1).*(interpOrder-2:-1:1)';
% To do Newton's method, need two derivatives

theta = 1/2;
% midpoint is a good initial guess
for newton = 1:2
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
  % Using filter is the same as polyval, but it is much
  % faster when only requiring a single polyval such as here.

  newtonNum = (zx-xTar)*Dzx + (zy-ytar)*Dzy;
  % numerator of Newton's method
  newtonDen = (zx-xTar)*D2zx + (zy-ytar)*D2zy + ...
      Dzx^2 + Dzy^2;
  % denominator of Newton's method
  theta = theta - newtonNum/newtonDen;
  % one step of Newton's method
end
% Do a few (no more than 3) Newton iterations

nearestx = filter(1,[1,-theta],px);
nearestx = nearestx(end);
nearesty = filter(1,[1,-theta],py);
nearesty = nearesty(end);
dist = sqrt((nearestx - xTar)^2 + (nearesty - ytar)^2);
% Compute nearest point and its distance from the target point

end % closestPnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnb,ynb] = nearbyCurves(geom)

oc = curve;
% extract x and y coordinates
[x,y] = oc.getXY(geom.X);
[xc,yc] = oc.getXY(geom.center);

x_nb = zeros(geom.N,geom.nb);
y_nb = zeros(geom.N,geom.nb);

for k = 1:geom.nb
  radii = 1 + 10*max(geom.sa(:,k))*2*pi/geom.N;
  % inflate the radii slightly
  xnb(:,k) = radii*(x(:,k) - xc(k)) + xc(k);
  ynb(:,k) = radii*(y(:,k) - yc(k)) + yc(k);
end


end % nearbyCurves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = yukawaRHS(geom)
% Build right-hand side for the Yukawa equaiton solver

rhs = zeros(geom.N,geom.nb);
oc = curve;
[xc,yc] = oc.getXY(geom.center);
[x,y] = oc.getXY(geom.X);

rhs = zeros(geom.N*geom.nb,1);

for i = 1:2
  rhs = rhs + besselk(0,sqrt((x(:)-xc(i)).^2 + (y(:)-yc(i)).^2)/geom.rho)/...
              besselk(0,geom.radii(i)/geom.rho);
end

rhs = rhs(:);

end % yukawaRHS


end % methods



end %capsules



