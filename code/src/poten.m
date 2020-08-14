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
  om;
  
end % properties

methods 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = poten(N)
% o = poten(N): constructor; N is the number of points per curve

o.interpMat = o.lagrangeInterp;
% load in the interpolation matrix which is precomputed with 7
% interpolation points

o.N = N;

end % poten: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingInt(o,geomSou,f,selfMat,Du,...
    NearStruct,kernel,kernelDirect,geomTar,tEqualS,idebug)
% LP = nearSingInt(geom,f,selfMat,NearStruct,kernel,kernelDirect,
% geomTar,tEqualS,idebug) computes a layer potential due to f at all
% points in geomTar.X.  If tEqualS==true, then the geomTar ==
% geomSou and the self-geom interaction is skipped.  selfMat is
% the diagonal of the potential needed to compute the layer potential of
% each geom indepenedent of all others.  kernel and kernelDirect are
% two (possibly the same) routines that compute the layer potential.
% kernelDirect always uses the direct method whereas kernel may use an
% FMM-accelerated method.  NearStruct is a structure containing the
% variables zone,dist,nearest,icp,argnear which are required by
% near-singular integration (they keep everything sorted and
% precomputed) Everything is in the 2*N x n format Can pass a final
% argument if desired so that plots of the near-singular integration
% algorithm are displayed

if (tEqualS && size(geomSou.X,2) == 1)
  LP = zeros(size(geomSou.X));
  return
end
% only a single geom, so velocity on all other geoms will always
% be zero

if (nargin == 10)
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

% upsample to N^(3/2).  
Xup = [interpft(Xsou(1:Nsou,:),Nup); interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];
fup = [interpft(f(1:Nsou,:),Nup); interpft(f(Nsou+1:2*Nsou,:),Nup)];

geomUp = capsules([],Xup);
% Build an object with the upsampled geom

interpOrder = size(o.interpMat,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right

if tEqualS % sources == targets
  if nsou > 1
    if (strfind(char(kernel),'fmm'))
      farField = kernel(geomUp,fup,Du);
      farField = farField(1:Nup/Ntar:end,:);
      % evaluate layer potential at all targets except ignore the
      % diagonal term
    else
      for k = 1:nsou
        K = [(1:k-1) (k+1:nsou)];
        [~,farField(:,k)] = kernelDirect(geomUp,fup,Du,Xtar(:,k),K);
      end
      % This is a huge savings if we are using a direct method rather
      % than the fmm to evaluate the layer potential.  The speedup is
      % more than N^{1/2}, where N is the resolution of the geoms
      % that we are computing with
    end
  else
    farField = zeros(2*Ntar,ntar);
  end

else % sources ~= targets
    [~,farField] = kernel(geomUp,fup,Du,Xtar,1:nsou);
    % evaluate layer potential due to all 'geoms' at all points in
    % Xtar
end
% Use upsampled trapezoid rule to compute layer potential

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
            % closest point on geom k1 to each point on geom k2
            % that is close to geom k1
            for j = 1:numel(J)
                pn = mod((indcp(j)-p+1:indcp(j)-p+interpOrder)' - 1,Nsou) + 1;
                % index of points to the left and right of the closest point
                v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
                    o.interpMat*vself(pn,k1));
                vel(J(j),k2,k1) = v(end);
                % x-component of the velocity at the closest point
                v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
                    o.interpMat*vself(pn+Nsou,k1));
                vel(J(j)+Ntar,k2,k1) = v(end);
                % y-component of the velocity at the closest point
            end
            %     compute values of velocity at required intermediate points
            %     using local interpolant
            
            if ((numel(J) + numel(fup)) >= 512 && numel(J) > 32)
                [~,potTar] = kernel(geomUp,fup,Du, [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
            else
                [~,potTar] = kernelDirect(geomUp,fup,Du,[Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
            end
            % Need to subtract off contribution due to geom k1 since its
            % layer potential will be evaulted using Lagrange interpolant of
            % nearby points
            nearField(J,k2) =  nearField(J,k2) - potTar(1:numel(J));
            nearField(J+Ntar,k2) =  nearField(J+Ntar,k2) - potTar(numel(J)+1:end);
            
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
                % Lagrange interpolation points coming off of geom k1 All
                % points are behind Xtar(J(i),k2) and are sufficiently far from
                % geom k1 so that the Nup-trapezoid rule gives sufficient
                % accuracy
            end
            
            if (numel(XLag)/2 > 100)
                [~,lagrangePts] = kernel(geomUp,fup,Du,XLag,k1);
            else
                [~,lagrangePts] = kernelDirect(geomUp,fup,Du,XLag,k1);
            end
            % evaluate velocity at the lagrange interpolation points
            
            for i = 1:numel(J)
                Px = o.interpMat*[vel(J(i),k2,k1) ...
                    lagrangePts(i,:)]';
                Py = o.interpMat*[vel(J(i)+Ntar,k2,k1) ...
                    lagrangePts(i+numel(J),:)]';
                % Build polynomial interpolant along the one-dimensional
                % points coming out of the geom
                dscaled = full(dist{k1}(J(i),k2)/(beta*h(k1)*(interpOrder-1)));
                % Point where interpolant needs to be evaluated
                
                v = filter(1,[1 -dscaled],Px);
                nearField(J(i),k2) = nearField(J(i),k2) +  v(end);
                
                v = filter(1,[1 -dscaled],Py);
                nearField(J(i)+Ntar,k2) = nearField(J(i)+Ntar,k2) + v(end);
                
                
                if idebug
                    figure(2); clf; hold on;
                    plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.','markersize',10)
                    plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.','markersize',10)
                    plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.','markersize',10)
                    plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx','markersize',10)
                    plot(XLag(i,:),XLag(numel(J)+i,:),'gx','markersize',10)
                    axis equal
                    
                    figure(1); clf; hold on
                    plot((0:interpOrder-1)*beta*h(k1),...
                        real([vel(J(i),k2,k1) lagrangePts(i,:)]),'g-o')
                    plot((0:interpOrder-1)*beta*h(k1),...
                        real([vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)]),'r--o')
                    
                    figure(3)
                    clf
                    hold on
                    plot(f(1:Nsou,k1));
                    plot(f(Nsou+1:2*Nsou,k1));
                    
                    drawnow;
                    pause
                    
                end
                % DEBUG: PASS IN idebug=true INTO THIS ROUTINE AND THEN YOU CAN SEE
                % THE INTERPOLATION POINTS AND CHECK THE SMOOTHNESS OF THE INTERPOLANT
                
            end % i
        end % numel(J) ~= 0
        % Evaluate layer potential at Lagrange interpolation
        % points if there are any
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

normal = [wall.xt(wall.N+1:2*wall.N,:);-wall.xt(1:wall.N,:)]; % Normal vector
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
function DLP = exactStokesDLdiag(~,geom,D,f)
% DLP = exactStokesDLdiag(geom,f,K) computes the diagonal term of
% the double-layer potential due to f around all geoms.  Source and
% target points are the same.  This uses trapezoid rule with the
% curvature at the diagonal in order to guarantee spectral accuracy.
% This routine can either compute the double-layer potential
% matrix-free, which may upsample the number of source points.  Or, if
% the matrix D is passed in and anti-aliasing is not requested, it will
% simply do the matrix-vector product with the precomputed matrix D.

DLP = zeros(2*geom.N,geom.nb);
for k = 1:geom.nb
  A = D(:,:,k);
  DLP(:,k) = A * f(:,k);
end
% 
% DLP = permute(sum(bsxfun(@times, D, permute(f,[3 1 2])),2), [1 3 2]);

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
function [stokesDLP,stokesDLPtar] = exactStokesDL(~, geom, f, ~, Xtar,K1)
% [stokesDLP,stokesDLPtar] = exactStokesDL(geom,f,Xtar,K1) computes the
% double-layer potential due to f around all parts of the geometry
% except itself.  Also can pass a set of target points Xtar and a
% collection of geom K1 and the double-layer potential due to components
% of the geometry in K1 will be evaluated at Xtar.  Everything but Xtar
% is in the 2*N x n format Xtar is in the 2*Ntar x ncol format

if nargin == 6
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

stokesDLP = zeros(2*geom.N,geom.n);
if (nargin == 4 && geom.n > 1)
  for k = 1:geom.n
    K = [(1:k-1) (k+1:geom.n)];
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
% END OF ROUTINES THAT EVALUATE LAYER-POTENTIALS
% WHEN SOURCES ~= TARGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
