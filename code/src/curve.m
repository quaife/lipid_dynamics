classdef curve 
% This class implements that basic calculus on the curve.
% The basic data structure is a matrix X in which the columns 
% represent periodic C^{\infty} closed curves with N points, 
% X(1:n,j) is the x-coordinate of the j_th curve and X(n+1:N,j) 
% is the y-coordinate of the j_th curve; here n=N/2
% X coordinates do not have to be periodic, but the curvature,
% normals, etc that they compute will be garbage.

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=getXY(o,X)
% [x,y] = getXY(X) get the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);
end % getXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=getXYperp(o,X)
% [x,y] = getXY(X) get the [x,y] component of curves X
N = size(X,1)/2;
x = X(N+1:end,:);
y = -X(1:N,:);
end % getXYperp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = setXY(o,x,y)
% V = setXY(x,y) set the [x,y] component of vector V on the curve
N = size(x,1);
V=zeros(2*N,size(x,2));
V(1:N,:) = x;
V(N+1:end,:) = y;
end % setXY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dx,Dy]=getDXY(o,X)
% [Dx,Dy]=getDXY(X), compute the derivatives of each component of X 
% these are the derivatives with respect the parameterization 
% not arclength
[x,y] = o.getXY(X);
N = size(x,1);
nv = size(x,2);
IK = o.modes(N,nv);
Dx = o.diffFT(x,IK);
Dy = o.diffFT(y,IK);
end % getDXY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jacobian,tangent,curvature] = diffProp(o,X)
% [jacobian,tangent,curvature] = diffProp(X) returns differential
% properties of the curve each column of the matrix X. Each column of 
% X should be a closed curve defined in plane. The tangent is the 
% normalized tangent.
%
N = size(X,1)/2;
nv = size(X,2);
IK = o.modes(N,nv);

% get the x y components
[Dx,Dy] = o.getDXY(X);

jacobian = sqrt( Dx.^2 + Dy.^2 ); 

if nargout>1  % if user requires tangent
  tangent = o.setXY( Dx./jacobian, Dy./jacobian);
end

if nargout>2  % if user requires curvature
  DDx = curve.arcDeriv(Dx,1,ones(N,1),IK);
  DDy = curve.arcDeriv(Dy,1,ones(N,1),IK);
  curvature = (Dx.*DDy - Dy.*DDx)./(jacobian.^3);
end
% curvature of the curve

end % diffProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [area,length] = geomProp(o,X)
% [area length] = geomProp(X) calculate the length, area and the reduced
% volume of domains inclose by columns of X. 

[x,y] = o.getXY(X);
N = size(x,1);
[Dx,Dy] = o.getDXY(X);

length = sum( sqrt(Dx.^2 + Dy.^2) )*2*pi/N;
area = sum( x.*Dy - y.*Dx)*pi/N;

end % geomProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = redistributeArcLength(o, X)
% X = resdistributeArcLength(o,X) resdistributes the geometry shape
% eqiuspaced in arclength

N = size(X,1)/2;
modes = [(0:N/2-1) (-N/2:-1)];

jac = o.diffProp(X);

tol = 1e-10;

if norm(jac - mean(jac),inf) > tol*mean(jac)
  theta = o.arcLengthParameter(X(1:end/2),...
      X(end/2+1:end));
  zX = X(1:end/2) + 1i*X(end/2+1:end);
  zXh = fft(zX)/N;
  zX = zeros(N,1);
  for j = 1:N
    zX = zX + zXh(j)*exp(1i*modes(j)*theta);
  end
  X = o.setXY(real(zX),imag(zX));
  % redistribute the geometry so that it is
  % equispaced in arclength
end

end % redistributeArcLength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta,arcLength] = arcLengthParameter(o,x,y)
% theta = arcLengthParamter(o,x,y) finds a discretization of parameter
% space theta so that the resulting geometry will be equispaced in
% arclength

uprate = 1;
N = numel(x);
Nup = uprate*N;
t = (0:Nup-1)'*2*pi/Nup; % this is not correct when you iterate
x = interpft(x,Nup);
y = interpft(y,Nup);
[~,len] = o.geomProp([x;y]);
% find total perimeter
[Dx,Dy] = o.getDXY([x;y]);
% find derivative
arc = sqrt(Dx.^2 + Dy.^2);
arch = fft(arc);
modes = -1i./[(0:Nup/2-1) 0 (-Nup/2+1:-1)]';
modes(1) = 0;
modes(Nup/2+1) = 0;
arcLength = real(ifft(modes.*arch) - sum(modes.*arch/Nup) + ...
    arch(1)*t/Nup);

z1 = [arcLength(end-6:end)-len;arcLength;arcLength(1:7)+len];
z2 = [t(end-6:end)-2*pi;t;t(1:7)+2*pi];
% put in some overlap to account for periodicity

theta = [interp1(z1,z2,(0:N-1)'*len/N,'spline')];

end % arcLengthParamter

end % methods


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)

function df = arcDeriv(f,m,sa,IK)
% df = arcDeriv(f,m,s,IK) is the arclength derivative of order m.
% f is a matrix of scalar functions (each function is a column)
% f is assumed to have an arbitrary parametrization
% sa = d s/ d a, where a is the aribtrary parameterization
% IK is the fourier modes which is saved and used to accelerate 
% this routine

if nargin<3, m=1; end;
if nargin<4, sa=ones(size(f,1),1); end;

col = size(f,2);
if col == 2
  sa = [sa sa];
elseif col == 3
  sa = [sa sa sa];
elseif col > 3
  sa = repmat(sa,1,col);
end
% This is much faster than always using repmat

df = f;
for j=1:m
  df = sa.*real(ifft(IK.*fft(df)));
end
% compute the mth order derivative

end % arcDeriv



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function df = diffFT(f,IK)
% df = diffFT(f,IK) Computes the first derivative of an array of 
% periodic functions f using fourier transform. The f(:,1) is the 
% first function, f(:,2) is the second function, etc.
% IK is used to speed up the code.  It is the index of the fourier
% modes so that fft and ifft can be used
% 

df = real(ifft(IK.*fft(f)));

end % end diffFT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IK = modes(N,nv)
% IK = modes(N) builds the order of the fourier modes required for using
% fft and ifft to do spectral differentiation

IK = 1i*[(0:N/2-1) 0 (-N/2+1:-1)]';
% diagonal term for Fourier differentiation with the -N/2 mode
% zeroed to avoid Nyquist frequency

if nv == 2
  IK = [IK IK];
elseif nv == 3
  IK = [IK IK IK];
elseif nv > 3
  IK = repmat(IK,1,nv);
end

end % modes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wn = wn_PnPoly(o, p1, p2, v1, v2, N)

% wn_PnPoly(): winding number test for a point in a polygon
%              useful for zeroing out interior
%      Input:   P = a point (p1, p2)
%               V = (v1, v2) = vertex points of a polygon with V(1) = V(N+1)
%      Return:  wn = the winding number (=0 only when P is outside)

    wn = 0*p1;

    for i = 1:N
    
        wn = wn + (v2(i) <= p2).*(v2(i+1)  > p2).*(o.isLeft(v1(i), v2(i), v1(i+1), v2(i+1), p1, p2) > 0);
        wn = wn - (v2(i)  > p2).*(v2(i+1) <= p2).*(o.isLeft(v1(i), v2(i), v1(i+1), v2(i+1), p1, p2) < 0);

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = isLeft(o, p1, p2, q1, q2, r1, r2)
    %    Input:  three points p, q, and r
    %    Return: >0 for r left of the line through p and q
    %            =0 for r  on the line
    %            <0 for r  right of the line     
    %http://geomalgorithms.com/a03-_inclusion.html

    out = (q1 - p1) .* (r2 - p2) - (r1 -  p1) .* (q2 - p2) ;

end 




end % methods (Static)

end % classdef


