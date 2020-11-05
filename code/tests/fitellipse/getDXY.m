function [Dx,Dy]=getDXY(X)
% [Dx,Dy]=getDXY(X), compute the derivatives of each component of X 
% these are the derivatives with respect the parameterization 
% not arclength
x = X(1:end/2,:);
y = X(end/2+1:end,:);
N = size(x,1);
nv = size(x,2);
IK = fft1.modes(N,nv);
Dx = fft1.diffFT(x,IK);
Dy = fft1.diffFT(y,IK);

end % getDXY