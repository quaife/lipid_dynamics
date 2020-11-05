function [x,y]=getXY(X)
% [x,y] = getXY(X) get the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);

end % getXY