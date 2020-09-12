classdef misc

%
% Some miscellaneous functions 
%
%
%

methods

function V = trapz2(o, x1, x2, z)

s = size(x1);
m = s(1); n = s(2);

I = 1:m-1;
J = 1:n-1;

%only one pair of substractions will be nonzero here. 

dx1 = x1(I+1,J) - x1(I,J) + x1(I,J+1) - x1(I,J); 
dx2 = x2(I+1,J) - x2(I,J) + x2(I,J+1) - x2(I,J);

V = 0.25*( z(I,J) + z(I+1,J) + z(I,J+1) + z(I+1,J+1) ).*dx1.*dx2;

V = sum(sum(V));

end

end %methods

end %classdef