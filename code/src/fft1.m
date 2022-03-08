classdef fft1 < handle
% class implements fft transformations.  This includes computing
% the fourier differentiation matrix, doing interplation required
% by Alpert's quadrature rules, and defining the Fourier frequencies

properties
N; % Number of points in the incoming periodic functions
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = diffFT(f,IK)
% f = diffFT(f,IK) Computes the first derivative of an array of 
% periodic functions f using fourier transform. The f(:,1) is the 
% first function, f(:,2) is the second function, etc.
% IK is used to speed up the code.  It is the index of the fourier
% modes so that fft and ifft can be used
% 
% EXAMPLE: 
%  N = 24; h = 2*pi/N; x = h*(1:N)';
%  IK = 1i*[(0:N/2-1) 0 (-N/2+1:-1)]';
%  v = exp(sin(x)); vprime = cos(x).*v;
%  subplot(2,2,1),plot(x,v,'.-','markersize',13);
%  axis([0 2*pi 0 3]);
%
%  w = fft1.diffFT(v,IK);
%  error = norm(w-vprime,inf);
%  subplot(2,2,2), plot(x,w,'.-','markersize',13);

f = real(ifft(IK.*fft(f)));

end % end diffFT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



end % methods (static)

end % classdef
