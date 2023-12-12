function [Xmid, Area, Len, Ra] = MidArcLen(x,y,tau,Nout,Nin)
% For a closed connected bilayer structure, this script generates the
% midplane configuration, area, arc length, and reduced area.
%
% Inputs:
% x = [x_in; x_out]; y = [y_in; y_out]; tau = [tau_in; tau_out];
% N = Nin + Nout;
%
% Outputs:
% Xmid = [xmid; ymid];
% Area: Enclosed area of the midplane curve
% Len : Midplane arc length
% Ra  : Reduced area


distance  = zeros(Nin,Nout);

% find the minimum distances from inner to outer particles.
for m = 1:Nin
    for n = 1:Nout
        distance(m,n) = norm([x(m) y(m)]-[x(n+Nin) y(n+Nin)]);    
    end
end

[val, ~] = min(distance);

% rough midplane positions 
xmid = x(1:Nin) + mean(val)/2*cos(tau(1:Nin));
ymid = y(1:Nin) + mean(val)/2*sin(tau(1:Nin));

[TH,R] = cart2pol(xmid,ymid);
THNEW = -pi:pi/32:(pi-pi/32);
RNEW = interp1(TH,R,THNEW,'spline');
XNEW = RNEW'.*cos(THNEW)';
YNEW = RNEW'.*sin(THNEW)';
XNEW = XNEW(~isnan(XNEW));
YNEW = YNEW(~isnan(YNEW));

% smooth midplane
z = fourierSmooth(XNEW+1i*YNEW,10);
xmid = real(z); ymid = imag(z);
Xmid = [xmid;ymid];

z = xmid'+1i*ymid';
modes = [0:32 -31:-1];
zh   = fft(z);
dzh  = 1i*modes.*zh;
dz   = ifft(dzh);
Len  = sum(abs(dz))/64*2*pi;
Area = 0.5*sum(real(z).*imag(dz)-imag(z).*real(dz))*2*pi/64;
Ra   = 4*pi*Area./Len.^2;


% test script: copy and run the following lines
%--------------------------------------------------------------------------
% data = load('test_vesicle.dat');
% x = data(:,1); y = data(:,2); tau = data(:,3);
% Nin = 26; Nout = 32;
% [Xmid, Area, Len, Ra] = MidArcLen(x,y,tau,Nout,Nin);
% plot(x,y,'o'); hold on
% plot(Xmid(1:64),Xmid(65:128),'o-')
%--------------------------------------------------------------------------

end