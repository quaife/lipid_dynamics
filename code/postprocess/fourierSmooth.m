function z = fourierSmooth(z,sig)

N = numel(z);
% sig = 5;

freq = (-N/2:N/2-1)';
gaussh = exp(-freq.^2/2/sig^2);

zh = fftshift(fft(z));
z = ifft(ifftshift(zh.*gaussh));