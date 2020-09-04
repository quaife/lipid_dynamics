addpath ../src

oc = curve;
N = 128;
theta = (0:N-1)'*2*pi/N;
X = [cos(theta); 2*sin(theta)];
sa = sqrt(sin(theta).^2 + 4*cos(theta).^2);
IK = oc.modes(N,1);

f = exp(cos(theta));
dfds = oc.arcDeriv(f,1,sa,IK);

plot(dfds)

