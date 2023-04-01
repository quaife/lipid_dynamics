addpath ../src

prams.N = 32;
prams.nb = 2;
prams.radii = [1 1];
prams.rho = 1;
prams.gam = 1;
prams.ar = [2 1.5];
prams.RepulLength = 1;
prams.RepulStrength = 1;
prams.bcShift = 0;
prams.bcType = 0;
prams.bcType        = 'cosine';
prams.NPeaks        = 1;
prams.shape         = 'star'; % options: 'circle'; 'star'
prams.petal         = 5;

xc = [[+1;0] [-1;3]];
tau = [-1 +1];
geom = capsules(prams,xc,tau);

op = poten(prams.N,prams.rho);

% define a somewhat arbitrary density function
theta = (0:prams.N-1)'*2*pi/prams.N;
f = [exp(cos(theta));sin(5*sin(theta))];
f = [f f];
force = 10*[1;2;3;-4];
torque = 10*[2;-3];

dx = 0.01;
[xtar,ytar] = meshgrid(2:dx:4,2:dx:4);
Ntar = numel(xtar);

% velocity due to the double-layer potential
[~,velocityDL] = op.exactStokesDL(geom,f,[xtar(:);ytar(:)],(1:2));
% velocity due to Stokeslets and Rotlets
[~,velocityRS] = op.StokesletRotlet(geom,force,torque,...
      [xtar(:);ytar(:)],(1:2));
% total velocity
velocity = velocityDL + velocityRS;

% pressure due to the double-layer potential
[~,pressureDL] = op.exactStokesDLpressure(geom,f,[xtar(:);ytar(:)],(1:2));
% pressure due to the Stokeslets and Rotlets
pressureRS = op.RSpressure(geom,force,torque,[xtar(:);ytar(:)],(1:2));
% total pressure
pressure = pressureDL + pressureRS;

% stress due to the double-layer potential
[~,stressDL] = op.exactStokesDLstress(geom,f,[xtar(:);ytar(:)],(1:2));
% stess due to the Stokeslets and Rotlets
stressRS = op.RSstress(geom,force,torque,[xtar(:);ytar(:)],(1:2));
% total stress
stress = stressDL + stressRS;

u = reshape(velocity(1:Ntar),size(xtar));
v = reshape(velocity(Ntar+1:end),size(xtar));
p = reshape(pressure,size(xtar));
sxx = reshape(stress(1:Ntar),size(xtar));
sxy = reshape(stress(Ntar+1:2*Ntar),size(xtar));
syy = reshape(stress(2*Ntar+1:end),size(xtar));

% gradient of the velocity
ux = (u(2:end-1,3:end) - u(2:end-1,1:end-2))/(2*dx);
uy = (u(3:end,2:end-1) - u(1:end-2,2:end-1))/(2*dx);
vx = (v(2:end-1,3:end) - v(2:end-1,1:end-2))/(2*dx);
vy = (v(3:end,2:end-1) - v(1:end-2,2:end-1))/(2*dx);

% Laplacian of velocity
Lu = (u(2:end-1,3:end) + u(2:end-1,1:end-2) + ...
      u(3:end,2:end-1) + u(1:end-2,2:end-1) - ...
      4*u(2:end-1,2:end-1))/dx^2;
Lv = (v(2:end-1,3:end) + v(2:end-1,1:end-2) + ...
      v(3:end,2:end-1) + v(1:end-2,2:end-1) - ...
      4*v(2:end-1,2:end-1))/dx^2;

% Gradient of pressure
px = (p(2:end-1,3:end) - p(2:end-1,1:end-2))/(2*dx);
py = (p(3:end,2:end-1) - p(1:end-2,2:end-1))/(2*dx);

% Laplacian of pressure
Lp = (p(2:end-1,3:end) + p(2:end-1,1:end-2) + ...
      p(3:end,2:end-1) + p(1:end-2,2:end-1) - ...
      4*p(2:end-1,2:end-1))/dx^2;

% remove boundary values
u = u(2:end-1,2:end-1);
v = v(2:end-1,2:end-1);
p = p(2:end-1,2:end-1);
sxx = sxx(2:end-1,2:end-1);
sxy = sxy(2:end-1,2:end-1);
syy = syy(2:end-1,2:end-1);
xtar = xtar(2:end-1,2:end-1);
ytar = ytar(2:end-1,2:end-1);

% stress tensor using finite differences
sxx2 = -p + 2*ux;
sxy2 = uy + vx;
syy2 = -p + 2*vy;

% compute the norms of the difference of the two methods to compute the
% stress
disp(norm(sxx(:) - sxx2(:),inf)/norm(sxx(:),inf))
disp(norm(sxy(:) - sxy2(:),inf)/norm(sxy(:),inf))
disp(norm(syy(:) - syy2(:),inf)/norm(syy(:),inf))
% compute the error in the momentum equation
disp(norm(Lu(:) - px(:),inf)/norm(Lu(:),inf));
disp(norm(Lv(:) - py(:),inf)/norm(Lv(:),inf));

clf; hold on;
plot(geom.X(1:end/2,:),geom.X(end/2+1:end,:));
plot(xtar,ytar,'r.');
axis equal;


