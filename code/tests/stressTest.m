addpath ../src

prams.N = 32;
prams.nb = 2;
prams.radii = [1 1];
prams.rho = 1;
prams.gam = 1;
prams.ar = [2 1.5];
prams.RepulLength = 1;
prams.RepulStrength = 1;

xc = [[+1;0] [-1;3]];
tau = [-1 +1];
geom = capsules(prams,xc,tau);

op = poten(prams.N,prams.rho);


% define a somewhat arbitrary density function
theta = (0:prams.N-1)'*2*pi/prams.N;
f = [exp(cos(theta));sin(5*sin(theta))];
f = [f f];
%f = [cos(theta);zeros(prams.N,1)];

dx = 0.01;
[xtar,ytar] = meshgrid(2:dx:4,2:dx:4);
Ntar = numel(xtar);
[~,velocity] = op.exactStokesDL(geom,f,[xtar(:);ytar(:)],(1:2));
[~,pressure] = op.exactStokesDLpressure(geom,f,[xtar(:);ytar(:)],(1:2));
%pressure = -pressure;
[~,stress] = op.exactStokesDLstress(geom,f,[xtar(:);ytar(:)],(1:2));

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

u = u(2:end-1,2:end-1);
v = v(2:end-1,2:end-1);
p = p(2:end-1,2:end-1);
sxx2 = -p + 2*ux;
sxy2 = uy + vx;
syy2 = -p + 2*vy;
sxx = sxx(2:end-1,2:end-1);
sxy = sxy(2:end-1,2:end-1);
syy = syy(2:end-1,2:end-1);
xtar = xtar(2:end-1,2:end-1);
ytar = ytar(2:end-1,2:end-1);

disp(norm(sxx(:) - sxx2(:),inf)/norm(sxx(:),inf))
disp(norm(sxy(:) - sxy2(:),inf)/norm(sxy(:),inf))
disp(norm(syy(:) - syy2(:),inf)/norm(syy(:),inf))
disp(norm(Lu(:) - px(:),inf)/norm(Lu(:),inf));
disp(norm(Lv(:) - py(:),inf)/norm(Lv(:),inf));

clf; hold on;
plot(geom.X(1:end/2,:),geom.X(end/2+1:end,:));
plot(xtar,ytar,'r.');
axis equal;


