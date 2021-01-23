clear
% clf
format long e
format compact

prams.N = 16;

dt      = 0.2;
prams.m = 500*10;      % number of time steps
prams.T = prams.m*dt; % time horizon

prams.order = 2;      % time-stepping order 

prams.rho = 2.0;      % screen length
prams.gam = 1.0;      % molecular tension

prams.RepulLength   = 0.5; % repulsion length
prams.RepulStrength = 4.0; % repulsion strength

options.farField  = 'shear';
options.shearRate = 0.001;
options.saveData  = false;
options.fileBase  = 'shear';
options.append    = false;
options.inear     = true;
options.usePreco  = false;
options.verbose   = true;
options.timeOrder = 2;
options.gmresTol  = 1e-10;
options.usePlot   = true; %false;
options.plotAxis  = [-20 20 -20 20];

% tracers 
axs = options.plotAxis;
N = 2000;
xx = (axs(2) - axs(1)).*rand(N,1) + axs(1);
yy = (axs(4) - axs(3)).*rand(N,1) + axs(3);
XX = [xx yy];
save("-ascii", "tracers.dat", "XX");

% initial centers
data = load('N52_0.dat');
x = data(:,1)';
y = data(:,2)';
xc = [x;y];

prams.nb = size(xc,2); % number of bodies

[options,prams] = initRigid2D(options,prams);

tau   = data(:,3)';
radii = 1.0*ones(1,prams.nb);
ar    = 1.0*ones(1,prams.nb);

prams.tau   = tau;
prams.radii = radii;
prams.ar    = ar;

[Xfinal, trajectory] = rigid2D(options,prams,xc,tau);
