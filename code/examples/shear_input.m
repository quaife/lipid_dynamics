clear
% clf
format long e
format compact

prams.N = 32;

dt      = 0.1;
prams.m = 50;      % number of time steps
prams.T = prams.m*dt; % time horizon

prams.order = 2;      % time-stepping order 

prams.rho = 2.0;      % screen length
prams.gam = 1.0;      % molecular tension

prams.RepulLength   = 0.5; % repulsion length
prams.RepulStrength = 4.0; % repulsion strength

options.farField  = 'shear';
options.shearRate = 0.0;
options.saveData  = false;
options.fileBase  = 'shear';
options.append    = false;
options.inear     = true;
options.usePreco  = false;
options.verbose   = true;
options.timeOrder = 1;
options.gmresTol  = 1e-5;
options.usePlot   = true; %false;
options.plotAxis  = 30*[-1 1 -1 1];

% initial centers
data = load('N3_0.dat');
x = 5*data(:,1)';
y = data(:,2)';
xc = [x;y];

prams.nb = size(xc,2); % number of bodies

[options,prams] = initRigid2D(options,prams);

tau   = data(:,3)';
radii = 1*ones(1,prams.nb);
ar    = 1*ones(1,prams.nb);

prams.tau   = tau;
prams.radii = radii;
prams.ar    = ar;

[Xfinal, trajectory] = rigid2D(options,prams,xc,tau);
