clear
% clf
format long e
format compact

prams.N = 32;

dt      = 0.2;
prams.m = 500*10;      % number of time steps
prams.T = prams.m*dt; % time horizon

prams.order = 2;      % time-stepping order 

prams.rho = 2.0;      % screen length
prams.gam = 1.0;      % molecular tension

prams.RepulLength   = 0.5; % repulsion length
prams.RepulStrength = 4.0; % repulsion strength

options.farField  = 'shear';
<<<<<<< HEAD
options.shearRate = 0.2;
=======
options.shearRate = 0.1;
>>>>>>> 9f7fe9d2af7a26f4ffebff473627289c70eed513
options.saveData  = false;
options.fileBase  = 'shear';
options.append    = false;
options.inear     = true;
options.usePreco  = false;
options.verbose   = true;
options.timeOrder = 1;
options.gmresTol  = 1e-10;
options.usePlot   = true; %false;
options.plotAxis  = 3*[-3 3 -3 3];

% initial centers
data = load('N18_0.dat');
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
