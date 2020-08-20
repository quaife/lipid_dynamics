% clc
%close all
clf
format long e
format compact
prams.N = 128; % points per body

prams.T = 50; % time horizon
prams.m = 500; % number of time steps

prams.order = 1;

prams.rho = 1 ; % screen length

options.farField = 'shear';
options.janusbc = 2;        % put power of function here
options.saveData = true;
options.fileBase = 'shear';
options.append = false;
options.inear = true;
options.usePreco = false;
options.verbose = true;
options.timeOrder = 1;
options.gmresTol = 1e-8;
options.usePlot = true;
options.plotAxis = 5*[-1 1 -1 1];

% initial center
% xc = [(20*rand(10,1)-10)' ;(20*rand(10,1)-10)' ]; 
%xc = [-1.3;-1.8];
 xc = [-3 3;1 -1];
%xc = [-1.5 1.5 1.5 -1.5; -1.5 -1.5 1.5 1.5];    
%xc = [0;0];

prams.nb = size(xc,2); % number of bodies

[options,prams] = initRigid2D(options,prams);

dtau = 2*pi/prams.nb;
tau  = (0:dtau:2*pi-dtau); % initial incliation angle
radii = 2*ones(1,prams.nb);
ar = ones(1,prams.nb);

prams.tau = tau;
prams.radii = radii;
prams.ar = ar;
Xfinal = rigid2D(options,prams,xc,tau);

