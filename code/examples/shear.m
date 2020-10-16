% clc
% close all

clear
clf
format long e
format compact

% for jj=0:3
% prams.N = 2^jj*32; % points per body

prams.N = 64;

prams.T = 10.0; % time horizon
prams.m = 1000; % number of time steps

prams.order = 2; % time stepping order

prams.rho = 4.0; % screen length
prams.RepulLength = 0.3; % repulsion length
prams.RepulStrength = 0.5; % repulsion strength

options.farField  = 'shear';
options.shearRate = 0;
options.janusbc   = 2;        % put power of function here
options.saveData  = true;
options.fileBase  = 'shear';
options.append    = false;
options.inear     = true;
options.usePreco  = false;
options.verbose   = true;
options.gmresTol  = 1e-8;
options.usePlot   = true;
options.plotAxis  = 5*[-1 1 -1 1];

% initial configuration
xc = [0 2.4 3.4;-0.5 0.0 0.5];
tau   = [0.6*pi 0.6*pi pi]; 
radii = [0.5 0.5 0.5];
ar    = [1 1 1];

%xc = [0 2.4;-2.5 0.0];
%tau   = [0.6*pi 0.6*pi]; 
%radii = [0.5 0.5];
%ar    = [1 1];

%tau = [pi/2];
%radii = [0.5]
%ar = [2];
%xc = [0;0];

prams.nb = size(xc,2); % number of bodies
[options,prams] = initRigid2D(options,prams);

prams.tau   = tau;
prams.radii = radii;
prams.ar    = ar;

[Xfinal, trajectory] = rigid2D(options,prams,xc,tau);
%end
