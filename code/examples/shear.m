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
prams.m = 400 ; % number of time steps

prams.order = 1;

prams.rho = 4.0; % screen length

options.farField  = 'shear';
options.janusbc   = 2;        % put power of function here
options.saveData  = true;
options.fileBase  = 'shear';
options.append    = false;
options.inear     = true;
options.usePreco  = false;
options.verbose   = true;
options.timeOrder = 1;
options.gmresTol  = 1e-8;
options.usePlot   = true;
options.plotAxis  = 5*[-1 1 -1 1];

% initial center
% xc = [(20*rand(10,1)-10)' ;(20*rand(10,1)-10)' ]; 
% xc = [-1.3;-1.8];
% xc = [-1 0.6; 0 -0.5];
% xc = [-1.5 1.5 1.5 -1.5; -1.5 -1.5 1.5 1.5];    
% xc = [0;0];

xc = [0 2.4 4.3; 0 0.0 2.5];

prams.nb = size(xc,2); % number of bodies

[options,prams] = initRigid2D(options,prams);

dtau  = 2*pi/prams.nb;
tau   = [0.6*pi 0.6*pi pi]; %(0:dtau:2*pi-dtau); % initial incliation angle
radii = [0.5 0.5 0.5];
ar    = [2 2 2];

prams.tau   = tau;
prams.radii = radii;
prams.ar    = ar;

[Xfinal, trajectory] = rigid2D(options,prams,xc,tau);
%end
