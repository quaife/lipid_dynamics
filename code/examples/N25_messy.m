% Self assembly for initially randomly placed particles
% with the von Mises boundary condition. 
%
% The boundary conditions use varied 

clear
% clf
format long e
format compact

prams.N = 16;      % number of point per body
prams.nb = 25;     % number of bodies
prams.sstep = 0;   % starting step

radii = [1.0342    1.0017    0.9652    0.9881    1.0321    0.9671    0.9830    1.0466    1.0306    0.9722    1.0500    0.9564    0.992 0.9904    0.9900    0.9612    0.9924    1.0114    1.0488    0.9720    0.9854    0.9766    0.9791    0.9688    0.9523];
ar    = radii([10:25 1:9]);

prams.rho = 2.0;      % screen length
prams.gam = 4.0;      % molecular tension

prams.RepulLength   = 0.5; % repulsion length
prams.RepulStrength = 4.0; % repulsion strength

prams.bcShift       = 1 + 10*(flip(radii) - 1); % shift constant for yukawaRHS
prams.bcType        = 'vonMises'; % options: 'cosine'; 'vonMises'

dt      = 0.3;
prams.m = 10000;      % number of time steps
prams.T = prams.m*dt; % time horizon
prams.order = 2;      % time-stepping order


EXTERNAL_OPTIONS = 1;
shear_inputEXT
