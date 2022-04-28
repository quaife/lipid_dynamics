clear
figure(1); clf
format long e
% prams.N = 2^jj*32; % points per body

prams.N = 32;
prams.NPeaks = 1;
prams.shape = 'circle';
prams.petal = [];

prams.T = 20; % time horizon
prams.m = 2000; % number of time steps

prams.order = 2; % time stepping order

%prams.rho = 4.0;  % screen length
prams.rho = 1.0;  % screen length
prams.gam = 1.0;  % HAP Strength

prams.RepulLength = 0.3; % repulsion length
prams.RepulStrength = 0.5; % repulsion strength

prams.Nwall = 4*128; % points on solid wall

options.farField  = 'channel';
options.shearRate = 10;
% power of janus bc function
% options.janusbc   = 2;
options.saveData  = true;
options.saveVel   = true;
options.fileBase  = 'confined_tube';
options.append    = false;
options.verbose   = true;
options.gmresTol  = 1e-8;
options.usePlot   = true;
options.plotAxis  = [-15 15 -3.1 3.1];
options.confined  = true;
options.walls     = 'choke';

if 0
  xc = [0 2.8;1.9 +1.9];
  tau   = [0 pi/2]; 
  radii = [1 1];
  ar    = [1.0 +1.0];
end
if 0
  xc = [-8;0.5];
  tau = +pi/4;
  radii = 0.4;
  ar = 1;
end
if 0
  xx = linspace(-14,-11,4); xx = [xx xx];
  yy = [0.8*rand(1,4)+0.5 -0.8*rand(1,4)-0.5];
%  yy = [1.4*ones(1,4) -1.4*ones(1,4)];
  xc = [xx;yy];
  tau = 2*pi*rand(size(xx));
  radii = 0.4*ones(size(xx));
  ar = ones(size(xx));
end
if 1
  xx = [-14 -13.1 -12 -11.1];
  yy = [0 0.1 -0.8 -0.9];
  xc = [xx;yy];
  tau = zeros(1,4);
  radii = 0.4*ones(size(xx));
  ar = ones(size(xx));
end

prams.nb = size(xc,2); % number of bodies
% don't think this is actually used for anything in the src codes
prams.sstep = []; 
[options,prams] = initRigid2D(options,prams);

%tau   = tau + rand(1,prams.nb);
prams.radii = radii;
prams.ar    = ar;

%oc = curve;
%if options.confined
%  Xwalls = oc.wallsGeom(prams.Nwall,options.walls);
%else
%  Xwalls = [];
%end


[Xfinal, trajectory] = rigid2D(options,prams,xc,tau);
