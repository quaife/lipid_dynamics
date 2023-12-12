% if ~EXTERNAL_OPTIONS

clear
% clf
format long e
format compact

% profile on

prams.N = 24;      % number of point per body
prams.nb = 5;      % number of bodies

dt      = 0.2;
prams.sstep = 0;   % starting step
prams.m = 1;       % number of time steps
prams.T = prams.m*dt; % time horizon
prams.order = 2;      % time-stepping order 

prams.rho = 2.0;      % screen length
prams.gam = 1.0;      % molecular tension

prams.RepulLength   = 0.5; % repulsion length
prams.RepulStrength = 4.0; % repulsion strength

% The energy is asymptotically equal to the \int_\sigma f^2 dS
% We normalize the boundary condition by a factor.
prams.bcShift       = 1.0*ones(1,prams.nb); % shift constant for yukawaRHS
prams.bcType        = 'cosine'; % options: 'cosine'; 'vonMises'
prams.NPeaks        = 1;   % for boundary condition 'cosine'

% end

options.farField  = 'shear'; 
% options: 'shear'; 'extensional'; 'parabolic'; 'taylorgreen'
options.shearRate = 0.00;

options.saveData  = false;
options.fileBase  = options.farField;

options.append    = false;
options.verbose   = true;
options.timeOrder = 2;
options.gmresTol  = 1e-10;

options.usePlot   = true;
options.tracer    = false;
options.plotAxis  = [-10 10 -10 10];


% shape parameters
prams.shape         = 'circle'; % options: 'circle'; 'star'
prams.petal         = 4;
% circle:  (x,y) = (ar*radii*cos(th), radii*sin(th))
% star  :  (x,y) = radii + ar*cos(petal*th)
radii = 1.0*ones(1,prams.nb);
ar    = 1.0*ones(1,prams.nb);


% 1. Please input the number of bodies and have the corresponding initial
% configuration file ready.
% 2. If a desired starting step is given, please modify prams.sstep and 
% prepare the corresponding configuration file w/ or w/o the tracer file.

if prams.sstep == 0
    data = load(['N' num2str(prams.nb) '_' num2str(prams.sstep) '.dat']);
    data = [data 0*zeros(prams.nb,3)];
    fileName = sprintf("../output/data/frames/N%d_%f_%d.dat", ...
                   prams.nb, options.shearRate, prams.sstep);
    save("-ascii", fileName, "data");
else
    fileName = sprintf("../output/data/frames/N%d_%f_%d.dat", ...
                   prams.nb, options.shearRate, prams.sstep);    
    data = load(fileName);
end

% initial centers and angles 
x = data(:,1)'; %x = x-mean(x);
y = data(:,2)'; %y = y-mean(y);
tau = data(:,3)';
xc = [x;y];


% tracers 
if options.tracer
    axs = options.plotAxis;
    
    if prams.sstep ==0
        N = 1000;
        xx = (axs(2) - axs(1)).*rand(N,1) + axs(1);
        yy = (axs(4) - axs(3)).*rand(N,1) + axs(3);
        XX = [xx yy];
        save("-ascii", "tracers.dat", "XX");
        fileName = sprintf("../output/data/frames/N%d_%f_%d.tracer", ...
                   prams.nb, options.shearRate, prams.sstep);
        save("-ascii", fileName, "XX");
    else
        fileName = sprintf("../output/data/frames/N%d_%f_%d.tracer", ...
                   prams.nb, options.shearRate, prams.sstep);
        XX = load(fileName);    
        save("-ascii", "tracers.dat", "XX");
    end
end

[options, prams] = initRigid2D(options, prams);

prams.tau   = tau;
prams.radii = radii;
prams.ar    = ar;

[Xfinal, trajectory] = rigid2D(options, prams, xc, tau);


% profile off;
% profile viewer;
% 
% 
% profsave
