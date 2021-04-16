if ~EXTERNAL_OPTIONS

clear
% clf
format long e
format compact

prams.N = 16;      % number of point per body
prams.nb = 3;      % number of bodies
prams.sstep = 0;   % starting step

radii = 1.0*ones(1,prams.nb);
ar    = 1.0*ones(1,prams.nb);

prams.rho = 4.0;      % screen length
prams.gam = 4.0;      % molecular tension

prams.RepulLength   = 0.5; % repulsion length
prams.RepulStrength = 4.0; % repulsion strength

prams.bcShift       = 0.0*ones(1,prams.nb); % shift constant for yukawaRHS
prams.bcType        = 'cosine'; % options: 'cosine'; 'vonMises'

dt      = 0.3;
prams.m = 500;        % number of time steps
prams.T = prams.m*dt; % time horizon
prams.order = 2;      % time-stepping order 

end


options.farField  = 'shear'; 
% options: 'shear'; 'extensional'; 'parabolic'; 'taylorgreen'
options.shearRate = 0.00;

options.saveData  = false;
options.fileBase  = options.farField;

options.append    = false;
options.inear     = true;
options.usePreco  = false;
options.verbose   = true;
options.timeOrder = 2;
options.gmresTol  = 1e-10;

options.usePlot   = true;
options.tracer    = true;
options.plotAxis  = [-20 20 -20 20];

% 1. Please input the number of bodies and have the corresponding initial
% configuration file ready.
% 2. If a desired starting step is given, please modify prams.sstep and 
% prepare the corresponding configuration file w/ or w/o the tracer file.

if prams.sstep == 0
    data = load(['N' num2str(prams.nb) '_' num2str(prams.sstep) '.dat']);
    fileName = sprintf("../output/data/frames/N%d_%f_%d.dat", ...
                   prams.nb, options.shearRate, prams.sstep);
    save("-ascii", fileName, "data");
else
    fileName = sprintf("../output/data/frames/N%d_%f_%d.dat", ...
                   prams.nb, options.shearRate, prams.sstep);    
    data = load(fileName);
end

% initial centers and angles 
x = data(:,1)';
y = data(:,2)';
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
