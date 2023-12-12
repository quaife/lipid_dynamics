clear
% clf
addpath ../src
addpath ../output/

format long e
format compact

% load options and parameters
load("../output/data/frames/options.mat");

%
A = options.plotAxis;

% grid for HAP Action field
xtmp = linspace(A(1),A(2),200);
ytmp = linspace(A(3),A(4),200);
[X,Y] = meshgrid(xtmp,ytmp);
X = X(:)'; Y = Y(:)'; Xtar = [X;Y];
Ntar = length(X);


%%% loop for loading data and post-processing
stepJump = 1;
for i = 0:stepJump:1 %prams.m

% load configurations, forces, torques, and densities
fileName1 = sprintf("../output/data/frames/N%d_%f_%d.dat", prams.nb, options.shearRate, i);    
fileName2 = sprintf("../output/data/frames/N%d_%f_%d.mat", prams.nb, options.shearRate, i);    
data = load(fileName1); load(fileName2);
x = data(:,1)'; y = data(:,2)'; tau   = data(:,3)';
xc = [x;y];

geom0 = capsules(prams,xc,tau);
op = poten(geom0.N,geom0.rho);

% Calculate fluid stress, pressure and velocity at target points Xtar in
% the computational domain
[stress, pressure, velocity] = fluidstress(geom0,etaS0,force,torque,Xtar);


% SPV = [X' Y' stress(1:Ntar) stress(Ntar+1:2*Ntar) stress(2*Ntar+1:end) ...
%        pressure velocity(1:Ntar) velocity(Ntar+1:end)];

% output pressure and velocity
pressure = reshape(pressure,length(xtmp),length(ytmp));
fileName = sprintf("../output/data/frames/N%d_%f_%d.pres", geom0.nb, options.shearRate, i);
save("-ascii", fileName, "pressure");
  
fileName = sprintf("../output/data/frames/N%d_%f_%d.vel", geom0.nb, options.shearRate, i);
save("-ascii", fileName, "velocity");



% Calculate HAP solution at target points Xtar
[yukawaDLP,yukawaDLPtar] = op.exactYukawaDL(geom0,etaY0,Xtar,(1:geom0.nb));

% output HAP solution
yukawaDLPtar = reshape(yukawaDLPtar,200,200);
fileName = sprintf("../output/data/frames/N%d_%f_%d.eta", geom0.nb, options.shearRate, i);
save("-ascii", fileName, "yukawaDLPtar");  
end