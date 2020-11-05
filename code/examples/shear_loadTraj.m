clear
clf
format long e
format compact

addpath ../src
addpath ../output/data
addpath ../output/velocity

file = 'shear.bin';

[yukawaRHS,posx,posy,xc,tau,time] = loadFile(file);

prams.N = size(posx,1);
prams.nb = size(posx,2);
ntime = size(posx,3);

dt = 0.05;
prams.m = ntime ; % number of time steps
prams.T = prams.m*dt; % time horizon

prams.order = 2;

prams.rho = 1.0; % screen length
prams.gam = 1.0;

prams.RepulLength = 2*0.3; % repulsion length
prams.RepulStrength = 2*0.5; % repulsion strength

options.farField  = 'shear';
options.shearRate = 0.005;
options.saveData  = true;
options.fileBase  = 'shear';
options.append    = false;
options.inear     = true;
options.usePreco  = false;
options.verbose   = true;
options.timeOrder = 1;
options.gmresTol  = 1e-5;
options.usePlot   = false; %true;
options.plotAxis  = 15*[-1 1 -1 1];

ratio = 1;
radii = 1/ratio*ones(1,prams.nb);
ar = ratio*ones(1,prams.nb);

prams.tau   = tau;
prams.radii = radii;
prams.ar    = ar;

om = monitor(options,prams);
tt = tstep(options,prams);

kstart = 0;
tstart = kstart*dt;
irate = 10;

if (irate ~= 1)
    kend = floor(ntime/irate)+1;
else
    kend = ntime;
end

Up = 0*xc;
wp = 0*tau;

for k = 1:irate:ntime
% main script
geom = capsules(prams,xc(:,:,k),tau(1,:,k));
[Up(:,:,k), wp(1,:,k),~,~,etaY,etaS] = tt.timeStep(geom,geom.X,geom.X);

% write current time to console
message = ['Completed time ', num2str(tstart+time(k),'%4.2e'), ...
      ' of time horizon ' , num2str(prams.T,'%4.2e')];
om.writeMessage(message);

end

for k = 1 : kend
ind = (k-1)*irate+1;
data = [xc(1,:,ind);xc(2,:,ind);tau(1,:,ind);Up(1,:,ind);Up(2,:,ind);wp(1,:,ind)];
filename = ['../output/velocity/', 'N' num2str(prams.nb) '_' num2str((k-1)*irate+kstart),'_vel.dat'];
fid = fopen(filename,'w');
fprintf(fid,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',data);
fclose(fid); 
end