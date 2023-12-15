function slipv = SlipVel(Nout,Nin,StepI,irate)
% This function uses all configuration data files to calculate the slip
% velocity of the closed tank-treading vesicle.
% Inputs: 
% Nout: number of particles in the outer leaflet
% Nin: number of particles in the inner leaflet
% StepI: Interval [index_1 index_2]
% irate: check data files every irate steps
%
% Output:
% slipv:  array of slip velocities

addpath ../src
oc = curve;

% load options and parameters
load("../output/data/frames/options.mat");

dt = prams.T/pram.m;

Tstart = StepI(1);
Tend = StepI(end);

% corresponding time t 
t = (Tstart:irate:Tend)*dt;

k=1;
slipv = [];
for i = Tstart:irate:Tend
    filename = sprintf("../output/data/frames/N%d_%f_%d.xcvel",...
                                        prams.nb, options.shearRate, i);    
    data=load(filename);
    x = data(:,1); 
    y = data(:,2);
    % 
    u = data(:,3);
    v = data(:,4);

    % x = [x_in; x_out];  y = [y_in; y_out];
    x_in  = x(1:Nin); y_in = y(1:Nin);
    x_out = x(Nin+1:end); y_out = y(Nin+1:end);
    Xi = [x_in;y_in];
    Xo = [x_out;y_out];

    % ****  notice **** Nin/2 has to be even 
    [jacI,tanI,kaI] = oc.diffProp(Xi);  
    [jacO,tanO,kaO] = oc.diffProp(Xo);
    
    % Normal vectors
    [ni1,ni2] = oc.getXYperp(tanI);
    [no1,no2] = oc.getXYperp(tanO);

    % Tangent vectors
    ti1 = -tanI(1:end/2);
    ti2 = -tanI(end/2+1:end);
    to1 = -tanO(1:end/2);
    to2 = -tanO(end/2+1:end);

    innerV = zeros(Nin,1);
    outerV = zeros(Nout,1);

    % tangential velocities
    for j = 1:Nin
        innerV(j) = dot([u(j) v(j)],[ti1(j) ti2(j)]);
    end

    for j = 1:Nout
        outerV(j) = dot([u(j+Nin) v(j+Nin)],[to1(j) to2(j)]);
    end
        
	slipv(k,:) = [ t(k)   mean(outerV)-mean(innerV)];
    k=k+1;
    end

    % save data
    save("SlipV.dat",'slipv','-ascii');
end