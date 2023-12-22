function EllipTar = EllipsesTargets(i, Nin, Nout,velind)
% This script generates some curves outside and inside the closed bilayers
% using fft. 
% Input:
% i: step number
% N = Nin+Nout
% velind: number of units from midplane
% Output: 
% EllipTar: matrix [X1 Y1 X2 Y2 ...] where X1 and Y1 are column vectors.

addpath ../src
addpath ../output/

% load options and parameters
load("../output/data/frames/options.mat");

EllipTar = [];

oc = curve;

if nargin~=4
% scalar multiples: "-" interior; "+" exterior
velind = [-2.0 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.35 -1.3 ... 
                        -0.05 0 0.05 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6];
end

fileName = sprintf("../output/data/frames/N%d_%f_%d.dat", prams.nb, options.shearRate, i);        
data = load(fileName); 
x = data(:,1)'; y = data(:,2)'; tau   = data(:,3)';

[X, ~, ~, ~] = MidArcLen(x,y,tau,Nout,Nin);

[~,tangent,~] = oc.diffProp(X);
[n1,n2] = oc.getXYperp(tangent);

% t1 = tangent(1:end/2);
% t2 = tangent(end/2+1:end);

for k = 1:length(velind)
% generate curves in normal directions
    Xn = X + velind(k)*[n1 ; n2]; 
    figure(111); hold on
    plot(Xn(1:end/2),Xn(end/2+1:end))
    EllipTar = [EllipTar Xn(1:end/2) Xn(end/2+1:end)];    
end

% for saving data
% fileName = sprintf("N%d_%f.midtar",nb, shearRate);
% save("-ascii", fileName, "EllipTar");


plot(x,y,'o');
plot(X(1:64),X(65:128),'o-');

% example 
% data = load('test_vesicle.dat');
% x = data(:,1); y = data(:,2); tau = data(:,3);
% Nin = 26; Nout = 32;
% EllipTar = EllipsesTargets(x, y, tau, Nin, Nout);

end