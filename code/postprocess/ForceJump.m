function dF = ForceJump(EllipTar,i,velind)
% This function uses the generated target points in ellipses done by
% EllipseTargets.m to calculate the force jump at the midplane.
% Input: 
% i: step number

if nargin~=3
    velind = [-2.0 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.35 -1.3 ...
              0 0.1 0.2 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6];  
end

addpath ../src
addpath ../output/

% load options and parameters
load("../output/data/frames/options.mat");

% load force, torque, and densities    
fileName1 = sprintf("../output/data/frames/N%d_%f_%d.dat", prams.nb, options.shearRate, i);    
fileName2 = sprintf("../output/data/frames/N%d_%f_%d.mat", prams.nb, options.shearRate, i);    
data = load(fileName1); 
load(fileName2);
x = data(:,1)'; y = data(:,2)'; tau   = data(:,3)';
xc = [x;y];

geom0 = capsules(prams,xc,tau);

xtar = EllipTar(:,1:2:end);
ytar = EllipTar(:,2:2:end);
xtar = xtar(:)';
ytar = ytar(:)';
Xtar = [xtar;ytar];
Ntar = length(xtar);

[stress, pressure, velocity] = fluidstress(geom0,etaS0,force,torque,Xtar);


SPV = [xtar' ytar' stress(1:Ntar) stress(Ntar+1:2*Ntar) stress(2*Ntar+1:end) ...
       pressure velocity(1:Ntar) velocity(Ntar+1:end)];

dF = CompForceOnly(SPV,velind);

end


function dF = CompForceOnly(SPV,velind)
    oc = curve;   
    
    x = SPV(:,1);
    y = SPV(:,2);
    sxx = SPV(:,3);
    sxy = SPV(:,4);
    syy = SPV(:,5);
    P = SPV(:,6);
    vx = SPV(:,7);
    vy = SPV(:,8);
    
    x = reshape(x,64,21);
    y = reshape(y,64,21);
    sxx = reshape(sxx,64,21);
    sxy = reshape(sxy,64,21);   
    syy = reshape(syy,64,21);
    P = reshape(P,64,21);
    vx = reshape(vx,64,21);
    vy = reshape(vy,64,21);

    mid = 11;
    [jacobian,tangent,ka] = oc.diffProp([x(:,mid);y(:,mid)]);

    cx = mean(x(:,mid));
    cy = mean(y(:,mid));
    
    [n1,n2] = oc.getXYperp(tangent);

    t1 = -tangent(1:end/2);   
    t2 = -tangent(end/2+1:end);
    
    
 
    FoutData = [];
    
    for out = 21:-1:13
    dSx =  sxx(:,out).*n1 + sxy(:,out).*n2;
    dSy =  sxy(:,out).*n1 + syy(:,out).*n2;
    dSo = dSx.*t1 + dSy.*t2;
    
    [jaout,tangent,ka] = oc.diffProp([x(:,out);y(:,out)]);
	FoutData = [FoutData; sum(dSo.*jaout)*2*pi/64];
    end
    
    FinData = [];
    for in = 1:5
    dSx =  sxx(:,in).*n1 + sxy(:,in).*n2;
    dSy =  sxy(:,in).*n1 + syy(:,in).*n2;
    dSi = dSx.*t1 + dSy.*t2;    
    
%   
    [jain,tangent,ka] = oc.diffProp([x(:,in);y(:,in)]);   
 	FinData = [FinData; sum(dSi.*jain)*2*pi/64];
    end
    
    
    d1 = velind(end:-1:13);
    d2 = velind(1:5);
    figure(2);
    hold on
    plot(d1,FoutData)
    p = polyfit(d1,FoutData,2);
    f = polyval(p,0:0.1:5);
    plot(d1,FoutData,'*',0:0.1:5,f)
    Fout = f(1);
    plot(d2,FinData)
    p = polyfit(d2,FinData,2);
    f = polyval(p,-2.0:0.1:0);
    plot(d2,FinData,'*',-2.0:0.1:0,f) 
    Fin= f(end);

    dF = Fout - Fin; 

end