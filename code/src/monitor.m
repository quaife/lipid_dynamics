classdef monitor
% Used for doing input/output of data, monitoring runs during execution

properties
verbose    % write data to console
save       % save data to the dat files and log file
dataFile   % name of data file containing fibre centres and orientations
velFile    % name of data file containing velocities
logFile    % name of log file
usePlot    % flag for plotting
plotAxis   % plotting axis
tracer     % flag for tracer
confined   % flag to determine if geometry is confined or not

OUTPUTPATH_DATA % folder in which to save data
OUTPUTPATH_LOG  % folder in which to save logs
OUTPUTPATH_VEL  % folder in which to save velocities
end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(options,prams)
% monitor(options,prams) saves options and parameters needed by the
% class.
% This is the constructor

o.OUTPUTPATH_DATA = '../output/data/';
o.OUTPUTPATH_VEL = '../output/velocity/';
o.OUTPUTPATH_LOG = '../output/logs/';

o.verbose = options.verbose;
% write data to console

o.save     = options.saveData;
o.dataFile = [o.OUTPUTPATH_DATA, options.fileBase, '.bin'];
o.velFile  = [o.OUTPUTPATH_VEL, options.fileBase, '_vel.bin'];
o.logFile  = [o.OUTPUTPATH_LOG, options.fileBase, '.log'];

o.confined = options.confined;
o.usePlot  = options.usePlot;
o.plotAxis = options.plotAxis;
o.tracer   = options.tracer;

% If saving output, reset data and log files and write the number of
% points on each body and number of bodies
if o.save
  o.clearFiles();
  o.welcomeMessage(prams.N,prams.nb);
end

end % constructor: monitor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearFiles(o)
% clearFiles() clears the previous log and data files so that there is 
% nothing from previous runs

fid1 = fopen(o.dataFile,'w');
fclose(fid1);

fid2 = fopen(o.logFile,'w');
fclose(fid2);

fid3 = fopen(o.velFile,'w');
fclose(fid3);

end % clearFiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function welcomeMessage(o,N,nb)
% welcomeMessage() writes specs from the simulation to the
% log file and console

o.writeStars;
message1 = ['RIGID FIBRE SIMULAION ', datestr(now)];
o.writeMessage(message1);
message2 = ['SIMULATING ' num2str(nb) ...
      ' BODIES DISCRETIZED WITH ' num2str(N) ' POINTS'];
o.writeMessage(message2);
o.writeStars;

% This file was being overwritten in initializeFiles
%fid = fopen(o.dataFile,'a');
%% write the number of points per body and number of bodies to the data
%% file
%fwrite(fid,[N;nb],'double');
%% write the solid walls to the data file
%if numel(Xwalls) ~= 0
%  fwrite(fid,Xwalls(:),'double');
%end
%fclose(fid);

end % welcomeMessage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStars(o)
% writeStars writes a message of stars to the console and the log file
% depending on verbose and saveData

stars = '***********************************************************';
o.writeMessage(stars)

end % writeStars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMessage(o,message)
% function writeMessage(message) appends message to o.fileName

% save to log file
if o.save
  fid = fopen(o.logFile,'a');
  fprintf(fid,'%s\n',message);
  fclose(fid);
end
% write to console if verbose==true
if o.verbose
  disp(message)
end

end % writeMessage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initializeFiles(o,geom,walls);

N = geom.N;
nb = geom.nb;
yukawaRHS = reshape(geom.yukawaRHS,geom.N,geom.nb);
oc = curve;
[x,y] = oc.getXY(geom.X);
if o.confined
  Nwall = walls.N;
  nbwall = walls.nb;
  Xwalls = walls.X;
else
  Nwall = 0;
  nbwall = 0;
  Xwalls = [];
end

if o.save;
  fid = fopen(o.dataFile,'w');
  % number of points per body, number of bodies, number of points on
  % solid wall
  fwrite(fid,[N;nb;Nwall;nbwall],'double'); 
  fwrite(fid,Xwalls(:),'double');
  fwrite(fid,yukawaRHS(:),'double'); % yukawa right hand side
  fclose(fid);

  fid = fopen(o.velFile,'w');
  fwrite(fid,[N;nb],'double');
  fclose(fid);  
end

end % initializeFiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(o,time,xc,tau,X)
% writeData(time,xc,tau,X) writes the current time, center, inclination
% angle, and geometry shape to a data file for postprocessing

oc = curve;
[x,y] = oc.getXY(X);
output = [x(:);y(:);xc(:);tau(:);time];

fid = fopen(o.dataFile,'a');
fwrite(fid,output,'double');
fclose(fid);

end % writeData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeVelData(o,time,Up,wp)
% writeVelData(time,xc,tau,X) writes the current time, center,
% inclination angle, and geometry shape to a data file for
% postprocessing

Upx = Up(1,:);
Upy = Up(2,:);

output = [Upx(:); Upy(:); wp(:); time];

fid = fopen(o.velFile,'a');
fwrite(fid,output,'double');
fclose(fid);

end % writeVelData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData(o,geom,Xwalls)

if o.usePlot
  oc = curve;
  figure(1); clf;
  for i = 1:geom.nb
    xx = geom.center(1,i);
    yy = geom.center(2,i);
    th = geom.tau(i);
    ra = geom.radii(i);
    quiver(xx,yy,ra*cos(th),ra*sin(th),'k','linewidth',3);
    hold on
  end

  rhs = reshape(geom.yukawaRHS,geom.N,geom.nb);
  [x,y] = oc.getXY(geom.X);
  for i = 1:geom.nb
    z = rhs(:,i); % color
    h = cline([x(:,i);x(1,i)],[y(:,i);y(1,i)],[z;z(1)]);
    set(h,'linewidth',2)
%    set(h,'EdgeColor','flat');
%    set(h,'FaceColor','none');
    colormap(jet);
  end

  if o.confined
    xx = Xwalls(1:end/2,:); xx = [xx;xx(1,:)];
    yy = Xwalls(end/2+1:end,:); yy = [yy;yy(1,:)];
    % solid walls always have 0 for the HAP model
    zz = zeros(size(Xwalls,1)/2+1,size(Xwalls,2));
%    plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k','linewidth',2)
    h = cline(xx,yy,zz);
    set(h,'linewidth',2);
  end

  axis equal
  axis(o.plotAxis);  
  hold on  
  colorbar
  if o.tracer
    XX = load("../examples/tracers.dat");
    plot(XX(:,1),XX(:,2),'k.');  
  end
  drawnow  
  pause(0.01);
  hold off
end


end % plotData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotField(o,geom,Unum,Xtest,Ytest)

if o.usePlot
  figure(1);
  surf(Xtest,Ytest,Unum,'edgecolor','none');
  colorbar
  caxis([-1 1])
    
  hold on
  for i=1:geom.nb
    xind = [geom.X(1:geom.N,i); geom.X(1,i)];
    yind = [geom.X(geom.N+1:2*geom.N,i); geom.X(geom.N+1,i)];
    xx = geom.center(1,i);
    yy = geom.center(2,i);
    th = geom.tau(i);
    ra = geom.radii(i);
    text(xx,yy,['a_' num2str(i)],'fontsize',16);
    hold on
    plot(xind,yind,'k','linewidth',3);
    hold on
    quiver(xx,yy,ra*cos(th),ra*sin(th),'k','linewidth',3,'MaxHeadSize',1);
    axis equal
    hold on
  end
  view(0,90)
  pause(0.01);
  hold off
end


end % plotField


end % methods


end % classdef
