
classdef monitor
% Used for doing input/output of data, monitoring runs during execution

properties
verbose    % write data to console
save       % save data to the dat files and log file
dataFile   % name of data file containing fibre centres and orientations
logFile    % name of log file
usePlot    % flag for plotting
plotAxis   % plotting axis

OUTPUTPATH_DATA % folder in which to save data
OUTPUTPATH_LOG  % folder in which to save logs

end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = monitor(options,prams)
% monitor(options,prams) saves options and parameters needed by the
% class.
% This is the constructor

o.OUTPUTPATH_DATA = '../output/data/';
o.OUTPUTPATH_LOG = '../output/logs/';

o.verbose = options.verbose;
% write data to console

o.save = options.saveData;
o.dataFile = [o.OUTPUTPATH_DATA, options.fileBase, '.bin'];
o.logFile = [o.OUTPUTPATH_LOG, options.fileBase, '.log'];

o.usePlot = options.usePlot;
o.plotAxis = options.plotAxis;

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

% write the number of points per body and number of bodies to the data
% file
fid = fopen(o.dataFile,'a');
fwrite(fid,[N;nb],'double');
fclose(fid);

end % welcomeMessage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeStars(o)
% writeStars writes a message of stars to the console
% and the log file depending on verbose and saveData

stars = '***********************************************************';
o.writeMessage(stars)

end % writeStars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeMessage(o,message)
% function writeMessage(message) appends message 
% to o.fileName

if o.save
  fid = fopen(o.logFile,'a');
  fprintf(fid,'%s\n',message);
  fclose(fid);
end
% save to log file
if o.verbose
  disp(message)
end
% write to console if verbose==true


end % writeMessage


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

end
% writeData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData(o,geom)
if o.usePlot
    figure(1);
%   fill(geom.X(1:end/2),geom.X(end/2+1:end),'k')
    for i=1:geom.nb
        xind = [geom.X(1:geom.N,i); geom.X(1,i)];
        yind = [geom.X(geom.N+1:2*geom.N,i); geom.X(geom.N+1,i)];
        xx = geom.center(1,i);
        yy = geom.center(2,i);
        th = geom.tau(i);
        ra = geom.radii(i);
        plot(xind,yind,'k');
        hold on
        quiver(xx,yy,ra*cos(th),ra*sin(th),'k','linewidth',3);
        axis equal
        axis(o.plotAxis);  
        hold on
    end
% pause(0.01);
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
%         axis(o.plotAxis);  
        hold on
    end
view(0,90)
pause(0.01);
hold off
end


end % plotData


end % methods


end % classdef
