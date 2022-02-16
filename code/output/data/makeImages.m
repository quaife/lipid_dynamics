addpath ../../src
addpath ../velocity
set(0,'DefaultAxesFontSize',22)
options.savefig = false;
options.saveVelData = true;

file = 'shearBryan.bin';
velfile = 'shear_vel.bin';

irate = 1; % controls the speed of the visualization
ax = 8*[-1 1 -1 1];
% ax = [-1 6 -1 4];

[yukawaRHS,posx,posy,xc,tau] = loadFile(file);
% load right hand side for Yukawa solve, positions, centers, inclination
% angles, and times
[velx,vely,torq] = loadVelFile(velfile);
% load x and y components of velocity profile, torque
N = size(posx,1);
nb = size(posx,2);
ntime = size(posx,3);
yukawaRHS = reshape(yukawaRHS,N,nb);

oc = curve;

figure(1); clf
count = 0;
for k = 1:irate:ntime;
  xx = posx(:,:,k);
  yy = posy(:,:,k);
  figure(1); clf;
  hold on
  for i = 1:nb
    xx = xc(1,i,k);
    yy = xc(2,i,k);
    th = tau(1,i,k);
    quiver(xx,yy,1*cos(th),1*sin(th),'k','linewidth',3);
  end

  x = posx(:,:,k);
  y = posy(:,:,k);
  for i = 1:nb
    z = yukawaRHS(:,i); % color
    h = cline([x(:,i);x(1,i)],[y(:,i);y(1,i)],[z;z(1)]);
    set(h,'linewidth',2)
    colormap(jet);
  end
  axis equal
  axis(ax);
  hold on
  colorbar

  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'xcolor','white')
  set(gca,'ycolor','white')
%  titleStr = ['t = ' num2str(time(k),'%4.2e')];
  title(titleStr)
  pause(0.01)
  if options.savefig
    filename = ['./frames/image', sprintf('%04d',count),'.pdf'];
    count = count+1;
    figure(1);
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

    print(gcf,'-dpdf','-r300',filename);
  end 
  pause(0.01)
  
  if options.saveVelData
    % centers, tau, velx, vely, torq
    data = [xc(1,:,k);xc(2,:,k);tau(1,:,k);velx(1,:,k);vely(1,:,k);torq(1,:,k);];
    count = count+1;
    filename = ['../velocity/', 'N' num2str(nb) '_' num2str(count),'_vel.dat'];
    fid = fopen(filename,'w');
    fprintf(fid,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',data);
    fclose(fid);      
  end
  pause(0.01) 
end

