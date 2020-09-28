addpath ../../src
set(0,'DefaultAxesFontSize',22)
options.savefig = false;

file = 'shear.bin';

irate = 1; % controls the speed of the visualization
ax = 5*[-1 1 -1 1];

[yukawaRHS,posx,posy,xc,tau,time] = loadFile(file);
% load right hand side for Yukawa solve, positions, centers, inclination
% angles, and times
N = size(posx,1);
nb = size(posx,2);
ntime = size(posx,3);
yukawaRHS = reshape(yukawaRHS,N,nb);

oc = curve;

figure(1); clf
for k = 1:ntime;
  xx = posx(:,:,k);
  yy = posy(:,:,k);
  figure(1); clf;
  hold on
  for i = 1:nb
    xx = xc(1,i,k);
    yy = xc(2,i,k);
    th = tau(1,i,k);
    quiver(xx,yy,0.5*cos(th),0.5*sin(th),'k','linewidth',3);
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
  titleStr = ['t = ' num2str(time(k),'%4.2e')];
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
end

