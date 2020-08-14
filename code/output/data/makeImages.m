addpath ../../src
set(0,'DefaultAxesFontSize',22)
options.savefig = false;

file = 'shear.bin';

irate = 1; % controls the speed of the visualization
ax = 5*[-1 1 -1 1];

[posx,posy,xc,tau,time] = loadFile(file);
% load positions, centers, inclination angles, and times

oc = curve;

figure(1); clf
for k = 1:numel(time);
  xx = posx(:,:,k);
  yy = posy(:,:,k);
  clf; hold on;
  fill(xx,yy,'k')
  axis equal
  axis(ax)
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

