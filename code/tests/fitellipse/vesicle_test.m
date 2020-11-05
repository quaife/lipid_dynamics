clear
% Please extract the zip file "vel_files.zip" first to have test data

figure('Renderer', 'painters', 'Position', [250 250 1500 600])
N = 16;
nb = 50;

inner = 1:20;
outer = 21:nb;

ar = 2;
radii = 0.5;

tl=[]; ta=[]; majorm=[]; minorm=[];

irate = 40;
T=15000;
dt = 0.05;

k=1;
for i = 0:irate:T
    filename = ['../vel_files/N50_'  num2str(i) '_vel.dat'];
    data=load(filename);
    x = data(:,1); 
    y = data(:,2);
    tau = data(:,3);

    
%%%
theta = (0:N-1)'*2*pi/N;
posx = zeros(N,nb); 
posy = zeros(N,nb); 

tt = 0:irate*dt:T*dt;

for j = 1:nb
  refX = kron([cos(tau(j)) -sin(tau(j)); +sin(tau(j)) cos(tau(j))],eye(N)) * ... 
      [ar*cos(theta);sin(theta)]*radii;  
  % shape of particle

  % rotated circle
  posx(:,j) = refX(1:N) + x(j);
  posy(:,j) = refX(N+1:end) + y(j);
  % shift to the correct center
end    
    
    
%%%
    
% figure(1);   
subplot(1,3,1)
    hold off
    [zo, ao, bo, alphao] = fitellipse([x((outer))';y(outer)'], 'linear', 'constraint', 'trace');
    [zi, ai, bi, alphai] = fitellipse([x(inner)';y(inner)'], 'linear', 'constraint', 'trace');
%%%    
npts = 100;
t = linspace(0, 2*pi, npts);

% Rotation matrix
Qo = [cos(alphao), -sin(alphao); sin(alphao) cos(alphao)];
Qi = [cos(alphai), -sin(alphai); sin(alphai) cos(alphai)];
% Ellipse points
Xo = Qo * [ao * cos(t); bo * sin(t)] + repmat(zo, 1, npts);  
Xi = Qi * [ai * cos(t); bi * sin(t)] + repmat(zi, 1, npts);  
%%%
% midplane    
    [zm, am, bm, alpham] = fitellipse([Xo Xi], 'linear', 'constraint', 'trace');
    Qm = [cos(alpham), -sin(alpham); sin(alpham) cos(alpham)];
    Xm = Qm * [am * cos(t); bm * sin(t)] + repmat(zm, 1, npts);  
    hm=plotellipse(zm, am, bm, alpham);
    set(hm,'linewidth',3)
    hold on    


    ho=plotellipse(zo, ao, bo, alphao);
    set(ho,'linewidth',3)
    hold on

    hi=plotellipse(zi, ai, bi, alphai);
    set(hi,'linewidth',3)
    hold on


  for j = 1:nb
    z = zeros(N,1); % color
    h = cline([posx(:,j);posx(1,j)],[posy(:,j);posy(1,j)],[z;z(1)]);
    set(h,'linewidth',2)
%     colormap(jet);
  end
    
    drawnow
    F = getframe;
%     imwrite(F.cdata, ['File' num2str(i) '.png']);
    
    
    sizeXm = size(Xm);
    nxm = sizeXm(2);
    XX = [Xm(1,:) Xm(2,:)]';   
    [x,y] = getXY(XX);
    [Dx,Dy] = getDXY(XX);
    tl(k) = sum(sqrt(Dx.^2+Dy.^2))*2*pi/nxm;
    ta(k) = sum(x.*Dy-y.*Dx)*pi/nxm;
    
    
    majorm(k) = am;
    minorm(k) = bm;

%%%    
subplot(1,3,2)
plot(tt(1:k), ta(1:k), 'r', 'linewidth',3)
xlabel('time t');
ylabel('total area');
set(gca,'fontsize',16);
axis([0 T*dt ta(1)-40 ta(1)])

%%%
subplot(1,3,3)
plot(tt(1:k), tl(1:k), 'r', 'linewidth',3)
xlabel('time t');
ylabel('total arc length');
set(gca,'fontsize',16); 
axis([0 T*dt tl(1)-10 tl(1)])

area=pi*majorm'.*minorm';

k=k+1;
end
