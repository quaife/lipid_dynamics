%eta = 0; eta2 = 0; % variables not used; 
%!!!! BLOCK COMMENTED OUT THROUGH LINE 381 !!!!
%{

rhs2 = o.janusbc(X,tau,center);  
% specify the boundary condition for Janus particles
rhs2 = rhs2(:);

op2 = poten_yukawa(N,rho);
% build double-layer potential matrix
geom.DLP2 = op2.yukawaDLmatrix(geom);

% max GMRES iterations
maxit2 = N*nb; 

% SOLVE SYSTEM USING GMRES
[Xn2,iflag2,res2,I2] = gmres(@(X) o.timeMatVecHalf(X,geom),rhs2,[],...
      o.gmresTol,maxit2);
  
iter2 = I2(2);

% REORGANIZE COLUMN VECTOR INTO MATRIX
% EXTRACT DENSITY FUNCTIONS ON FIBERS AND WALLS
% each column of eta corresponds to the density function of a rigid body
eta2 = zeros(N,nb);
for k = 1:nb
  eta2(:,k) = Xn2((k-1)*N+1:k*N);
end


%%%%%%%%%%%%%%%%%
% plot solution field of screen laplace problem
NX = 500; NY = 500;
% 50 x 50 grid

xmin = o.plotAxis(1);xmax = o.plotAxis(2);
ymin = o.plotAxis(3);ymax = o.plotAxis(4);

xx = linspace(xmin,xmax,NX);
yy = linspace(ymin,ymax,NY);

% test 07/31/2020 
% rr = linspace(1.1, 2, NX);
% tt = linspace(0,2*pi,NY);
% 
% [RR, TT] = meshgrid(rr,tt);
% Uxact = 1/2*(besselk(0,RR/rho)/besselk(0,1/rho)+...
%             besselk(1,RR/rho)/besselk(1,1/rho).*cos(TT));
% Xtest = RR.*cos(TT);
% Ytest = RR.*sin(TT);

[Xtest, Ytest] = meshgrid(xx,yy);
Ztest = Xtest+1i*Ytest;
ind_int = []; %indices of grid points inside particles
% we only want data outside particles
zc = center(1,:) + 1i*center(2,:);
for j = 1:nb
ind_int = [ind_int; find(abs(Ztest-zc(j))<radii(j))];
end

Unum = Xtest; 

[xsou,ysou] = oc.getXY(geom.X);
sa = geom.sa(:,:);
sa = sa(:);
% Jacobian
zt = xsou(:)+1i*ysou(:);
normalx = geom.xt(N+1:2*N,:);
normaly = -geom.xt(1:N,:);


for k = 1:NX
  for j = 1:NY       
    indtmp = (k-1)*NY+j;
    xtest = Xtest(indtmp);
    ytest = Ytest(indtmp);
    ztest = xtest + 1i*ytest;
    % calculate the term  (xsou-xtar)\cdot normal(y)
    const = -((xsou(:)-xtest).*normalx(:) + (ysou(:)-ytest).*normaly(:));
    rr = abs(zt-ztest);
    tmp = 1i*rho*rr;
    kernel = besselh(1,tmp);       
    % boundary integration
    Unum(indtmp)= Xn2'*(-1i/4*tmp.*kernel(:).*const./rr.^2.*sa)*2*pi/N;
  end
end
Unum(ind_int)=0;

% test 07/31/2020 
% figure(4);
% surf(Xtest,Ytest,Uxact-Unum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bug is still there
% compute the body force and torque

force = zeros(2*nb,1);  
torque = zeros(nb,1);   

GradUmat = op2.yukawaGradDLmatrix(geom);
% spy(GradUmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           |           |           | %
%           |x component|     O     | %
%           |           |           | %
%  GradUmat=|-----------|-----------| %
%           |           |           | %
%           |     O     |y component| %
%           |           |           | %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

sa = geom.sa(:,:);
% sa2 = [sa; sa];
% jacobian

% density = [eta2; eta2];

grad = zeros(2*N,nb);
sol = zeros(N,nb);
stress = zeros(2,2,N,nb);
TdotNorm = zeros(2,N,nb);
r0 = zeros(2,N,nb);

% everything is on the boundary
% First, obtain the solution and gradient on the boundary
% by using double layer potential
for k = 1:nb
  sol(:,k) = eta2(:,k)'*(geom.DLP2(:,:,k).*sa(:,k))*2*pi/N;
  grad(1:N,k) = eta2(:,k)'*(GradUmat(1:N,1:N,k).*sa(:,k))*2*pi/N;
  grad(N+1:2*N,k) = eta2(:,k)'*(GradUmat(N+1:2*N,N+1:2*N,k).*sa(:,k))*2*pi/N;
end

figure(3)   % test numerical sol on boundary
% truesol = 0.5*(1+cos(0:2*pi/512:2*pi-2*pi/512));
% plot(rhs2-truesol')
% truegradx = (0.5*sin(1+cos(0:2*pi/512:2*pi-2*pi/512)))'...
%      ./(1+X(N+1:2*N).^2./X(1:N).^2).*(X(N+1:2*N)./X(1:N).^2);
%  size(truegradx)
% plot(truegradx)
plot(grad(1:N))

normal = [geom.xt(N+1:2*N,:); -geom.xt(1:N,:)]; 
  
% Second, calculate the hydrophobic stress, force and, torque.  
for k=1:nb   
    for j=1:N
        gradvec = [grad(j,k); grad(j+N,k)];

% calculuate point stress
        stress(:,:,j,k) = 1/rho*sol(j,k)^2*eye(2) + ...
            2*rho*(0.5*norm(gradvec)^2*eye(2)-gradvec*gradvec');
       
        normali = [normal(j,k); normal(j+N,k)];
        
%         TdotNorm(:,j,k) = (rho*norm(gradvec)^2 + 1/rho*sol(j,k)^2)*normali ...
%              - 2*rho*gradvec*(gradvec'*normali);

        TdotNorm(:,j,k) = stress(:,:,j,k)'*normali;           
        TnuVec = TdotNorm(:,j,k);

% we need the vector from center to the correpsonding boundary
        r0(:,j,k) = [xsou(j,k)-center(1,k); ysou(j,k)-center(2,k)];
        r0vec = r0(:,j,k);
        
        % cross product r0 x TdotN
        cpval = cross ([r0vec;0], [TnuVec;0]);

% force and torque calculated by using integrations
        force(k) = force(k) - TnuVec(1)*sa(j,k)*2*pi/N;
        force(k+nb) = force(k+nb) - TnuVec(2)*sa(j,k)*2*pi/N;
        torque(k) = torque(k) - cpval(3)*sa(j,k)*2*pi/N;
    end
    
end
%}
