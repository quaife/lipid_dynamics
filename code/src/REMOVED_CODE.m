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


% Code removed on Jan 23, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F1,F2,Tq] = evalForcesQBXOld(o, Nb, N, x1, x2, nu1, nu2, dS, rho, h)

% Uses * Tpq + Tqp identity, 
%      * Jpq jump value to evaluate on curve itself 
%      * QBX expansion to derive uq and uq_n on curve p, p ~= q
%      * FFT for tangential derivatives     

F1 = zeros(Nb,1);
F2 = zeros(Nb,1);
Tq = zeros(Nb,1);

rad   = 0.3;
m_max = 6;
tol   = rad;
B_m   = zeros(N,2*m_max + 1);

%errors = zeros(Nb, 3);

IK = 1i*[(0:N/2-1) 0 (-N/2+1:-1)]'; % Fourier modes
dT = 2*pi/N;     % to convert dS to sa

for p = 1:Nb
  for q = [1:p-1, p+1:Nb]

    x1p = x1(:,p);
    x2p = x2(:,p);        
    dSp = dS(:,p);
    nu1p = nu1(:,p);
    nu2p = nu2(:,p);
    tau1p = -nu2p;
    tau2p = +nu1p;

    % setting argin Nb = 1 tricks evaluators to using only one geometry
    % column
    hp = h(:, p);

    % decide which strategy to use : expansion or standard evaluation 
              
    D = pdist2([x1(:,p) x2(:,p)],[x1(:,q), x2(:,q)]);
    D = min(D, [], "all");

    if D < tol
      uq = zeros(N,1);
      uq_x1 = zeros(N,1);
      uq_x2 = zeros(N,1);

      c1 = x1p - rad*nu1p;
      c2 = x2p - rad*nu2p;

      for m = -m_max:m_max
        B_m(:,m+m_max+1) = o.QBX_coeff(c1,c2,m,x1(:,q),x2(:,q), ...
              nu1(:,q),nu2(:,q),dS(:,q),rho,h(:,q));
      end                

      [uq,uq_x1,uq_x2] = o.QBX_exp(x1p,x2p,c1,c2,B_m,m_max,rho); 

    else
      [uq,uq_x1,uq_x2] = o.evalDL(x1p,x2p,1,N,x1(:,q),x2(:,q), ...
            nu1(:,q),nu2(:,q),dS(:,q),rho, h(:,q));
    end % if D < tol                                   

    hp_t = o.tanDeriv(hp,dSp/dT,IK);
    uq_t = o.tanDeriv(uq,dSp/dT,IK);
    uq_n = uq_x1.*nu1p + uq_x2.*nu2p;

    Jpq1 = 2.0/rho*hp.*uq.*nu1p + 2.0*rho*hp_t.*uq_t.*nu1p - ...
        2.0*rho*hp_t.*uq_n.*tau1p;
    Jpq2 = 2.0/rho*hp.*uq.*nu2p + 2.0*rho*hp_t.*uq_t.*nu2p - ...
        2.0*rho*hp_t.*uq_n.*tau2p;        

    F1(p) = F1(p) + sum(Jpq1.*dSp);
    F2(p) = F2(p) + sum(Jpq2.*dSp);
    Tq(p) = Tq(p) + sum(( + (x1p - x1pc).*Jpq2 - (x2p - x2pc).*Jpq1 ) .* dSp);

  end    

end     
    
end % evalForcesQBXOld


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F1, F2, Tq] = evalForcesTaylorOld(o,Nb,N,x1,x2,nu1,nu2,dS,rho,h)
% NOTE: THIS CODE IS CURRENTLY NOT BEING CALLED, SO NO NEED TO CLEAN IT
% UP. IDEA IS TO USE TAYLOR EXPANSIONS RATHER THAN QBX

% Uses * Tpq + Tqp identity, 
%      * Jpq jump value to evaluate on curve itself 
%      * Taylor expansion to derive uq and uq_n on curve p, p ~= q
%      * FFT for tangential derivatives     

F1 = zeros(Nb,1);
F2 = zeros(Nb,1);
Tq = zeros(Nb,1);

tol   = 1.0;

p_max = 8;
r_max = 0.2;
r_m   = linspace(0, r_max, p_max+1);
r_m   = r_m(2:end)';
r_0   = r_m(p_max);
pw    = (0:p_max-1);
AA    = (r_m - r_0).^pw;
BB    = (0 - r_0).^pw;

%errors = zeros(Nb, 3);

IK = 1i*[(0:N/2-1) 0 (-N/2+1:-1)]'; % Fourier modes    
dT = 2*pi/N;     % to convert dS to sa

for p = 1:Nb
  for q = [1:p-1, p+1:Nb]

    x1p   = x1(:,p);
    x2p   = x2(:,p);        
    dSp   = dS(:,p);
    nu1p  = nu1(:,p);
    nu2p  = nu2(:,p);
    tau1p = -nu2p;
    tau2p =  nu1p;
    hp    = h(:, p);

    % decide which strategy to use : expansion or standard evaluation 
              
    D = pdist2([x1(:,p) x2(:,p)],[x1(:,q), x2(:,q)]);
    D = min(D, [], "all");
               
    if D < tol
      [Uq,Uq_x1,Uq_x2] = o.evalDL(x1p - nu1p*r_m',x2p - nu2p*r_m', ...
         1,N,x1(:,q),x2(:,q),nu1(:,q),nu2(:,q),dS(:,q),rho,h(:,q));

      coef    = AA\Uq.';
      coef_x1 = AA\Uq_x1.';
      coef_x2 = AA\Uq_x2.';                

      uq      = (BB*coef)';
      uq_x1   = (BB*coef_x1)';
      uq_x2   = (BB*coef_x2)';

    else

      [uq,uq_x1,uq_x2] = o.evalDL(x1p,x2p,1,N,x1(:,q),x2(:,q),...
          nu1(:,q),nu2(:,q),dS(:,q),rho,h(:,q));

    end % if D < tol                                               

    hp_t = o.tanDeriv(hp,dSp/dT,IK);
    uq_t = o.tanDeriv(uq,dSp/dT,IK);
    uq_n = uq_x1.*nu1p + uq_x2.*nu2p;

    Jpq1 = 2.0/rho*hp.*uq.*nu1p + 2.0*rho*hp_t.*uq_t.*nu1p - ...
        2.0*rho*hp_t.*uq_n.*tau1p;
    Jpq2 = 2.0/rho*hp.*uq.*nu2p + 2.0*rho*hp_t.*uq_t.*nu2p - ...
        2.0*rho*hp_t.*uq_n.*tau2p;        

    F1(p) = F1(p) + sum(Jpq1.*dSp);

    F2(p) = F2(p) + sum(Jpq2.*dSp);

    Tq(p) = Tq(p) + sum(( + (x1p - x1pc).*Jpq2 - (x2p - x2pc).*Jpq1 ) .* dSp);

  end    
end     

end % evalForcesTaylorOld


function [F1, F2, Tq] = evalForcesTaylor(o,geom,eta)
% NOTE: THIS CODE IS CURRENTLY NOT BEING CALLED, SO NO NEED TO CLEAN IT
% UP. IDEA IS TO USE TAYLOR EXPANSIONS RATHER THAN QBX

% Uses * Tpq + Tqp identity, 
%      * Jpq jump value to evaluate on curve itself 
%      * Taylor expansion to derive uq and uq_n on curve p, p ~= q
%      * FFT for tangential derivatives     

oc = curve;
N = geom.N;
[x1,x2] = oc.getXY(geom.X);
[tau1,tau2] = oc.getXY(geom.xt);
nu1 = +tau2;
nu2 = -tau1;
dS = geom.sa*2*pi/N;

F1 = zeros(geom.nb,1);
F2 = zeros(geom.nb,1);
Tq = zeros(geom.nb,1);

tol   = 1.0;

p_max = 8;
r_max = 0.2;
r_m   = linspace(0, r_max, p_max+1);
r_m   = r_m(2:end)';
r_0   = r_m(p_max);
pw    = (0:p_max-1);
AA    = (r_m - r_0).^pw;
BB    = (0 - r_0).^pw;

IK = oc.modes(geom.N,1); % Fourier modes

for p = 1:geom.nb
  for q = [1:p-1, p+1:geom.nb]
    x1p   = x1(:,p);
    x2p   = x2(:,p);        
    dSp   = dS(:,p);
    nu1p  = nu1(:,p);
    nu2p  = nu2(:,p);
    tau1p = -nu2p;
    tau2p =  nu1p;
    etap    = eta(:, p);

    % decide which strategy to use : expansion or standard evaluation 
    D = min(pdist2([x1(:,p) x2(:,p)],[x1(:,q), x2(:,q)]),[],"all");
               
    if D < tol
      [Uq,Uq_x1,Uq_x2] = o.evalDL(x1p - nu1p*r_m',x2p - nu2p*r_m', ...
         1,N,x1(:,q),x2(:,q),nu1(:,q),nu2(:,q),dS(:,q),geom.rho,eta(:,q));

      coef    = AA\Uq.';
      coef_x1 = AA\Uq_x1.';
      coef_x2 = AA\Uq_x2.';                

      uq      = (BB*coef)';
      uq_x1   = (BB*coef_x1)';
      uq_x2   = (BB*coef_x2)';

    else

      [uq,uq_x1,uq_x2] = o.evalDL(x1p,x2p,1,N,x1(:,q),x2(:,q),...
          nu1(:,q),nu2(:,q),dS(:,q),geom.rho,eta(:,q));

    end % if D < tol                                               

%    etap_t = o.tanDeriv(etap,geom.sa(:,p),IK);
%    uq_t = o.tanDeriv(uq,geom.sa(:,p),IK);

    % compute tangent derivatives of etap and uq.
    etap_t = oc.diffFT(etap,IK)./geom.sa(:,p);
    uq_t = oc.diffFT(uq,IK)./geom.sa(:,p);
    uq_n = uq_x1.*nu1p + uq_x2.*nu2p;

    Jpq1 = 2.0/geom.rho*etap.*uq.*nu1p + ...
        2.0*geom.rho*etap_t.*uq_t.*nu1p - ...
        2.0*geom.rho*etap_t.*uq_n.*tau1p;
    Jpq2 = 2.0/geom.rho*etap.*uq.*nu2p + ...
        2.0*geom.rho*etap_t.*uq_t.*nu2p - ...
        2.0*geom.rho*etap_t.*uq_n.*tau2p;        

    F1(p) = F1(p) + sum(Jpq1.*dSp);

    F2(p) = F2(p) + sum(Jpq2.*dSp);

    Tq(p) = Tq(p) + sum( ( x1p.*Jpq2 - x2p.*Jpq1 ) .* dSp );        
             
  end    
end     

end % evalForcesTaylor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [force,torque] = bodyForceTorque(o)

% will direct the calculated force, torque here
force = zeros(2*o.nb,1);  % original [1;2]
torque = zeros(o.nb,1);   % original -10


end % bodyForceTorque


