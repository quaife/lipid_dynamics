% Calculate force and torque for multiple bodies

% Here we translate Bryan's geometry class into extant variables 

N = prams.N;                 % number of points per componenet
Nb = prams.nb;               % number of rigid bodies
x1 = geom.X(1:N,:);         % grid points on curves 
x2 = geom.X(N+1:2*N,:);            
pc = geom.center;           % center of each rigid body
orn = geom.tau;             % orientation of each rigid body
a = geom.radii;             % radius of each rigid body
rho = geom.rho;             % screen length of particles
xt = geom.xt;          % tangent unit vector
tau1 = xt(1:N,:);      
tau2 = xt(N+1:2*N,:);
nu1  = - tau2;         % outward normal : (nu, tau) is right-handed
nu2  = + tau1;

dS = geom.sa*2*pi/N;   % Jacobian
cur = geom.cur;        % curvature


% Now we design some nearby curves which are dilated versions of the primal curves.
% This leaves the unit tangent and normal unchanged. The Jacobian scales directly. 
%
% At the same time we assign the Dirichlet data uD.

x1_nbr  = 0*x1;
x2_nbr  = 0*x2;
nu1_nbr = nu1;
nu2_nbr = nu2;
dS_nbr  = dS;

uD      = 0*x1;
for p = 1:Nb

    fac = (a(p) + 10*max(dS(:,p)))/a(p);
    x1_nbr(:,p) = fac*(x1(:,p) - pc(1,p)) + pc(1,p);
    x2_nbr(:,p) = fac*(x2(:,p) - pc(2,p)) + pc(2,p);        
    dS_nbr(:,p) = fac*dS(:,p);    

    tmp     = atan2(x2(:,p) - pc(2,p), x1(:,p) - pc(1,p)); 
    uD(:,p) = 0.5*(1 + cos( tmp - orn(p) ));
    uD(:,p) = ones(size(tmp));

end 


% visualizaion 
%clf
[X1, X2] = meshgrid(linspace(min(x1,[],'all')-2, max(x1,[],'all')+2, 201), linspace(min(x2,[],'all')-2, max(x2,[],'all')+2, 201));
%{
plot(x1, x2, 'k');
hold on
plot(x1_nbr, x2_nbr, 'm');
plot3(x1,x2,uD,'r','linewidth',2)
hold off
axis equal
view(3)
%}

%///////////////////////////////////////////////
%Setting up linear system

K = zeros(N*Nb,N*Nb);

% vectorized version of the above 
for q = 1:Nb
    for j = 1:N

        r1 = x1(j,q) - x1(:);
        r2 = x2(j,q) - x2(:);
        r  = sqrt( r1.^2 + r2.^2 );

        Lnu = ( r1.*nu1(j,q) + r2.*nu2(j,q) )./r.^2;

        J = (q-1)*N + j;
        K(:, J) = 1/(2*pi)*(r/rho).*besselk(1, r/rho).*Lnu*dS(j,q);

        K(J, J) = -cur(j,q)*dS(j,q)/(4*pi); %removable singularity 
        
    end
end
disp('here')
pause

RHS = reshape(uD,N*Nb,1);
h   = (1/2*eye(N*Nb) + K)\RHS;
h   = reshape(h, N, Nb);

% next two commented lines needed only for solution visualization in bulk
 Dh = evalDL(X1, X2, Nb, N, x1, x2, nu1, nu2, dS, rho, h);
% clf
% surf(X1,X2, Dh,'edgecolor','none')
% pause

disp('The system of forces and the net force [F1 F2 Tq; sum(F1) sum(F2) sum(Tq)]')
[F1, F2, Tq] = evalForces(Nb, N, x1_nbr, x2_nbr, nu1_nbr, nu2_nbr, dS_nbr, x1, x2, nu1, nu2, dS, rho, h);

[[F1, F2, Tq];
sum([F1 F2 Tq])]

%END of yukawa_force.m


function [F1, F2, Tq] = evalForces(Nb, N, x1_nbr, x2_nbr, nu1_nbr, nu2_nbr, dS_nbr, x1, x2, nu1, nu2, dS, rho, h)

% evaluates the force and torque on the particles 

% x1, x2 : parametrization of particle curves
% x1_nbr, x2_nbr : parametrization of nearby outer curve; the stress tensor T is divergence free and so we
%            can evaluate force and torque on nearby curves, versus directly on particle curve 

F1 = zeros(Nb,1);
F2 = zeros(Nb,1);
Tq = zeros(Nb,1);

u            = evalDL(    x1_nbr, x2_nbr, Nb, N, x1, x2, nu1, nu2, dS, rho, h);
[u_x1, u_x2] = evalGradDL(x1_nbr, x2_nbr, Nb, N, x1, x2, nu1, nu2, dS, rho, h);

%{
subplot(111) %reset single plot
hold off
plot(x1, x2, 'k', 'linewidth', 2);
hold on
plot(x1_nbr, x2_nbr, 'r');
quiver(x1_nbr, x2_nbr, u_x1, u_x2,'b');

hold off
axis equal;
%}

T11 = 1/rho*u.^2 + 2*rho*( 0.5*(u_x1.^2 + u_x2.^2) - u_x1.*u_x1 );
T12 =              2*rho*(                         - u_x1.*u_x2 );
T21 = T12;
T22 = 1/rho*u.^2 + 2*rho*( 0.5*(u_x1.^2 + u_x2.^2) - u_x2.*u_x2 );

for p = 1:Nb

    for i = 1:N

        F1(p) = F1(p) + (T11(i,p)*nu1_nbr(i,p) + T12(i,p)*nu2_nbr(i,p))*dS_nbr(i,p);
        
        F2(p) = F2(p) + (T21(i,p)*nu1_nbr(i,p) + T22(i,p)*nu2_nbr(i,p))*dS_nbr(i,p);
        
        Tq(p) = Tq(p) + ( x1_nbr(i,p)*(T21(i,p)*nu1_nbr(i,p) + T22(i,p)*nu2_nbr(i,p)) ...
                        - x2_nbr(i,p)*(T11(i,p)*nu1_nbr(i,p) + T12(i,p)*nu2_nbr(i,p)) ) * dS_nbr(i,p);
                        
    end

end

end

function Kh = selfEvalDL(Nb, N, x1, x2, nu1, nu2, cur, dS, rho, h)

    % evaluates double layer potential on the boundary curves x1, x2

    Kh = zeros(N*Nb,1);
    for p = 1:Nb
        for i = 1:N

            r1 = x1(:) - x1(i,p);
            r2 = x2(:) - x2(i,p);
            r = sqrt( r1.^2 + r2.^2 );
            
            rdotnu = r1.*nu1(:) + r2.*nu2(:);
            
            K1 = besselk(1, r/rho);

            I      = (p-1)*N + i;
            
            KI = -1/(2*pi)*(r/rho).*K1.*rdotnu./r.^2.*dS(:);

            KI(I) = cur(i,p)*dS(i,p)/(4*pi);
            
            Kh(I)  = dot(KI(:),h(:));
            
        end
    end

end

function Dh = evalDL(X1, X2, Nb, N, x1, x2, nu1, nu2, dS, rho, h)

    % evaluates double layer potential at (X1, X2)

    Dh = 0*X1;
    for q = 1:Nb
        for j = 1:N

            r1 = X1 - x1(j,q); r2 = X2 - x2(j,q);
            r = sqrt( r1.^2 + r2.^2 );                        
            rdotnu = r1.*nu1(j,q) + r2.*nu2(j,q);            
            K1 = besselk(1, r/rho);
            %using r = x - y; hence lack of minus below
            Dh = Dh + 1/(2*pi)*(r/rho).*K1.*rdotnu./r.^2.*dS(j,q).*h(j,q);
            
        end
    end

    uexact = besselk(0,sqrt(X1.^2 + X2.^2))/besselk(0,1);
    for j = 1:numel(X1)
      if X1(j).^2 + X2(j).^2 < 1.4
        Dh(j) = 0;
        uexact(j) = 0;
      end
    end

    clf;
    surf(X1,X2,Dh + 0*uexact)
    shading interp;
    pause 

end


function [Dh_X1, Dh_X2] = evalGradDL(X1, X2, Nb, N, x1, x2, nu1, nu2, dS, rho, h)
    % evaluates gradient of double layer potential at (X1, X2)
    
    Dh_X1 = 0*X1;
    Dh_X2 = 0*X1;    

    for q = 1:Nb
        for j = 1:N

            r1 = X1 - x1(j,q);
            r2 = X2 - x2(j,q);

            rdotnu = r1.*nu1(j,q) + r2.*nu2(j,q);

            r  = sqrt( r1.^2 + r2.^2 );
            
            K0 = besselk(0, r/rho);
            K1 = besselk(1, r/rho);
            K2 = besselk(2, r/rho);
            dK1 = -0.5*(K0 + K2); %see identity https://functions.wolfram.com/Bessel-TypeFunctions/BesselK/20/01/02/
            
%           Dh = Dh + 1/(2*pi)*(r/rho).*K1.*rdotnu./r.^2*dsD(q)*h(q,j);

            Dh_X1 = Dh_X1 + 1/(2*pi)*( ...
                  + r1./(r*rho).*K1.*rdotnu./r.^2 ...
                  + (r/rho).*dK1.*r1./(r*rho).*rdotnu./r.^2 ...
                  + (r/rho).*K1.*(1./r.^2).*(nu1(j,q) - 2*r1.*rdotnu./r.^2) ...
                  )*dS(j,q)*h(j,q);

            Dh_X2 = Dh_X2 + 1/(2*pi)*( ...
                  + r2./(r*rho).*K1.*rdotnu./r.^2 ...
                  + (r/rho).*dK1.*r2./(r*rho).*rdotnu./r.^2 ...
                  + (r/rho).*K1.*(1./r.^2).*(nu2(j,q) - 2*r2.*rdotnu./r.^2) ...
                  )*dS(j,q)*h(j,q);

        end
    end

end

%////////////////////////////////////////////////////////
% Following not used but useful for future simulations 

function [x1 x2 tau1 tau2 nu1 nu2 cur dS] = particle_curve(a, pc1, pc2, orn, param_fun, N)

% a   : characteristic radius
% pc  : particle center
% orn : orientation 
% typ : type, eg circle, ellipse, star, etc
% N   : numper of points per particle curve

T  = linspace(0, 2*pi, N+1)';
T  = T(1:end-1);
dT = T(2) - T(1);

[X1,  X2,  dX1,  dX2, ddX1, ddX2]  = param_fun(a, T  );
V   = sqrt( dX1.^2 + dX2.^2 );
%dT1 = 1./V.*( ddX1 - (ddX1.*dX1 + ddX2.*dX2).*dX1./V.^2)
%dT2 = 1./V.*( ddX1 - (ddX1.*dX1 + ddX2.*dX2).*dX1./V.^2)

cur = (ddX1.*dX2 - ddX2.*dX1)./V.^3;

% Perform final assignments along with rigid body motion 
x1 = pc1 + cos(orn)*X1 - sin(orn)*X2;
x2 = pc2 + sin(orn)*X1 + cos(orn)*X2;

dx1 = cos(orn)*dX1 - sin(orn)*dX2;
dx2 = sin(orn)*dX1 + cos(orn)*dX2;

tau1 = dx1./V; tau2 = dx2./V; 
nu1  = tau2;   nu2  = -tau1; 

dS   = V*dT; %dS*ones(size(x1));

% confirmation for ellipse only
% b = 0.05*a;
% e = sqrt(1-b^2/a^2);
% [K,E] = ellipke(e^2);
% [s(end) - 4*a*E]

%plot(t, s, T, S,'.');

%hold off
%plot(x1, x2);
%hold on
%plot3(x1, x2, cur);
%quiver(x1, x2, dx1./V, dx2./V);
%quiver(x1, x2, dT1, dT2);
%axis equal
%hold off
%pause
%return

end



function [x1, x2, dx1, dx2, ddx1, ddx2] = star_shape(a, t)

fac = 0.3;
P   = 5;

r   = (1 + fac*cos(P*t));
dr  = -P*fac*sin(P*t);
ddr = -P^2*fac*cos(P*t);

x1  =   a*cos(t).*r;
x2  =   a*sin(t).*r;

dx1 = - a*sin(t).*r + a*cos(t).*dr;
dx2 = + a*cos(t).*r + a*sin(t).*dr;

ddx1 = - a*cos(t).*r - 2*a*sin(t).*dr + a*cos(t).*ddr;
ddx2 = - a*sin(t).*r + 2*a*cos(t).*dr + a*sin(t).*ddr;

end

function [x1, x2, dx1, dx2, ddx1, ddx2] = bean_shape(a, t)

s = (t-pi)/2;

r   =     cos(s).^3 + sin(s).^3;
dr  =  -3*cos(s).^2.*sin(s)/2 + 3*sin(s).^2.*cos(s)/2;
ddr =  +6*cos(s).*sin(s).^2/4 - 3*cos(s).^3/4 + 6*sin(s).*cos(s).^2/4 - 3*sin(s).^3/4; 

x1  =   a*cos(s).*r - a*0.25;
x2  =   a*sin(s).*r - a*0.25;

dx1 = - a*sin(s)/2.*r + a*cos(s).*dr;
dx2 = + a*cos(s)/2.*r + a*sin(s).*dr;

ddx1 = - a*cos(s)/4.*r - 2*a*sin(s)/2.*dr + a*cos(s).*ddr;
ddx2 = - a*sin(s)/4.*r + 2*a*cos(s)/2.*dr + a*sin(s).*ddr;

end


function [z z_soft] = CUT_OFF(X1, X2, pc1, pc2, a, orn, param_fun, tol)

% decide whether or not a point X1, X2 lies in the curve or not
% transform point back to primal coordinates

% w is just a cutoff if too close to the curve set by tol

tmp1 = cos(-orn)*(X1-pc1) - sin(-orn)*(X2-pc2);
tmp2 = sin(-orn)*(X1-pc1) + cos(-orn)*(X2-pc2);

X1 = tmp1; X2 = tmp2;

t = linspace(0, 2*pi, 200);
dt = t(2) - t(1);
[y1 y2 dy1 dy2] = param_fun(a, t);

ds  = sqrt(dy1.^2 + dy2.^2);
nu1 = dy2./ds; nu2 = -dy1./ds;

z = 0*X1;
D = sqrt(  (X1 - y1(10)).^2 + (X2 - y2(10)).^2 );

for i = 1:length(t)
    r1 = X1 - y1(i); r2 = X2 - y2(i);
    z = z + (r1.*nu1(i) + r2.*nu2(i))./(r1.^2 + r2.^2).*ds(i)*dt/(-2*pi);
    D = min(D, sqrt( r1.^2 + r2.^2 ) );
end 

z = z > 0.5;

z_soft = D < tol | z;

end



