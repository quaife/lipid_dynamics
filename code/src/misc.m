classdef misc

%
% Some miscellaneous functions 
%
%
%

methods

function wn = wn_PnPoly(o, p1, p2, v1, v2, N)

% wn_PnPoly(): winding number test for a point in a polygon
%              useful for zeroing out interior
%      Input:   P = a point (p1, p2)
%               V = (v1, v2) = vertex points of a polygon with V(1) = V(N+1)
%      Return:  wn = the winding number (=0 only when P is outside)

    wn = 0*p1;

    for i = 1:N
    
        wn = wn + (v2(i) <= p2).*(v2(i+1)  > p2).*(o.isLeft(v1(i), v2(i), v1(i+1), v2(i+1), p1, p2) > 0);
        wn = wn - (v2(i)  > p2).*(v2(i+1) <= p2).*(o.isLeft(v1(i), v2(i), v1(i+1), v2(i+1), p1, p2) < 0);

    end
    
end

function out = isLeft(o, p1, p2, q1, q2, r1, r2)
    %    Input:  three points p, q, and r
    %    Return: >0 for r left of the line through p and q
    %            =0 for r  on the line
    %            <0 for r  right of the line     
    %http://geomalgorithms.com/a03-_inclusion.html

    out = (q1 - p1) .* (r2 - p2) - (r1 -  p1) .* (q2 - p2) ;

end 

function V = trapz2(o, x1, x2, z)

s = size(x1);
m = s(1); n = s(2);

I = 1:m-1;
J = 1:n-1;

%only one pair of substractions will be nonzero here. 

dx1 = x1(I+1,J) - x1(I,J) + x1(I,J+1) - x1(I,J); 
dx2 = x2(I+1,J) - x2(I,J) + x2(I,J+1) - x2(I,J);

V = 0.25*( z(I,J) + z(I+1,J) + z(I,J+1) + z(I+1,J+1) ).*dx1.*dx2;

V = sum(sum(V));

end


function [F1, F2, Tq] = evalForcesAlt2(o, Nb, N, x1, x2, nu1, nu2, cur, dS, rho, RHS, u)

    K = zeros(N*Nb,N*Nb);

    for q = 1:Nb
        for j = 1:N

            r1 = x1(j,q) - x1(:);
            r2 = x2(j,q) - x2(:);
            r  = sqrt( r1.^2 + r2.^2 );

            Lnu = ( r1.*nu1(j,q) + r2.*nu2(j,q) )./r.^2;

            J = (q-1)*N + j;
            K(:, J) = -1/(2*pi)*(r/rho).*besselk(1, r/rho).*Lnu*dS(j,q);

            K(J, J) = cur(j,q)*dS(j,q)/(4*pi); %removable singularity 

        end
    end
    %RHS = reshape(uD,N*Nb,1);
    size(RHS)
    h   = (1/2*eye(N*Nb) + K)\RHS;
    h   = reshape(h, N, Nb);
    
    % evaluates the force and torque on the particles 

    % x1, x2 : parametrization of Nb particle curves with normal (nu1, nu2)
    % and arc length dS
    % h is the surface density 
    % This 3rd version (Alternative 2) uses a specialized Tpq + Tqp identity, 
    % and jump identities to evaluate forces on the curves themselves 

    F1 = zeros(Nb,1);
    F2 = zeros(Nb,1);
    Tq = zeros(Nb,1);

    %periodic central difference indexing 
    Ir = [2:N, 1]';  Il = [N, 1:N-1]';
        
    for p = 1:Nb
        for q = [1:p-1, p+1:Nb]

            x1p   = x1(:,p);
            x2p   = x2(:,p);        
            dSp   = dS(:,p);
            nu1p  = nu1(:,p);
            nu2p  = nu2(:,p);
            tau1p = -nu2p;
            tau2p =  nu1p;

            %setting argin Nb = 1 tricks evaluators to using only one geometry column
            hp             = h(:, p);
            uq             = o.evalDL(     x1p, x2p, 1, N, x1(:,q), x2(:,q), nu1(:,q), nu2(:,q), dS(:,q), rho, h(:,q));
            [uq_x1, uq_x2] = o.evalGradDL( x1p, x2p, 1, N, x1(:,q), x2(:,q), nu1(:,q), nu2(:,q), dS(:,q), rho, h(:,q));

            hpt            = (hp(Ir) - hp(Il))./(2*dS(:,p));
            uqt            = uq_x1.*tau1p + uq_x2.*tau2p;   %could also use arclength derivative here ..., but already have gradient 
            uqn            = uq_x1.*nu1p  + uq_x2.*nu2p;

            Jpq1           = 2.0/rho*hp.*uq.*nu1p + 2.0*rho*hpt.*uqt.*nu1p - 2.0*rho*hpt.*uqn.*tau1p;
            Jpq2           = 2.0/rho*hp.*uq.*nu2p + 2.0*rho*hpt.*uqt.*nu2p - 2.0*rho*hpt.*uqn.*tau2p;        

            F1(p) = F1(p) + sum( Jpq1.*dSp );

            F2(p) = F2(p) + sum( Jpq2.*dSp );

            Tq(p) = Tq(p) + sum( ( x1p.*Jpq2 - x2p.*Jpq1 ) .* dSp );        

        end    
    end     

end

function [F1, F2, Tq] = evalForcesExact(o, u, u_x1, u_x2, Nb, N, x1, x2, nu1, nu2, dS, rho)

    % evaluates the force and torque on the particles 
    % using the "exact", manufactured solution 
    % x1, x2 : parametrization of particle curves

    F1 = zeros(Nb,1);
    F2 = zeros(Nb,1);
    Tq = zeros(Nb,1);
    
    T11 = 1/rho*u.^2 + 2*rho*( 0.5*(u_x1.^2 + u_x2.^2) - u_x1.*u_x1 );
    T12 =              2*rho*(                         - u_x1.*u_x2 );
    T21 = T12;
    T22 = 1/rho*u.^2 + 2*rho*( 0.5*(u_x1.^2 + u_x2.^2) - u_x2.*u_x2 );

    for p = 1:Nb

        for i = 1:N

            F1(p) = F1(p) + (T11(i,p)*nu1(i,p) + T12(i,p)*nu2(i,p))*dS(i,p);

            F2(p) = F2(p) + (T21(i,p)*nu1(i,p) + T22(i,p)*nu2(i,p))*dS(i,p);

            Tq(p) = Tq(p) + ( x1(i,p)*(T21(i,p)*nu1(i,p) + T22(i,p)*nu2(i,p)) ...
                            - x2(i,p)*(T11(i,p)*nu1(i,p) + T12(i,p)*nu2(i,p)) ) * dS(i,p);

        end

    end

end

function Dh = evalDL(o, X1, X2, Nb, N, x1, x2, nu1, nu2, dS, rho, h)

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

end


function [Dh_X1, Dh_X2] = evalGradDL(o, X1, X2, Nb, N, x1, x2, nu1, nu2, dS, rho, h)
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

function Kh = selfEvalDL(o, Nb, N, x1, x2, nu1, nu2, cur, dS, rho, h)

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




end %methods

end %classdef