%solves shallow water equations with 2nd order FD scheme that conserves
%energy and enstrophy. Momentum equations are written using velocity-form
%(non-conservative), but nonlinear terms written in energy-vorticity
%formulation.

%compare this to Fourier in ~/work/SWNH2d/test/SWNHkarn.m

clear;
close all;

Lx=4e3;
Ly=4e3;
%Ly=4.15e3;
Nx=128;   %this sets resolution
g = 9.81;
H0 = 10;
%f=1e-1;
f=0;
FINTIME = 3600;


d = Lx/Nx;
dinv = 1/d;
Ny = ceil(Ly/d);
Ly = Ny*d;   %Ly may be slightly adjusted from input

H=H0*ones(Ny-1,Nx-1);
c0 = sqrt(g*H);

dt = 0.15*min(d./c0(:));

%make q grid (cell vertices)
xq = linspace(0,Lx,Nx);
yq = linspace(0,Ly,Ny);
[xxq,yyq] = meshgrid(xq,yq);

%make h grid (cell centers)
xh = linspace(d/2,Lx-d/2,Nx-1);
yh = linspace(d/2,Ly-d/2,Ny-1);
[xxh,yyh] = meshgrid(xh,yh);

%make u grid (center of vertical cell edges)
xu = linspace(0,Lx,Nx);
yu = linspace(d/2,Ly-d/2,Ny-1);
[xxu,yyu] = meshgrid(xu,yu);

%make v grid (center of horizontal cell edges)
xv = linspace(d/2,Lx-d/2,Nx-1);
yv = linspace(0,Ly,Ny);
[xxv,yyv] = meshgrid(xv,yv);

%%begin user input
eta0 = .05*exp(-((xxh-0.55*Lx)/3e2).^2 - ((yyh-0.5*Ly)/3e2).^2);
%eta0 = 2*exp(-((xxh-0.55*Lx)/3e2).^2);

%interpolate h to u grid
h0 = H+eta0;
h0_ghostx = [h0(:,1) h0 h0(:,end)];
h0_u = zeros(Ny-1,Nx);
for i =1:Nx
    h0_u(:,i) = 0.5*(h0_ghostx(:,i)+h0_ghostx(:,i+1));
end

eta0_u = h0_u-H0;

u0 = -sqrt(g/H0)*eta0_u;

pcolor(xxu,yyu,u0); shading flat; colorbar; drawnow;

%u0 = xxu*0+10;
v0 = xxv*0;
%end user input


NUMSTEPS = ceil(FINTIME/dt);

t=0;
h1=h0;
u1=u0;
v1=v0;
eta1=eta0;
for k = 1:NUMSTEPS
    %interpolate h to u grid
    h1_ghostx = [h1(:,1) h1 h1(:,end)];
    h1_u = zeros(Ny-1,Nx);
    for i =1:Nx
        h1_u(:,i) = 0.5*(h1_ghostx(:,i)+h1_ghostx(:,i+1));
    end

    %interpolate h to v grid
    h1_ghosty = [h1(1,:);
                 h1;
                 h1(end,:)];
    h1_v = zeros(Ny,Nx-1);
    for j=1:Ny
        h1_v(j,:) = 0.5*(h1_ghosty(j,:)+h1_ghosty(j+1,:));    
    end

    %interpolate h to q grid
    h1_ghostq = [h1(1,1)   h1(1,:)     h1(1,end);
                 h1(:,1)   h1          h1(:,end);
                 h1(end,1) h1(end,:) h1(end,end)];

    h1_q = zeros(Ny,Nx);
    for j=1:Ny
        for i=1:Nx
            h1_q(j,i) = 0.25*( h1_ghostq(j,i)   + h1_ghostq(j,i+1)  + ...
                               h1_ghostq(j+1,i) + h1_ghostq(j+1,i+1) );
        end
    end

    zeta = zeros(Ny,Nx);

    u1ghosty = [-u1(1,:); u1; -u1(end,:)];
    v1ghostx = [-v1(:,1) v1 -v1(:,end)];

    %compute relative vorticity
    for j=1:Ny
        for i=1:Nx
            zeta(j,i) = dinv*(u1ghosty(j,i)-u1ghosty(j+1,i) + v1ghostx(j,i+1)-v1ghostx(j,i));
        end
    end
    %compute potential vorticity
    q = (f+zeta)./h1_q;

    %form transports
    ustar = h1_u.*u1;
    vstar = h1_v.*v1;
    
    %put ghost strips around ustar
    ustarghost = [-ustar(1,1)   ustar(1,:)     -ustar(1,end);
                 -ustar(:,1)   ustar          -ustar(:,end);
                 -ustar(end,1) ustar(end,:) ustar(end,end)];
             
   %put ghost strips around vstar
   vstarghost = [-vstar(1,1)   -vstar(1,:)     -vstar(1,end);
     vstar(:,1)   vstar          vstar(:,end);
     -vstar(end,1) -vstar(end,:) -vstar(end,end)];
             
    %RHSmass lives at cell centres, can calculate it vectorized
    %-div(vstar)_{i+.5,j+.5} = (1/d)*(ustar_{i+1,j+.5} - ustar_{i,j+.5} + 
    %                                 vstar_{i+.5,j+1} - vstar_{i+.5,j})
    %don't need to use ghost variables here, since RHSmass lives at
    %cell centers
    RHSmass = -dinv*(ustar(:,2:end) - ustar(:,1:end-1) + ...
                     vstar(2:end,:) - vstar(1:end-1,:));



    %calculate coefficients of u in "u cross omega terms:"
    epsilon = zeros(Ny+1,Nx+1);
    phi = zeros(Ny+1,Nx+1);
    alpha = zeros(Ny+1,Nx+1);
    beta = zeros(Ny+1,Nx+1);
    gamma= zeros(Ny+1,Nx+1);
    delta= zeros(Ny+1,Nx+1);

    %may need to put ghost strips around q for this to work
    
    %put ghost strips around q
    qghost = [q(1,1)   q(1,:)     q(1,end);
                 q(:,1)   q          q(:,end);
                 q(end,1) q(end,:) q(end,end)];
    
    for j=1:Ny+1
        for i=1:Nx
            epsilon(j,i) = (1/24)*(qghost(j+1,i+1) + qghost(j+1,i) - qghost(j,i) - qghost(j,i+1));    
            phi(j,i) = (1/24)*(-qghost(j+1,i+1) + qghost(j+1,i) + qghost(j,i) - qghost(j,i+1));
            alpha(j,i) = (1/24)*(2*qghost(j+1,i+1) + qghost(j+1,i) + 2*qghost(j,i) + qghost(j,i+1));
            beta(j,i+1) = (1/24)*(qghost(j+1,i+1) + 2*qghost(j+1,i) + qghost(j,i) + 2*qghost(j,i+1));
            gamma(j,i+1) = (1/24)*(2*qghost(j+1,i+1) + qghost(j+1,i) + 2*qghost(j,i) + qghost(j,i+1));
            delta(j,i)   = (1/24)*(qghost(j+1,i+1) + 2*qghost(j+1,i) + qghost(j,i) + 2*qghost(j,i+1));
        end
    end

    %Phi is geo-potential
    %K is kinetic energy, lives at cell centers
    Phi = g*eta1;

    %interpolate v^2 to h grid
    v2_h = zeros(Ny-1,Nx-1);
    for j=1:Ny-1
        v2_h(j,:) = 0.5*(v1(j,:).^2+v1(j+1,:).^2);    
    end

    %interpolate u^2 to h grid
    u2_h = zeros(Ny-1,Nx-1);
    for i=1:Nx-1
       u2_h(:,i) = 0.5*(u1(:,i).^2+u1(:,i+1).^2);  
    end

    %form kinetic energy on h-grid
    K = 0.5*(u2_h + v2_h);

    %define total energy
    E = K + Phi;
    
       %put ghost strips around E
    Eghost = [E(1,1)   E(1,:)     E(1,end);
     E(:,1)   E          E(:,end);
     E(end,1) E(end,:) E(end,end)];

    RHSmomu = zeros(Ny-1,Nx);
    RHSmomv = zeros(Ny,Nx-1);

    %calculate FD approximation to RHS of momentum equations:
    for j=1:Ny-1
        for i = 1:Nx
            RHSmomu(j,i) = alpha(j+1,i+1)*vstarghost(j+2,i+1) + beta(j+1,i+1)*vstarghost(j+2,i) ...
                         + gamma(j+1,i+1)*vstarghost(j+1,i) - delta(j+1,i+1)*vstarghost(j+1,i+1) ...
                         - epsilon(j+1,i+1)*ustarghost(j+1,i+2) + epsilon(j+1,i)*ustarghost(j+1,i)  ...
                         - dinv*(Eghost(j+1,i+1)-Eghost(j+1,i)); %why j+1?
        end
    end
    
    for j=1:Ny
        for i = 1:Nx-1
            RHSmomv(j,i) =-gamma(j+1,i+2)*ustarghost(j+1,i+2) - delta(j+1,i+1)*ustarghost(j+1,i+1) ...
                      -alpha(j,i+1)*ustarghost(j,i+1) - beta(j,i+2)*ustarghost(j,i+2) ...
                      -phi(j+1,i+1)*vstarghost(j+2,i+1) + phi(j,i+1)*vstarghost(j,i+1) ...
                      - dinv*(Eghost(j+1,i+1)-Eghost(j,i+1)); %why i+1?
        end
    end

    %time-step all equations
    if k == 1
        %forward euler:
        h1 = h0 + dt*RHSmass;
        u1 = u0 + dt*RHSmomu;
        v1 = v0 + dt*RHSmomv;
        %update eta
        eta1 = h1-H;
        %cycle fields for next time-step.
        RHSmass0=RHSmass;
        RHSmomu0=RHSmomu;
        RHSmomv0=RHSmomv;
        eta0=eta1;
        h0=h1;
        u0=u1;
        v0=v1;
    else
        %Euler is unstable for long times, need something else for
        %non-initial time-steps...
        %Leapfrog:
        %h2 = h0 + 2*dt*RHSmass;
        %u2 = u0 + 2*dt*RHSmomu;
        %v2 = v0 + 2*dt*RHSmomv;
        %AB2:
        RHSmass1=RHSmass;
        RHSmomu1=RHSmomu;
        RHSmomv1=RHSmomv;
        h2 = h1 + dt*(1.5*RHSmass1-0.5*RHSmass0);
        u2 = u1 + dt*(1.5*RHSmomu1-0.5*RHSmomu0);
        v2 = v1 + dt*(1.5*RHSmomv1-0.5*RHSmomv0);
        %update eta
        eta2 = h2-H;
        %cycle fields for next time-step.
        eta0=eta1;
        eta1=eta2;
        h0=h1;
        h1=h2;
        u0=u1;
        u1=u2;
        v0=v1;
        v1=v2;
        RHSmass0=RHSmass1;
        RHSmomu0=RHSmomu1;
        RHSmomv0=RHSmomv1;
    end
    t=t+dt;

    
    
    
    figure(2);
    pcolor(xxh,yyh,eta1); title(['t=' num2str(t) ' s']); shading flat; colorbar; axis equal; axis tight; drawnow; 

end
%figure(1);
%subplot(2,1,1);
%plot(xxq(:),yyq(:),'.k',xxh(:),yyh(:),'.b',xxu(:),yyu(:),'.y',xxv(:),yyv(:),'.g')
%subplot(2,1,2);
%pcolor(xxh,yyh,eta0); colorbar; axis equal; axis tight;
