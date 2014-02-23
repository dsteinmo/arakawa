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

dt = 0.25*min(d./c0(:));

%make q grid ('cell vertices)
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
eta0 = 4*exp(-((xxh-0.55*Lx)/3e2).^2 - ((yyh-0.5*Ly)/3e2).^2);
%eta0 = 2*exp(-((xxh-0.55*Lx)/3e2).^2);

%interpolate h to u grid
h0 = H+eta0;
h0_ghostx = [h0(:,1) h0 h0(:,end)];
h0_u = zeros(Ny-1,Nx);
for i =1:Nx
    h0_u(:,i) = 0.5*(h0_ghostx(:,i)+h0_ghostx(:,i+1));
end

eta0_u = h0_u-H0;

u0 = sqrt(g/H0)*eta0_u;

pcolor(xxu,yyu,u0); shading flat; colorbar; drawnow;

%u0 = xxu*0+10;
v0 = xxv*0;
%end user input


NUMSTEPS = ceil(FINTIME/dt);

t=0;
for k = 1:NUMSTEPS
    %interpolate h to u grid
    h0_ghostx = [h0(:,1) h0 h0(:,end)];
    h0_u = zeros(Ny-1,Nx);
    for i =1:Nx
        h0_u(:,i) = 0.5*(h0_ghostx(:,i)+h0_ghostx(:,i+1));
    end

    %interpolate h to v grid
    h0_ghosty = [h0(1,:);
                 h0;
                 h0(end,:)];
    h0_v = zeros(Ny,Nx-1);
    for j=1:Ny
        h0_v(j,:) = 0.5*(h0_ghosty(j,:)+h0_ghosty(j+1,:));    
    end

    %interpolate h to q grid
    h0_ghostq = [h0(1,1)   h0(1,:)     h0(1,end);
                 h0(:,1)   h0          h0(:,end);
                 h0(end,1) h0(end,:) h0(end,end)];

    h0_q = zeros(Ny,Nx);
    for j=1:Ny
        for i=1:Nx
            h0_q(j,i) = 0.25*( h0_ghostq(j,i)   + h0_ghostq(j,i+1)  + ...
                               h0_ghostq(j+1,i) + h0_ghostq(j+1,i+1) );
        end
    end

    zeta = zeros(Ny,Nx);

    u0ghosty = [u0(1,:); u0; u0(end,:)];
    v0ghostx = [v0(:,1) v0 v0(:,end)];

    %compute relative vorticity
    for j=1:Ny
        for i=1:Nx
            zeta(j,i) = dinv*(u0ghosty(j,i)-u0ghosty(j+1,i) + v0ghostx(j,i+1)-v0ghostx(j,i));
        end
    end
    %compute potential vorticity
    q = (f+zeta)./h0_q;

    %form transports
    ustar = h0_u.*u0;
    vstar = h0_v.*v0;

    %RHSmass lives at cell centres, can calculate it vectorized
    %-div(vstar)_{i+.5,j+.5} = (1/d)*(ustar_{i+1,j+.5} - ustar_{i,j+.5} + 
    %                                 vstar_{i+.5,j+1} - vstar_{i+.5,j})
    RHSmass = -dinv*(ustar(:,2:end) - ustar(:,1:end-1) + ...
                     vstar(2:end,:) - vstar(1:end-1,:));



    %calculate coefficients of u in "u cross omega terms:"
    epsilon = zeros(Ny-1,Nx-1);
    phi = zeros(Ny-1,Nx-1);
    alpha = zeros(Ny-1,Nx-1);
    beta = zeros(Ny-1,Nx-1);
    gamma= zeros(Ny-1,Nx-1);
    delta= zeros(Ny-1,Nx-1);

    %may need to put ghost strips around q for this to work
    for j=1:Ny-1
        for i=1:Nx-1
            epsilon(j,i) = (1/24)*(q(j+1,i+1) + q(j+1,i) - q(j,i) - q(j,i+1));    
            phi(j,i) = (1/24)*(-q(j+1,i+1) + q(j+1,i) + q(j,i) - q(j,i+1));
            alpha(j,i) = (1/24)*(2*q(j+1,i+1) + q(j+1,i) + 2*q(j,i) + q(j,i+1));
            beta(j,i+1) = (1/24)*(q(j+1,i+1) + 2*q(j+1,i) + q(i,j) + 2*q(j,i+1));
            gamma(j,i+1) = (1/24)*(2*q(j+1,i+1) + q(j+1,i) + 2*q(i,j) + q(j,i+1));
            delta(j,i)   = (1/24)*(q(j+1,i+1) + 2*q(j+1,i) + q(j,i) + 2*q(j,i+1));
        end
    end

    %Phi is geo-potential
    %K is kinetic energy, lives at cell centers
    Phi = g*eta0;
    K = zeros(Ny-1,Nx-1);

    %interpolate v^2 to h grid
    v2_h = zeros(Ny-1,Nx-1);
    for j=1:Ny-1
        v2_h(j,:) = 0.5*(v0(j,:).^2+v0(j+1,:).^2);    
    end

    %interpolate u^2 to h grid
    u2_h = zeros(Ny-1,Nx-1);
    for i=1:Nx-1
       u2_h(:,i) = 0.5*(u0(:,i).^2+u0(:,i+1).^2);  
    end

    %form kinetic energy on h-grid
    K = 0.5*(u2_h + v2_h);

    %define total energy
    E = K + Phi;

    RHSmomu = zeros(Ny-1,Nx);
    RHSmomv = zeros(Ny,Nx-1);

    %calculate FD approximation to RHS of momentum mequations:
    for j=2:Ny-1
        for i = 2:Nx-1
            RHSmomu(j,i) = alpha(j,i)*vstar(j+1,i) + beta(j,i)*vstar(j+1,i-1) ...
                         + gamma(j,i)*vstar(j,i-1) - delta(j,i)*vstar(j,i) ...
                         - epsilon(j,i)*ustar(j,i+1) + epsilon(j,i-1)*ustar(j,i-1)  ...
                         - dinv*(E(j,i)-E(j,i-1));

            RHSmomv(j,i) =-gamma(j,i+1)*ustar(j,i+1) - delta(j,i)*ustar(j,i) ...
                          -alpha(j-1,i)*ustar(j-1,i) - beta(j-1,i+1)*ustar(j-1,i+1) ...
                          -phi(j,i)*vstar(j+1,i) + phi(j-1,i)*vstar(j-1,i) ...
                          - dinv*(E(j,i)-E(j-1,i));
        end
    end

    %time-step all equations
    h1 = h0 + dt*RHSmass;
    u1 = u0 + dt*RHSmomu;
    v1 = v0 + dt*RHSmomv;
    t=t+dt;

    %update eta
    eta1 = h1-H;

    %cycle fields for next time-step.
    eta0=eta1;
    h0=h1;
    u0=u1;
    v0=v1;
    
    figure(2);
    pcolor(xxh,yyh,eta0); title(['t=' num2str(t) ' s']); shading flat; colorbar; axis equal; axis tight; drawnow; 

end
%figure(1);
%subplot(2,1,1);
%plot(xxq(:),yyq(:),'.k',xxh(:),yyh(:),'.b',xxu(:),yyu(:),'.y',xxv(:),yyv(:),'.g')
%subplot(2,1,2);
%pcolor(xxh,yyh,eta0); colorbar; axis equal; axis tight;
