function [A,ce,cw,cn,cs,cen,cwn,ces,cws,co,h,xr,zr,vr,zw] = set_nhlap_ND(nx,nz)

L = 10e3;
H = 4e3;

dx = L/nx;
dz = H/nz;
xr = repmat([1/2:1:nx-1/2]*dx,[nz 1])';
zr = repmat(-[nz-1/2:-1:1/2]*dz,[nx 1]);

xu = repmat([0:1:nx]*dx,[nz 1])';
zw = repmat(-[nz:-1:0]*dz,[nx 1]);
vr = ones(nx,nz);

h = -H*ones(1,nx);

%% masks
pmask = ones(nx,nz);
umask = pmask(1:end-1,:).*pmask(2:end,:);
wmask = pmask(:,1:end-1).*pmask(:,2:end);

%% coefficient matrix - size (nx,nz,5)
%  hom Neumann boundary conditions at the bottom and sides, hom Dirichlet at the top
clear cA; 
cA = zeros(nx,nz,5);
iz = 1; % bottom level
    ix = 1;
    cA(ix,iz,1) = wmask(ix,iz)/dz^2;
    cA(ix,iz,2) = umask(ix,iz)/dx^2;
    cA(ix,iz,5) = -(umask(ix,iz))/dx^2-(wmask(ix,iz))/dz^2;
    for ix = 2:nx-1
        cA(ix,iz,1) = wmask(ix,iz)/dz^2; 
        cA(ix,iz,2) = umask(ix,iz)/dx^2;
        cA(ix,iz,4) = umask(ix-1,iz)/dx^2;
        cA(ix,iz,5) = -(umask(ix,iz)+umask(ix-1,iz))/dx^2-(wmask(ix,iz))/dz^2;
    end
    ix = nx;
    cA(ix,iz,1) = wmask(ix,iz)/dz^2;
    cA(ix,iz,4) = umask(ix-1,iz)/dx^2;
    cA(ix,iz,5) = -(umask(ix-1,iz))/dx^2-(wmask(ix,iz))/dz^2;
for iz = 2:nz-1 % levels 2 to nz-1
    ix = 1;
    cA(ix,iz,1) = wmask(ix,iz)/dz^2;
    cA(ix,iz,2) = umask(ix,iz)/dx^2;
    cA(ix,iz,3) = wmask(ix,iz-1)/dz^2;
    cA(ix,iz,5) = -(umask(ix,iz))/dx^2-(wmask(ix,iz)+wmask(ix,iz-1))/dz^2;
    for ix = 2:nx-1
        cA(ix,iz,1) = wmask(ix,iz)/dz^2; 
        cA(ix,iz,2) = umask(ix,iz)/dx^2;
        cA(ix,iz,3) = wmask(ix,iz-1)/dz^2;
        cA(ix,iz,4) = umask(ix-1,iz)/dx^2;
        cA(ix,iz,5) = -(umask(ix,iz)+umask(ix-1,iz))/dx^2-(wmask(ix,iz)+wmask(ix,iz-1))/dz^2;
    end
    ix = nx;
    cA(ix,iz,1) = wmask(ix,iz)/dz^2; 
    cA(ix,iz,3) = wmask(ix,iz-1)/dz^2;
    cA(ix,iz,4) = umask(ix-1,iz)/dx^2;
    cA(ix,iz,5) = -(umask(ix-1,iz))/dx^2-(wmask(ix,iz)+wmask(ix,iz-1))/dz^2;
end
iz = nz; % top level
    ix = 1;
    cA(ix,iz,2) = umask(ix,iz)/dx^2;
    cA(ix,iz,3) = wmask(ix,iz-1)/dz^2;
    cA(ix,iz,5) = -(umask(ix,iz))/dx^2-(2+wmask(ix,iz-1))/dz^2; % hom Dirichlet
%    cA(ix,iz,5) = -(umask(ix,iz))/dx^2-(wmask(ix,iz-1))/dz^2; % hom Neumann
    for ix = 2:nx-1
        cA(ix,iz,2) = umask(ix,iz)/dx^2;
        cA(ix,iz,3) = wmask(ix,iz-1)/dz^2;
        cA(ix,iz,4) = umask(ix-1,iz)/dx^2;
        cA(ix,iz,5) = -(umask(ix,iz)+umask(ix-1,iz))/dx^2-(2+wmask(ix,iz-1))/dz^2; % hom Dirichlet
%        cA(ix,iz,5) = -(umask(ix,iz)+umask(ix-1,iz))/dx^2-(wmask(ix,iz-1))/dz^2; % hom Neumann
    end
    ix = nx;
    cA(ix,iz,3) = wmask(ix,iz-1)/dz^2;
    cA(ix,iz,4) = umask(ix-1,iz)/dx^2;
    cA(ix,iz,5) = -(umask(ix-1,iz))/dx^2-(2+wmask(ix,iz-1))/dz^2; % hom Dirichlet
%    cA(ix,iz,5) = -(umask(ix-1,iz))/dx^2-(wmask(ix,iz-1))/dz^2; % hom Neumann

ce = cA(:,:,2);
cw = cA(:,:,4);
cn = cA(:,:,1);
cs = cA(:,:,3);
cen = zeros(nx,nz);
cwn = zeros(nx,nz);
ces = zeros(nx,nz);
cws = zeros(nx,nz);
co = cA(:,:,5);

% %% rhs vector - size(nx,nz)
% bet = 600/L^2;
% x1 = L*0.65;
% z1 = H*(0.75-1);
% r1 = (x-x1).^2 + (z-z1).^2;
% x2 = L*0.75;
% z2 = H*(0.65-1);
% r2 = (x-x2).^2 + (z-z2).^2;
% rhs = exp(-bet*r1) - exp(-bet*r2);
% rhs1d = reshape(rhs,[nx*nz 1]);

%% solve after arranging coefficient in a sparse matrix 
clear ii isp zsp ssp A;
ii = 0;
% main diagonal 
for iz = 1:nz
    for ix = 1:nx
        ii = ii+1; 
        isp(ii) = (iz-1)*nx+ix;
        zsp(ii) = (iz-1)*nx+ix;
        ssp(ii) = cA(ix,iz,5);
    end
end
% up diagonal
for iz = 1:nz
    for ix = 1:nx-1
        ii = ii+1;
        isp(ii) = (iz-1)*nx+ix;
        zsp(ii) = (iz-1)*nx+ix+1;
        ssp(ii) = cA(ix,iz,2);
    end
end
% down diagonal
for iz = 1:nz
    for ix = 2:nx
        ii = ii+1;
        isp(ii) = (iz-1)*nx+ix;
        zsp(ii) = (iz-1)*nx+ix-1;
        ssp(ii) = cA(ix,iz,4);
    end
end
% up up diagonal
for iz = 1:nz-1
    for ix = 1:nx
        ii = ii+1;
        isp(ii) = (iz-1)*nx+ix;
        zsp(ii) = (iz-1+1)*nx+ix;
        ssp(ii) = cA(ix,iz,1);
    end
end
% down down diagonal
for iz = 2:nz
    for ix = 1:nx
        ii = ii+1;
        isp(ii) = (iz-1)*nx+ix;
        zsp(ii) = (iz-1-1)*nx+ix;
        ssp(ii) = cA(ix,iz,3);
    end
end
%
A = sparse(isp,zsp,ssp,nx*nz,nx*nz); % test the propertie of the matrix. ellipticity?
%
% p1d = A\rhs1d;
% p = reshape(p1d,[nx nz]);
% p(pmask==0)=NaN;
% 
% %% plot
% figure; 
% pcolor(x-dx/2,z-dz/2,rhs); shading flat; colorbar;
% figure; 
% pcolor(x-dx/2,z-dz/2,p); shading flat; colorbar;
% hold on;
% pmx = max(max(p)); contour(x,z,p,[-1:0.1:1]*pmx,'k')
% 
% figure;
% plot(z(96,:),rhs(96,:),'k*-');
% figure;
% plot(z(96,:),p(96,:),'k*-');
% 
% %% restrict
% 
% x_l2   = 0.25*(x(1:2:end-1,1:2:end-1)+x(1:2:end-1,2:2:end)+x(2:2:end,1:2:end-1)+x(2:2:end,2:2:end));
% z_l2   = 0.25*(z(1:2:end-1,1:2:end-1)+z(1:2:end-1,2:2:end)+z(2:2:end,1:2:end-1)+z(2:2:end,2:2:end));
% rhs_l2 = 0.25*(rhs(1:2:end-1,1:2:end-1)+rhs(1:2:end-1,2:2:end)+rhs(2:2:end,1:2:end-1)+rhs(2:2:end,2:2:end));
% 
% plot(z_l2(48,:),rhs_l2(48,:),'k*-');
% 
% x_l3   = 0.25*(x_l2(1:2:end-1,1:2:end-1)+x_l2(1:2:end-1,2:2:end)+x_l2(2:2:end,1:2:end-1)+x_l2(2:2:end,2:2:end));
% z_l3   = 0.25*(z_l2(1:2:end-1,1:2:end-1)+z_l2(1:2:end-1,2:2:end)+z_l2(2:2:end,1:2:end-1)+z_l2(2:2:end,2:2:end));
% rhs_l3 = 0.25*(rhs_l2(1:2:end-1,1:2:end-1)+rhs_l2(1:2:end-1,2:2:end)+rhs_l2(2:2:end,1:2:end-1)+rhs_l2(2:2:end,2:2:end));
% 
% plot(z_l3(24,:),rhs_l3(24,:),'k*-');
% 
% x_l4   = 0.25*(x_l3(1:2:end-1,1:2:end-1)+x_l3(1:2:end-1,2:2:end)+x_l3(2:2:end,1:2:end-1)+x_l3(2:2:end,2:2:end));
% z_l4   = 0.25*(z_l3(1:2:end-1,1:2:end-1)+z_l3(1:2:end-1,2:2:end)+z_l3(2:2:end,1:2:end-1)+z_l3(2:2:end,2:2:end));
% rhs_l4 = 0.25*(rhs_l3(1:2:end-1,1:2:end-1)+rhs_l3(1:2:end-1,2:2:end)+rhs_l3(2:2:end,1:2:end-1)+rhs_l3(2:2:end,2:2:end));
% 
% plot(z_l4(12,:),rhs_l4(12,:),'k*-');
% 
% x_l5   = 0.25*(x_l4(1:2:end-1,1:2:end-1)+x_l4(1:2:end-1,2:2:end)+x_l4(2:2:end,1:2:end-1)+x_l4(2:2:end,2:2:end));
% z_l5   = 0.25*(z_l4(1:2:end-1,1:2:end-1)+z_l4(1:2:end-1,2:2:end)+z_l4(2:2:end,1:2:end-1)+z_l4(2:2:end,2:2:end));
% rhs_l5 = 0.25*(rhs_l4(1:2:end-1,1:2:end-1)+rhs_l4(1:2:end-1,2:2:end)+rhs_l4(2:2:end,1:2:end-1)+rhs_l4(2:2:end,2:2:end));
% 
% plot(z_l5(6,:),rhs_l5(6,:),'k*-');
% 
% x_l6   = 0.25*(x_l5(1:2:end-1,1:2:end-1)+x_l5(1:2:end-1,2:2:end)+x_l5(2:2:end,1:2:end-1)+x_l5(2:2:end,2:2:end));
% z_l6   = 0.25*(z_l5(1:2:end-1,1:2:end-1)+z_l5(1:2:end-1,2:2:end)+z_l5(2:2:end,1:2:end-1)+z_l5(2:2:end,2:2:end));
% rhs_l6 = 0.25*(rhs_l5(1:2:end-1,1:2:end-1)+rhs_l5(1:2:end-1,2:2:end)+rhs_l5(2:2:end,1:2:end-1)+rhs_l5(2:2:end,2:2:end));
% 
% plot(z_l6(3,:),rhs_l6(3,:),'k*-');