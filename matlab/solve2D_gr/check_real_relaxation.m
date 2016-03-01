%% check sigma coord relaxation

% param
nx = 128;
ny = 128
nz = 128;
Lx = 1e4;
Ly = 1e4;
H = 4e3;

% matlab results
rhs_mat = rhs;
p_mat = p;

% mgroms results
zr_000 = ncread('../../src/zr_zr_000.nc','zr');
zr_001 = ncread('../../src/zr_zr_001.nc','zr');
zr_002 = ncread('../../src/zr_zr_002.nc','zr');
zr_003 = ncread('../../src/zr_zr_003.nc','zr');
zw_000 = ncread('../../src/zw_zw_000.nc','zw');
zw_001 = ncread('../../src/zw_zw_001.nc','zw');
zw_002 = ncread('../../src/zw_zw_002.nc','zw');
zw_003 = ncread('../../src/zw_zw_003.nc','zw');
rhs_000 = ncread('../../src/rhs_rhs_000.nc','rhs');
rhs_001 = ncread('../../src/rhs_rhs_001.nc','rhs');
rhs_002 = ncread('../../src/rhs_rhs_002.nc','rhs');
rhs_003 = ncread('../../src/rhs_rhs_003.nc','rhs');
p_000 = ncread('../../src/p_p_000.nc','p');
p_001 = ncread('../../src/p_p_001.nc','p');
p_002 = ncread('../../src/p_p_002.nc','p');
p_003 = ncread('../../src/p_p_003.nc','p');
r_000 = ncread('../../src/r_r_000.nc','r');
r_001 = ncread('../../src/r_r_001.nc','r');
r_002 = ncread('../../src/r_r_002.nc','r');
r_003 = ncread('../../src/r_r_003.nc','r');

% pcolor plot sections
i = 60;

zr = cat(2,squeeze(zr_000(:,2:end-1,i+1)),squeeze(zr_002(:,2:end-1,i+1)));
zw = cat(2,squeeze(zw_000(:,2:end-1,i+1)),squeeze(zw_002(:,2:end-1,i+1)));
rhs = cat(2,squeeze(rhs_000(:,2:end-1,i+1)),squeeze(rhs_002(:,2:end-1,i+1)));
p = cat(2,squeeze(p_000(:,2:end-1,i+1)),squeeze(p_002(:,2:end-1,i+1)));
r = cat(2,squeeze(r_000(:,2:end-1,i+1)),squeeze(r_002(:,2:end-1,i+1)));

yr = repmat(((0.5:1:ny)*Ly/ny),[nz 1]);
yv_pc = cat(2,zeros(nz,1),0.5*(yr(:,1:end-1)+yr(:,2:end)),Ly*ones(nz,1));
yv_pc = cat(1,yv_pc,yv_pc(end,:));
zwv_pc = cat(2,zw(:,1),0.5*(zw(:,1:end-1)+zw(:,2:end)),zw(:,end));
rhs_pc = cat(2,rhs,rhs(:,end));
rhs_pc = cat(1,rhs_pc,rhs_pc(end,:));
p_pc = cat(2,p,p(:,end));
p_pc = cat(1,p_pc,p_pc(end,:));
r_pc = cat(2,r,r(:,end));
r_pc = cat(1,r_pc,r_pc(end,:));

figure;
pcolor(yv_pc,zwv_pc,p_pc); shading flat;
hold on
contour(yr,zr,p,[-14:1:14]*1e4,'k')
plot(yr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 Ly])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p(:));
cmax = max(p(:));
caxis([cmin cmax]);
title(['pressure  mgroms i=' num2str(i) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

% pcolor plot sections
j = 60;

zr = cat(2,squeeze(zr_000(:,j+1,2:end-1)),squeeze(zr_001(:,j+1,2:end-1)));
zw = cat(2,squeeze(zw_000(:,j+1,2:end-1)),squeeze(zw_001(:,j+1,2:end-1)));
rhs = cat(2,squeeze(rhs_000(:,j+1,2:end-1)),squeeze(rhs_001(:,j+1,2:end-1)));
p = cat(2,squeeze(p_000(:,j+1,2:end-1)),squeeze(p_001(:,j+1,2:end-1)));
r = cat(2,squeeze(r_000(:,j+1,2:end-1)),squeeze(r_001(:,j+1,2:end-1)));

xr = repmat(((0.5:1:nx)*Lx/nx),[nz 1]);
xu_pc = cat(2,zeros(nz,1),0.5*(xr(:,1:end-1)+xr(:,2:end)),Lx*ones(nz,1));
xu_pc = cat(1,xu_pc,xu_pc(end,:));
zwu_pc = cat(2,zw(:,1),0.5*(zw(:,1:end-1)+zw(:,2:end)),zw(:,end));
rhs_pc = cat(2,rhs,rhs(:,end));
rhs_pc = cat(1,rhs_pc,rhs_pc(end,:));
p_pc = cat(2,p,p(:,end));
p_pc = cat(1,p_pc,p_pc(end,:));
r_pc = cat(2,r,r(:,end));
r_pc = cat(1,r_pc,r_pc(end,:));

figure;
pcolor(xu_pc,zwu_pc,rhs_pc); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 Lx])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(rhs(:));
cmax = max(rhs(:));
caxis([cmin cmax]);
title(['rhs  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

figure;
pcolor(xu_pc,zwu_pc,p_pc); shading flat;
hold on
contour(xr,zr,p,[-14:1:14]*1e4,'k')
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 Lx])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p(:));
cmax = max(p(:));
caxis([cmin cmax]);
title(['pressure mgroms j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

figure;
pcolor(xu_pc,zwu_pc,r_pc); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 Lx])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(r(:));
cmax = max(r(:));
caxis([cmin cmax]);
title(['residual mgroms j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

p_mat_pc = cat(2,p_mat,p_mat(:,end));
p_mat_pc = cat(1,p_mat_pc,p_mat_pc(end,:));

figure;
pcolor(xu_pc,zwu_pc,p_pc-p_mat_pc); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 Lx])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p(:)-p_mat(:));
cmax = max(p(:)-p_mat(:));
caxis([cmin cmax]);
title(['pressure error  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

% imagesc plot maps
k = 1;

p_1 = cat(2,squeeze(p_000(k,2:end-1,2:end-1)),squeeze(p_001(k,2:end-1,2:end-1)));
p_2 = cat(2,squeeze(p_002(k,2:end-1,2:end-1)),squeeze(p_003(k,2:end-1,2:end-1)));
p = cat(1,p_1,p_2); clear p_1 p_2;
r_1 = cat(2,squeeze(r_000(k,2:end-1,2:end-1)),squeeze(r_001(k,2:end-1,2:end-1)));
r_2 = cat(2,squeeze(r_002(k,2:end-1,2:end-1)),squeeze(r_003(k,2:end-1,2:end-1)));
r = cat(1,r_1,r_2); clear r_1 r_2;

figure;
imagesc(r); axis xy
colorbar
cmin = min(r(:));
cmax = max(r(:));
caxis([cmin cmax]);
title(['residual mgroms k=' num2str(k) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

figure;
imagesc(p-repmat(p_mat(k,:),[nx 1])); axis xy
colorbar
cmin = min(min(p-repmat(p_mat(k,:),[nx 1])));
cmax = max(max(p-repmat(p_mat(k,:),[nx 1])));
caxis([cmin cmax]);
title(['pressure error  k=' num2str(k) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
