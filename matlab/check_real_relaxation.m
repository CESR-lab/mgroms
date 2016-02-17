clear all

% %% Nico's solution (only flat bottom)
% solve_pressure_ND;
% rhs_ND = rhs';
% p_ND = p';
% clear rhs p;

%% Jeroen's solution
solve_pressure_JM_ND;
rhs_JM = rhs';
p_JM = p';
clear rhs p;

%% Guillaume's solution
solve_pressure_GR_ND;
rhs_GR = rhs;
p_GR = p;
clear rhs p;

%% mgroms solution
umask_000 = ncread('../src/umask_umask_000.nc','umask');
umask_001 = ncread('../src/umask_umask_001.nc','umask');
umask_002 = ncread('../src/umask_umask_002.nc','umask');
umask_003 = ncread('../src/umask_umask_003.nc','umask');
vmask_000 = ncread('../src/vmask_vmask_000.nc','vmask');
vmask_001 = ncread('../src/vmask_vmask_001.nc','vmask');
vmask_002 = ncread('../src/vmask_vmask_002.nc','vmask');
vmask_003 = ncread('../src/vmask_vmask_003.nc','vmask');
zr_000 = ncread('../src/zr_zr_000.nc','zr');
zr_001 = ncread('../src/zr_zr_001.nc','zr');
zr_002 = ncread('../src/zr_zr_002.nc','zr');
zr_003 = ncread('../src/zr_zr_003.nc','zr');
zw_000 = ncread('../src/zw_zw_000.nc','zw');
zw_001 = ncread('../src/zw_zw_001.nc','zw');
zw_002 = ncread('../src/zw_zw_002.nc','zw');
zw_003 = ncread('../src/zw_zw_003.nc','zw');
cA_000 = ncread('../src/cA_cA_000.nc','cA');
cA_001 = ncread('../src/cA_cA_001.nc','cA');
cA_002 = ncread('../src/cA_cA_002.nc','cA');
cA_003 = ncread('../src/cA_cA_003.nc','cA');
rhs_000 = ncread('../src/rhs_rhs_000.nc','rhs');
rhs_001 = ncread('../src/rhs_rhs_001.nc','rhs');
rhs_002 = ncread('../src/rhs_rhs_002.nc','rhs');
rhs_003 = ncread('../src/rhs_rhs_003.nc','rhs');
p_000 = ncread('../src/p_p_000.nc','p');
p_001 = ncread('../src/p_p_001.nc','p');
p_002 = ncread('../src/p_p_002.nc','p');
p_003 = ncread('../src/p_p_003.nc','p');
r_000 = ncread('../src/r_r_000.nc','r');
r_001 = ncread('../src/r_r_001.nc','r');
r_002 = ncread('../src/r_r_002.nc','r');
r_003 = ncread('../src/r_r_003.nc','r');

%% imagesc plot matrix
cA1_1 = cat(3,squeeze(cA_000(1,:,1:end-1,1:end-1)),squeeze(cA_001(1,:,1:end-1,2:end)));
cA1_2 = cat(3,squeeze(cA_002(1,:,2:end,1:end-1)),squeeze(cA_003(1,:,2:end,2:end)));
cA1 = cat(2,cA1_1,cA1_2); clear cA1_1 cA1_2;
cA2_1 = cat(3,squeeze(cA_000(2,:,1:end-1,1:end-1)),squeeze(cA_001(2,:,1:end-1,2:end)));
cA2_2 = cat(3,squeeze(cA_002(2,:,2:end,1:end-1)),squeeze(cA_003(2,:,2:end,2:end)));
cA2 = cat(2,cA2_1,cA2_2); clear cA2_1 cA2_2;
cA3_1 = cat(3,squeeze(cA_000(3,:,1:end-1,1:end-1)),squeeze(cA_001(3,:,1:end-1,2:end)));
cA3_2 = cat(3,squeeze(cA_002(3,:,2:end,1:end-1)),squeeze(cA_003(3,:,2:end,2:end)));
cA3 = cat(2,cA3_1,cA3_2); clear cA3_1 cA3_2;
cA4_1 = cat(3,squeeze(cA_000(4,:,1:end-1,1:end-1)),squeeze(cA_001(4,:,1:end-1,2:end)));
cA4_2 = cat(3,squeeze(cA_002(4,:,2:end,1:end-1)),squeeze(cA_003(4,:,2:end,2:end)));
cA4 = cat(2,cA4_1,cA4_2); clear cA4_1 cA4_2;
cA5_1 = cat(3,squeeze(cA_000(5,:,1:end-1,1:end-1)),squeeze(cA_001(5,:,1:end-1,2:end)));
cA5_2 = cat(3,squeeze(cA_002(5,:,2:end,1:end-1)),squeeze(cA_003(5,:,2:end,2:end)));
cA5 = cat(2,cA5_1,cA5_2); clear cA5_1 cA5_2;
cA6_1 = cat(3,squeeze(cA_000(6,:,1:end-1,1:end-1)),squeeze(cA_001(6,:,1:end-1,2:end)));
cA6_2 = cat(3,squeeze(cA_002(6,:,2:end,1:end-1)),squeeze(cA_003(6,:,2:end,2:end)));
cA6 = cat(2,cA6_1,cA6_2); clear cA6_1 cA6_2;
cA7_1 = cat(3,squeeze(cA_000(7,:,1:end-1,1:end-1)),squeeze(cA_001(7,:,1:end-1,2:end)));
cA7_2 = cat(3,squeeze(cA_002(7,:,2:end,1:end-1)),squeeze(cA_003(7,:,2:end,2:end)));
cA7 = cat(2,cA7_1,cA7_2); clear cA7_1 cA7_2;
cA8_1 = cat(3,squeeze(cA_000(8,:,1:end-1,1:end-1)),squeeze(cA_001(8,:,1:end-1,2:end)));
cA8_2 = cat(3,squeeze(cA_002(8,:,2:end,1:end-1)),squeeze(cA_003(8,:,2:end,2:end)));
cA8 = cat(2,cA8_1,cA8_2); clear cA8_1 cA8_2;

%
i_p = 60;
figure; 
subplot(3,3,1); imagesc(squeeze(cA6(:,:,i_p+1))); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA7(:,:,i_p+1))); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA1(:,:,i_p+1))); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA8(:,:,i_p+1))); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA2(:,:,i_p+1))); axis xy; colorbar
figure; 
subplot(3,3,1); imagesc(squeeze(cA3(:,:,i_p+1))); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA4(:,:,i_p+1))); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA1(:,:,i_p+1))); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA5(:,:,i_p+1))); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA2(:,:,i_p+1))); axis xy; colorbar

%
j_p = 60;
figure; 
subplot(3,3,1); imagesc(squeeze(cA6(:,j_p+1,:))); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA7(:,j_p+1,:))); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA1(:,j_p+1,:))); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA8(:,j_p+1,:))); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA2(:,j_p+1,:))); axis xy; colorbar
figure; 
subplot(3,3,1); imagesc(squeeze(cA3(:,j_p+1,:))); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA4(:,j_p+1,:))); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA1(:,j_p+1,:))); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA5(:,j_p+1,:))); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA2(:,j_p+1,:))); axis xy; colorbar

%
k_p = 1;
figure; 
subplot(3,3,1); imagesc(squeeze(cA6(k_p,:,:))); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA7(k_p,:,:))); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA1(k_p,:,:))); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA8(k_p,:,:))); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA2(k_p,:,:))); axis xy; colorbar
figure; 
subplot(3,3,1); imagesc(squeeze(cA3(k_p,:,:))); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA4(k_p,:,:))); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA1(k_p,:,:))); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA5(k_p,:,:))); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA2(k_p,:,:))); axis xy; colorbar

%% pcolor sections
j = 60;

zr = cat(2,squeeze(zr_000(:,j+1,2:end-1)),squeeze(zr_001(:,j+1,2:end-1)));
zw = cat(2,squeeze(zw_000(:,j+1,2:end-1)),squeeze(zw_001(:,j+1,2:end-1)));
rhs = cat(2,squeeze(rhs_000(:,j+1,2:end-1)),squeeze(rhs_001(:,j+1,2:end-1)));
p = cat(2,squeeze(p_000(:,j+1,2:end-1)),squeeze(p_001(:,j+1,2:end-1)));
r = cat(2,squeeze(r_000(:,j+1,2:end-1)),squeeze(r_001(:,j+1,2:end-1)));

xr = repmat(((0.5:1:nx)*L/nx),[nz 1]);
xu_pc = cat(2,zeros(nz,1),0.5*(xr(:,1:end-1)+xr(:,2:end)),L*ones(nz,1));
xu_pc = cat(1,xu_pc,xu_pc(end,:));
zwu_pc = cat(2,zw(:,1),0.5*(zw(:,1:end-1)+zw(:,2:end)),zw(:,end));
rhs_pc = cat(2,rhs,rhs(:,end));
rhs_pc = cat(1,rhs_pc,rhs_pc(end,:));
p_pc = cat(2,p,p(:,end));
p_pc = cat(1,p_pc,p_pc(end,:));
r_pc = cat(2,r,r(:,end));
r_pc = cat(1,r_pc,r_pc(end,:));

% figure;
% pcolor(xu_pc,zwu_pc,rhs_pc); shading flat;
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(rhs(:));
% cmax = max(rhs(:));
% caxis([cmin cmax]);
% title(['rhs  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

figure;
subplot(2,1,1)
pcolor(xu_pc,zwu_pc,p_pc); shading flat;
hold on
contour(xr,zr,p,[-14:1:14]*1e4,'k')
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p(:));
cmax = max(p(:));
caxis([cmin cmax]);
title(['pressure mgroms j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
subplot(2,1,2)
pcolor(xu_pc,zwu_pc,r_pc); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(r(:));
cmax = max(r(:));
caxis([cmin cmax]);
title(['residual mgroms j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

% p_ND_pc = cat(2,p_ND,p_ND(:,end));
% p_ND_pc = cat(1,p_ND_pc,p_ND_pc(end,:));
% 
% figure;
% subplot(2,1,1)
% pcolor(xu_pc,zwu_pc,p_pc-p_ND_pc); shading flat;
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(p(:)-p_ND(:));
% cmax = max(p(:)-p_ND(:));
% caxis([cmin cmax]);
% title(['pressure error mgroms-ND  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
% subplot(2,1,2)
% pcolor(xu_pc,zwu_pc,log10(abs(p_pc-p_ND_pc)./abs(p_ND_pc))); shading flat;
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(log10(abs(p_pc(:)-p_ND_pc(:))./abs(p_ND_pc(:))));
% cmax = max(log10(abs(p_pc(:)-p_ND_pc(:))./abs(p_ND_pc(:))));
% caxis([cmin cmax]);
% title(['relative pressure error mgroms-ND  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

p_GR_pc = cat(2,p_GR,p_GR(:,end));
p_GR_pc = cat(1,p_GR_pc,p_GR_pc(end,:));

figure;
subplot(2,1,1)
pcolor(xu_pc,zwu_pc,p_pc-p_GR_pc); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p(:)-p_GR(:));
cmax = max(p(:)-p_GR(:));
caxis([cmin cmax]);
title(['pressure error mgroms-GR  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
subplot(2,1,2)
pcolor(xu_pc,zwu_pc,log10(abs(p_pc-p_GR_pc)./abs(p_GR_pc))); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(log10(abs(p_pc(:)-p_GR_pc(:))./abs(p_GR_pc(:))));
cmax = max(log10(abs(p_pc(:)-p_GR_pc(:))./abs(p_GR_pc(:))));
caxis([cmin cmax]);
title(['relative pressure error mgroms-GR  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

p_JM_pc = cat(2,p_JM,p_JM(:,end));
p_JM_pc = cat(1,p_JM_pc,p_JM_pc(end,:));

figure;
subplot(2,1,1)
pcolor(xu_pc,zwu_pc,p_pc-p_JM_pc); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p(:)-p_JM(:));
cmax = max(p(:)-p_JM(:));
caxis([cmin cmax]);
title(['pressure error mgroms-JM  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
subplot(2,1,2)
pcolor(xu_pc,zwu_pc,log10(abs(p_pc-p_JM_pc)./abs(p_JM_pc))); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(log10(abs(p_pc(:)-p_JM_pc(:))./abs(p_JM_pc(:))));
cmax = max(log10(abs(p_pc(:)-p_JM_pc(:))./abs(p_JM_pc(:))));
caxis([cmin cmax]);
title(['relative pressure error mgroms-JM  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

%
% figure;
% subplot(2,1,1)
% pcolor(xu_pc,zwu_pc,p_ND_pc-p_GR_pc); shading flat;
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(p_ND(:)-p_GR(:));
% cmax = max(p_ND(:)-p_GR(:));
% caxis([cmin cmax]);
% title(['pressure error ND-GR  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
% 
% subplot(2,1,2)
% pcolor(xu_pc,zwu_pc,log10(abs(p_ND_pc-p_GR_pc)./abs(p_GR_pc))); shading flat;
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(log10(abs(p_ND_pc(:)-p_GR_pc(:))./abs(p_GR_pc(:))));
% cmax = max(log10(abs(p_ND_pc(:)-p_GR_pc(:))./abs(p_GR_pc(:))));
% caxis([cmin cmax]);
% title(['relative pressure error ND-GR  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

% figure;
% subplot(2,1,1)
% pcolor(xu_pc,zwu_pc,p_ND_pc-p_JM_pc); shading flat;
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(p_ND(:)-p_JM(:));
% cmax = max(p_ND(:)-p_JM(:));
% caxis([cmin cmax]);
% title(['pressure error ND-JM  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
% 
% subplot(2,1,2)
% pcolor(xu_pc,zwu_pc,log10(abs(p_ND_pc-p_JM_pc)./abs(p_JM_pc))); shading flat;
% log10(abs(p_ND_pc-p_GR_pc)./abs(p_GR_pc))
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(log10(abs(p_ND_pc(:)-p_JM_pc(:))./abs(p_JM_pc(:)))); cmin(isinf(cmin)) = -15;
% cmax = max(log10(abs(p_ND_pc(:)-p_JM_pc(:))./abs(p_JM_pc(:))));
% caxis([cmin cmax]);
% title(['relative pressure error ND-JM  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

% figure;
% subplot(2,1,1)
% pcolor(xu_pc,zwu_pc,p_GR_pc-p_JM_pc); shading flat;
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(p_GR(:)-p_JM(:));
% cmax = max(p_GR(:)-p_JM(:));
% caxis([cmin cmax]);
% title(['pressure error GR-JM  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
% subplot(2,1,2)
% pcolor(xu_pc,zwu_pc,log10(abs(p_GR_pc-p_JM_pc)./abs(p_JM_pc))); shading flat;
% hold on
% plot(xr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 L])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(log10(abs(p_GR_pc(:)-p_JM_pc(:))./abs(p_JM_pc(:))));
% cmax = max(log10(abs(p_GR_pc(:)-p_JM_pc(:))./abs(p_JM_pc(:))));
% caxis([cmin cmax]);
% title(['relative pressure error GR-JM  j=' num2str(j) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

% % pcolor sections
% i = 60;
% 
% zr = cat(2,squeeze(zr_000(:,2:end-1,i+1)),squeeze(zr_002(:,2:end-1,i+1)));
% zw = cat(2,squeeze(zw_000(:,2:end-1,i+1)),squeeze(zw_002(:,2:end-1,i+1)));
% rhs = cat(2,squeeze(rhs_000(:,2:end-1,i+1)),squeeze(rhs_002(:,2:end-1,i+1)));
% p = cat(2,squeeze(p_000(:,2:end-1,i+1)),squeeze(p_002(:,2:end-1,i+1)));
% r = cat(2,squeeze(r_000(:,2:end-1,i+1)),squeeze(r_002(:,2:end-1,i+1)));
% 
% yr = repmat(((0.5:1:ny)*Ly/ny),[nz 1]);
% yv_pc = cat(2,zeros(nz,1),0.5*(yr(:,1:end-1)+yr(:,2:end)),Ly*ones(nz,1));
% yv_pc = cat(1,yv_pc,yv_pc(end,:));
% zwv_pc = cat(2,zw(:,1),0.5*(zw(:,1:end-1)+zw(:,2:end)),zw(:,end));
% rhs_pc = cat(2,rhs,rhs(:,end));
% rhs_pc = cat(1,rhs_pc,rhs_pc(end,:));
% p_pc = cat(2,p,p(:,end));
% p_pc = cat(1,p_pc,p_pc(end,:));
% r_pc = cat(2,r,r(:,end));
% r_pc = cat(1,r_pc,r_pc(end,:));
% 
% figure;
% pcolor(yv_pc,zwv_pc,p_pc); shading flat;
% hold on
% contour(yr,zr,p,[-14:1:14]*1e4,'k')
% plot(yr(1,:),zw(1,:),'k','linewidth',1)
% xlim([0 Ly])
% ylim([-H 0])
% axis equal tight
% colorbar
% cmin = min(p(:));
% cmax = max(p(:));
% caxis([cmin cmax]);
% title(['pressure  mgroms i=' num2str(i) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

% imagesc maps
k = 1;

p_1 = cat(2,squeeze(p_000(k,2:end-1,2:end-1)),squeeze(p_001(k,2:end-1,2:end-1)));
p_2 = cat(2,squeeze(p_002(k,2:end-1,2:end-1)),squeeze(p_003(k,2:end-1,2:end-1)));
p = cat(1,p_1,p_2); clear p_1 p_2;
r_1 = cat(2,squeeze(r_000(k,2:end-1,2:end-1)),squeeze(r_001(k,2:end-1,2:end-1)));
r_2 = cat(2,squeeze(r_002(k,2:end-1,2:end-1)),squeeze(r_003(k,2:end-1,2:end-1)));
r = cat(1,r_1,r_2); clear r_1 r_2;

figure;
subplot(2,1,1)
imagesc(p); axis xy
axis equal tight
colorbar
cmin = min(min(p(:)));
cmax = max(max(p(:)));
caxis([cmin cmax]);
title(['pressure mgroms  k=' num2str(k) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
subplot(2,1,2)
imagesc(r); axis xy
axis equal tight
colorbar
cmin = min(r(:));
cmax = max(r(:));
caxis([cmin cmax]);
title(['residual mgroms k=' num2str(k) ' [' num2str(cmin) ' ' num2str(cmax) ']']);

% 
% figure;
% imagesc(p-repmat(p_GR(k,:),[nx 1])); axis xy
% colorbar
% cmin = min(min(p-repmat(p_GR(k,:),[nx 1])));
% cmax = max(max(p-repmat(p_GR(k,:),[nx 1])));
% caxis([cmin cmax]);
% title(['pressure error  k=' num2str(k) ' [' num2str(cmin) ' ' num2str(cmax) ']']);
