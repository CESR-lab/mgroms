clear all

nx = 128; 
nz = 128;

%% Jeroen's matrix
[A_JM,cec_JM,cwc_JM,cct_JM,ccb_JM,cet_JM,cwt_JM,ceb_JM,cwb_JM,ccc_JM,h_JM,xr_JM,zr_JM,vr_JM,zw_JM] ...
    = set_nhlap_JM_ND(nx,nz);

%% Guillaume's matrix
[A_GR,ce_GR,cw_GR,cn_GR,cs_GR,cen_GR,cwn_GR,ces_GR,cws_GR,co_GR,h_GR,xr_GR,zr_GR,vr_GR,zw_GR] ...
    = set_nhlap_GR_ND(nx,nz);

%% mgroms matrix
cA_000 = ncread('../src/cA_cA_000.nc','cA');
cA_001 = ncread('../src/cA_cA_001.nc','cA');
cA_002 = ncread('../src/cA_cA_002.nc','cA');
cA_003 = ncread('../src/cA_cA_003.nc','cA');
cA = cat(3,squeeze(cA_000(:,:,10,2:end-1)),squeeze(cA_001(:,:,10,2:end-1)));
%cA = cat(3,squeeze(cA_000(:,:,10,1:end-1)),squeeze(cA_001(:,:,10,2:end)));

figure; 
subplot(3,3,1); imagesc(squeeze(cA(6,:,:))); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA(7,:,:))); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA(1,:,:))); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA(8,:,:))); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA(2,:,:))); axis xy; colorbar

% figure; 
% subplot(3,3,1); imagesc(squeeze(cA(6,:,:))-cwn_GR,[-1.5 1.5]*1e-3); axis xy; colorbar
% subplot(3,3,4); imagesc(squeeze(cA(7,:,:))-cw_GR,[-1.5 1.5]*1e-3); axis xy; colorbar
% subplot(3,3,5); imagesc(squeeze(cA(1,:,:))-co_GR,[-1.5 1.5]*1e-3); axis xy; colorbar
% subplot(3,3,7); imagesc(squeeze(cA(8,:,:))-cws_GR,[-1.5 1.5]*1e-3); axis xy; colorbar
% subplot(3,3,8); imagesc(squeeze(cA(2,:,:))-cs_GR,[-1.5 1.5]*1e-3); axis xy; colorbar
% figure; 
% subplot(3,3,1); imagesc(squeeze(cA(6,:,:))-cwt_JM',[-1.5 1.5]*1e-3); axis xy; colorbar
% subplot(3,3,4); imagesc(squeeze(cA(7,:,:))-cwc_JM',[-1.5 1.5]*1e-3); axis xy; colorbar
% subplot(3,3,5); imagesc(squeeze(cA(1,:,:))-ccc_JM',[-1.5 1.5]*1e-3); axis xy; colorbar
% subplot(3,3,7); imagesc(squeeze(cA(8,:,:))-cwb_JM',[-1.5 1.5]*1e-3); axis xy; colorbar
% subplot(3,3,8); imagesc(squeeze(cA(2,:,:))-ccb_JM',[-1.5 1.5]*1e-3); axis xy; colorbar

% figure; 
% k_p = 10;
% subplot(3,3,1); plot(squeeze(cA(6,k_p,:)),'*'); hold on; plot(squeeze(cwt_JM(:,k_p)),'o');
% subplot(3,3,4); plot(squeeze(cA(7,k_p,:)),'*'); hold on; plot(squeeze(cwc_JM(:,k_p)),'o');
% subplot(3,3,5); plot(squeeze(cA(1,k_p,:)),'*'); hold on; plot(squeeze(ccc_JM(:,k_p)),'o');
% subplot(3,3,7); plot(squeeze(cA(8,k_p,:)),'*'); hold on; plot(squeeze(cwb_JM(:,k_p)),'o');
% subplot(3,3,8); plot(squeeze(cA(2,k_p,:)),'*'); hold on; plot(squeeze(ccb_JM(:,k_p)),'o');
% figure; 
% subplot(3,3,1); plot(squeeze(cA(6,k_p,:))-squeeze(cwt_JM(:,k_p)),'*'); 
% subplot(3,3,4); plot(squeeze(cA(7,k_p,:))-squeeze(cwc_JM(:,k_p)),'*'); 
% subplot(3,3,5); plot(squeeze(cA(1,k_p,:))-squeeze(ccc_JM(:,k_p)),'*'); 
% subplot(3,3,7); plot(squeeze(cA(8,k_p,:))-squeeze(cwb_JM(:,k_p)),'*'); 
% subplot(3,3,8); plot(squeeze(cA(2,k_p,:))-squeeze(ccb_JM(:,k_p)),'*');

% figure; 
% i_p = 48;
% subplot(3,3,1); plot(squeeze(cA(6,:,i_p)),1:nz,'*'); hold on; plot(squeeze(cwt_JM(i_p,:)),1:nz,'o');
% subplot(3,3,4); plot(squeeze(cA(7,:,i_p)),1:nz,'*'); hold on; plot(squeeze(cwc_JM(i_p,:)),1:nz,'o');
% subplot(3,3,5); plot(squeeze(cA(1,:,i_p)),1:nz,'*'); hold on; plot(squeeze(ccc_JM(i_p,:)),1:nz,'o');
% subplot(3,3,7); plot(squeeze(cA(8,:,i_p)),1:nz,'*'); hold on; plot(squeeze(cwb_JM(i_p,:)),1:nz,'o');
% subplot(3,3,8); plot(squeeze(cA(2,:,i_p)),1:nz,'*'); hold on; plot(squeeze(ccb_JM(i_p,:)),1:nz,'o');
% figure; 
% subplot(3,3,1); plot(squeeze(cA(6,:,i_p))-squeeze(cwt_JM(i_p,:)),1:nz,'*'); 
% subplot(3,3,4); plot(squeeze(cA(7,:,i_p))-squeeze(cwc_JM(i_p,:)),1:nz,'*'); 
% subplot(3,3,5); plot(squeeze(cA(1,:,i_p))-squeeze(ccc_JM(i_p,:)),1:nz,'*'); 
% subplot(3,3,7); plot(squeeze(cA(8,:,i_p))-squeeze(cwb_JM(i_p,:)),1:nz,'*'); 
% subplot(3,3,8); plot(squeeze(cA(2,:,i_p))-squeeze(ccb_JM(i_p,:)),1:nz,'*'); 

% translate to Guillaume's notation
ce_mgroms = zeros(nz,nx);
cw_mgroms = zeros(nz,nx);
cn_mgroms = zeros(nz,nx);
cs_mgroms = zeros(nz,nx);
cen_mgroms = zeros(nz,nx);
cwn_mgroms = zeros(nz,nx);
ces_mgroms = zeros(nz,nx);
cws_mgroms = zeros(nz,nx);
co_mgroms = zeros(nz,nx);

i = 1:nx; im = 1:nx-1; ip = 2:nx;
k = 1:nz; km = 1:nz-1; kp = 2:nz;
siz = [nz,nx];

ce_mgroms(k,im) = ce_mgroms(k,im) + squeeze(cA(7,k,ip));
cw_mgroms(k,ip) = cw_mgroms(k,ip) + squeeze(cA(7,k,ip));
cn_mgroms(km,i) = cn_mgroms(km,i) + squeeze(cA(2,kp,i));
cs_mgroms(kp,i) = cs_mgroms(kp,i) + squeeze(cA(2,kp,i));
cen_mgroms(km,im) = cen_mgroms(km,im) + squeeze(cA(8,kp,ip));
cwn_mgroms(km,ip) = cwn_mgroms(km,ip) + squeeze(cA(6,km,ip));
cws_mgroms(kp,ip) = cws_mgroms(kp,ip) + squeeze(cA(8,kp,ip));
ces_mgroms(kp,im) = ces_mgroms(kp,im) + squeeze(cA(6,km,ip));
co_mgroms(k,i) = co_mgroms(k,i) + squeeze(cA(1,k,i));

figure; 
subplot(3,3,1); imagesc(cwn_mgroms); axis xy; colorbar
subplot(3,3,2); imagesc(cn_mgroms); axis xy; colorbar
subplot(3,3,3); imagesc(cen_mgroms); axis xy; colorbar
subplot(3,3,4); imagesc(cw_mgroms); axis xy; colorbar
subplot(3,3,5); imagesc(co_mgroms); axis xy; colorbar
subplot(3,3,6); imagesc(ce_mgroms); axis xy; colorbar
subplot(3,3,7); imagesc(cws_mgroms); axis xy; colorbar
subplot(3,3,8); imagesc(cs_mgroms); axis xy; colorbar
subplot(3,3,9); imagesc(ces_mgroms); axis xy; colorbar
figure; 
subplot(3,3,1); imagesc(cwn_mgroms-cwn_GR); axis xy; colorbar
subplot(3,3,2); imagesc(cn_mgroms-cn_GR); axis xy; colorbar
subplot(3,3,3); imagesc(cen_mgroms-cen_GR); axis xy; colorbar
subplot(3,3,4); imagesc(cw_mgroms-cw_GR); axis xy; colorbar
subplot(3,3,5); imagesc(co_mgroms-co_GR); axis xy; colorbar
subplot(3,3,6); imagesc(ce_mgroms-ce_GR); axis xy; colorbar
subplot(3,3,7); imagesc(cws_mgroms-cws_GR); axis xy; colorbar
subplot(3,3,8); imagesc(cs_mgroms-cs_GR); axis xy; colorbar
subplot(3,3,9); imagesc(ces_mgroms-ces_GR); axis xy; colorbar

% figure; 
% k_p = 1;
% subplot(3,3,1); plot(squeeze(cwn_mgroms(k_p,:)),'*'); hold on; plot(squeeze(cwn_GR(k_p,:)),'o');
% subplot(3,3,2); plot(squeeze(cn_mgroms(k_p,:)),'*'); hold on; plot(squeeze(cn_GR(k_p,:)),'o');
% subplot(3,3,3); plot(squeeze(cen_mgroms(k_p,:)),'*'); hold on; plot(squeeze(cen_GR(k_p,:)),'o');
% subplot(3,3,4); plot(squeeze(cw_mgroms(k_p,:)),'*'); hold on; plot(squeeze(cw_GR(k_p,:)),'o');
% subplot(3,3,5); plot(squeeze(co_mgroms(k_p,:)),'*'); hold on; plot(squeeze(co_GR(k_p,:)),'o');
% subplot(3,3,6); plot(squeeze(ce_mgroms(k_p,:)),'*'); hold on; plot(squeeze(ce_GR(k_p,:)),'o');
% subplot(3,3,7); plot(squeeze(cws_mgroms(k_p,:)),'*'); hold on; plot(squeeze(cws_GR(k_p,:)),'*');
% subplot(3,3,8); plot(squeeze(cs_mgroms(k_p,:)),'*'); hold on; plot(squeeze(cs_GR(k_p,:)),'o');
% subplot(3,3,9); plot(squeeze(ces_mgroms(k_p,:)),'*'); hold on; plot(squeeze(ces_GR(k_p,:)),'o');
% figure; 
% subplot(3,3,1); plot(squeeze(cwn_mgroms(k_p,:))-squeeze(cwn_GR(k_p,:)),'*'); 
% subplot(3,3,2); plot(squeeze(cn_mgroms(k_p,:))-squeeze(cn_GR(k_p,:)),'*'); 
% subplot(3,3,3); plot(squeeze(cen_mgroms(k_p,:))-squeeze(cen_GR(k_p,:)),'*'); 
% subplot(3,3,4); plot(squeeze(cw_mgroms(k_p,:))-squeeze(cw_GR(k_p,:)),'*'); 
% subplot(3,3,5); plot(squeeze(co_mgroms(k_p,:))-squeeze(co_GR(k_p,:)),'*'); 
% subplot(3,3,6); plot(squeeze(ce_mgroms(k_p,:))-squeeze(ce_GR(k_p,:)),'*'); 
% subplot(3,3,7); plot(squeeze(cws_mgroms(k_p,:))-squeeze(cws_GR(k_p,:)),'*'); 
% subplot(3,3,8); plot(squeeze(cs_mgroms(k_p,:))-squeeze(cs_GR(k_p,:)),'*');
% subplot(3,3,9); plot(squeeze(ces_mgroms(k_p,:))-squeeze(ces_GR(k_p,:)),'*');

% arranged in a sparse matrix
I=[];
J=[];
s=[];
% east
c=ce_mgroms(k,im);
[ii,kk]=meshgrid(ip,k);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% west
c=cw_mgroms(k,ip);
[ii,kk]=meshgrid(ip,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,k);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% north
c=cn_mgroms(km,i);
[ii,kk]=meshgrid(i,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(i,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% south
c=cs_mgroms(kp,i);
[ii,kk]=meshgrid(i,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(i,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% central
c=co_mgroms(k,i);
[ii,kk]=meshgrid(i,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;I0];
s=[s;c(:)];
% northeast
c=cen_mgroms(km,im);
[ii,kk]=meshgrid(ip,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% southeast
c=ces_mgroms(kp,im);
[ii,kk]=meshgrid(ip,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% northwest
c=cwn_mgroms(km,ip);
[ii,kk]=meshgrid(im,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(ip,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% southwest
c=cws_mgroms(kp,ip);
[ii,kk]=meshgrid(im,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(ip,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% A
A_mgroms=sparse(I,J,s,nx*nz,nx*nz);

%% define rhs
L = 10e3;
H = 4e3;
bet = 600 / (L*L);
x1 = L * 0.65;
z1 = H * (0.75 - 1);
x2 = L * 0.75;
z2 = H * (0.65 - 1);

rhs_GR = vr_GR .* (exp(-bet * ((xr_GR-x1).^2 + (zr_GR-z1).^2)) ...
             -  exp(-bet * ((xr_GR-x2).^2 + (zr_GR-z2).^2)));
rhs_JM = vr_JM .* (exp(-bet * ((xr_JM-x1).^2 + (zr_JM-z1).^2)) ...
             -  exp(-bet * ((xr_JM-x2).^2 + (zr_JM-z2).^2)));
         
%% solve for p
p_GR =  A_GR \ rhs_GR(:);
p_GR = reshape(p_GR,nz,nx);

rhs2_JM = cat(2,zeros(nx,1),rhs_JM);
p_JM = A_JM \ rhs2_JM(:);
p_JM = reshape(p_JM,[nx nz+1]);
p_JM = p_JM(:,2:end);
p_JM = p_JM';

p_mgroms =  A_mgroms \ rhs_GR(:);
p_mgroms = reshape(p_mgroms,nz,nx);

%% pcolor plot
xu_GR_pc = cat(2,zeros(nz,1),0.5*(xr_GR(:,1:end-1)+xr_GR(:,2:end)),L*ones(nz,1));
xu_GR_pc = cat(1,xu_GR_pc,xu_GR_pc(end,:));
zwu_GR_pc = cat(2,zw_GR(:,1),0.5*(zw_GR(:,1:end-1)+zw_GR(:,2:end)),zw_GR(:,end));

xu_JM_pc = cat(2,zeros(nz,1),0.5*(xr_JM(:,1:end-1)+xr_JM(:,2:end)),L*ones(nz,1));
xu_JM_pc = cat(1,xu_JM_pc,xu_JM_pc(end,:));
zwu_JM_pc = cat(2,zw_JM(:,1),0.5*(zw_JM(:,1:end-1)+zw_JM(:,2:end)),zw_JM(:,end));

p_GR_pc = cat(2,p_GR,p_GR(:,end));
p_GR_pc = cat(1,p_GR_pc,p_GR_pc(end,:));
p_JM_pc = cat(2,p_JM,p_JM(:,end));
p_JM_pc = cat(1,p_JM_pc,p_JM_pc(end,:));
p_mgroms_pc = cat(2,p_mgroms,p_mgroms(:,end));
p_mgroms_pc = cat(1,p_mgroms_pc,p_mgroms_pc(end,:));

figure;
pcolor(xu_GR_pc,zwu_GR_pc,p_GR_pc); shading flat;
hold on
contour(xr_GR,zr_GR,p_GR,[-14:1:14]*1e4,'k')
plot(xr_GR(1,:),zw_GR(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p_GR(:));
cmax = max(p_GR(:));
caxis([cmin cmax]);
title(['pressure GR  [' num2str(cmin) ' ' num2str(cmax) ']']);

figure;
pcolor(xu_GR_pc,zwu_GR_pc,p_JM_pc); shading flat;
hold on
contour(xr_GR,zr_GR,p_JM,[-14:1:14]*1e4,'k')
plot(xr_GR(1,:),zw_GR(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p_JM(:));
cmax = max(p_JM(:));
caxis([cmin cmax]);
title(['pressure JM  [' num2str(cmin) ' ' num2str(cmax) ']']);

figure;
pcolor(xu_GR_pc,zwu_GR_pc,p_mgroms_pc); shading flat;
hold on
contour(xr_GR,zr_GR,p_mgroms,[-14:1:14]*1e4,'k')
plot(xr_GR(1,:),zw_GR(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p_mgroms(:));
cmax = max(p_mgroms(:));
caxis([cmin cmax]);
title(['pressure mgroms  [' num2str(cmin) ' ' num2str(cmax) ']']);

figure
subplot(2,1,1)
pcolor(xu_GR_pc,zwu_GR_pc,p_mgroms_pc-p_GR_pc); shading flat;
hold on
plot(xr_GR(1,:),zw_GR(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p_mgroms(:)-p_GR(:));
cmax = max(p_mgroms(:)-p_GR(:));
caxis([cmin cmax]);
title(['pressure error mgroms-GR  [' num2str(cmin) ' ' num2str(cmax) ']']);

subplot(2,1,2)
pcolor(xu_GR_pc,zwu_GR_pc,log10(abs(p_mgroms_pc-p_GR_pc)./abs(p_GR_pc))); shading flat;
hold on
plot(xr_GR(1,:),zw_GR(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(log10(abs(p_mgroms_pc(:)-p_GR_pc(:))./abs(p_GR_pc(:))));
cmax = max(log10(abs(p_mgroms_pc(:)-p_GR_pc(:))./abs(p_GR_pc(:))));
caxis([cmin cmax]);
title(['relative pressure error mgroms-GR  [' num2str(cmin) ' ' num2str(cmax) ']']);

figure;
subplot(2,1,1)
pcolor(xu_GR_pc,zwu_GR_pc,p_mgroms_pc-p_JM_pc); shading flat;
hold on
plot(xr_GR(1,:),zw_GR(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(p_mgroms_pc(:)-p_JM_pc(:));
cmax = max(p_mgroms_pc(:)-p_JM_pc(:));
caxis([cmin cmax]);
title(['pressure error mgroms-JM  [' num2str(cmin) ' ' num2str(cmax) ']']);

subplot(2,1,2)
pcolor(xu_GR_pc,zwu_GR_pc,log10(abs(p_mgroms_pc-p_JM_pc)./abs(p_JM_pc))); shading flat;
hold on
plot(xr_GR(1,:),zw_GR(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(log10(abs(p_mgroms_pc(:)-p_JM_pc(:))./abs(p_JM_pc(:))));
cmax = max(log10(abs(p_mgroms_pc(:)-p_JM_pc(:))./abs(p_JM_pc(:))));
caxis([cmin cmax]);
title(['relative pressure error mgroms-JM  [' num2str(cmin) ' ' num2str(cmax) ']'])

