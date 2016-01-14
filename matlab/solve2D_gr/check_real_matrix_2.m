clear all

%% 

nx = 128; 
nz = 128;

% Guillaume's matrix

[A_GR,ce_GR,cw_GR,cn_GR,cs_GR,cen_GR,cwn_GR,ces_GR,cws_GR,co_GR,h,xr,zr,vr,zw,dzr,dzu,dzw,dxu,dxw,alphauw,alphaw] = set_nhlap_ND(nx,nz);

figure; 
subplot(3,3,1); imagesc(cwn_GR); axis xy; colorbar
subplot(3,3,2); imagesc(cn_GR); axis xy; colorbar
subplot(3,3,3); imagesc(cen_GR); axis xy; colorbar
subplot(3,3,4); imagesc(cw_GR); axis xy; colorbar
subplot(3,3,5); imagesc(co_GR); axis xy; colorbar
subplot(3,3,6); imagesc(ce_GR); axis xy; colorbar
subplot(3,3,7); imagesc(cws_GR); axis xy; colorbar
subplot(3,3,8); imagesc(cs_GR); axis xy; colorbar
subplot(3,3,9); imagesc(ces_GR); axis xy; colorbar
figure;
imagesc(A_GR-A_GR'); colorbar

% Jeroen's matrix coefficients
cA_000 = ncread('../../src/cA_cA_000.nc','cA');
cA_001 = ncread('../../src/cA_cA_001.nc','cA');
cA_002 = ncread('../../src/cA_cA_002.nc','cA');
cA_003 = ncread('../../src/cA_cA_003.nc','cA');

j = 10;
cA = cat(3,squeeze(cA_000(:,:,j+1,2:end-1)),squeeze(cA_001(:,:,j+1,2:end-1)));

figure; 
subplot(3,3,1); imagesc(squeeze(cA(6,:,:))); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA(7,:,:))); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA(1,:,:))); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA(8,:,:))); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA(2,:,:))); axis xy; colorbar
figure; 
subplot(3,3,1); imagesc(squeeze(cA(6,:,:))-cwn_GR); axis xy; colorbar
subplot(3,3,4); imagesc(squeeze(cA(7,:,:))-cw_GR); axis xy; colorbar
subplot(3,3,5); imagesc(squeeze(cA(1,:,:))-co_GR); axis xy; colorbar
subplot(3,3,7); imagesc(squeeze(cA(8,:,:))-cws_GR); axis xy; colorbar
subplot(3,3,8); imagesc(squeeze(cA(2,:,:))-cs_GR); axis xy; colorbar

figure; 
k = 1;
subplot(3,3,1); plot(squeeze(cA(6,k,:)),'*'); hold on; plot(squeeze(cwn_GR(k,:)),'o');
subplot(3,3,4); plot(squeeze(cA(7,k,:)),'*'); hold on; plot(squeeze(cw_GR(k,:)),'o');
subplot(3,3,5); plot(squeeze(cA(1,k,:)),'*'); hold on; plot(squeeze(co_GR(k,:)),'o');
subplot(3,3,7); plot(squeeze(cA(8,k,:)),'*'); hold on; plot(squeeze(cws_GR(k,:)),'o');
subplot(3,3,8); plot(squeeze(cA(2,k,:)),'*'); hold on; plot(squeeze(cs_GR(k,:)),'o');
figure; 
subplot(3,3,1); plot(squeeze(cA(6,k,:))-squeeze(cwn_GR(k,:))','*'); 
subplot(3,3,4); plot(squeeze(cA(7,k,:))-squeeze(cw_GR(k,:))','*'); 
subplot(3,3,5); plot(squeeze(cA(1,k,:))-squeeze(co_GR(k,:))','*'); 
subplot(3,3,7); plot(squeeze(cA(8,k,:))-squeeze(cws_GR(k,:))','*'); 
subplot(3,3,8); plot(squeeze(cA(2,k,:))-squeeze(cs_GR(k,:))','*'); 
figure; 
i = 48;
subplot(3,3,1); plot(squeeze(cA(6,:,i)),1:nz,'*'); hold on; plot(squeeze(cwn_GR(:,i)),1:nz,'o');
subplot(3,3,4); plot(squeeze(cA(7,:,i)),1:nz,'*'); hold on; plot(squeeze(cw_GR(:,i)),1:nz,'o');
subplot(3,3,5); plot(squeeze(cA(1,:,i)),1:nz,'*'); hold on; plot(squeeze(co_GR(:,i)),1:nz,'o');
subplot(3,3,7); plot(squeeze(cA(8,:,i)),1:nz,'*'); hold on; plot(squeeze(cws_GR(:,i)),1:nz,'o');
subplot(3,3,8); plot(squeeze(cA(2,:,i)),1:nz,'*'); hold on; plot(squeeze(cs_GR(:,i)),1:nz,'o');
figure; 
subplot(3,3,1); plot(squeeze(cA(6,:,i))-squeeze(cwn_GR(:,i))',1:nz,'*'); 
subplot(3,3,4); plot(squeeze(cA(7,:,i))-squeeze(cw_GR(:,i))',1:nz,'*'); 
subplot(3,3,5); plot(squeeze(cA(1,:,i))-squeeze(co_GR(:,i))',1:nz,'*'); 
subplot(3,3,7); plot(squeeze(cA(8,:,i))-squeeze(cws_GR(:,i))',1:nz,'*'); 
subplot(3,3,8); plot(squeeze(cA(2,:,i))-squeeze(cs_GR(:,i))',1:nz,'*'); 

% translate to Guillaume's notation
ce_JM = zeros(nz,nx);
cw_JM = zeros(nz,nx);
cn_JM = zeros(nz,nx);
cs_JM = zeros(nz,nx);
cen_JM = zeros(nz,nx);
cwn_JM = zeros(nz,nx);
ces_JM = zeros(nz,nx);
cws_JM = zeros(nz,nx);
co_JM = zeros(nz,nx);

i=1:nx;im=1:nx-1;ip=2:nx;
k=1:nz;km=1:nz-1;kp=2:nz;
siz=[nz,nx];

ce_JM(k,im) = ce_JM(k,im) + squeeze(cA(7,k,ip));
cw_JM(k,ip) = cw_JM(k,ip) + squeeze(cA(7,k,ip));
cn_JM(km,i) = cn_JM(km,i) + squeeze(cA(2,kp,i));
cs_JM(kp,i) = cs_JM(kp,i) + squeeze(cA(2,kp,i));
cen_JM(km,im) = cen_JM(km,im) + squeeze(cA(8,kp,ip));
cwn_JM(km,ip) = cwn_JM(km,ip) + squeeze(cA(6,km,ip));
cws_JM(kp,ip) = cws_JM(kp,ip) + squeeze(cA(8,kp,ip));
ces_JM(kp,im) = ces_JM(kp,im) + squeeze(cA(6,km,ip));
%co_JM(k,i) = co_JM(k,i) + squeeze(cA(1,k,i));
co_JM(km,i) = -ce_JM(km,i)-cw_JM(km,i)-cn_JM(km,i)-cs_JM(km,i);
co_JM(nz,i) = -ce_JM(nz,i)-cw_JM(nz,i)-3*cs_JM(nz,i);

figure; 
subplot(3,3,1); imagesc(cwn_JM); axis xy; colorbar
subplot(3,3,2); imagesc(cn_JM); axis xy; colorbar
subplot(3,3,3); imagesc(cen_JM); axis xy; colorbar
subplot(3,3,4); imagesc(cw_JM); axis xy; colorbar
subplot(3,3,5); imagesc(co_JM); axis xy; colorbar
subplot(3,3,6); imagesc(ce_JM); axis xy; colorbar
subplot(3,3,7); imagesc(cws_JM); axis xy; colorbar
subplot(3,3,8); imagesc(cs_JM); axis xy; colorbar
subplot(3,3,9); imagesc(ces_JM); axis xy; colorbar
figure; 
subplot(3,3,1); imagesc(cwn_JM-cwn_GR); axis xy; colorbar
subplot(3,3,2); imagesc(cn_JM-cn_GR); axis xy; colorbar
subplot(3,3,3); imagesc(cen_JM-cen_GR); axis xy; colorbar
subplot(3,3,4); imagesc(cw_JM-cw_GR); axis xy; colorbar
subplot(3,3,5); imagesc(co_JM-co_GR); axis xy; colorbar
subplot(3,3,6); imagesc(ce_JM-ce_GR); axis xy; colorbar
subplot(3,3,7); imagesc(cws_JM-cws_GR); axis xy; colorbar
subplot(3,3,8); imagesc(cs_JM-cs_GR); axis xy; colorbar
subplot(3,3,9); imagesc(ces_JM-ces_GR); axis xy; colorbar

figure; 
k = 1;
subplot(3,3,1); plot(squeeze(cwn_JM(k,:)),'*'); hold on; plot(squeeze(cwn_GR(k,:)),'o');
subplot(3,3,2); plot(squeeze(cn_JM(k,:)),'*'); hold on; plot(squeeze(cn_GR(k,:)),'o');
subplot(3,3,3); plot(squeeze(cen_JM(k,:)),'*'); hold on; plot(squeeze(cen_GR(k,:)),'o');
subplot(3,3,4); plot(squeeze(cw_JM(k,:)),'*'); hold on; plot(squeeze(cw_GR(k,:)),'o');
subplot(3,3,5); plot(squeeze(co_JM(k,:)),'*'); hold on; plot(squeeze(co_GR(k,:)),'o');
subplot(3,3,6); plot(squeeze(ce_JM(k,:)),'*'); hold on; plot(squeeze(ce_GR(k,:)),'o');
subplot(3,3,7); plot(squeeze(cws_JM(k,:)),'*'); hold on; plot(squeeze(cws_GR(k,:)),'o');
subplot(3,3,8); plot(squeeze(cs_JM(k,:)),'*'); hold on; plot(squeeze(cs_GR(k,:)),'o');
subplot(3,3,9); plot(squeeze(ces_JM(k,:)),'*'); hold on; plot(squeeze(ces_GR(k,:)),'o');
figure; 
subplot(3,3,1); plot(squeeze(cwn_JM(k,:))-squeeze(cwn_GR(k,:)),'*'); 
subplot(3,3,2); plot(squeeze(cn_JM(k,:))-squeeze(cn_GR(k,:)),'*'); 
subplot(3,3,3); plot(squeeze(cen_JM(k,:))-squeeze(cen_GR(k,:)),'*'); 
subplot(3,3,4); plot(squeeze(cw_JM(k,:))-squeeze(cw_GR(k,:)),'*'); 
subplot(3,3,5); plot(squeeze(co_JM(k,:))-squeeze(co_GR(k,:)),'*'); 
subplot(3,3,6); plot(squeeze(ce_JM(k,:))-squeeze(ce_GR(k,:)),'*'); 
subplot(3,3,7); plot(squeeze(cws_JM(k,:))-squeeze(cws_GR(k,:)),'*'); 
subplot(3,3,8); plot(squeeze(cs_JM(k,:))-squeeze(cs_GR(k,:)),'*');
subplot(3,3,9); plot(squeeze(ces_JM(k,:))-squeeze(ces_GR(k,:)),'*');

% arranged in a sparse matrix
I=[];
J=[];
s=[];
% east
c=ce_JM(k,im);
[ii,kk]=meshgrid(ip,k);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% west
c=cw_JM(k,ip);
[ii,kk]=meshgrid(ip,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,k);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% north
c=cn_JM(km,i);
[ii,kk]=meshgrid(i,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(i,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% south
c=cs_JM(kp,i);
[ii,kk]=meshgrid(i,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(i,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% central
c=co_JM(k,i);
[ii,kk]=meshgrid(i,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;I0];
s=[s;c(:)];
% northeast
c=cen_JM(km,im);
[ii,kk]=meshgrid(ip,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% southeast
c=ces_JM(kp,im);
[ii,kk]=meshgrid(ip,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% northwest
c=cwn_JM(km,ip);
[ii,kk]=meshgrid(im,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(ip,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% southwest
c=cws_JM(kp,ip);
[ii,kk]=meshgrid(im,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(ip,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% A
A_JM=sparse(I,J,s,nx*nz,nx*nz);

figure;
imagesc(A_JM-A_JM'); colorbar

% define rhs
L = 10e3;
H = 4e3;
alpha = 25./(L^2);

i2x = inline('ix','ix');

ix = (0.5:nx)/nx;
x = L*i2x(ix);
dxu = diff(x);

ix = (0:nx)/nx;
xu = L*i2x(ix);
dxw = diff(xu);

xr = ones(nz,1)*x;

topo = inline('(H - 0.5*H*exp(-alpha*(x-0.5*L).^2))','x','H','L','alpha');

h = topo(x,H,L,alpha).*ones(1,nx);

zw = zlevs_ND(h,nz,'w');
zw = reshape(zw,nz+1,nx);
dzr = zw(2:end,:)-zw(1:end-1,:);
dxr = ones(nz,1)*diff(xu,1,2);
vr = (zw(2:end,:)-zw(1:end-1,:)).*dxr;

zr = zlevs_ND(h,nz,'r');
zr = reshape(zr,nz,nx);

bet = 600 / (L*L);
x1 = L * 0.65;
z1 = H * (0.75 - 1);
x2 = L * 0.75;
z2 = H * (0.65 - 1);

rhs = vr .* (exp(-bet * ((xr-x1).^2 + (zr-z1).^2)) ...
          -  exp(-bet * ((xr-x2).^2 + (zr-z2).^2)));

% solve for p
p =  A_JM \ rhs(:);
p = reshape(p,nz,nx);

% pcolor plot
xu_pc = cat(2,zeros(nz,1),0.5*(xr(:,1:end-1)+xr(:,2:end)),L*ones(nz,1));
xu_pc = cat(1,xu_pc,xu_pc(end,:));
zwu_pc = cat(2,zw(:,1),0.5*(zw(:,1:end-1)+zw(:,2:end)),zw(:,end));
rhs_pc = cat(2,rhs,rhs(:,end));
rhs_pc = cat(1,rhs_pc,rhs_pc(end,:));
p_pc = cat(2,p,p(:,end));
p_pc = cat(1,p_pc,p_pc(end,:));

figure;
pcolor(xu_pc,zwu_pc,rhs_pc); shading flat;
hold on
plot(xr(1,:),zw(1,:),'k','linewidth',1)
xlim([0 L])
ylim([-H 0])
axis equal tight
colorbar
cmin = min(rhs(:));
cmax = max(rhs(:));
caxis([cmin cmax]);
title(['rhs [' num2str(cmin) ' ' num2str(cmax) ']']);

figure;
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
title(['pressure mgroms coeff [' num2str(cmin) ' ' num2str(cmax) ']']);
