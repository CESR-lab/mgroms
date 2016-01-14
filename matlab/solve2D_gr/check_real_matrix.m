%% check sigma coord metrics, slopes and matrix 

% Guillaume's matlab code
nx = 128; 
nz = 128;
[A,ce,cw,cn,cs,cen,cwn,ces,cws,co,h,xr,zr,vr,zw,dzr,dzu,dzw,dxu,dxw,alphauw,alphaw] = set_nhlap_ND(nx,nz);

% Jeroen's matlab code

% fortran code
h_000 = ncread('../../src/h_h_000.nc','h');
h_001 = ncread('../../src/h_h_001.nc','h');
h_002 = ncread('../../src/h_h_002.nc','h');
h_003 = ncread('../../src/h_h_003.nc','h');

zr_000 = ncread('../../src/zr_zr_000.nc','zr');
zr_001 = ncread('../../src/zr_zr_001.nc','zr');
zr_002 = ncread('../../src/zr_zr_002.nc','zr');
zr_003 = ncread('../../src/zr_zr_003.nc','zr');

zw_000 = ncread('../../src/zw_zw_000.nc','zw');
zw_001 = ncread('../../src/zw_zw_001.nc','zw');
zw_002 = ncread('../../src/zw_zw_002.nc','zw');
zw_003 = ncread('../../src/zw_zw_003.nc','zw');

dzr_000 = ncread('../../src/dz_dz_000.nc','dz');
dzr_001 = ncread('../../src/dz_dz_001.nc','dz');
dzr_002 = ncread('../../src/dz_dz_002.nc','dz');
dzr_003 = ncread('../../src/dz_dz_003.nc','dz');

dzw_000 = ncread('../../src/dzw_dzw_000.nc','dzw');
dzw_001 = ncread('../../src/dzw_dzw_001.nc','dzw');
dzw_002 = ncread('../../src/dzw_dzw_002.nc','dzw');
dzw_003 = ncread('../../src/dzw_dzw_003.nc','dzw');

dxu_000 = ncread('../../src/dxu_dxu_000.nc','dxu');
dxu_001 = ncread('../../src/dxu_dxu_001.nc','dxu');
dxu_002 = ncread('../../src/dxu_dxu_002.nc','dxu');
dxu_003 = ncread('../../src/dxu_dxu_003.nc','dxu');

cA_000 = ncread('../../src/cA_cA_000.nc','cA');
cA_001 = ncread('../../src/cA_cA_001.nc','cA');
cA_002 = ncread('../../src/cA_cA_002.nc','cA');
cA_003 = ncread('../../src/cA_cA_003.nc','cA');

% cmp h : okay
figure; 
cmin = min(h_000(:)); 
cmax = max(h_000(:));
subplot(2,2,1); imagesc(h_002,[cmin cmax]); axis xy;
subplot(2,2,2); imagesc(h_003,[cmin cmax]); axis xy;
subplot(2,2,3); imagesc(h_000,[cmin cmax]); axis xy;
subplot(2,2,4); imagesc(h_001,[cmin cmax]); axis xy;
figure;
plot(cat(2,h_000(2:end-1,2:end-1),h_001(2:end-1,2:end-1))','*'); hold on
plot(cat(2,h_002(2:end-1,2:end-1),h_003(2:end-1,2:end-1))','*');
plot(h,'o')

% cmp zr : okay
figure;
j = 2
cmin = min(min(squeeze(zr_000(:,j,:))));
cmax = max(max(squeeze(zr_000(:,j,:))));
subplot(1,2,1); imagesc(cat(2,squeeze(zr_000(:,j,2:end-1)),squeeze(zr_001(:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(zr,[cmin cmax]); axis xy;
figure;
k = 1
plot(cat(1,squeeze(zr_000(k,j,2:end-1)),squeeze(zr_001(k,j,2:end-1))),'*'); hold on
plot(zr(k,:),'o')

% cmp zw : okay
figure
j = 2
cmin = min(min(squeeze(zw_000(:,j,:))));
cmax = max(max(squeeze(zw_000(:,j,:))));
subplot(1,2,1); imagesc(cat(2,squeeze(zw_000(:,j,2:end-1)),squeeze(zw_001(:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(zw,[cmin cmax]); axis xy;
figure;
k = 1
plot(cat(1,squeeze(zw_000(k,j,2:end-1)),squeeze(zw_001(k,j,2:end-1))),'*'); hold on
plot(zw(k,:),'o')

% cmp dzr : okay
figure
j = 2
cmin = min(min(squeeze(dzr_000(:,j,:))));
cmax = max(max(squeeze(dzr_000(:,j,:))));
subplot(1,2,1); imagesc(cat(2,squeeze(dzr_000(:,j,2:end-1)),squeeze(dzr_001(:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(dzr,[cmin cmax]); axis xy;
figure;
k = 10
plot(cat(1,squeeze(dzr_000(k,j,2:end-1)),squeeze(dzr_001(k,j,2:end-1))),'*'); hold on
plot(dzr(k,:),'o')

% cmp dzw : problem at top level -> solved by adding a level and recomputing
figure
j = 2
cmin = min(min(cat(2,squeeze(dzw_000(:,j,2:end-1)),squeeze(dzw_001(:,j,2:end-1)))));
cmax = max(max(cat(2,squeeze(dzw_000(:,j,2:end-1)),squeeze(dzw_001(:,j,2:end-1)))));
subplot(1,2,1); imagesc(cat(2,squeeze(dzw_000(:,j,2:end-1)),squeeze(dzw_001(:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(dzw,[cmin cmax]); axis xy;
figure;
k = 1
plot(cat(1,squeeze(dzw_000(k,j,2:end-1)),squeeze(dzw_001(k,j,2:end-1))),'*'); hold on
plot(dzw(k,:),'o')

% cmp dxu : okay
figure;
plot(cat(2,dxu_000(2:end-1,2:end-1),dxu_001(2:end-1,2:end-1))','*'); hold on
plot(dxu,'o')

% cmp cA with i-1 : problem at the bottom -> solved by redefining zx and zx
figure;
c = 7;
j = 2;
cmin = min(min(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))))); 
cmax = max(max(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1)))));
subplot(1,2,1); imagesc(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(ce,[cmin cmax]); axis xy;
figure;
k = 1;
plot(cat(1,squeeze(cA_000(c,k,j,2:end-1)),squeeze(cA_001(c,k,j,2:end-1))),'*'); hold on
plot(ce(k,:),'o')

% cmp cA with j-1 : not to be taken into account when comparing with the 2D version
figure;
c = 4;
j = 2;
cmin = min(min(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))))); 
cmax = max(max(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1)))));
subplot(1,2,1); imagesc(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))),[cmin cmax]); axis xy;

% cmp cA with k-1 : problem at the top -> solved by redefining dzw
figure;
c = 2;
j = 2;
cmin = min(min(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))))); 
cmax = max(max(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1)))));
subplot(1,2,1); imagesc(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(cs,[cmin cmax]); axis xy;
figure;
k = 128;
plot(cat(1,squeeze(cA_000(c,k,j,2:end-1)),squeeze(cA_001(c,k,j,2:end-1))),'*'); hold on
plot(cs(k,:),'o')

% cmp cA with i-1, k+1 : problem with sign?
figure;
c = 6;
j = 2;
cmin = min(min(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))))); 
cmax = max(max(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1)))));
subplot(1,2,1); imagesc(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(cen,[cmin cmax]); axis xy;
figure;
k = 10;
plot(cat(1,squeeze(cA_000(c,k,j,2:end-1)),squeeze(cA_001(c,k,j,2:end-1))),'*'); hold on
plot(cen(k,:),'o')

% cmp cA with j-1, k+1 : not to be taken into account when comparing with the 2D version
figure;
c = 3
j = 2

% cmp cA with i-1, k-1 : problem with sign?
%                        bottom value in fortran code represent coupling
%                        with i-1, j-1 : 0 in that case of no meridional slope
figure;
c = 8;
j = 2;
cmin = min(min(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))))); 
cmax = max(max(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1)))));
subplot(1,2,1); imagesc(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(ces,[cmin cmax]); axis xy;
figure;
k = 1;
plot(cat(1,squeeze(cA_000(c,k,j,2:end-1)),squeeze(cA_001(c,k,j,2:end-1))),'*'); hold on
plot(ces(k,:),'o')

% cmp cA with j-1, k-1 : not to be taken into account when comparing with the 2D version
figure;
c = 5
j = 2

% cmp cA with itself : interior boundaries problem -> solved with fill_halo, 
%                      problem at the top -> solved by changing from 2 to 3
figure;
c = 1;
j = 2;
cmin = min(min(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))))); 
cmax = max(max(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1)))));
subplot(1,2,1); imagesc(cat(2,squeeze(cA_000(c,:,j,2:end-1)),squeeze(cA_001(c,:,j,2:end-1))),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(co,[cmin cmax]); axis xy;
figure;
k = 128;
plot(cat(1,squeeze(cA_000(c,k,j,2:end-1)),squeeze(cA_001(c,k,j,2:end-1))),'*'); hold on
plot(co(k,:),'o')

