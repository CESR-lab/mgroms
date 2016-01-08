%% check sigma coord relaxation

% coord
zr_000 = ncread('../../src/zr_zr_000.nc','zr');
zr_001 = ncread('../../src/zr_zr_001.nc','zr');
zr_002 = ncread('../../src/zr_zr_002.nc','zr');
zr_003 = ncread('../../src/zr_zr_003.nc','zr');

zw_000 = ncread('../../src/zw_zw_000.nc','zw');
zw_001 = ncread('../../src/zw_zw_001.nc','zw');
zw_002 = ncread('../../src/zw_zw_002.nc','zw');
zw_003 = ncread('../../src/zw_zw_003.nc','zw');

% rhs
rhs_000 = ncread('../../src/rhs_rhs_000.nc','rhs');
rhs_001 = ncread('../../src/rhs_rhs_001.nc','rhs');
rhs_002 = ncread('../../src/rhs_rhs_002.nc','rhs');
rhs_003 = ncread('../../src/rhs_rhs_003.nc','rhs');

figure; 
k = 120;
cmin = min(min(cat(1,squeeze(rhs_002(k,:,:)),squeeze(rhs_003(k,:,:)),squeeze(rhs_000(k,:,:)),squeeze(rhs_001(k,:,:)))));
cmax = max(max(cat(1,squeeze(rhs_002(k,:,:)),squeeze(rhs_003(k,:,:)),squeeze(rhs_000(k,:,:)),squeeze(rhs_001(k,:,:)))));              
subplot(2,2,1); imagesc(squeeze(rhs_002(k,:,:)),[cmin cmax]); axis xy;
subplot(2,2,2); imagesc(squeeze(rhs_003(k,:,:)),[cmin cmax]); axis xy;
subplot(2,2,3); imagesc(squeeze(rhs_000(k,:,:)),[cmin cmax]); axis xy;
subplot(2,2,4); imagesc(squeeze(rhs_001(k,:,:)),[cmin cmax]); axis xy;
figure;
j = 20;
cmin = min(min(cat(2,squeeze(rhs_000(:,j,:)),squeeze(rhs_001(:,j,:)))));
cmax = max(max(cat(2,squeeze(rhs_000(:,j,:)),squeeze(rhs_001(:,j,:)))));
subplot(1,2,1); imagesc(squeeze(rhs_000(:,j,:)),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(squeeze(rhs_001(:,j,:)),[cmin cmax]); axis xy;

% p
p_000 = ncread('../../src/p_p_000.nc','p');
p_001 = ncread('../../src/p_p_001.nc','p');
p_002 = ncread('../../src/p_p_002.nc','p');
p_003 = ncread('../../src/p_p_003.nc','p');

figure; 
k = 120;
cmin = min(min(cat(1,squeeze(p_002(k,:,:)),squeeze(p_003(k,:,:)),squeeze(p_000(k,:,:)),squeeze(p_001(k,:,:)))));
cmax = max(max(cat(1,squeeze(p_002(k,:,:)),squeeze(p_003(k,:,:)),squeeze(p_000(k,:,:)),squeeze(p_001(k,:,:)))));              
subplot(2,2,1); imagesc(squeeze(p_002(k,:,:)),[cmin cmax]); axis xy;
subplot(2,2,2); imagesc(squeeze(p_003(k,:,:)),[cmin cmax]); axis xy;
subplot(2,2,3); imagesc(squeeze(p_000(k,:,:)),[cmin cmax]); axis xy;
subplot(2,2,4); imagesc(squeeze(p_001(k,:,:)),[cmin cmax]); axis xy;
figure;
j = 20;
cmin = min(min(cat(2,squeeze(p_000(:,j,:)),squeeze(p_001(:,j,:)))));
cmax = max(max(cat(2,squeeze(p_000(:,j,:)),squeeze(p_001(:,j,:)))));
subplot(1,2,1); imagesc(squeeze(p_000(:,j,:)),[cmin cmax]); axis xy;
subplot(1,2,2); imagesc(squeeze(p_001(:,j,:)),[cmin cmax]); axis xy;
figure;
cmin = min(min(cat(2,squeeze(p_000(:,j,2:end-1)),squeeze(p_001(:,j,2:end-1)))));
cmax = max(max(cat(2,squeeze(p_000(:,j,2:end-1)),squeeze(p_001(:,j,2:end-1)))));
imagesc(cat(2,squeeze(p_000(:,j,2:end-1)),squeeze(p_001(:,j,2:end-1))),[cmin cmax]); axis xy;

