%% solve for pressure

nx = 128; 
nz = 128;

[A,ce,cw,cn,cs,cen,cwn,ces,cws,co,h,xr,zr,vr,zw] = set_nhlap_GR_ND(nx,nz);

% define rhs
L = 10e3;
H = 4e3;
bet = 600 / (L*L);
x1 = L * 0.65;
z1 = H * (0.75 - 1);
x2 = L * 0.75;
z2 = H * (0.65 - 1);

rhs = vr .* (exp(-bet * ((xr-x1).^2 + (zr-z1).^2)) ...
          -  exp(-bet * ((xr-x2).^2 + (zr-z2).^2)));

% solve for p
p =  A \ rhs(:);
p = reshape(p,nz,nx);

% pcolor plot
xu_pc = cat(2,zeros(nz,1),0.5*(xr(:,1:end-1)+xr(:,2:end)),L*ones(nz,1));
xu_pc = cat(1,xu_pc,xu_pc(end,:));
zwu_pc = cat(2,zw(:,1),0.5*(zw(:,1:end-1)+zw(:,2:end)),zw(:,end));
rhs_pc = cat(2,rhs,rhs(:,end));
rhs_pc = cat(1,rhs_pc,rhs_pc(end,:));
p_pc = cat(2,p,p(:,end));
p_pc = cat(1,p_pc,p_pc(end,:));

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
% title(['rhs GR [' num2str(cmin) ' ' num2str(cmax) ']']);

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
title(['pressure GR [' num2str(cmin) ' ' num2str(cmax) ']']);
