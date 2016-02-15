%% retrieve the divergence free velocity field that satisfies the boundary conditions
% u=w==0 on all sides of the mask

nx = 128; 
nz = 128;

[A,ce,cw,cn,cs,cen,cwn,ces,cws,co,h,xr,zr,vr,zw,dzr,dzu,dzw,dxu,dxw,alphauw,alphaw] = set_nhlap_ND(nx,nz);

%% define coord
xu=0.5*(xr(:,2:end)+xr(:,1:end-1));
zu=0.5*(zr(:,2:end)+zr(:,1:end-1));
xw=ones(nz+1,1)*xr(1,:);

%% define (ur,wr) and compute (u,w)
testcase = 1; 
switch testcase
    case 1
        % to make comparisons 
        sizex=1e4
        wr=zw;
        wr=cos(xw*pi/sizex);
        ur=xu*0;
        [u,w]=real2momentum(ur,wr,alphauw);
    case 2
        % uniform stratification
        b=zr;
        w=zeros(nz+1,nx);
        u=zeros(nz,nx-1);
        w(2:nz,:)=0.5*(b(2:nz,:)+b(1:nz-1,:));
        w(nz+1,:)=2*w(nz,:)-w(nz-1,:);
        w(1,:)=2*w(2,:)-w(3,:);
        u=0.5*(b(:,2:end)+b(:,1:end-1)).*(zr(:,2:end)-zr(:,1:end-1));
        u=u./(ones(nz,1)*dxu);
end

%% compute (U,W);
[U,W]=momentum2flux(u,w,alphauw,alphaw,dzu,dzw,dxu,dxw);

%% compute div(U,W)
div = flux2div(U,W);

%% plot
figure;
cmin = min([ur(:);wr(:)]);
cmax = max([ur(:);wr(:)]);
subplot(3,2,1); imagesc(ur,[cmin cmax]); axis xy
subplot(3,2,2); imagesc(wr,[cmin cmax]); axis xy
subplot(3,2,3); imagesc(u,[cmin cmax]); axis xy
subplot(3,2,4); imagesc(w,[cmin cmax]); axis xy
%subplot(3,2,5); imagesc(U./dzu,[cmin cmax]); axis xy
%subplot(3,2,6); imagesc(W./(ones(129,1)*dxw),[cmin cmax]); axis xy
figure;
cmin = min([U(:);W(:)]);
cmax = max([U(:);W(:)]);
subplot(1,3,1); imagesc(U,[cmin cmax]); axis xy
subplot(1,3,2); imagesc(W,[cmin cmax]); axis xy
subplot(1,3,3); imagesc(div); axis xy; colorbar

%% solve for p
p =  A \ (div(:));
p = reshape(p,nz,nx);

%% plot
figure;
imagesc(p); axis xy; colorbar

%% correct (u,w)
du = (p(:,2:end)-p(:,1:end-1))./(ones(nz,1)*dxu);
dw = w*0.;
dw(2:nz,:) = (p(2:end,:)-p(1:end-1,:))./dzw(2:nz,:);
dw(1,:)    = 0;
dw(nz+1,:) = -2*p(end,:)./dzw(nz+1,:);

u = u-du;
w = w-dw;

%% compute (ur,wr)
[ur,wr]=momentum2real(u,w,alphauw);

%% compute (U,W)
[U,W]=momentum2flux(u,w,alphauw,alphaw,dzu,dzw,dxu,dxw);

%% compute div(U,W)
div = flux2div(U,W);

%% plot
figure;
cmin = min([ur(:);wr(:)]);
cmax = max([ur(:);wr(:)]);
subplot(3,2,1); imagesc(ur,[cmin cmax]); axis xy
subplot(3,2,2); imagesc(wr,[cmin cmax]); axis xy
subplot(3,2,3); imagesc(u,[cmin cmax]); axis xy
subplot(3,2,4); imagesc(w,[cmin cmax]); axis xy
%subplot(3,2,5); imagesc(U./dzu,[cmin cmax]); axis xy
%subplot(3,2,6); imagesc(W./(ones(129,1)*dxw)); axis xy
figure;
cmin = min([U(:);W(:)]);
cmax = max([U(:);W(:)]);
subplot(1,3,1); imagesc(U,[cmin cmax]); axis xy
subplot(1,3,2); imagesc(W,[cmin cmax]); axis xy
subplot(1,3,3); imagesc(div); axis xy; colorbar
