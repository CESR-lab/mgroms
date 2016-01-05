%% retrieve the divergence free velocity field that satisfies the boundary conditions
% u=w==0 on all sides of the mask
%
% from w = cos(x*pi/sizex).

if ~exist('nx','var')
nz=80*2;nx=200*2;
end

fprintf('--------- retrieve (u,w) with nx=%i / nz=%i\n',nx,nz)

[A,xr,zr,vr,h,alphauw,alphaw,dzu,dzw,dxu,dxw] = set_nhlap(nx,nz);
%A=0.5*(A+A');

figure(4)
d=sum(abs(A-A'));
d=sum(A);

B=A-diag(diag(A));
d=sum(B);

d=reshape(d,nz,nx);
imagesc(d(1:nz-1,:))
axis xy

%return

xu=0.5*(xr(:,2:end)+xr(:,1:end-1));


sizex=1e4;
%dx=sizex/nx;

zw=zeros(nz+1,nx);
zw(1,:)=-h;
zw(2:nz,:)=0.5*(zr(1:end-1,:)+zr(2:end,:));
zw(nz+1,:)=0.;
zu=0.5*(zr(:,2:end)+zr(:,1:end-1));
xw=ones(nz+1,1)*xr(1,:);

%%

testcase =1; % testcase=1 is the good one to make comparisons 


% first testcase
switch testcase
    case 1
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

% ---end of testcase definition ----

% compute the flux (U,W);

[U,W]=momentum2flux(u,w,alphauw,alphaw,dzu,dzw,dxu,dxw);

% compute div(U,W)

div = flux2div(U,W);

% solve for p

p =  A \ (div(:));
%p = p ./ vr(:);

p= reshape(p,nz,nx);

% correct (u,w)

du = (p(:,2:end)-p(:,1:end-1))./(ones(nz,1)*dxu);

dw = w*0.;
dw(2:nz,:) = (p(2:end,:)-p(1:end-1,:))./dzw(2:nz,:);
dw(1,:)    = 0;%2*dw(2,:)-dw(3,:);
%dw(end,:) = 2*dw(end-1,:)-dw(end-2,:);
dw(nz+1,:) = -2*p(end,:)./dzw(nz+1,:);

figure(4)
imagesc(div(2:end,:));axis xy
%return
%%
pmx=max(div(:));
%figure(3)
%clf
%pcolor(xr,zr,reshape(div,nz,nx));shading flat
%hold on
%contour(xr,zr,reshape(p,nz,nx),[-1:0.1:1]*pmx,'k');
%caxis([-1 1]*pmx)
%plot(xr,-h,'k','linewidth',2)


div0=div;
U0=U;
W0=W;
u0=u;
w0=w;
u=u-du;
w=w-dw;

[ur,wr]=momentum2real(u,w,alphauw);


% compute (U,W)
[U,W]=momentum2flux(u,w,alphauw,alphaw,dzu,dzw,dxu,dxw);

format long
KE0 = integrated_ke(u0,w0,U0,W0,vr,dxu,dzw)
KE  = integrated_ke(u,w,U,W,vr,dxu,dzw)
format short
%%

figure(1)
clf
subplot(2,1,1)
pcolor(xu,zu,reshape(u,nz,nx-1));shading flat
axis equal tight
%pcolor(xw,zw,reshape(W-W0,nz-1,nx));shading flat
hold on
plot(xr,-h,'k','linewidth',2)
%caxis([-0.95 0])
caxis([-1 1])

colorbar
title('u')
%
%figure(2)
subplot(2,1,2)
pcolor(xw,zw,reshape(w,nz+1,nx));shading flat
axis equal tight
hold on
plot(xr,-h,'k','linewidth',2)
%caxis([-1 1]*.55)
caxis([-1 1])
colorbar
title('w')

%print('-dpdf','UWdivergencefree.pdf')
%%
div = flux2div(U,W);
figure(3)
clf
pcolor(xr,zr,reshape(-p,nz,nx));shading flat
hold on
contourf(xr,zr,reshape(-p,nz,nx),20);
axis equal tight
colorbar

figure(2)
i=nx/2;
plot(zu(:,i),ur(:,i),'r')
xlabel('z [m]')
ylabel('u [m/s]')
title('ur profile at the seamount top')

format long
utop=interp1(zu(:,i),ur(:,i),zw(1,i),'linear','extrap')
format short
%%
figure(6)
imagesc(div(:,2:nx-1));axis xy
colorbar
title('divergence after projection')
%%
figure(5)
plot(W(2:nz+1,nx/4),'+-')

%%
%save('res_gr_400x160.mat','ur','wr','xu','zu','xw','zw','KE')