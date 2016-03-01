
%--------------------------------------------
sizex = 1e4;
sizez = 4e3;

nx = 128; %100/4;
nz = 128; %40/4;


x1 = sizex*[0.5:1:nx-0.5]'/nx;
x1u= sizex*[0.0:1:nx-0.0]'/nx;
% if 1
%   x1 = sizex*(x1/sizex+0.1*sin(2*x1*pi/sizex) );
%   x1u= sizex*(x1u/sizex+0.1*sin(2*x1u*pi/sizex) );
% end

dx = x1u(2:end)-x1u(1:end-1);
dxi = 1.0./dx;
dxu = zeros(nx+1,1);
dxu(2:nx) = 0.5*(dx(1:nx-1)+dx(2:nx)); dxu(1) = dxu(2);dxu(nx+1) = dxu(nx);

dxu(2:nx) = x1(2:nx)-x1(1:nx-1);  dxu(1) = dxu(2);dxu(nx+1) = dxu(nx);
dxiu = 1.0./dxu;

alpha = 25./(sizex^2);
h1 = (sizez - 0.5*sizez*exp(-alpha*(x1-0.5*sizex).^2));
h1u = (sizez - 0.5*sizez*exp(-alpha*(x1u-0.5*sizex).^2));

xr = zeros(nx,nz);
zr = zeros(nx,nz);
zw = zeros(nx,nz+1);
xw = zeros(nx,nz+1);
xu = zeros(nx+1,nz);
zu = zeros(nx+1,nz);
zp = zeros(nx,nz+1);
xp = zeros(nx,nz+1);

for k=1:nz
   xr(:,k) = x1;
   xw(:,k) = x1;
   xp(:,k) = x1;
   xu(:,k) = x1u;
end
xw(:,nz+1) = x1;
xp(:,nz+1) = x1;

% if 0
   for i=1:nx
      zr(i,:) = h1(i)*[0.5:1:nz-0.5]/nz -h1(i);
      zw(i,:) = h1(i)*[0.0:1:nz]/nz -h1(i);
   end
   for k=1:nz
      zu(:,k) = h1u*(k-0.5)/nz -h1u;
   end
% else
%    zu = zlevs(h1u,0*h1u,6.,1.,200,nz,'r')';
%    zr = zlevs(h1 ,0*h1 ,6.,1.,200,nz,'r')';
%    zw = zlevs(h1 ,0*h1 ,6.,1.,200,nz,'w')';
% end

dz = zw(:,2:nz+1)-zw(:,1:nz);

dzu = zeros(nx+1,nz);
dzu(2:nx,:) = 0.5*(dz(1:nx-1,:)+dz(2:nx,:));
dzu(1,:) = dzu(2,:);dzu(nx+1,:) = dzu(nx,:);

dzw = zeros(nx,nz+1);
dzw(:,2:nz) = zr(:,2:nz) - zr(:,1:nz-1);
dzw(:,1   ) = 2*(zr(:,1)-zw(:,1));
dzw(:,nz+1) = 2*(zw(:,nz+1)-zr(:,nz));

dziw = 1./dzw;

zr_bot = zr(:,1)-dzw(:,1);
zp(:,2:nz+1) = zr;
zp(:,1) = zr_bot;


zxu = zeros(nx+1,nz);
zxu(2:nx,:) = (zr(2:nx,:)-zr(1:nx-1,:));

zx = zeros(nx,nz);
for k=1:nz
 zx(:,k) = 0.5*(zxu(2:nx+1,k)+zxu(1:nx,k)).*dxi;
end
for k=1:nz
 zxu(:,k) = zxu(:,k).*dxiu;
end

delzw = zeros(nx+1,nz+1);
delzw(2:nx,:) = zw(2:nx,:)-zw(1:nx-1,:);
zxw = zeros(nx,nz+1);
for k=1:nz+1
   zxw(:,k) = 0.5*(delzw(2:nx+1,k)+delzw(1:nx,k)).*dxi;
end


%G = -makeG(nx,nz,dxiu,dziw);
%disp('G is done')
D = makeD(nx,nz,dxiu,dziw);
disp('D is done')
T = makeTm(nx,nz,zx,zxw,dx,dzu,dzw);
Tr2m = makeR2M(nx,nz,zx,dzu,dzw);
disp('matrices are done')

%P = D*T*G;
P = D*T*D';
return

ur = zeros(nx+1,nz);
wr = cos(xw*pi/sizex); %:,1) = w(:,1) = 1;
wr = zw/sizez;
Xr = uw2x(ur,wr);  %% physical velocity
Xm = Tr2m*Xr;      %% physical to momentum 
[um,wm] = x2uw(Xm,nx,nz);

%wm(:,1) = 100;
Xm = uw2x(um,wm);

%% check Ek
  Xf = T*Xm;
  [uf,wf] = x2uw(Xf,nx,nz);
  [um,wm] = x2uw(Xm,nx,nz);
  for i=1:nx
    uf(i,:) = uf(i,:)/dxu(i)./dzu(i,:);
    wf(i,:) = wf(i,:)/dx(i)./dzw(i,:);
  end
  u2 = 0.5*uf.*um;
  w2 = 0.5*wf.*wm;
  Ek = u2(2:nx+1,:) + u2(1:nx,:)+w2(:,2:nz+1) + w2(:,1:nz);

  


rhs = D*T*Xm;

p = P\rhs; p2 = reshape(p,[nx nz+1]);pnh = -p2(:,2:nz+1);

%dX = G*p;
dX = D'*p;
Xm = Xm-dX;

ndu = (nx+1)*nz;
for i = 1:nx
   udx = i;
   wdx = i + ndu;
%   Xm(wdx) = 0.5*zx(i,1)*(dxu(i)*Xm(udx)+dxu(i+1)*Xm(udx+1))/dx(i)/(1+zxw(i,1)*zxw(i,1));
end


Ek = dX'*T*Xm

res = max(abs(D*T*Xm))

Xf = T*Xm;
Xr = Tr2m\Xm;
[ur,wr] = x2uw(Xr,nx,nz);
[uf,wf] = x2uw(Xf,nx,nz); %% these have units of m/s times m*m
for i=1:nx
   uf(i,:) = uf(i,:)/dxu(i)./dzu(i,:);
   wf(i,:) = wf(i,:)/dx(i)./dzw(i,:);
end
uf(nx+1,:)= uf(nx+1,:)./dxu(nx+1)./dzu(nx+1,:);
[um,wm] = x2uw(Xm,nx,nz);





figure(1)
pmx = max(max(abs(pnh)));
pcolor(xr,zr,pnh);colorbar;shading flat;
title('Sigma coord')
hold on;
plot(x1,h1,'k','linewidth',2)
set(gca,'linewidth',2)
set(gca,'fontsize',15)
contour(xr,zr,pnh,[-1:0.1:1]*pmx,'k');
hold off