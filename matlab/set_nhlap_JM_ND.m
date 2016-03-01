function [A,cec,cwc,cct,ccb,cet,cwt,ceb,cwb,ccc,h1,xr,zr,vr,zw] = set_nhlap_JM_ND(nx,nz)

sizex = 1e4;
sizez = 4e3;

%% horiz metrics

x1 = sizex*[0.5:1:nx-0.5]'/nx;
x1u= sizex*[0.0:1:nx-0.0]'/nx;

dx = x1u(2:end)-x1u(1:end-1);
dxi = 1.0./dx;

dxu = zeros(nx+1,1);
dxu(2:nx) = x1(2:nx)-x1(1:nx-1);  dxu(1) = dxu(2);dxu(nx+1) = dxu(nx);
dxiu = 1.0./dxu;

xr = zeros(nx,nz);
xw = zeros(nx,nz+1);
xu = zeros(nx+1,nz);
xp = zeros(nx,nz+1);

for k=1:nz
   xr(:,k) = x1;
   xw(:,k) = x1;
   xp(:,k) = x1;
   xu(:,k) = x1u;
end
xw(:,nz+1) = x1;
xp(:,nz+1) = x1;

%% vert metrics

alpha = 25./(sizex^2);
h1 = (sizez - 0.5*sizez*exp(-alpha*(x1-0.5*sizex).^2));
h1u = (sizez - 0.5*sizez*exp(-alpha*(x1u-0.5*sizex).^2));
% h1 = sizez*ones(nx,1);
% h1u = sizez*ones(nx+1,1);

zr = zeros(nx,nz);
zw = zeros(nx,nz+1);
zu = zeros(nx+1,nz);
zp = zeros(nx,nz+1);

for i=1:nx
    zr(i,:) = h1(i)*[0.5:1:nz-0.5]/nz -h1(i);
    zw(i,:) = h1(i)*[0.0:1:nz]/nz -h1(i);
end
for k=1:nz
    zu(:,k) = h1u*(k-0.5)/nz -h1u;
end
   
dz = zw(:,2:nz+1)-zw(:,1:nz);

vr = dz.*repmat(dx,[1 nz]);

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

%% slopes

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
   
% Slope in the bottom ghost cells

dz_bot = 2*(zr(:,1)-zw(:,1));
zr_bot = zr(:,1)-dz_bot;

zxu_bot = zeros(nx+1,1);
zxu_bot(2:nx) = (zr_bot(2:nx)-zr_bot(1:nx-1));
zx_bot = 0.5*(zxu_bot(2:nx+1)+zxu_bot(1:nx)).*dxi;

% Slope in the surface ghost cells

dz_sur = 2*(zw(:,nz+1)-zr(:,nz));
zr_sur = zr(:,nz)+dz_sur;

zxu_sur = zeros(nx+1,1);
zxu_sur(2:nx) = (zr_sur(2:nx)-zr_sur(1:nx-1));
zx_sur = 0.5*(zxu_sur(2:nx+1)+zxu_sur(1:nx)).*dxi;

%% coeff

zx4 = 0.25 * zx;

% central, up, down, right and left terms
ccc = zeros(nx,nz);
for k = 1:nz
   ccc(:,k) = - dzu(2:nx+1,k).*dxiu(2:nx+1)-dzu(1:nx,k).*dxiu(1:nx) ...
                - dx.*dziw(:,k+1).*(1+zxw(:,k+1).^2) ...
                - dx.*dziw(:,k  ).*(1+zxw(:,k).^2);
end

ccb = zeros(nx,nz+1);       
for k = 1:nz+1
   ccb(:,k) = dx.*dziw(:,k).*(1+zxw(:,k).^2); 
end
cct = ccb(:,2:nz+1);
ccb = ccb(:,1:nz);

cwc = zeros(nx+1,nz);
for k = 1:nz
   cwc(:,k) =  dzu(:,k).*dxiu;
end
cec = cwc(2:nx+1,:);
cwc = cwc(1:nx,:);

% cross terms
cwt = zeros(nx,nz);
cwt(2:nx,:) = 0.25*zx(1:nx-1,:);
cwt(:,1:nz-1) = cwt(:,1:nz-1) + 0.25*zx(:,2:nz);
cwt(:,nz) = cwt(:,nz) + 0.25*zx_sur; % ultimately used in cwc

ceb = zeros(nx,nz);
ceb(1:nx-1,:) = 0.25*zx(2:nx,:);
ceb(:,2:nz) = ceb(:,2:nz) + 0.25*zx(:,1:nz-1);
ceb(:,1) = ceb(:,1) + 0.25*zx_bot;

cwb = zeros(nx,nz);
cwb(2:nx,:) = -0.25*zx(1:nx-1,:);
cwb(:,2:nz) = cwb(:,2:nz) - 0.25*zx(:,1:nz-1);
cwb(:,1) = cwb(:,1) - 0.25*zx_bot;

cet = zeros(nx,nz);
cet(1:nx-1,:) = -0.25*zx(2:nx,:);
cet(:,1:nz-1) = cet(:,1:nz-1) - 0.25*zx(:,2:nz);
cet(:,nz) = cet(:,nz) - 0.25*zx_sur; % ultimately used in cec

% bc's at the lef:  no flux and no px, slightly subtle!
ccc(1,:) = ccc(1,:) + cwc(1,:);
cct(1,1:nz-1) = cct(1,1:nz-1) + zx4(1,2:nz) - zx4(1,1:nz-1);
ccb(1,2:nz  ) = ccb(1,2:nz  ) + zx4(1,2:nz) - zx4(1,1:nz-1);

% bc's at the right
ccc(nx,:) = ccc(nx,:) + cec(nx,:);
cct(nx,1:nz-1) = cct(nx,1:nz-1) - zx4(nx,2:nz) + zx4(nx,1:nz-1);
ccb(nx,2:nz  ) = ccb(nx,2:nz  ) - zx4(nx,2:nz) + zx4(nx,1:nz-1);

% bc's at the top
ccc(:,nz) = ccc(:,nz) - cct(:,nz);
cwc(:,nz) = cwc(:,nz) - cwt(:,nz);
cec(:,nz) = cec(:,nz) - cet(:,nz);

%% matrix

nd = nx*(nz+1);
is = 1;
ks = nx;
A = sparse(nd,nd,9*nd);         
for i = 1:nx
   for k = 1:nz-1
       idx = (k-1)*nx + i-1 +nx;
       if i>1
         A(idx+1,idx+1-is)   = cwc(i,k);
         A(idx+1,idx+1-is+ks)= cwt(i,k);
         A(idx+1,idx+1-is-ks)= cwb(i,k);
       end
       if i<nx
         A(idx+1,idx+1+is)   = cec(i,k);
         A(idx+1,idx+1+is+ks)= cet(i,k);
         A(idx+1,idx+1+is-ks)= ceb(i,k);
       end
       A(idx+1,idx+1-ks)= ccb(i,k);
       A(idx+1,idx+1+ks)= cct(i,k);
       A(idx+1,idx+1   )= ccc(i,k); 
   end
   k = nz;
   idx = (k-1)*nx + i-1 + nx;
   if i>1
      A(idx+1,idx+1-is)   = cwc(i,k);
      A(idx+1,idx+1-is-ks)= cwb(i,k);
   end
   if i<nx
      A(idx+1,idx+1+is)   = cec(i,k);
      A(idx+1,idx+1+is-ks)= ceb(i,k);
   end
   A(idx+1,idx+1-ks)= ccb(i,k);
   A(idx+1,idx+1   )= ccc(i,k);    
end
% Bottom boundary condition
for i = 1:nx
   idx =  i-1;
   if i>1
     A(idx+1,idx+1-is+ks) =  0.25*zx(i,1);
   end
   if i<nx
      A(idx+1,idx+1+is+ks) = -0.25*zx(i,1);
   end
   A(idx+1,idx+1      ) = -0.5*dx(i)*dziw(i,1)*(1+zxw(i,1)*zxw(i,1));
   A(idx+1,idx+1   +ks) = +0.5*dx(i)*dziw(i,1)*(1+zxw(i,1)*zxw(i,1));
end
