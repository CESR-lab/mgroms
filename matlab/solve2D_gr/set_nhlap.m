function [A,xr,zr,vr,h,alphauw,alphaw,dzu,dzw,dxu,dxw] = set_nhlap(nx,nz)
%% setup the 9 points pressure operator in sigma coordinates

L=10e3;
H=4e3;
hb=2500;
thetas=6;
thetab=1;
hc=200;

% regular grid
dx=L/nx;
x=(0.5:nx)*dx;

i2x=inline('ix+0.1*sin(2*ix*pi)','ix');
topo = inline('(H - 0.5*H*exp(-alpha*(x-0.5*L).^2))','x','H','L','alpha');
%i2x=inline('ix','ix');

% stretched grid
ix=(0.5:nx)/nx;
x=L*i2x(ix);
dxu=diff(x);

ix=(0:nx)/nx;
xu=L*i2x(ix);%( ix+0.1*sin(2*ix*pi) );
dxw=diff(xu);

xr=ones(nz,1)*x;

% topo
%h=H-hb*exp( -(x-L*.5).^2/(2*(0.05*L)^2));

% Jeroen's topo alpha = 25./(L^2);
alpha = 25./(L^2);

%alpha = alpha*2;

h = topo(x,H,L,alpha);

zw=zlevs(h,h*0,thetas,thetab,hc,nz,'w');
zw=reshape(zw,nz+1,nx);

dxr=ones(nz,1)*diff(xu,1,2);
%whos xu dxr
vr=(zw(2:end,:)-zw(1:end-1,:)).*dxr;

hu=topo(xu,H,L,alpha);%0.5*(h(:,2:end)+h(:,1:end-1));

zwu=zlevs(hu,hu*0,thetas,thetab,hc,nz,'w');
zwu=reshape(zwu,(nz+1),nx+1);
dzu = zwu(2:end,2:end-1)-zwu(1:end-1,2:end-1);

dzr=zw(2:end,:)-zw(1:end-1,:);

zr=zlevs(h,h*0,thetas,thetab,hc,nz,'r');
zr=reshape(zr,nz,nx);
dzw=zeros(nz+1,nx);
dzw(2:end-1,:) = zr(2:end,:)-zr(1:end-1,:);

% 2nd order interpolation to get 
%dzw(nz+1,:) = 2*(zw(nz+1,:)-zr(nz,:));
dzw(nz+1,:) = ( 3*zw(nz+1,:)-3*zr(nz,:)+zw(nz,:))-zr(nz,:);



%dzw(1,:) = 2*dzw(2,:)-dzw(3,:);
%dzw(1,:)=-(3*zw(1,:)-2*zr(1,:));
dzw(1,:)=zr(1,:)-( 3*zw(1,:)-3*zr(1,:)+zw(2,:));
%dzw(1,:)=3*dzw(2,:)-dzr(1,:)-dzw(3,:);

whos dzw zwu xu
alphaw=dzw*0;
alphaw=(zwu(:,2:end)-zwu(:,1:end-1))./(ones(nz+1,1)*diff(xu,1,2));

%whos dxu dzu

dxu=diff(x,1,2);
alphauw = zeros(nz+1,nx+1);
alphauw(1:nz+1,2:nx) =(zw(1:nz+1,2:end)-zw(1:nz+1,1:end-1))./(ones(nz+1,1)*dxu);
%alphauw(1,:) = 2*alphauw(2,:)-alphauw(3,:);
%alphauw(nz+1,:) = 0;% free surface is flat
alphauw(:,1) = 2*alphauw(:,2)-alphauw(:,3);
alphauw(:,nx+1) = 2*alphauw(:,nx)-alphauw(:,nx-1);
%alphau=(zr(:,2:end)-zr(:,1:end-1))./(ones(nz,1)*dxu);

%% compute - d_i U = - U_i+1/2 + U_i-1/2
% by sweeping over columns at i+1/2:
% - determine the flux U_i+1/2 = (p_i+1 - p_i)*dz/dx
% - ce(:,i  )= +dz/dx
% - co(:,i  )= -dz/dx
%
% - cw(:,i+1)= +dz/dx
% - co(:,i+1)= -dz/dx

ce=zeros(nz,nx);
cw=zeros(nz,nx);
cn=zeros(nz,nx);
cs=zeros(nz,nx);
cen=zeros(nz,nx);
cwn=zeros(nz,nx);
ces=zeros(nz,nx);
cws=zeros(nz,nx);
co=zeros(nz,nx);





% sweep over columns
cff=zeros(nz,1);
for i=1:nx-1
    % compute the flux U
   cff(:)   = dzu(:,i)/dxu(i);
   ce(:,i)  = ce(:,i)  +cff;
   co(:,i)  = co(:,i)  -cff;
   cw(:,i+1)= cw(:,i+1)+cff;
   co(:,i+1)= co(:,i+1)-cff;
   
   % bottom flux correction
  k=1;
  cff2 = -0.5*alphauw(k,i+1)^2 / (1+alphauw(k,i+1)^2) *dzw(k,i)/dxu(i);
  ce(k,i)  = ce(k,i)  +cff2;
  co(k,i)  = co(k,i)  -cff2;
  cw(k,i+1)= cw(k,i+1)+cff2;
  co(k,i+1)= co(k,i+1)-cff2;
   
end

%sweep over rows
cff=zeros(1,nx);
for k=1:nz-1
    cff(:)   = dxw(:)'./dzw(k+1,:).*(1+alphaw(k+1,:).^2);
    cn(k,:)  = cn(k,:)  +cff;
    co(k,:)  = co(k,:)  -cff;
    cs(k+1,:)= cs(k+1,:)+cff;
    co(k+1,:)= co(k+1,:)-cff;
end
k=nz+1;
cff=2*dxw(:)'./dzw(k,:);
co(nz,:)=co(nz,:)-cff; % Dirichlet cond at the upper boundary (2nd order)
 


%% cross term
if 1==1
im=1:nx-1;
ip=2:nx;
for k=1:nz-1
    cff=-0.5*alphauw(k+1,2:nx);
%    c=0;
    s=1;
    co(k ,im)  = co(k ,im) - s*cff;
%    cn(k ,im)  = cn(k ,im) - s*cff*c;
    cen(k,im)  = cen(k,im) + s*cff;
%    ce(k ,im)  = ce(k ,im) - s*cff*c;

    s=-1;
%    cw(k ,ip)  = cw(k ,ip) - s*cff*c;
    cwn(k,ip)  = cwn(k,ip) + s*cff;
%    cn(k ,ip)  = cn(k ,ip) - s*cff*c;
    co(k ,ip)  = co(k ,ip) - s*cff;
    
    s=1;
    cws(k+1,ip)= cws(k+1,ip)+ cff*s;
%    cw(k+1,ip) = cw(k+1,ip) - cff*s*c;
    co(k+1,ip) = co(k+1,ip) - cff*s;
%    cs(k+1,ip) = cs(k+1,ip) - cff*s*c;
    
    s=-1;
%    cs(k+1 ,im)  = cs(k+1 ,im) - s*cff*c;
    co(k+1 ,im)  = co(k+1 ,im) - s*cff;
%    ce(k+1 ,im)  = ce(k+1 ,im) - s*cff*c;
    ces(k+1,im)  = ces(k+1,im) + s*cff;
    
    
end
end
%% define sparse matrix
i=1:nx;im=1:nx-1;ip=2:nx;
k=1:nz;km=1:nz-1;kp=2:nz;

siz=[nz,nx];
%
I=[];
J=[];
s=[];
% east
c=ce(k,im);
[ii,kk]=meshgrid(ip,k);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% west
c=cw(k,ip);
[ii,kk]=meshgrid(ip,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,k);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% north
c=cn(km,i);
[ii,kk]=meshgrid(i,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(i,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% south
c=cs(kp,i);
[ii,kk]=meshgrid(i,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(i,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];
% central
c=co(k,i);
[ii,kk]=meshgrid(i,k);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;I0];
s=[s;c(:)];

% northeast
c=cen(km,im);
[ii,kk]=meshgrid(ip,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];

% southeast
c=ces(kp,im);
[ii,kk]=meshgrid(ip,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(im,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];

% northwest
c=cwn(km,ip);
[ii,kk]=meshgrid(im,kp);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(ip,km);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];

% southwest
c=cws(kp,ip);
[ii,kk]=meshgrid(im,km);ii=ii(:);kk=kk(:);
J0=sub2ind(siz,kk,ii);
[ii,kk]=meshgrid(ip,kp);ii=ii(:);kk=kk(:);
I0=sub2ind(siz,kk,ii);
I=[I;I0];
J=[J;J0];
s=[s;c(:)];



A=sparse(I,J,s,nx*nz,nx*nz);



