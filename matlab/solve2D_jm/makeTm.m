function [T,T2] = makeTm(nx,nz,zx,zxw,dx,dzu,dzw)

%% Flux representation: (u1,w1)

%% Fom Flux to real
%% ur = uf
%% wr = uf*z_x + wf

%% Ek in terms of uf,wf
%%
%%  Ek = 0.5*uf^2  + 0.5*(uf*z_x + wf)^2

%% find (u2,w2) such that Ek = u1u2 + w1w2

%%  um = Ek_uf = uf*(1+zx^2) + wf*z_x
%%  wm = Ek_wf = uf*z_x       + wf

%%  Invert this to find:
%%
%%  uf =  um     - wm*zx
%%  wf = -um*z_x + wm*(1+zx^2)


ndu = (nx+1)*nz;
ndw = nx*(nz+1);



dxu = zeros(nx+1,1);
dxu(2:nx) = 0.5*(dx(1:nx-1)+dx(2:nx));
dxu(1) = dxu(2); dxu(nx+1) = dxu(nx);

if 0
    tic
T = sparse(ndu+ndw,ndu+ndw,10*ndu);
is = 1;
ks = nx;
for i = 1:nx+1
   k=1;
   udx = (k-1)*(nx+1) + i-1;
   wdx = (k-1)*(nx+0) + i-1 + ndu;
   T(udx+1,udx+1) = dxu(i)*dzu(i,k)*(1-0.5*zx^2/(1+zx^2));   %%  uf =  um - wm*z_x
   if i>1&i<nx+1 %%
     %% use:  wf = -um*z_x + wm*(1+zx^2) = 0 -> wm = um *zx/(1+zx^2)
     T(udx+1,wdx+1   +ks) = -0.25*zx(i  ,k)*dxu(i)*dzw(i  ,k+1);
     T(udx+1,wdx+1-is+ks) = -0.25*zx(i-1,k)*dxu(i)*dzw(i-1,k+1);
   end
   for k = 2:nz
      udx = (k-1)*(nx+1) + i-1;
      wdx = (k-1)*(nx+0) + i-1 + ndu;
      
      T(udx+1,udx+1) = dxu(i)*dzu(i,k);       %%  uf =  um - wm*z_x
      if i>1&i<nx+1 %%
        T(udx+1,wdx+1      ) = -0.25*zx(i  ,k)*dxu(i)*dzw(i  ,k);
        T(udx+1,wdx+1-is   ) = -0.25*zx(i-1,k)*dxu(i)*dzw(i-1,k);
        T(udx+1,wdx+1   +ks) = -0.25*zx(i  ,k)*dxu(i)*dzw(i  ,k+1);
        T(udx+1,wdx+1-is+ks) = -0.25*zx(i-1,k)*dxu(i)*dzw(i-1,k+1);
      end
      if i==1   %%  If the flux u(1,k) == 0 these guys are irrelevant!!!
        T(udx+1,wdx+1   ) = -0.25*zx(i,k)*dxu(i)*dzw(i,k);
        T(udx+1,wdx+1+ks) = -0.25*zx(i,k)*dxu(i)*dzw(i,k+1);   
      end
      if i==nx+1
        T(udx+1,wdx+1-is   ) = -0.25*zx(i-1,k)*dxu(i)*dzw(i-1,k);
        T(udx+1,wdx+1-is+ks) = -0.25*zx(i-1,k)*dxu(i)*dzw(i-1,k+1);
      end
   end
end

ks = nx+1;
for i = 1:nx
   for k = 1:nz+1
      udx = (k-1)*(nx+1) + i-1;
      wdx = (k-1)*(nx+0) + i-1 + ndu;
      
      if k<nz+1
       T(wdx+1,wdx+1) = dx(i)*dzw(i,k)*(1.0+zxw(i,k)*zxw(i,k)); %% wf = -zx*um+(1+zx^2)*wm
      else
        T(wdx+1,wdx+1) = 0.5*dx(i)*dzw(i,k)*(1.0+zxw(i,k)*zxw(i,k));
      end
      if k == 1
         T(wdx+1,wdx+1) = 0.5*dx(i)*dzw(i,k)*(1.0+zxw(i,k)*zxw(i,k));
      end
      if k>1&&k<nz+1
        T(wdx+1,udx+1   -ks)=  -0.25*zx(i,k-1)*dzw(i,k)*dxu(i  );
        T(wdx+1,udx+1      )=  -0.25*zx(i,k  )*dzw(i,k)*dxu(i  );
        T(wdx+1,udx+1+is-ks)=  -0.25*zx(i,k-1)*dzw(i,k)*dxu(i+1);
        T(wdx+1,udx+1+is   )=  -0.25*zx(i,k  )*dzw(i,k)*dxu(i+1);
      end
      if k==1  
        T(wdx+1,udx+1   )   =  -0.25*zx(i,k  )*dzw(i,k)*dxu(i  );
        T(wdx+1,udx+1+is   )=  -0.25*zx(i,k  )*dzw(i,k)*dxu(i+1);    
      end
      if k==nz+1
        T(wdx+1,udx+1   -ks)=  -0.25*zx(i,k-1)*dzw(i,k)*dxu(i  );
        T(wdx+1,udx+1+is-ks)=  -0.25*zx(i,k-1)*dzw(i,k)*dxu(i+1);
      end
      
   end
end
toc
else
tic

   cuucc = zeros(nx+1,nz);
   cuwcc = zeros(nx+1,nz);
   cuwwc = zeros(nx+1,nz);
   cuwcu = zeros(nx+1,nz);
   cuwwu = zeros(nx+1,nz);
   
   cuucc = repmat(dxu,1,nz).*dzu;   
   cuwcc(1:nx,:) =  -0.25.*zx.*repmat(dxu(1:nx  ),1,nz).*dzw(:,1:nz);
   cuwcu(1:nx,:) =  -0.25.*zx.*repmat(dxu(1:nx  ),1,nz).*dzw(:,2:nz+1);
   cuwwc(2:nx+1,:)= -0.25.*zx.*repmat(dxu(2:nx+1),1,nz).*dzw(:,1:nz);
   cuwwu(2:nx+1,:)= -0.25.*zx.*repmat(dxu(2:nx+1),1,nz).*dzw(:,2:nz+1);
   
   cwwcc = zeros(nx,nz+1);
   cwucc = zeros(nx,nz+1);
   cwuec = zeros(nx,nz+1);
   cwucd = zeros(nx,nz+1);
   cwued = zeros(nx,nz+1);
   
   cwwcc = repmat(dx,1,nz+1).*dzw.*(1+zxw.*zxw); 
   cwwcc(:,1   ) = 0.5*cwwcc(:,1);
   cwwcc(:,nz+1) = 0.5*cwwcc(:,nz+1);
   
   cwucc(:,1:nz  ) = -0.25.*zx.*repmat(dxu(1:nx  ),1,nz).*dzw(:,1:nz);
   cwuec(:,1:nz  ) = -0.25.*zx.*repmat(dxu(2:nx+1),1,nz).*dzw(:,1:nz);
   cwucd(:,2:nz+1) = -0.25.*zx.*repmat(dxu(1:nx  ),1,nz).*dzw(:,2:nz+1);
   cwued(:,2:nz+1) = -0.25.*zx.*repmat(dxu(2:nx+1),1,nz).*dzw(:,2:nz+1);
   
   nnz = 5*(ndu + ndw);
   iT = zeros(nnz,1);
   jT = zeros(nnz,1);
   cT = zeros(nnz,1); 

   is = 1; ks = nx;
   
   [ii,kk] = meshgrid([1:nx+1],[1:nz]);
   ii1 = reshape(ii',(nx+1)*nz,1);
   kk1 = reshape(kk',(nx+1)*nz,1);
   udx = (kk1-1)*(nx+1) + ii1;
   wdx = (kk1-1)*(nx+0) + ii1 + ndu; 

   iT(1+0*ndu:1*ndu) = udx;
   jT(1+0*ndu:1*ndu) = udx;
   cT(1+0*ndu:1*ndu) = reshape(cuucc,ndu,1);
   
   iT(1+1*ndu:2*ndu) = udx;
   jT(1+1*ndu:2*ndu) = wdx;
   cT(1+1*ndu:2*ndu) = reshape(cuwcc,ndu,1);
   
   iT(1+2*ndu:3*ndu) = udx;
   jT(1+2*ndu:3*ndu) = wdx +ks;
   cT(1+2*ndu:3*ndu) = reshape(cuwcu,ndu,1);
   
   iT(1+3*ndu:4*ndu) = udx;
   jT(1+3*ndu:4*ndu) = wdx -is;
   cT(1+3*ndu:4*ndu) = reshape(cuwwc,ndu,1);
   
   iT(1+4*ndu:5*ndu) = udx;
   jT(1+4*ndu:5*ndu) = wdx -is+ks;
   cT(1+4*ndu:5*ndu) = reshape(cuwwu,ndu,1);
   
   is = 1; ks = nx+1;
   [ii,kk] = meshgrid([1:nx],[1:nz+1]);
   ii1 = reshape(ii',nx*(nz+1),1);
   kk1 = reshape(kk',nx*(nz+1),1);
   udx = (kk1-1)*(nx+1) + ii1;
   wdx = (kk1-1)*(nx+0) + ii1 + ndu; 
   
   ib = 5*ndu+1+0*ndw; ie = 5*ndu+1*ndw;
   iT(ib:ie) = wdx;
   jT(ib:ie) = wdx;
   cT(ib:ie) = reshape(cwwcc,ndw,1);
      
   ib = 5*ndu+1+1*ndw; ie = 5*ndu+2*ndw;
   iT(ib:ie) = wdx;
   jT(ib:ie) = udx;
   cT(ib:ie) = reshape(cwucc,ndw,1);
      
   ib = 5*ndu+1+2*ndw; ie = 5*ndu+3*ndw;
   iT(ib:ie) = wdx;
   jT(ib:ie) = udx +is;
   cT(ib:ie) = reshape(cwuec,ndw,1);
   
   ib = 5*ndu+1+3*ndw; ie = 5*ndu+4*ndw;  
   iT(ib:ie) = wdx;
   jT(ib:ie) = udx -ks;
   cT(ib:ie) = reshape(cwucd,ndw,1);
   
   ib = 5*ndu+1+4*ndw; ie = 5*ndu+5*ndw;  
   iT(ib:ie) = wdx;
   jT(ib:ie) = udx +is-ks;
   cT(ib:ie) = reshape(cwued,ndw,1);
   
   iT(cT==0) = [];
   jT(cT==0) = [];
   cT(cT==0) = [];
   
   T = sparse(iT,jT,cT,ndu+ndw,ndu+ndw,nnz);
toc
end