function D = makeD(nx,nz,dxiu,dziw)



nd = nx*nz;
ndu = (nx+1)*nz;
ndw = nx*(nz+1);



if 0
tic
D = sparse(nd+nx,ndu+ndw,4*nd+nx);
  
for i = 1:nx
   idx = i-1;
   wdx = i-1 + ndu;
   D(idx+1,wdx+1) = dziw(i,1);
end
is = 1;
ks = nx+1;
for i = 1:nx
   for k = 1:nz
      idx = (k-1)*nx + i-1 +nx;
      udx = (k-1)*(nx+1) + i-1;
      if i>1
        D(idx+1,udx+1)    = -dxiu(i);
      end
      if i<nx
        D(idx+1,udx+1+is) =  dxiu(i+1);
      end
   end
end

ks = nx;
for i = 1:nx
   for k = 1:nz
      idx = (k-1)*nx + i-1 +nx;     
      wdx = (k-1)*(nx+0) + i-1 + ndu;

      D(idx+1,wdx+1)    = -dziw(i,k);

      if k<nz
       D(idx+1,wdx+1+ks) =  dziw(i,k+1);
      else
       D(idx+1,wdx+1+ks) =  2*dziw(i,k+1);
      end
          
   end
end
toc

else
tic
nnz = 4*nx*nz + nx;
iD = zeros(nnz,1);
jD = zeros(nnz,1);
cD = zeros(nnz,1);

[ii,kk] = meshgrid([1:nx],[1:nz]);
ii1 = reshape(ii',nx*nz,1);
kk1 = reshape(kk',nx*nz,1);
idx = (kk1-1)*(nx  ) + ii1 +nx;
udx = (kk1-1)*(nx+1) + ii1;
wdx = (kk1-1)*(nx+0) + ii1 + ndu;

%% Bottom Flux is zero
 iD(1:nx) = [1:nx];
 jD(1:nx) = [1:nx] + ndu;
 cD(1:nx) = dziw(:,1);
 
%% Divergence
 cuc = -repmat(dxiu(1:nx),1,nz);cuc(1,:) = 0;
 cue = repmat(dxiu(2:nx+1),1,nz);cue(nx,:) = 0;
 cwc = -dziw(:,1:nz);
 cwu = dziw(:,2:nz+1);cwu(:,nz) = 2*dziw(:,nz+1);

 is = 1;
 iD(nx+0*nd+1:nx+1*nd) = idx;
 jD(nx+0*nd+1:nx+1*nd) = udx;
 cD(nx+0*nd+1:nx+1*nd) = reshape(cuc,nd,1);
 iD(nx+1*nd+1:nx+2*nd) = idx;
 jD(nx+1*nd+1:nx+2*nd) = udx + is;
 cD(nx+1*nd+1:nx+2*nd) = reshape(cue,nd,1);
 
 ks = nx;
 iD(nx+2*nd+1:nx+3*nd) = idx;
 jD(nx+2*nd+1:nx+3*nd) = wdx;
 cD(nx+2*nd+1:nx+3*nd) = reshape(cwc,nd,1);
 iD(nx+3*nd+1:nx+4*nd) = idx;
 jD(nx+3*nd+1:nx+4*nd) = wdx+ks;
 cD(nx+3*nd+1:nx+4*nd) = reshape(cwu,nd,1);
 iD(cD==0) = [];
 jD(cD==0) = [];
 cD(cD==0) = [];
 D = sparse(iD,jD,cD,nd+nx,ndu+ndw,nnz);
 
 toc
end


