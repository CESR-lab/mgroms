function [U,W]=momentum2flux(u,w,alphauw,alphaw,dzu,dzw,dxu,dxw)

nx=size(w,2);
nz=size(u,1);

U=u*0;
W=w*0;

u=u.*(ones(nz,1)*dxu);
w=w.*dzw;

for i=1:nx-1
   % compute the flux U
   cff   = dzu(:,i)/dxu(i);
   U(:,i) = u(:,i) .* cff(:);
      
end


for i=1:nx-1
    cff = -0.25*alphauw(:,i+1);
%    U(2:nz,i) = U(2:nz,i) + cff(2:nz).*(w(1:end-1,i)+w(1:end-1,i+1));
%    U(1:nz,i) = U(1:nz,i) + cff(1:nz).*(w(1:end,i)+w(1:end,i+1));
    U(1:nz,i) = U(1:nz,i) + cff(2:nz+1).*( w(2:nz+1,i)+w(2:nz+1,i+1));
    U(2:nz,i) = U(2:nz,i) + cff(2:nz  ).*( w(2:nz  ,i)+w(2:nz  ,i+1));
    
    k=1;
    cff2=-0.5*alphauw(k,i+1).^2 / (1+alphauw(k,i+1)^2) *dzw(k,i)/dxu(i);
    U(1,i)    = U(1,i) + cff2*u(1,i);
    
%    k=1;
%    cff1 = -0.5*alphau(k,i)^2 / (1+alphau(k,i)^2) *dzu(k,i)/dxu(i);
%    U(k,i) = U(k,i) + u(k,i) *cff1;

%    k=nz;
%    cff1= -0.25*alphau(k,i);
%    U(k,i) = U(k,i) + u(k,i) *cff1*(w(k,i)+w(k,i+1));
end

for k=2:nz+1
    cff   = dxw(:)'./dzw(k,:).*(1+alphaw(k,:).^2);
    W(k,:) = w(k,:) .* cff;
end

%uu = u.*alphau;
for k=2:nz
    cff=-0.25*alphauw(k,:);
    W(k,2:nx-1) = W(k,2:nx-1) +cff(2:nx-1).*(u(k-1,1:nx-2)+u(k,1:nx-2));
    W(k,2:nx-1) = W(k,2:nx-1) +cff(3:nx  ).*(u(k-1,2:nx-1)+u(k,2:nx-1));
%    W(k,2:nx-1) = W(k,2:nx-1) +cff(2:nx-1).*u(k,1:nx-2)+cff(3:nx).*u(k,2:nx-1);
%    if k<nz
%    W(k,2:nx-1) = W(k,2:nx-1) +cff(2:nx-1).*u(k+1,1:nx-2)+cff(3:nx).*u(k+1,2:nx-1);
%    else
%    W(k,2:nx-1) = W(k,2:nx-1) +cff(2:nx-1).*u(k,1:nx-2)+cff(3:nx).*u(k,2:nx-1);
%    end
end
k=nz+1;
cff=-0.5*alphauw(k,:);
W(k,2:nx-1) = W(k,2:nx-1) +cff(2:nx-1).*u(k-1,1:nx-2)+cff(3:nx).*u(k-1,2:nx-1);

