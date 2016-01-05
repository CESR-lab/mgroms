function [ur,wr]=momentum2real(u,w,alphauw)

% compute the components of the "real" velocity (=horiz/vert) from the momentum from 
nx=size(w,2);
nz=size(u,1);

ur=u;
wr=w;


for i=1:nx-1
    cff = -0.25*alphauw(:,i+1);
    ur(:,i) = ur(:,i) + cff(1:nz)  .*( w(1:nz,i)+w(1:nz,i+1));
    ur(:,i) = ur(:,i) + cff(2:nz+1).*( w(2:nz+1,i)+w(2:nz+1,i+1));
end

