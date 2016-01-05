function [u,w]=real2momentum(ur,wr,alphauw)

% compute the components of the momentum from the "real" velocity (=horiz/vert)

nx=size(wr,2);
nz=size(ur,1);

u=ur;
w=wr;


for i=1:nx-1
    cff = 0.25*alphauw(:,i+1);
    u(:,i) = u(:,i) + cff(1:nz)  .*( w(1:nz,i)+w(1:nz,i+1));
    u(:,i) = u(:,i) + cff(2:nz+1).*( w(2:nz+1,i)+w(2:nz+1,i+1));
end


