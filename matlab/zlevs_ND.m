function z = zlevs_ND(h,N,type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,L]=size(h);
%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
%
if type=='w'
  sc=((0:N)-N)/N;
  N=N+1;
else
  sc=((1:N)-N-0.5)/N;
end
%
% Create S-coordinate system: based on model topography h(i,j)
%
z=zeros(N,M,L);
for k=1:N
  z0=sc(k)*h;
  z(k,:,:)=z0;
end
return

