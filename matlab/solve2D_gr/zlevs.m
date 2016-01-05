function z = zlevs(h,zeta,theta_s,theta_b,hc,N,type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  pierrick 2002                                                            %
%                                                                           %
% function z = zlevs(h,zeta,theta_s,theta_b,hc,N,type);                     %
%                                                                           %
% this function compute the depth of rho or w points for ROMS               %
%                                                                           %
% On Input:                                                                 %
%                                                                           %
%    type    'r': rho point 'w': w point                                    %
%                                                                           %
% On Output:                                                                %
%                                                                           %
%    z       Depths (m) of RHO- or W-points (3D matrix).                    %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M,L]=size(h);
%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
%
cff1=1./sinh(theta_s);
cff2=0.5/tanh(0.5*theta_s);
if type=='w'
  sc=((0:N)-N)/N;
  N=N+1;
else
  sc=((1:N)-N-0.5)/N;
end
Cs=(1.-theta_b)*cff1*sinh(theta_s*sc)...
    +theta_b*(cff2*tanh(theta_s*(sc+0.5))-0.5);
%
% Create S-coordinate system: based on model topography h(i,j),
% fast-time-averaged free-surface field and vertical coordinate
% transformation metrics compute evolving depths of of the three-
% dimensional model grid.
%    
hinv=1./h;
cff=hc*(sc-Cs);
cff1=Cs;
cff2=sc+1;
z=zeros(N,M,L);
for k=1:N
  z0=cff(k)+cff1(k)*h;
  z(k,:,:)=z0+zeta.*(1.+z0.*hinv);
end
return

