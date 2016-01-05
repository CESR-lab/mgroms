function KE = integrated_ke(u,w,U,W,vr,dxu,dzw)

% compute the integrated KE = sum (u_r^2 + w_r^2 )*vol 

nz=size(u,1);

u=u.*(ones(nz,1)*dxu);
w=w.*dzw;
w(nz+1,:)=w(nz+1,:)*0.5;

vol = sum(vr(:));

KE = sum( u(:).*U(:) ) + sum( w(:).*W(:) ) ;

KE = 0.5 * KE /(1e4*4e3);%/ vol;

disp(sprintf('Jeroen''s value   : %12.8f m^2/s^2',0.02933817))
disp(sprintf('Guillaume''s value: %12.8f m^2/s^2',KE))