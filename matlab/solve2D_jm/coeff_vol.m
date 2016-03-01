%% construct P directly



%                           % u1 =  u2 - w2*z_x
%                           % w1 = -zx*u2+(1+zx^2)*w2
%
%    w  u     D*T*G          u+  - u-  + w+ - w-
%
%           u+: ( p(i+1,k  ) - p(i  ,k  ))*dzu(i+1,k)*dxiu(i+1) ...
%              -( p(i  ,k  ) - p(i  ,k-1))*zx4(i  ,k) ...
%              -( p(i  ,k+1) - p(i  ,k  ))*zx4(i  ,k) ...
%              -( p(i+1,k  ) - p(i+1,k-1))*zx4(i+1,k) ... 
%              -( p(i+1,k+1) - p(i+1,k  ))*zx4(i+1,k);
%
%           u-:-( p(i  ,k  ) - p(i-1,k  ))*dzu(i  ,k)*dxui(i  ) ...
%              +( p(i-1,k  ) - p(i-1,k-1))*zx4(i-1,k) ...
%              +( p(i-1,k+1) - p(i-1,k)  )*zx4(i-1,k) ...
%              +( p(i  ,k  ) - p(i  ,k-1))*zx4(i  ,k) ...
%              +( p(i  ,k+1) - p(i  ,k)  )*zx4(i  ,k);
%
%           w+: ( p(i  ,k+1) - p(i  ,k  ))*dx(i)*dzwi(i,k+1)*(1+zxw^2) ...
%              -( p(i  ,k  ) - p(i-1,k  ))*zx4(i,k  ) ...
%              -( p(i  ,k+1) - p(i-1,k+1))*zx4(i,k+1) ...
%              -( p(i+1,k  ) - p(i  ,k  ))*zx4(i,k  ) ...
%              -( p(i+1,k+1) - p(i  ,k+1))*zx4(i,k+1);
%    
%           wi:-( p(i  ,k  ) - p(i  ,k-1))*dx(i)*dzwi(i,k  )*(1+zxw^2) ...
%              +( p(i  ,k-1) - p(i-1,k-1))*zx4(i,k-1) ...
%              +( p(i  ,k  ) - p(i-1,k  ))*zx4(i,k  ) ...
%              +( p(i+1,k-1) - p(i  ,k-1))*zx4(i,k-1) ...
%              +( p(i+1,k  ) - p(i  ,k  ))*zx4(i,k  );

%%% Upper Boundary:

%           u+: ( p(i+1,k  ) - p(i  ,k  ))*dzu(i+1,k)*dxiu(i+1) ...
%              -( p(i  ,k  ) - p(i  ,k-1))*zx4(i  ,k) ...
%              -( -2*p(i  ,k  ))*zx4(i  ,k) ...
%              -( p(i+1,k  ) - p(i+1,k-1))*zx4(i+1,k) ... 
%              -( -2*p(i+1,k  ))*zx4(i+1,k);
%
%           u-:-( p(i  ,k  ) - p(i-1,k  ))*dzu(i  ,k)*dxui(i  ) ...
%              +( p(i-1,k  ) - p(i-1,k-1))*zx4(i-1,k) ...
%              +( -2*p(i-1,k)  )*zx4(i-1,k) ...
%              +( p(i  ,k  ) - p(i  ,k-1))*zx4(i  ,k) ...
%              +( -2*p(i  ,k)  )*zx4(i  ,k);
%
%           w+:2*( p(i  ,k+1) - p(i  ,k  ))*dx(i)*dzwi(i,k+1)*(1+zxw^2) ...
%              -( p(i  ,k  ) - p(i-1,k  ))*zx4(i,k  ) ...
%              -( p(i+1,k  ) - p(i  ,k  ))*zx4(i,k  ) ...
%    
%           wi:-( p(i  ,k  ) - p(i  ,k-1))*dx(i)*dzwi(i,k  )*(1+zxw^2) ...
%              +( p(i  ,k-1) - p(i-1,k-1))*zx(i,k-1) ...
%              +( p(i  ,k  ) - p(i-1,k  ))*zx(i,k  ) ...
%              +( p(i+1,k-1) - p(i  ,k-1))*zx(i,k-1) ...
%              +( p(i+1,k  ) - p(i  ,k  ))*zx(i,k  );



%                       
%%% West Boundary:

%           u+: ( p(i+1,k  ) - p(i  ,k  ))*dzu(i+1,k)*dxiu(i+1) ...
%              -( p(i  ,k  ) - p(i  ,k-1))*zx4(i  ,k) ...
%              -( p(i  ,k+1) - p(i  ,k  ))*zx4(i  ,k) ...
%              -( p(i+1,k  ) - p(i+1,k-1))*zx4(i+1,k) ... 
%              -( p(i+1,k+1) - p(i+1,k  ))*zx4(i+1,k);
%
%
%           w+: ( p(i  ,k+1) - p(i  ,k  ))*dx(i)*dzwi(i,k+1)*(1+zxw^2) ...
%              -( p(i+1,k  ) - p(i  ,k  ))*zx4(i,k  ) ...
%              -( p(i+1,k+1) - p(i  ,k+1))*zx4(i,k+1);
%    
%           wi:-( p(i  ,k  ) - p(i  ,k-1))*dx(i)*dzwi(i,k  )*(1+zxw^2) ...
%              +( p(i+1,k-1) - p(i  ,k-1))*zx(i,k-1) ...
%              +( p(i+1,k  ) - p(i  ,k  ))*zx(i,k  );

%% look at cct and ccb (p(i,k+1), p(i,k-1).
%%  p(i,k+1):    dx(i)*dziw(i,k+1)(1+zxw^2) 
%                + zx4(i,k+1) - zx4(i,k)

%%  p(i,k-1):    dx(i)*dziw(i,k)(1+zxw^2) 
%                - zx4(i,k-1) + zx4(i,k)


%%% East Boundary:


%           u-:-( p(i  ,k  ) - p(i-1,k  ))*dzu(i  ,k)*dxui(i  ) ...
%              +( p(i-1,k  ) - p(i-1,k-1))*zx4(i-1,k) ...
%              +( p(i-1,k+1) - p(i-1,k)  )*zx4(i-1,k) ...
%              +( p(i  ,k  ) - p(i  ,k-1))*zx4(i  ,k) ...
%              +( p(i  ,k+1) - p(i  ,k)  )*zx4(i  ,k);
%
%           w+: ( p(i  ,k+1) - p(i  ,k  ))*dx(i)*dzwi(i,k+1)*(1+zxw^2) ...
%              -( p(i  ,k  ) - p(i-1,k  ))*zx4(i,k  ) ...
%              -( p(i  ,k+1) - p(i-1,k+1))*zx4(i,k+1) ...

%    
%           wi:-( p(i  ,k  ) - p(i  ,k-1))*dx(i)*dzwi(i,k  )*(1+zxw^2) ...
%              +( p(i  ,k-1) - p(i-1,k-1))*zx4(i,k-1) ...
%              +( p(i  ,k  ) - p(i-1,k  ))*zx4(i,k  );


%% look at cct and ccb (p(i,k+1), p(i,k-1).
%%  p(i,k+1):    dx(i)*dziw(i,k+1)(1+zxw^2) 
%                - zx4(i,k+1) + zx4(i,k)

%%  p(i,k-1):    dx(i)*dziw(i,k)(1+zxw^2) 
%                + zx4(i,k-1) - zx4(i,k)


%%% Bottom Flux:

%    
%           wi:-0.5*( p(i  ,k  ) - p(i  ,k-1))*dx(i)*dzwi(i,k  )*(1+zxw^2) ...
%              +( p(i  ,k  ) - p(i-1,k  ))*zx(i,k  ) ...
%              +( p(i+1,k  ) - p(i  ,k  ))*zx(i,k  );



%% Slope in the bottom ghost cells
dz_bot = 2*(zr(:,1)-zw(:,1));
zr_bot = zr(:,1)-dz_bot;
zxu_bot = zeros(nx+1,1);
zxu_bot(2:nx) = (zr_bot(2:nx)-zr_bot(1:nx-1));
zx_bot = 0.5*(zxu_bot(2:nx+1)+zxu_bot(1:nx)).*dxi;

%% Slope in the surface ghost cells
dz_sur = 2*(zw(:,nz+1)-zr(:,nz));
zr_sur = zr(:,nz)+dz_sur;
zxu_sur = zeros(nx+1,1);
zxu_sur(2:nx) = (zr_sur(2:nx)-zr_sur(1:nx-1));
zx_sur = 0.5*(zxu_sur(2:nx+1)+zxu_sur(1:nx)).*dxi;

zx4 = 0.25*zx;

ccc = zeros(nx,nz);
for k =1:nz
   ccc(:,k) = - dzu(2:nx+1,k).*dxiu(2:nx+1)-dzu(1:nx,k).*dxiu(1:nx) ...
                - dx.*dziw(:,k+1).*(1+zxw(:,k+1).^2) ...
                - dx.*dziw(:,k  ).*(1+zxw(:,k).^2);
end

ccb = zeros(nx,nz+1);       
for k =1:nz+1
   ccb(:,k) = dx.*dziw(:,k).*(1+zxw(:,k).^2); 
end
cct = ccb(:,2:nz+1);
ccb = ccb(:,1:nz);

cwc = zeros(nx+1,nz);  % for p(i-1,j,k)
for k = 1:nz
   cwc(:,k) =  dzu(:,k).*dxiu;
end
cec = cwc(2:nx+1,:);
cwc = cwc(1:nx,:);

cwt = zeros(nx,nz);
cwt(2:nx,:) = 0.25*zx(1:nx-1,:);
cwt(:,1:nz-1) = cwt(:,1:nz-1) + 0.25*zx(:,2:nz);
cwt(:,nz) = cwt(:,nz) + 0.25*zx_sur;

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
cet(:,nz) = cet(:,nz) - 0.25*zx_sur;

          
%% bc's at the lef:  no flux and no px, slightly subtle!
ccc(1,:) = ccc(1,:) + cwc(1,:);
cct(1,1:nz-1) = cct(1,1:nz-1) + zx4(1,2:nz) - zx4(1,1:nz-1);
ccb(1,2:nz  ) = ccb(1,2:nz  ) + zx4(1,2:nz) - zx4(1,1:nz-1);

%% bc's at the right
ccc(nx,:) = ccc(nx,:) + cec(nx,:);
cct(nx,1:nz-1) = cct(nx,1:nz-1) - zx4(nx,2:nz) + zx4(nx,1:nz-1);
ccb(nx,2:nz  ) = ccb(nx,2:nz  ) - zx4(nx,2:nz) + zx4(nx,1:nz-1);

%% bc's at the top
ccc(:,nz) = ccc(:,nz) - cct(:,nz);
cwc(:,nz) = cwc(:,nz) - cwt(:,nz);
cec(:,nz) = cec(:,nz) - cet(:,nz);

figure; 
subplot(3,3,1); imagesc(cwt'); axis xy; colorbar
subplot(3,3,2); imagesc(cct'); axis xy; colorbar
subplot(3,3,3); imagesc(cet'); axis xy; colorbar
subplot(3,3,4); imagesc(cwc'); axis xy; colorbar
subplot(3,3,5); imagesc(ccc'); axis xy; colorbar
subplot(3,3,6); imagesc(cec'); axis xy; colorbar
subplot(3,3,7); imagesc(cwb'); axis xy; colorbar
subplot(3,3,8); imagesc(ccb'); axis xy; colorbar
subplot(3,3,9); imagesc(ceb'); axis xy; colorbar

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

%% Bottom boundary condition
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

%% Get a coefficient out of the (P4) matrix for comparison
cf = zeros(nx,nz);
is = 1;
ks = nx;
for i = 2:nx-1
   for k = 2:nz-1
      idx = (k-1)*nx + i-1 +nx;
      cf(i,k) = P(idx+1,idx+1+ks);
%      cf(i,k) = P4(idx+1,idx+1+is);   %% couples with p(i,k+1)
   end
end

u = zeros(nx+1,nz);
w = cos(xw*pi/sizex); w(:,1) = 0;
Xr = uw2x(u,w);  %% physical velocity
%Xm = Trm*Xr;    %% physical to momentum (not neccesary here because u==0)
div = D*T*Xr;
Fbot = zeros(nx,1);%
%rhs = [Fbot' div']';
rhs = div;
p = A\rhs;
p2 = reshape(p,[nx nz+1]);


xp = zeros(nx,nz+1);
zp = zeros(nx,nz+1);
xp(:,2:nz+1) = xr;
zp(:,2:nz+1) = zr;
xp(:,1) = xp(:,2);
zp(:,1) = zr_bot;

