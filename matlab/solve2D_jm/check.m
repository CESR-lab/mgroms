%% Check 2 different approaches to bounday condition implementation.
%%
%% Change on line 24 to switch. We run 100 relaxations each time with a 
%% random rhs and random initial p. We do this 100 times to get statistics.

sol = [];
res = [];

nx = 40;
ny = 40;


res1 = zeros(nx,ny);
res2 = zeros(nx,ny);
nit = 100;
onethird = 1./3.;

cc = -4*ones(nx,ny);
ce = ones(nx,ny);
cw = ones(nx,ny);
cs = ones(nx,ny);
cn = ones(nx,ny);

if 0 %% Set to 1 to include bc's in matrix coefficients
  cc(2,:) = cc(2,:)+1;
  cc(:,2) = cc(:,2)+1;
  cc(nx-1,:) = cc(nx-1,:) +1;
  cc(:,ny-1) = cc(:,ny-1) +1;
  ce(nx-1,:) = 0;
  cw(2,:) = 0;
  cs(:,2) = 0;
  cn(:,ny-1) = 0;
end

for k = 1:100
 p = rand(nx,ny);
 b = 1.*rand(nx,ny);
 b = b - mean(mean(b(2:nx-1,2:ny-1)));
 for it = 1:nit
   %% Fill halo
   p(1,:)  = p(2,:);
   p(nx,:) = p(nx-1,:);
   p(:,1)  = p(:,2);
   p(:,ny) = p(:,ny-1);
   
   %% Compute residual
   for i = 2:nx-1
      for j=2:ny-1
         res1(i,j) = b(i,j) - cc(i,j)*p(i,j) ...
                     - cw(i,j)*p(i-1,j) - ce(i,j)*p(i+1,j) ...
                     - cs(i,j)*p(i,j-1) - cn(i,j)*p(i,j+1);
      end
   end
   
   % relax
   for i = 2:nx-1
     for j=2:ny-1
       p(i,j) =-(cw(i,j)*p(i-1,j)+ce(i,j)*p(i+1,j)+...
                 cs(i,j)*p(i,j-1)+cn(i,j)*p(i,j+1) - b(i,j))/cc(i,j);        
     end
   end
   if 1
   %% Compute residual
   for i = 2:nx-1
      for j=2:ny-1
          res2(i,j) = b(i,j) - cc(i,j)*p(i,j) ...
                      - cw(i,j)*p(i-1,j) - ce(i,j)*p(i+1,j) ...
                      - cs(i,j)*p(i,j-1) - cn(i,j)*p(i,j+1);
      end
   end
   end
   if mod(it,nit)==0
      % it
    figure(1)
    subplot(1,3,1)
    imagesc(res1(2:nx-1,2:ny-1)');colorbar;
    title('residual after halo update')
    subplot(1,3,2)
    imagesc(res2(2:nx-1,2:ny-1)');colorbar;
    title('residual before halo update')
    subplot(1,3,3)
    imagesc(p(2:nx-1,2:ny-1)');colorbar;
    title('solution')
    drawnow
    
   end
   
 end
 pm = mean(mean(p));
 res = [res max(max(abs(res2))) ];
 sol = [sol mean(mean( (p-pm).^2))];
end