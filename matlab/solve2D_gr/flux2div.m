function div = flux2div(U,W)

div = zeros(size(U,1),size(W,2));

div(:,2:end-1) =                 U(:,2:end)-U(:,1:end-1);
div(:,1)       =div(:,1)        +U(:,1);
div(:,end)     =div(:,end)      -U(:,end);

div(2:end,:) = div(2:end,:) + W(3:end,:)-W(2:end-1,:);
div(1,:)       = div(1,:)   + W(2,:);
%div(end,:)     = div(end,:)     - W(end,:);
