%% test convergence rate in sigma coord
clear conv
conv=[];
pps=0:3;
for pp=pps
    nz=40*2^pp;
    nx=100*2^pp;
    enforce_divergencefree;
    conv(pp+1)=utop; % get the horizontal gradient at the seamount top
end

%% monitor how utop is converging with resolution (slope should be -2 because of second order accuracy)
figure(5)
loglog(2.^(pps(2:end)),diff(conv),'+-')