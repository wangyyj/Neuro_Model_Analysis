%------------------------------------------------------------------------
% raster.m
%------------------------------------------------------------------------
% function h = raster(x,Fs); x are spiketimes(samples)
% Makes raster plot of matrix x(m,n) with m trials
% YY Wang
%------------------------------------------------------------------------
function h = raster(x,Fs);
Fs=Fs/1000; % plot in ms
% figure('Color',[1 1 1]);
hold on
counter=1;
[L r]=size(x);
for rast_i=1:L
    for rast_m=1:r
        line([x(rast_i,rast_m) x(rast_i,rast_m)]/Fs,[counter-1 counter],'LineWidth',1,'Color','k')
    end
    counter=counter+1;
end
% xlabel('Time (ms)')
% ylabel('Trials')
