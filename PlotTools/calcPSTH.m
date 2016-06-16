%------------------------------------------------------------------------
% CalcPSTH.m
%------------------------------------------------------------------------
% function h = plotpsth(x,bin,start,end); x are spiketimes
% Plots psth with bin, tracelength all units must be the same (samples)
% YY Wang
%------------------------------------------------------------------------

function [psth,time] = plotpsth(x,bin,start_t,end_t);
maxLength=end_t-start_t;
bins=ceil(maxLength/bin);
time=start_t:bin:end_t;
for n=1:length(time)-1
    psth(n)=length(find(x>=time(n) & x<=time(n+1)));
end

