function y= strf(data,Fs,cf)
%------------------------------------------------------------------------
% STRF.m
%------------------------------------------------------------------------
% Plots spectraltemporal power of stimulus 
%------------------------------------------------------------------------

fcoefs=MakeERBFilters(Fs,cf); %cf: vector containing center frequencies

strf=ERBFilterBank(data,fcoefs);

for nn=1:size(strf,1)
    strf(nn,:)=abs(hilbert(strf(nn,:)));
end


for t=2:size(data,1)
    tm=ERBFilterBank(data(t,:),fcoefs);
    for nn=1:size(strf,1)
        tm(nn,:)=abs(hilbert(tm(nn,:)));
    end
    strf=strf+tm;

end
figure;
imagesc(strf);







