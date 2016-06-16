%%%%%%%%%%%%%%%%%%ISI.mat%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short function to plot histogram for interstimulus intervals
isi_n=1;
for i=1:22 %do up to 22 only, since no data for 23
    for n=1:length(StimulusPulseTimes_trimmed)-1
        isi(1,isi_n)= StimulusPulseTimes_trimmed(1,n+1)-StimulusPulseTimes_trimmed(1,n);
        isi_n=isi_n+1;
    end
end
isi(2,:)=isi(1,:)*1000/Stimulus.Fs;
