%Calculate Variance for a given interval of ITD Curve 
%JW 11/6/09

ITD
lower= input('Input lower limit   ');
upper= input('Input upper limit   ');

for v=1:length(ITD)
    if lower==ITD(v)
        lower=v;
    end
    if upper==ITD(v)
        upper=v;
    end
    v=v+1;
end

variance=mean(var(spikecount_sorted2(:,lower:upper)))/max(max(spikecount_sorted2(:,lower:upper)))