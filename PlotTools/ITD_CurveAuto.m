% ITD_CurveAuto
function [norm] = ITD_CurveAuto(file,thres)
% For loading ITD data and plotting ITD curves obtained from HPSearch Program 
file=[file '.dat']
[data, datainfo] = readHPData(file);
rep=datainfo.curve.nreps;
trials=datainfo.curve.nTrials;  %number of ITDs presented

%Calculate spikes
spiketimes=zeros(length(data),50);
spikecounts=zeros(length(data),1);
n=1;
for n=1:length(data)
    spikes = spikeschmitt2(data{n}.datatrace(1220:4150), thres, 1,datainfo.indev.Fs); %set analysis window here
    spikecounts(n,1)=length(spikes);
    spiketimes(n,1:length(spikes))=spikes;
    n=n+1;    
end

%sort spikecount into reps
i=1; %sorted rows
j=1; %unsorted index
for m=1:rep
    spikecount_sorted(i,1:trials)=spikecounts(j:j+trials-1);
    i=i+1;
    j=j+trials;
end

%sort spikecount according to ITD
order = datainfo.curve.trialRandomSequence;
row=1;
col=1;
spikecount_sorted2=zeros(rep,trials);
for row=1:rep
    for col=1:trials
        spikecount_sorted2(row,order(row,col))= spikecount_sorted(row,col);
        col=col+1;
    end
    row=row+1;
end

%Normalize to Peak=1
norm=scaledata(mean(spikecount_sorted2),0,1);

%plot ITD curve
figure('Color',[1 1 1]);
ITD = datainfo.curve.ITDrange;
errorbar(ITD,mean(spikecount_sorted2),std(spikecount_sorted2)/sqrt(rep)); % standard error of the mean
title('ITD Curve');
xlabel('ITD (us)');
ylabel('# Spike');
