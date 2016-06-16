% For loading ITD data and plotting ITD curves obtained from HPSearch Program 
% Created 7/7/09    Y Wang 

clear;
%Query for file to be analyzed, load data
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file', 'cd');
full_filename=[PathName FileName];
[data, datainfo] = readHPData(full_filename);
rep=datainfo.curve.nreps;
trials=datainfo.curve.nTrials;  %number of ITDs presented

%Demultiplex data (4 channels recorded, only chann2 is real data)
%m=1;    %rows in de-multiplexed data matrix
%for  m=1:trials*rep
%data{m}.datatrace= mcDeMux(data{m}.datatrace, datainfo.tdt.nChannels);
%data{m}.datatrace=data{m}.datatrace(:,datainfo.tdt.MonitorChannel);
%end

%Draw Sample traces
figure('Color',[1 1 1]);
subplot(211)
plot([(1:length(data{1}.datatrace))/datainfo.indev.Fs],data{1}.datatrace, 'b');
subplot(212)
hold on
for n=1:6
plot([(1:length(data{1}.datatrace))/datainfo.indev.Fs],data{n}.datatrace, 'b');
end
title('Raw Data- Trace 1,2,3,4. Please Determine Threshold');
xlabel('Time (s)');
ylabel('Resp (mV)');
grid;
hold off

%Query for Threshold
thres= input('Enter Threshold Value   ');

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

%plot ITD curve
figure('Color',[1 1 1]);
ITD = datainfo.curve.ITDrange;
errorbar(ITD,mean(spikecount_sorted2),std(spikecount_sorted2)/sqrt(rep)); % standard error of the mean
title('ITD Curve');
xlabel('ITD (us)');
ylabel('# Spike');

%Calculate Response Latency
bin_size=0.2; %ms
bin=round((bin_size/1000)*datainfo.indev.Fs); %sampling rate factor
spiketimes(spiketimes==0)=1;
for row=1:length(spiketimes)
    for col=1:50
    psth(row, spiketimes(row,col))=1;
    col=col+1;
    end
    row=row+1;
end

%take out ones, sum
psth_sum=sum(psth);
psth_sum(1)=0;
%bin
psth_binned(1,1)= sum(psth_sum(1:bin));
for psth_n=1:floor(length(psth_sum)/bin)-2
    psth_binned(1,psth_n+1)=sum(psth_sum(psth_n*bin+1:(psth_n+1)*bin));
    psth_n=psth_n+1;
end

%plot psth
psth_time=(1:length(psth_binned))*bin_size;
figure('Color',[1 1 1]);
bar(psth_time, psth_binned);
%Normalize to Peak=1
norm=mean(spikecount_sorted2);
norm=norm/max(norm);

% %Calculate Latency
% [maxfiring,indexmax]= max(psth_binned(datainfo.stim.Delay/bin_size:2*datainfo.stim.Delay/bin_size));
% disp('Latency(ms)   = ');
% latency = indexmax(1)*0.2
