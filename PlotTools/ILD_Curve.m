% For loading ILD data and plotting curves 7/7/09



clear;
%Query for file to be analyzed, load data
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file','cd');
full_filename=[PathName FileName];
[data, datainfo] = readHPData(full_filename);
rep=datainfo.curve.nreps;
trials=datainfo.curve.nTrials;  %number of ILDs presented

% Demultiplex data (4 channels recorded, only chann2 is real data)
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
    spikes = spikeschmitt2(data{n}.datatrace, thres, 1,datainfo.indev.Fs);
    spikecounts(n,1)=length(spikes);
    spiketimes(n,1:length(spikes))=spikes;
    n=n+1;
end

%sort spikecount into reps
rep=datainfo.curve.nreps;
trials=datainfo.curve.nTrials;
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
figure('Color',[1 1 1]);;
ILD = datainfo.curve.ILDrange;
errorbar(ILD,mean(spikecount_sorted2),std(spikecount_sorted2)/sqrt(rep)); % standard error of the mean
title('ILD Curve');
xlabel('ILD (dB)');
ylabel('# Spike');
