% For loading ITD data and plotting ITD curves
% 9/30/09    YY Wang
% Revision - 1/5/09 added filtering option

clear;
%%%%%%%%%%%%%%%Query for file to be analyzed, load data%%%%%%%%%%%%%%%%%%%%
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file','Y:\Jennifer\Dichotic_Osc\13Oct09');
full_filename=[PathName FileName];
[data, datainfo] = readOwscData(full_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%Constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rep=datainfo.curve.nreps;
trials =datainfo.curve.nTrials;
phase=4;    %total number of phases
Fs= datainfo.indev.Fs;
%%%%%%%%%%%%%%%%%%%%%%%Find Max trace length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lengths=zeros(trials,1);
length_n=1;
for length_n=1:trials
    lengths(length_n)=length(data{length_n}.datatrace);
    lengthn=length_n+1;
end
L_trace=max(lengths);

%%%%%%%%%%%%%%%%%%%%%%%Put raw data in a matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_n=1;
data_raw=zeros(trials,L_trace);
for raw_n=1:trials
data_raw(raw_n,1:length(data{raw_n}.datatrace))= data{raw_n}.datatrace(1:end);
raw_n=raw_n+1;
end
% Rows are trials, columns are datapoints
disp('...Data Loaded');

%%%%%%%%%%Draw Traces, Determine if Filtering needed%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(211)
plot(data{1}.datatrace, 'b');
hold on
subplot(212)
plot(data{1}.datatrace, 'b');
plot(data{2}.datatrace, 'g');
plot(data{3}.datatrace, 'r');
plot(data{4}.datatrace, 'k');
title('Raw Data- Trace 1,2,3,4. Please Determine Threshold');
xlabel('Time');
ylabel('Resp (mV)');
grid;
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Filtering%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filt=0;   %tells if data was filtered and how many times
filt_q=input('Filter Data Before Analysis? (y/n)        ', 's');
while(filt_q=='y')
    HiPass=input('Enter Hi-Pass Value          ');
    LoPass=input('Enter Lo-Pass Value          ');
    filt_n=1;
        for n=1:trials
        trace=data_raw(filt_n,:);
        data_filt(filt_n,:)=filterdata(trace,Fs,HiPass,LoPass);
        filt_n=filt_n+1;
        end
    %% Draw Sample traces again
    figure;
    hold on    
    subplot 411
    plot(data{1}.datatrace, 'b'); 
    title('Trace 1 raw')
    subplot 412    
    plot(data_filt(1,:), 'b');
    title('Trace 1 FILTERED')
    subplot 413    
    plot(data{2}.datatrace, 'g');
    title('Trace 2 raw')
    subplot 414    
    plot(data_filt(2,:), 'g');
    title('Trace 2 FILTERED')
    xlabel('Time');
    ylabel('Resp (mV)');
    grid;
    hold off
    %ask if want to refilter
    filt_q=input('Refilter? (y/n)        ', 's');
    filt=filt+1;
end

if (filt==0)           %%% data_ITD will be used to calculate ITD Curve
    data_ITD=data_raw;
else
    data_ITD=data_filt;
end

%Query for Threshold
    figure;
    plot(data{1}.datatrace, 'b');
    hold on
    plot(data{2}.datatrace, 'g');
    plot(data{3}.datatrace, 'r');
    plot(data{4}.datatrace, 'k');
    title('Determine Threshold');
    xlabel('Time');
    ylabel('Resp (mV)');
    grid;
    hold off
thres= input('Enter Threshold Value   ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate spikes%%%%%%%%%%%%%%%%%%%%%%%%%
spikeTimes=zeros(length(data),50); %%spikes calculated, unsorted
spikeCounts=zeros(length(data),1); %%unsorted spike counts
n=1;
for n=1:length(data)
    spikes= spikeschmitt2(data{n}.datatrace(:), thres, 1,datainfo.indev.Fs);  
    spikeTimes(n,1:length(spikes))=spikes;
    spikeCounts(n)=length(spikes); 
    n=n+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sorting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
order=datainfo.curve.StimSequence;

%Spike Count Sorting
%Sort into Phases and Reps and ITD
spikecountSorted = zeros(phase*rep+rep,length(datainfo.curve.ITDrange)); %actual data starts from row 6 
i=1;
for i=1:trials   
    %sorted matrix: rows= phase*4+rep, order(m,1)=ITD order; order(m,2=phase)   
    spikecountSorted((order(i,2)*rep+order(i,4)), order(i,1))=spikeCounts(i);    
    i=i+1;
end

%Trim data set so phase1 rep1 starts at row 1
spikecountSorted=spikecountSorted(rep+1:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%plot ITD curve - Single panel%%%%%%%%%%%%%%%%%%%%%
figure
hold on
ITD = datainfo.curve.ITDrange;
%Phase 1 Blue
errorbar(ITD,mean(spikecountSorted(1:rep,:)),std(spikecountSorted(1:rep,:))/sqrt(rep), 'b'); % standard error of the mean
%Phase 2 Green
errorbar(ITD,mean(spikecountSorted(rep+1:rep*2,:)),std(spikecountSorted(rep+1:rep*2,:))/sqrt(rep),'g'); 
%Phase 3 Red
errorbar(ITD,mean(spikecountSorted(rep*2+1:rep*3,:)),std(spikecountSorted(rep*2+1:rep*3,:))/sqrt(rep), 'r'); 
%Phase 4 Black
errorbar(ITD,mean(spikecountSorted(rep*3+1:rep*4,:)),std(spikecountSorted(rep*3+1:rep*4,:))/sqrt(rep), 'k'); 
hold off
title('ITD Curve');
xlabel('ITD (us)');
ylabel('# Spike');
legend('Phase 1','Phase 2','Phase 3', 'Phase 4');
grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot ITD Curve%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure;
%Phase 1 Blue
subplot(4,1,1) 
hold on
errorbar(ITD,mean(spikecountSorted(1:rep,:)),std(spikecountSorted(1:rep,:))/sqrt(rep), 'b'); % standard error of the mean
title('Phase 1'); 
ylabel('# Spikes');
grid;
%Phase 2 green
subplot(4,1,2)
errorbar(ITD,mean(spikecountSorted(rep+1:rep*2,:)),std(spikecountSorted(rep+1:rep*2,:))/sqrt(rep),'g'); 
title('Phase 2'); 
ylabel('# Spikes');
grid;
%Phase 3 red
subplot(4,1,3)
errorbar(ITD,mean(spikecountSorted(rep*2+1:rep*3,:)),std(spikecountSorted(rep*2+1:rep*3,:))/sqrt(rep), 'r'); 
title('Phase 3'); 
ylabel('# Spikes');
grid;
%Phase 4 black
subplot(4,1,4)
errorbar(ITD,mean(spikecountSorted(rep*3+1:rep*4,:)),std(spikecountSorted(rep*3+1:rep*4,:))/sqrt(rep), 'k'); 
title('Phase 4'); 
xlabel('ITD (us)');
ylabel('# Spikes');
grid;
hold off


%%%%%%%%%%%%%Calculate Variance - normalized by highest value%%%%%%%%%%%%%%
variance(1,1)=mean(var(spikecountSorted(1:rep,:)))/max(max(spikecountSorted(1:rep,:)));
variance(1,2)=mean(var(spikecountSorted(rep+1:rep*2,:)))/max(max(spikecountSorted(rep+1:rep*2,:)));
variance(1,3)=mean(var(spikecountSorted(rep*2+1:rep*3,:)))/max(max(spikecountSorted(rep*2+1:rep*3,:)));
variance(1,4)=mean(var(spikecountSorted(rep*3+1:rep*4,:)))/max(max(spikecountSorted(rep*3+1:rep*4,:)));
variance
mean(variance)


%%%%%%%%%%%%%Spike Times Sorting - sort into individual phases(1-4)%%%%%%%% 

Phase1_spiketimes(1,:)=spikeTimes(1,:);
Phase2_spiketimes(1,:)=spikeTimes(2,:);
Phase3_spiketimes(1,:)=spikeTimes(3,:);
Phase4_spiketimes(1,:)=spikeTimes(4,:);
place=2;
for m=1:(trials/phase)-1
    Phase1_spiketimes(place,:) = spikeTimes(m*phase+1,:);
    Phase2_spiketimes(place,:) = spikeTimes(m*phase+2,:);
    Phase3_spiketimes(place,:) = spikeTimes(m*phase+3,:);
    Phase4_spiketimes(place,:) = spikeTimes(m*phase+4,:);
    m=m+1;
    place=place+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Raster%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure
hold on
subplot(4,1,1)
title('Phase 1');
counter=1;
for rast_i=1:trials/4
for rast_m=1:length(Phase1_spiketimes(1,:))
	line([Phase1_spiketimes(rast_i,rast_m) Phase1_spiketimes(rast_i,rast_m)],[counter-1 counter])
end
counter=counter+1;
end
ylim([0 trials/4])
xlim([1 length(data{1}.datatrace)-2])

subplot(4,1,2)
title('Phase 2');
counter=1;
for rast_i=1:trials/4
for rast_m=1:length(Phase2_spiketimes(1,:))
	line([Phase2_spiketimes(rast_i,rast_m) Phase2_spiketimes(rast_i,rast_m)],[counter-1 counter])
end
counter=counter+1;
end
ylim([0 trials/4])
xlim([1 length(data{1}.datatrace)-2])

subplot(4,1,3)
title('Phase 3');
counter=1;
for rast_i=1:trials/4
for rast_m=1:length(Phase3_spiketimes(1,:))
	line([Phase3_spiketimes(rast_i,rast_m) Phase3_spiketimes(rast_i,rast_m)],[counter-1 counter])
end
counter=counter+1;
end
ylim([0 trials/4])
xlim([1 length(data{1}.datatrace)-2])

subplot(4,1,4)
title('Phase 4');
counter=1;
for rast_i=1:trials/4
for rast_m=1:length(Phase4_spiketimes(1,:))
	line([Phase4_spiketimes(rast_i,rast_m) Phase4_spiketimes(rast_i,rast_m)],[counter-1 counter])
end
counter=counter+1;
end
ylim([0 trials/4])
xlim([1 length(data{1}.datatrace)-2])

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%More Processing for PSTH plotting%%%%%%%%%%%%%%%%%
%must replace zeros with ones, so we don't have indexing problems
Phase1_spiketimes(Phase1_spiketimes==0)=1;
Phase2_spiketimes(Phase2_spiketimes==0)=1;
Phase3_spiketimes(Phase3_spiketimes==0)=1;
Phase4_spiketimes(Phase4_spiketimes==0)=1;

for row=1:trials/phase
    for col=1:50
    psth1(row, Phase1_spiketimes(row,col))=1;
    psth2(row, Phase2_spiketimes(row,col))=1;
    psth3(row, Phase3_spiketimes(row,col))=1;
    psth4(row, Phase4_spiketimes(row,col))=1;
    col=col+1;
    end
    row=row+1;
end

%take out all the ones placed at the first sample, sum
psth1=sum(psth1);
psth1(1)=0;
psth2=sum(psth2);
psth2(1)=0;
psth3=sum(psth3);
psth3(1)=0;
psth4=sum(psth4);
psth4(1)=0;

%bin
bin_size=0.2; %ms
bin=round((bin_size/1000)*datainfo.indev.Fs); %sampling rate factor
psth_binned1(1,1)= sum(psth1(1:bin));
psth_binned2(1,1)= sum(psth2(1:bin));
psth_binned3(1,1)= sum(psth3(1:bin));
psth_binned4(1,1)= sum(psth4(1:bin));
for psth_n=1:floor(length(psth1)/bin)-phase
    psth_binned1(1,psth_n+1)=sum(psth1(psth_n*bin+1:(psth_n+1)*bin));
    psth_n=psth_n+1;
end
for psth_n=1:floor(length(psth2)/bin)-phase
    psth_binned2(1,psth_n+1)=sum(psth2(psth_n*bin+1:(psth_n+1)*bin));
    psth_n=psth_n+1;
end

for psth_n=1:floor(length(psth3)/bin)-phase
    psth_binned3(1,psth_n+1)=sum(psth3(psth_n*bin+1:(psth_n+1)*bin));
    psth_n=psth_n+1;
end
for psth_n=1:floor(length(psth4)/bin)-phase
    psth_binned4(1,psth_n+1)=sum(psth4(psth_n*bin+1:(psth_n+1)*bin));
    psth_n=psth_n+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot psths%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psth_time1=(1:length(psth_binned1))*0.2;
psth_time2=(1:length(psth_binned2))*0.2;
psth_time3=(1:length(psth_binned3))*0.2;
psth_time4=(1:length(psth_binned4))*0.2;
figure
hold on
subplot(411)
bar(psth_time1, psth_binned1);
title('Phase1');
subplot(412)
bar(psth_time2, psth_binned2);
title('Phase2');
subplot(413)
bar(psth_time3, psth_binned3);
title('Phase3');
subplot(414)
bar(psth_time4, psth_binned4);
title('Phase4');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate Latency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[maxfiring1,indexmax1]= max(psth_binned1(datainfo.stim.Delay/bin_size:2*datainfo.stim.Delay/bin_size));
[maxfiring2,indexmax2]= max(psth_binned2(datainfo.stim.Delay/bin_size:2*datainfo.stim.Delay/bin_size));
[maxfiring3,indexmax3]= max(psth_binned3(datainfo.stim.Delay/bin_size:2*datainfo.stim.Delay/bin_size));
[maxfiring4,indexmax4]= max(psth_binned4(datainfo.stim.Delay/bin_size:2*datainfo.stim.Delay/bin_size));

latency(1) = indexmax1(1)*0.2;
latency(2) = indexmax2(1)*0.2;
latency(3) = indexmax3(1)*0.2;
latency(4) = indexmax4(1)*0.2;

latency
