% For analyzing control data - to determine if neuron is motion sensitive
% Modified from Osc_ITd_Curve
% 1/5/10   YY Wang

clear all;
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
subplot(212)
plot(data{1}.datatrace, 'b');
hold on
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
%Spike Count Sorting
%Sort into Phases

%spikecountSorted(1,1)=spikeCounts(1,:); %Phase 1 trial 1 is all zeros
%because a sound is played here
spikecountSorted(1,2)=spikeCounts(2);
spikecountSorted(1,3)=spikeCounts(3);
spikecountSorted(1,4)=spikeCounts(4);
place=2;
for m=1:(trials/phase)-1
    spikecountSorted(place,1) = spikeCounts(m*phase+1);
    spikecountSorted(place,2) = spikeCounts(m*phase+2);
    spikecountSorted(place,3) = spikeCounts(m*phase+3);
    spikecountSorted(place,4) = spikeCounts(m*phase+4);
    m=m+1;
    place=place+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Mean Spikes for Each Phase%%%%%%%%%%%%%%%%%%
figure
hold on
%Phase 1 blue
errorbar(mean(spikecountSorted(:,1)),std(spikecountSorted(:,1)/sqrt(trials/4)),'b');
%Phase 2 green
errorbar(mean(spikecountSorted(:,2)),std(spikecountSorted(:,2)/sqrt(trials/4)),'g');
%Phase 3 red
errorbar(mean(spikecountSorted(:,3)),std(spikecountSorted(:,3)/sqrt(trials/4)),'r');
%Phase 4 black
errorbar(mean(spikecountSorted(:,4)),std(spikecountSorted(:,4)/sqrt(trials/4)),'k');
xlim([0.85 1.15]);
title('Average Spikes for Each Phase');
ylabel('Avg # Spike');
legend('Phase 1','Phase 2','Phase 3', 'Phase 4');
hold off


%%%%%%%%%%%%%Spike Times Sorting - sort into individual phases(1-4)%%%%%%%% 

%Phase1_spiketimes(1,:)=spikeTimes(1,:); %Phase 1 trial 1 is all zeros
%because a sound is played here
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
xlim([1 length(data{5}.datatrace)-2])  % sometimes trace1 is not complete, due to buffer limitations

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
xlim([1 length(data{5}.datatrace)-2])

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
xlim([1 length(data{5}.datatrace)-2])

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
xlim([1 length(data{5}.datatrace)-2])

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
%%%%%%%%%%%%%%%%%%%%%%Kruskall-Wallis Test btwn Phases%%%%%%%%%%%%%%%%%%%%%
kw=zeros(1,6);
group=['s';'s']; %argument for krusalwallis so plots don't show
%Phase 1,3
kw(1)=kruskalwallis([spikecountSorted(:,1) spikecountSorted(:,3)],group,'off');
%Phase 2,4
kw(2)=kruskalwallis([spikecountSorted(:,2) spikecountSorted(:,4)],group,'off');
%Phase 1,2
kw(3)=kruskalwallis([spikecountSorted(:,1) spikecountSorted(:,2)],group,'off');
%Phase 3,4
kw(4)=kruskalwallis([spikecountSorted(:,3) spikecountSorted(:,4)],group,'off');
%Phase 1,4
kw(5)=kruskalwallis([spikecountSorted(:,1) spikecountSorted(:,4)],group,'off');
%Phase 2,3
kw(6)=kruskalwallis([spikecountSorted(:,2) spikecountSorted(:,3)],group,'off');
disp('Phase 1,3 	2,4	   1,2	 3,4	  1,4	 2,3');
kw


