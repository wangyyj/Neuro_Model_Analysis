% For analyzing control data - to determine if neuron is motion sensitive
% Modified from Osc_ITd_Curve
% 1/5/10   YY Wang

clear all;
%%%%%%%%%%%%%%%Query for file to be analyzed, load data%%%%%%%%%%%%%%%%%%%%
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file',cd);
full_filename=[PathName FileName];
disp(full_filename);
[posFileName,posPathName] = uigetfile('*.mat','Select the matching Position File(_posData.mat)',cd);
pos_full_filename=[posPathName posFileName ];
load (pos_full_filename);
disp(pos_full_filename);
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

spikecountSorted(1,1)=spikeCounts(1,:); %Phase 1 trial 1 is all zeros
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
%sorting to get the right phases 8/16/10
spikecountSorted(:,5)= spikecountSorted(:,1);
spikecountSorted(:,1)= spikecountSorted(:,2);
spikecountSorted(:,2)= spikecountSorted(:,3);
spikecountSorted(:,3)= spikecountSorted(:,4);
spikecountSorted(:,4)= spikecountSorted(:,5);
spikecountSorted(:,5)= 0; %clear extra column

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
%sorting to get the right phase 8/16/10
temp=Phase1_spiketimes;
Phase1_spiketimes=Phase2_spiketimes;
Phase2_spiketimes=Phase3_spiketimes;
Phase3_spiketimes=Phase4_spiketimes;
Phase4_spiketimes=temp;

%Concatenate
Phase2_spiketimes=Phase2_spiketimes+L_trace;
Phase2_spiketimes(Phase2_spiketimes==L_trace)=0;
Phase4_spiketimes=Phase4_spiketimes+L_trace;
Phase4_spiketimes(Phase4_spiketimes==L_trace)=0;
right=[Phase1_spiketimes Phase2_spiketimes];
left= [Phase3_spiketimes Phase4_spiketimes];

%Time manipulation using position data
posdata_step=floor(length(posData)/(trials/4));
maxmin=[];
start=1;
for n=1:trials/4
 maxmin(1,n)=max(posData(start:start+posdata_step)); %max
 maxmin(2,n)=min(posData(start:start+posdata_step)); %min
start=start+posdata_step;
end

maxmin(3,:)= (maxmin(1,:)-max(maxmin(1,:)))/(maxmin(1,1)-maxmin(2,1)); % percent difference from inception
maxmin(4,:)= floor(-1*maxmin(3,:)*L_trace*2/1.5);  % stagger amount in samples (/1.5 - time-position equivalent factor)
maxmin(5,1)= (maxmin(1,1)-maxmin(1,trials/4))/(trials/4);   %Stagger Slope
%stagger right
for row_n=1:trials/4
	for col_n=1:length(right)
		if right(row_n,col_n)<=1
			right_stag(row_n,col_n)=0;
		else
		right_stag(row_n,col_n)=right(row_n,col_n)+maxmin(4,row_n);
		end
	end
end

%stagger left
for row_n=1:trials/4
	for col_n=1:length(left)
		if left(row_n,col_n)<=1
			left_stag(row_n,col_n)=0;
		else
		left_stag(row_n,col_n)=left(row_n,col_n)-maxmin(4,row_n);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Raster%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot raster - original data
figure
hold on
subplot(2,1,1)
title('Right (Phase 1,2)     (Original Data)');
counter=1;
for rast_i=1:trials/4
    for rast_m=1:length(right(1,:))
        line([right(rast_i,rast_m) right(rast_i,rast_m)],[counter-1 counter])
    end
    counter=counter+1;
end
ylim([0 trials/4])

subplot(2,1,2)
title('Left (Phase 3,4)     (Original Data)');
counter=1;
for rast_i=1:trials/4
for rast_m=1:length(left(1,:))
	line([left(rast_i,rast_m) left(rast_i,rast_m)],[counter-1 counter])
end
counter=counter+1;
end
ylim([0 trials/4])
%xlim([1 length(data{5}.datatrace)-2])

hold off

% Plot raster - position-compensated data
figure
hold on
subplot(2,1,1)
title('Right (Phase 1,2)     (Position Corrected)');
counter=1;
for rast_i=1:trials/4
    for rast_m=1:length(right_stag(1,:))
        line([right_stag(rast_i,rast_m) right_stag(rast_i,rast_m)],[counter-1 counter])
    end
    counter=counter+1;
end
ylim([0 trials/4])
%xlim([max(maxmin(4,:)) L_trace*2-max(maxmin(4,:))])  

subplot(2,1,2)
title('Left (Phase 3,4)     (Position Corrected)');
counter=1;
for rast_i=1:trials/4
for rast_m=1:length(left_stag(1,:))
	line([left_stag(rast_i,rast_m) left_stag(rast_i,rast_m)],[counter-1 counter])
end
counter=counter+1;
end
ylim([0 trials/4])


hold off

%%%%%%%%%%%%%%%%%%%%%%%%%More Processing for PSTH plotting%%%%%%%%%%%%%%%%%
%must replace zeros with ones, so we don't have indexing problems
right_stag(right_stag==0)=1;
left_stag(left_stag==0)=1;
psth_right=zeros(trials/4, max(max(right_stag)));
for row_n=1:trials/phase
    for col_n=1:length(right_stag(1,:))
        psth_right(row_n, right_stag(row_n,col_n))=psth_right(row_n, right_stag(row_n,col_n))+1;
    end
end

%Make all L values positive
for row_n=1:trials/4
    for col_n=1:length(left_stag)
        if left_stag(row_n,col_n) == 1
            left_stag_pos(row_n,col_n)=1;
        else
            left_stag_pos(row_n,col_n)=left_stag(row_n,col_n)+max(maxmin(4,:));
        end 
    end
end
psth_left=zeros(trials/4,max(max(left_stag_pos)));
for row_n=1:trials/phase
    for col_n=1:length(left_stag_pos(1,:))
        psth_left(row_n, left_stag_pos(row_n,col_n))=psth_left(row_n, left_stag_pos(row_n,col_n))+1;
    end
end
%take out all the ones placed at the first sample, sum
right_sum=sum(psth_right);
right_sum(1)=0;
left_sum=sum(psth_left);
left_sum(1)=0;

%bin
bin_size=10; %ms
bin=round((bin_size/1000)*datainfo.indev.Fs); %sampling rate factor
index=1;
for bin_n=1:floor((length(right_sum)/bin)-1)
    right_binned(bin_n)=sum(right_sum(index:index+bin-1));
    index=index+bin;
end
index=1;
for bin_n=1:floor((length(left_sum)/bin)-1)
    left_binned(bin_n)=sum(left_sum(index:index+bin-1));
    index=index+bin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot psths%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psth_time_left=(1:length(left_binned))*bin_size;
psth_time_right=(1:length(right_binned))*bin_size;
y_max=max([left_binned right_binned]);
figure
hold on
subplot(211)
bar(psth_time_right, right_binned);
title('Right (Phase 1,2)');
ylim([0 y_max+2]);
%xlim([max(maxmin(4,:))*1000/Fs L_trace*2*1000/Fs-max(maxmin(4,:))*1000/Fs]); 

xlabel('Latency (ms)');
ylabel('# Spikes');

subplot(212)
bar(psth_time_left, left_binned);
title('Left (Phase 3,4)');
ylim([0 y_max+2])
%xlim([max(maxmin(4,:))*1000/Fs L_trace*2*1000/Fs-max(maxmin(4,:))*1000/Fs]); 
xlabel('Latency (ms)');
ylabel('# Spikes');
hold off

%Plot Position Data
figure;
title('Position Data');
plot(posData);
%%%%%%%%%%%%%%%%%%%%%%Kruskall-Wallis Test btwn Phases%%%%%%%%%%%%%%%%%%%%%
%kw=zeros(1,6);
%group=['s';'s']; %argument for krusalwallis so plots don't show
%Phase 1,3
%kw(1)=kruskalwallis([spikecountSorted(:,1) spikecountSorted(:,3)],group,'off');
%Phase 2,4
%kw(2)=kruskalwallis([spikecountSorted(:,2) spikecountSorted(:,4)],group,'off');
%Phase 1,2
%kw(3)=kruskalwallis([spikecountSorted(:,1) spikecountSorted(:,2)],group,'off');
%Phase 3,4
%kw(4)=kruskalwallis([spikecountSorted(:,3) spikecountSorted(:,4)],group,'off');
%Phase 1,4
%kw(5)=kruskalwallis([spikecountSorted(:,1) spikecountSorted(:,4)],group,'off');
%Phase 2,3
%kw(6)=kruskalwallis([spikecountSorted(:,2) spikecountSorted(:,3)],group,'off');
%disp('Phase 1,3 	2,4	   1,2	 3,4	  1,4	 2,3');
%kw


