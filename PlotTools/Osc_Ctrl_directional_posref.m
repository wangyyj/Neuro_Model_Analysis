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
w_Fs=round(Fs);
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

% Get Reference Points - Max, Min, Ref
posdata_step=floor(length(posData)/(trials/4)); 
maxmin=[];
start=1;
for n=1:(trials/4)
    [maxmin(1,n) maxmin(2,n)]=max(posData(start:start+posdata_step)); % max
    maxmin(2,n)=maxmin(2,n)+posdata_step*(n-1); %MAX index in pos Fs
    [maxmin(4,n) maxmin(5,n)]=min(posData(start:start+posdata_step)); % min
    maxmin(5,n)=maxmin(5,n)+posdata_step*(n-1); % MIN  index in pos Fs 
    start=start+posdata_step;
end
maxmin(3,:)=round(maxmin(2,:)/500*Fs); %Max index in neuron Fs
maxmin(6,:)=round(maxmin(5,:)/500*Fs); %Min index in neuron Fs
% middle of phase
index=1;
for n=1:length(maxmin)-1
    mid(1,index)= (maxmin(1,n)+maxmin(4,n))/2; % middle V
    mid(2,index)= (maxmin(5,n)+maxmin(2,n))/2; % middle time (pos)
    mid(3,index)= (maxmin(6,n)+maxmin(3,n))/2; % middle time (neuron)
    mid(1,index+1)= (maxmin(1,n)+maxmin(4,n+1))/2; % middle V
    mid(2,index+1)= (maxmin(5,n)+maxmin(2,n+1))/2; % middle time (pos)
    mid(3,index+1)= (maxmin(6,n)+maxmin(3,n+1))/2; % middle time (neuron)
    index=index+2;
end

%Get Reference Point
pos_thres=(maxmin(1,18)+maxmin(4,18))/2;
rounded=floor(posData*10);
intercepts=find(rounded==round(pos_thres*10));
ref(1,1)=intercepts(1,1);
for i=2:length(intercepts)
   if intercepts(i)>= intercepts(i-1)+500
	ref(i)=intercepts(i);
    end
end
ref=ref(ref~=0);
ref(2,:)=pos_thres;
%Plot Reference points for sanity check
figure;
plot(posData(1,:));
hold on
plot(ref(1,:),ref(2,:),'k*'); 
plot(maxmin(2,:),maxmin(1,:),'go'); %max
plot(maxmin(5,:),maxmin(4,:),'rx'); %min
plot(mid(2,:),mid(1,:),'cx'); %middle
hold off
   