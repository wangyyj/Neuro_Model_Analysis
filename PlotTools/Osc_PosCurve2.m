% For analyzing control data - to determine if neuron is motion sensitive
% Modified from Osc_ITd_Curve
% 1/5/10   YY Wang

%clear all;
%%%%%%%%%%%%%%%Query for file to be analyzed, load data%%%%%%%%%%%%%%%%%%%%
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file','cd');
full_filename=[PathName FileName];
disp(full_filename);
[posFileName,posPathName] = uigetfile('*.mat','Select the matching Position File(_posData.mat)','cd');
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
end
% Rows are trials, columns are datapoints
disp('...Data Loaded');

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
    spikes= spikeschmitt2(data{n}.datatrace(:), thres, 1,Fs);  
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
spikecountSorted(:,1:3)= spikecountSorted(:,2:4);
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
phase_len= 1.8*w_Fs; %length of each phase=2s
Phase2_spiketimes=Phase2_spiketimes+phase_len;
Phase2_spiketimes(Phase2_spiketimes==phase_len)=0;
Phase4_spiketimes=Phase4_spiketimes+phase_len;
Phase4_spiketimes(Phase4_spiketimes==phase_len)=0;
%reduce to 35 trials since 1st phase is phase4
Phase4_spiketimes=Phase4_spiketimes(2:end,:);
Phase1_spiketimes=Phase1_spiketimes(1:end-1,:);
Phase2_spiketimes=Phase2_spiketimes(1:end-1,:);
Phase3_spiketimes=Phase3_spiketimes(1:end-1,:);
right=[Phase1_spiketimes Phase2_spiketimes];
left= [Phase3_spiketimes Phase4_spiketimes];

%Position Alignment
%Get reference points from PosData - Max, Min
posdata_step=floor(length(posData)/(trials/4)); 
maxmin=[];
start=1;
for n=1:(trials/4)
    [maxmin(1,n) maxmin(2,n)]=max(posData(start:start+posdata_step)); % max
    maxmin(2,n)=maxmin(2,n)+posdata_step*(n-1); %MAX index in pos Fs
    [maxmin(6,n) maxmin(7,n)]=min(posData(start:start+posdata_step)); % min
    maxmin(7,n)=maxmin(7,n)+posdata_step*(n-1); % MIN  index in pos Fs 
    start=start+posdata_step;
end
maxmin(3,:)=round(maxmin(2,:)/500*Fs); %Max index in neuron Fs
maxmin(8,:)=round(maxmin(7,:)/500*Fs); %Min index in neuron Fs
%Plot PosData with Reference Points
figure;
plot(posData(1,:));
hold on
plot(maxmin(2,:),maxmin(1,:),'go'); %max
plot(maxmin(7,:),maxmin(6,:),'rx'); %min
hold off
% Put spiketimes in position sampling rate
right_Fs500= right*500/Fs;
left_Fs500= left*500/Fs;
% Find position of spike, using reference points
% Right starts at max, Left starts at min
for row=1:length(right(:,1))
    for col=1:length(right)
        if right_Fs500(row,col)==0
            right_v(row,col)=0;
        else
            right_v(row,col)= interp1(posData,maxmin(2,row)+right_Fs500(row,col));
        end
        if left_Fs500(row,col)==0
            left_v(row,col)=0;
        else
            left_v(row,col)= interp1(posData,maxmin(7,row)+left_Fs500(row,col));
        end
    end
end
%Transform into degrees
left_v(left_v==0)=NaN;
right_v(right_v==0)=NaN;
for row_i=1:length(left_v(:,1))
    for col_i=1:length(left_v(1,:))
        left_deg(row_i,col_i)=left_v(row_i,col_i)*-201.6 + 152.36;
        right_deg(row_i,col_i)=right_v(row_i,col_i)*-201.6 + 152.36;
    end 
end
left_deg=left_deg*-1; %Deg relative to RF of NEURON, not position in space
right_deg=right_deg*-1;
%plot position
%upper limit
xlim_lower=-90; 
xlim_upper=90;
%plot entire raster
figure;
subplot(2,1,1)
title('Movement to Right (Phase 1,2)');
counter=1;
for rast_i=1:length(right_deg(:,1))
    for rast_m=1:length(right_deg(1,:))
        line([right_deg(rast_i,rast_m) right_deg(rast_i,rast_m)],[counter-1 counter])
    end
    counter=counter+1;
end
xlim([xlim_lower xlim_upper])
ylim([0 length(right_deg(:,1))])
set(gca,'XTick',(xlim_lower:5:xlim_upper))
subplot(2,1,2)
title('Movement to Left (Phase 3,4))');
counter=1;
for rast_i=1:1:length(left_deg(:,1))
    for rast_m=1:length(left_deg(1,:))
        line([left_deg(rast_i,rast_m) left_deg(rast_i,rast_m)],[counter-1 counter])
    end
    counter=counter+1;
end
xlim([xlim_lower xlim_upper])
ylim([0 length(left_deg(:,1))])
set(gca,'XTick',(xlim_lower:5:xlim_upper))

%Binning for PSTH
xlim_lower= input('Enter Lower Binning Limit     ');
xlim_upper= input('Enter Upper Binning Limit     ');
bin_size= 2;  %deg
bins=floor(xlim_upper-xlim_lower)/bin_size;
bin_index=xlim_lower:bin_size:xlim_upper;
for n=1:bins-1
    left_bins(n)=length(find((left_deg>=bin_index(n))&(left_deg<=bin_index(n+1))));
    right_bins(n)=length(find((right_deg>=bin_index(n))&(right_deg<=bin_index(n+1))));
end
%plot raster and PSTH
%Right raster
figure;
subplot(4,1,1)
title('Leftward Motion'); %owl moving to right
counter=1;
for rast_i=1:length(right_deg(:,1))
    for rast_m=1:length(right_deg(1,:))
        line([right_deg(rast_i,rast_m) right_deg(rast_i,rast_m)],[counter-1 counter])
    end
    counter=counter+1;
end
xlim([xlim_lower xlim_upper])
ylim([0 length(right_deg(:,1))])
%Right PSTH
subplot(4,1,2)
xlim_lower:bin_size:xlim_upper;
bar(right_bins);
ylim([0 max([left_bins right_bins])+5]);
xlim([0.5 length(right_bins)+0.5]) %add 0.5 on both sides to get whole bar widths
set(gca,'XTick',(0))
%Left Raster
subplot(4,1,3)
title('Rightward Motion'); %owl Moving to left
counter=1;
for rast_i=1:length(left_deg(:,1))
    for rast_m=1:length(left_deg(1,:))
        line([left_deg(rast_i,rast_m) left_deg(rast_i,rast_m)],[counter-1 counter])
    end
    counter=counter+1;
end
xlim([xlim_lower xlim_upper])
ylim([0 length(right_deg(:,1))])
%Left PSTH
subplot(4,1,4)
xlim_lower:bin_size:xlim_upper;
bar(left_bins);
ylim([0 max([left_bins right_bins])+5]);
xlim([0.5 length(left_bins)+0.5]) %add 0.5 on both sides to get whole bar widths
set(gca,'XTick',(0))
hold off;

%Calculate Tot Spikes for each Direction
left_spikes=sum(left_bins);
right_spikes=sum(right_bins);
figure;
title('Total Spikes In Each Direction')
hold on
bar(1,left_spikes,'b')
bar(2,right_spikes,'r')
legend('Left','Right')
xlim([0 3])
set(gca,'XTick',(0))
hold off;
%Signed Directional Index (Wagner 1994)CCW=L->R
%if left_spikes>right_spikes
%    DI= 100*(1-left_spikes/right_spikes);
%else
%    DI= 100*(right_spikes/left_spikes-1);
%end
%disp('  Left  Right')
%disp([left_spikes right_spikes])
%disp('Signed Directional Index = ')
%disp(DI)
