% For loading data, filtering (if necessary) then drawing ITD curve

%%%%%%%%%%%%%%%%%%%%%Query for file to be analyzed, load data%%%%%%%%%%%%%%
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file','B:\Jennifer\');
full_filename=[PathName FileName];
[data, datainfo] = readOwscData(full_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trials=length(data);
rep=datainfo.curve.nreps;
trials=datainfo.curve.nTrials;
Fs=datainfo.indev.Fs;
phase=4;    %total number of phases

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Draw Sample Traces%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot 311
plot(data_raw(1,:), 'r');
hold on;
subplot 312
plot(data_raw(2,:), 'g');
subplot 313
plot(data_raw(3,:), 'b');
title('Raw Data- Trace 1,2,3');
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
    subplot 311
    plot(data_filt(1,:), 'r');
    hold on;
    subplot 312
    plot(data_filt(2,:), 'g');
    subplot 313
    plot(data_filt(3,:), 'b');
    title('Filtered Data');
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
    plot(data_ITD(1,:), 'r');
    hold on;
    plot(data_ITD(2,:), 'g');
    plot(data_ITD(3,:), 'b');
    title('Determine Threshold');
    xlabel('Time');
    ylabel('Resp (mV)');
    grid;
    hold off
thres= input('Enter Threshold Value   ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get Spikes%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate spikes
spikeTimes=zeros(trials,50); %%spikes calculated, unsorted
spikeCounts=zeros(trials,1); %%unsorted spike counts
n=1;
for n=1:trials
    spikes= spikeschmitt2(data_ITD(n,:), thres, 1,Fs);  
    spikeTimes(n,1:length(spikes))=spikes;
    spikeCounts(n)=length(spikes); 
    n=n+1;
end

%Sorting
%Sort into Phases and Reps and ITD
phase=4;    %total number of phases
spikecountSorted = zeros(phase*rep+rep,length(datainfo.curve.ITDrange)); %actual data starts from row 6 
order=datainfo.curve.StimSequence;
for m=1:trials      
    spikecountSorted((order(m,2)*rep+order(m,4)), order(m,1))=spikeCounts(m); %sorted matrix- rows= phase*4+rep, order(m,1)=ITD order
    m=m+1;
end

%Trim data set so phase1 rep1 starts at row 1
spikecountSorted=spikecountSorted(rep+1:end,:);

%plot ITD curve - Single panel
figure;
ITD = datainfo.curve.ITDrange;
%Phase 1 Blue
errorbar(ITD,mean(spikecountSorted(1:rep,:)),std(spikecountSorted(1:rep,:))/sqrt(rep), 'b'); % standard error of the mean
hold on
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

%Plot ITD Curve - Subplot
figure;
%Phase 1 Blue
subplot(4,1,1) 
hold on
errorbar(ITD,mean(spikecountSorted(1:rep,:)),std(spikecountSorted(1:rep,:))/sqrt(rep), 'b'); % standard error of the mean
title('Phase 1'); 
ylabel('# Spikes');
grid;
subplot(4,1,2)
errorbar(ITD,mean(spikecountSorted(rep+1:rep*2,:)),std(spikecountSorted(rep+1:rep*2,:))/sqrt(rep),'g'); 
title('Phase 2'); 
ylabel('# Spikes');
grid;
subplot(4,1,3)
errorbar(ITD,mean(spikecountSorted(rep*2+1:rep*3,:)),std(spikecountSorted(rep*2+1:rep*3,:))/sqrt(rep), 'r'); 
title('Phase 3'); 
ylabel('# Spikes');
grid;
subplot(4,1,4)
errorbar(ITD,mean(spikecountSorted(rep*3+1:rep*4,:)),std(spikecountSorted(rep*3+1:rep*4,:))/sqrt(rep), 'k'); 
title('Phase 4'); 
hold off
xlabel('ITD (us)');
ylabel('# Spikes');
grid;



