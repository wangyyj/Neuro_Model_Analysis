clear all;
%%%%%%%%%%%%%%%Query for file to be analyzed, load data%%%%%%%%%%%%%%%%%%%%
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file','Y:\Jennifer\Dichotic_Osc\13Oct09');
full_filename=[PathName FileName];
[data, datainfo] = readOwscData(full_filename);

rep=datainfo.curve.nreps;
trials =40;
Fs= datainfo.indev.Fs;

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
for n=1:trials
    spikes= spikeschmitt2(data{n}.datatrace(:), thres, 1,datainfo.indev.Fs);  
    spikeTimes(n,1:length(spikes))=spikes;
    spikeCounts(n)=length(spikes); 
    n=n+1;
end

%left side
index=2;
for n=1:trials/2
	left(n,1)=spikeCounts(index);
	index=index+2;
end
%right side
index=1;
for n=1:trials/2
	right(n,1)=spikeCounts(index);
	index=index+2;
end

%summary  col1 = left; col 2=right; row 1=mean; row 2=std
output(1,1)= mean(left(:));
output(2,1)= std(left(:));
output(1,2)= mean(right(:));
output(2,2)=std(right(:));
output
