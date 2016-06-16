clear;



%Query for file to be analyzed, load data
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file','Y:\Louisa');
full_filename=[PathName FileName];
[data, datainfo] = readHPData(full_filename);
rep=datainfo.curve.nreps;
trials=datainfo.curve.nTrials;  

Fs=datainfo.indev.Fs;



%Draw Sample traces
figure
subplot(211)
plot([(1:length(data{1}.datatrace))/datainfo.indev.Fs],data{1}.datatrace, 'b');
subplot(212)
plot([(1:length(data{1}.datatrace))/datainfo.indev.Fs],data{1}.datatrace, 'b');
hold on
plot([(1:length(data{2}.datatrace))/datainfo.indev.Fs],data{2}.datatrace, 'g');
plot([(1:length(data{3}.datatrace))/datainfo.indev.Fs],data{3}.datatrace, 'r');
plot([(1:length(data{4}.datatrace))/datainfo.indev.Fs],data{4}.datatrace, 'k');
title('Raw Data- Trace 1,2,3,4. Please Determine Threshold');
xlabel('Time (s)');
ylabel('Resp (mV)');
grid;
hold off

%Query for Threshold
thres= input('Enter Threshold Value   ');


%Calculate spikes
spiketimes=zeros(length(data),100);
spikecounts=zeros(length(data),1);
datastart=round((datainfo.stim.Delay*0.001)*Fs);
dataend=round(((datainfo.tdt.StimDuration+datainfo.stim.Delay+5)*0.001)*Fs);

n=1;
for n=1:length(data)
    spikes = spikeschmitt2(data{n}.datatrace(datastart:dataend), thres, 1,datainfo.indev.Fs); %set analysis window here
    spikecounts(n,1)=length(spikes);
    spiketimes(n,1:length(spikes))=spikes;
    n=n+1;    
end



%extract epochs

alleps=zeros(length(data),trials);

for i=1:length(data);
	
	alleps(i,data{i}.trialNumber)=1;
	
end

epochs1=find(alleps(:,1)==1)';
epochs2=find(alleps(:,2)==1)';


%sort spiketimes according to epochs

spikes1=spiketimes(epochs1,:)./Fs;
spikes2=spiketimes(epochs2,:)./Fs;


a=size(spikes1);
newspikes1=cell(a(1),1);

for i=1:a(1);
	
	temp=find(spikes1(i,:)~=0);
	newspikes1{i}=spikes1(i, 1:length(temp));
	
	clear temp
	
end


b=size(spikes2);
newspikes2=cell(b(1),1);

for i=1:b(1);
	
	temp=find(spikes2(i,:)~=0);
	newspikes2{i}=spikes2(i, 1:length(temp));
	
	clear temp
	
end



SAMFREQ1=int2str(data{epochs1(1)}.dataID);
SAMFREQ2=int2str(data{epochs2(1)}.dataID);


%Query for filenum
filenum= input('Enter filenumber  ');


NUM=int2str(filenum);


filename1=['Spikes' NUM 'SAMFRQ' SAMFREQ1 '.mat'];
save(filename1, 'newspikes1');

filename2=['Spikes' NUM 'SAMFRQ' SAMFREQ2 '.mat'];
save(filename2, 'newspikes2');


