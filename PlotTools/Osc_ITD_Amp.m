%% For Determining changes in amplitude for OSC ITD Curve files
%% 12/1/09 JW
clear;
[FileName,PathName] = uigetfile('*.dat','Select the .dat-file','Y:\Jennifer\Dichotic_Osc\13Oct09');
full_filename=[PathName FileName];
[data, datainfo] = readOwscData(full_filename);
figure;
plot(data{1}.datatrace, 'r');
hold on;
plot(data{2}.datatrace, 'g');
plot(data{3}.datatrace, 'b');
title('Raw Data- Trace 1,2,3. Please Determine Threshold');
xlabel('Time');
ylabel('Resp (mV)');
grid;
hold off
amps(1,1)= max(data{1}.datatrace(1600:end,1));
amps(1,2)= max(data{2}.datatrace(1600:end,1));
amps(1,3)= max(data{3}.datatrace(1600:end,1));
amps(1,4)= max(data{4}.datatrace(1600:end,1));
amps(1,5)= max(data{5}.datatrace(1600:end,1));
amps(1,6)= max(data{6}.datatrace(1600:end,1));
amps(1,7)= max(data{7}.datatrace(1600:end,1));
amps(1,8)= max(data{8}.datatrace(1600:end,1));
amps(1,9)= max(data{9}.datatrace(1600:end,1));
amps(1,10)= max(data{10}.datatrace(1600:end,1));
max_amp= mean(amps)
var_amp= var(amps)