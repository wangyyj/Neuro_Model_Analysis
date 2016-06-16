% script to analyze owlscillator sounds from B&K calibration mic
%
% Created 30 September, 2009

Flo = 500;
Fhi = 15000;

% load data file
[FileName,PathName] = uigetfile('*.mat','Select the .mat-file for Analysis',pwd);
full_filename=[PathName FileName];
load (full_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create time vector for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample interval
dt = 1/Fs;

% # points in bkData (microphone data vector)
nBK = length(bkData);

tVec = (0:(nBK-1)) * dt;

% convert bkData to Pascals instead of volts

bkDataPa = bkData  / BKsense;

% plot
figure
subplot(311)
plot(tVec, bkDataPa)

% position sample rate
posFs = Fs / PosDecimateFactor;
posdt = 1/posFs;
nPos = length(posData);
posT = (0:(nPos-1))*posdt;

subplot(312)
plot(posT, posData)


% plot spectrogram of data
% figure(2)
% spectrogram(bkData);

% period of owlscillation
Fosc = 0.05;
Period = 1/Fosc;

% 1 second binsize
windowtime = 0.5;
binsize = floor(windowtime * Fs);
nbins = floor(length(bkData) / binsize);

rmsdata = zeros(nbins, 1);

bkDataPaFilt = filterdata(bkData, Fs, Flo, Fhi);

for n = 1:nbins
	startbin = 1 + (n-1) * binsize;
	endbin = n*binsize;
	
	rmsdata(n) = rms(bkDataPaFilt(startbin:endbin));
end

subplot(313)
plot(dbspl(rmsdata))

[FileName_bg,PathName_bg] = uigetfile('*.mat','Select Background .mat-file',pwd);
full_filename_bg=[PathName_bg FileName_bg];
load (full_filename_bg);
nbins = floor(length(bkData) / binsize);
bg_rmsdata = zeros(nbins, 1);

bg_Filt = filterdata(bkData, Fs, Flo, Fhi);

for n = 1:nbins
	startbin = 1 + (n-1) * binsize;
	endbin = n*binsize;
	
	bg_rmsdata(n) = rms(bg_Filt(startbin:endbin));
end

subplot(313)
hold on
plot(dbspl(bg_rmsdata), 'r:')
hold off
	



