% oscsoundmeasure.m
% 
% Script to measure sound output with B&K while 
% delivering vestibular stimuli
% 
% Procedure: 
% 
% 	1) collect baseline (background) noise
% 	
% 	2) collect noise and rotation during owlscillation
% 
% Uses RX8_2_Owlscillator_ContinuousAcq.rco
% 
% 
% See Also:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad Shanbhag
% sharad.shanbhag@einstein.yu.edu
%------------------------------------------------------------------------
% Created: 29 September, 2009 (SJS) from osc_ITD.m 
% 
% Input Chan - Analog IN 4(potentiometer)
% Output Chan - Analog Out 19 (motor command) 
% Output Chan - Analog Out 20 (masking noise)
%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data files
[filename, pathname] = uiputfile('*.mat', 'Save File As');
    dataFile = [pathname filename];
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some universal constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	L = 1;
	R = 2;
	HI = 1;
	LO = 0;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Owlscillation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% frequency
	OscFreq = 1;
	OscPeriod = 1/OscFreq;
	
	% # of owscillations to run (periods)
	Cycles = 60*6;    %ea cycle=1 second
    
    % Masking noise
    maskAmp= 0.085;
    mask=input('Play Masking Noise During Run? (y/n)      ', 's');
    if mask=='y'
      mask=1;
    else 
        mask=0;
    end  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TDT settings (Recording Settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	PosDecimateFactor = 100;
	BKDecimateFactor = 1;

	% set the feedback low pass filter
	FBFiltFc = 20;
	
	% set the input high pass filter
	BK_HPFc = 25;
	% set the input low pass filter
	BK_LPFc = 20000;
	
	% BK sensitivity (V/Pa)
	BKsense = 10;
	
	% BKscaling factor (Pa/V)
	BKscale = 1/BKsense;	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RX8 settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	RX8.Fs = 50000;
	% set this to wherever the circuits are stored
	RX8.Circuit_Path = 'H:\Code\TytoLogy\toolbox\TDT\Circuits\RX8_2\50KHz';
	RX8.Circuit_Name = 'RX8_2_Owlscillator_ContinuousAcq_mask';
	% Dnum = device number - this is for RX8 (2)
	RX8.Dnum=2;
	RX8.C = [];
	RX8.status = 0;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize some cells and arrays for storing data and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TDT Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Zbus and Devices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('...starting TDT hardware...');

	disp('...starting zBUS...')
	zBUS = zBUSinit('GB');

	disp('...starting RX8 for owlscillator control and BK signal acq...')
	% Initialize RX8
	tmpdev = RX8init('GB', 2);
	RX8.C = tmpdev.C;
	RX8.handle = tmpdev.handle;

	% Loads circuits
	RX8.status = RPload(RX8);

	% Starts Circuits
	RPrun(RX8)

	% Get circuit information
	RX8.Fs = RPsamplefreq(RX8);
	Fs = RX8.Fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, setup the circuits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% RX8 %%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% set number of Owlscillator cycles (periods)
	RPsettag(RX8, 'NOscCycles', Cycles);
	NOscCycles = RPgettag(RX8, 'NOscCycles')
	stim.NOscCycles = NOscCycles;

	% set owlscillator oscillation frequency
	RPsettag(RX8, 'OscFreq', OscFreq);
	stim.OscFreq = OscFreq;
	
	% read the period back from the circuit
	oscperiod = RPgettag(RX8, 'OscPeriod');
	disp(sprintf('Owlscillator Period is %.4f seconds', oscperiod));
	stim.OscPeriod = oscperiod;
	
	% set the decimation factor for the position (potentiometer voltage) data
	RPsettag(RX8, 'FBDeciFactor', PosDecimateFactor);
	% set the decimation factor for the input (BKvoltage) data
	RPsettag(RX8, 'INDeciFactor', BKDecimateFactor);

	% set the feedback low pass filter
	RPsettag(RX8, 'FBFiltFc', FBFiltFc);
	
	% set the input high pass filter
	RPsettag(RX8, 'INHPFc', BK_HPFc);
	% set the input low pass filter
	RPsettag(RX8, 'INLPFc', BK_LPFc);
	
	% enable coefficient computation
	RPsettag(RX8, 'FBFiltEnable', 1);
	RPsettag(RX8, 'INHPEnable', 1);
	RPsettag(RX8, 'INLPEnable', 1);
	pause(0.1);
	% disable coefficient computation ( to save cycles - set this to 1 again
	% if Fc value is changed!!!!)
	RPsettag(RX8, 'FBFiltEnable', 0);
	RPsettag(RX8, 'INHPEnable', 0);
	RPsettag(RX8, 'INLPEnable', 0);
	
	% get the buffer sizes for the input channels
	posBufSize = RPgettag(RX8, 'FBBufsize');
	bkBufSize = RPgettag(RX8, 'INBufsize');
	
	% set the input scaling factor
	RPsettag(RX8, 'INScale', BKscale);
	
	% compute the threshold at which buffer data will be read
	posBufThreshold = floor(ms2samples(OscPeriod*1000, RX8.Fs)*0.1/PosDecimateFactor)
	bkBufThreshold = floor(ms2samples(OscPeriod*1000, RX8.Fs)*0.1/BKDecimateFactor)
    
    % masking noise
    RPsettag(RX8,'Mask',mask);
    RPsettag(RX8,'MaskAmp',maskAmp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% set zBUS trigger A to high to synchronize and enable circuits
	zBUStrigA(zBUS, 0, HI, 10);

	% reset the buffers software triggers
	RPtrig(RX8, 2);
	RPtrig(RX8, 4);
	RPtrig(RX8, 5);
	
	ncycles = RPgettag(RX8, 'Ncycles');

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% setup Plot window
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	F = figure(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Give user a chance to abort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	cFlag = query_user('Continue', 1);

	if ~cFlag
		osc_tdtexit;
		return;
	else
		disp('press key to start');
		pause;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	posData = [];
	bkData = [];
	nPosData = 0;
	nBKData = 0;
	
	% send soft trig 1 to RX8 to start sine output
	RPtrig(RX8, 1);
	oscRunFlag = RPgettag(RX8, 'OscRun')

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Main experimental loop!
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% loop while oscillator output is ON
	posIndexStart = RPgettag(RX8, 'FBIndex');

	while RPgettag(RX8, 'OscRun')
		
		% check if there's data in the input buffer
		bkPoints = RPgettag(RX8, 'INIndex');
		if bkPoints > bkBufThreshold
			%reads from the buffer
			tmpresp = RPreadV(RX8, 'INData', bkPoints, 0);
			% reset the posn acq buffer
			RPtrig(RX8, 5);
			% ndata counts # of data events captured
			nBKData = nBKData + 1;
			disp(['reading BK data ' num2str(nBKData)])
			% store data
			bkData = [bkData tmpresp];
			subplot(211)
			plot(bkData);
			drawnow
		end
		
		% check if there's data in the position buffer
		posPoints = RPgettag(RX8, 'FBIndex');		
		if posPoints > posBufThreshold
			%reads from the buffer
			tmpresp = RPreadV(RX8, 'FBData', posPoints, 0);
			% reset the posn acq buffer
			RPtrig(RX8, 4);
			% ndata counts # of data events captured
			nPosData = nPosData + 1;
			disp(['reading position data ' num2str(nPosData)])
			% store data
			posData = [posData tmpresp];
			subplot(212)
			plot(posData);
			drawnow
		end
	end

	pause(0.5);

	% check if there's data in the input buffer
	bkPoints = RPgettag(RX8, 'INIndex');
	if bkPoints
		%reads from the buffer
		tmpresp = RPreadV(RX8, 'INData', bkPoints, 0);
		% reset the posn acq buffer
		RPtrig(RX8, 5);
		% ndata counts # of data events captured
		nBKData = nBKData + 1;
		disp(['reading BK data ' num2str(nBKData)])
		% store data
		bkData = [bkData tmpresp];
	end

	% check if there's data in the position buffer
	posPoints = RPgettag(RX8, 'FBIndex');		
	if posPoints
		%reads from the buffer
		tmpresp = RPreadV(RX8, 'FBData', posPoints, 0);
		% reset the posn acq buffer
		RPtrig(RX8, 4);
		% ndata counts # of data events captured
		nPosData = nPosData + 1;
		disp(['reading position data ' num2str(nPosData)])
		% store data
		posData = [posData tmpresp];
		subplot(212)
		plot(posData);
		drawnow
	end
	
	ncycles = RPgettag(RX8, 'Ncycles')

	disp('done');
	
	% get time stamp
	time_end = now;

	disp(sprintf('feedback max (volts): %f', RPgettag(RX8, 'FeedbackHV')));
	disp(sprintf('feedback min (volts): %f', RPgettag(RX8, 'FeedbackLV')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in full position feedback data and save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shutdown TDT things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% set zBUS trigger A to LO to disable circuits
    zBUStrigA(zBUS, 0, LO, 10);
    RPclose(RX8);
    zBUSclose(zBUS);

	% rms the bkData
	bk_rms = rms(bkData)
	% convert to dB SPL
	bk_dbspl = dbspl(bk_rms/BKsense)
	
	bk_filtered_dbspl = dbspl(rms(filterdata(bkData, RX8.Fs, 60, 15000))/BKsense)

	save(dataFile, 'posData', 'bkData', 'Fs', ...
						'PosDecimateFactor', 'BKDecimateFactor', ...
						'nBKData', 'nPosData', ...
						'BKsense', 'BKscale', ...
						'-MAT');
	
	
