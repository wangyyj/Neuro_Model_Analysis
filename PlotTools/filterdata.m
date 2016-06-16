function out = filterdata(in, inFs, hipassFc, lopassFc)
%out = filterdata(in, hipassFc, lopassFc, inFs)
% 
% 
% 
% See Also:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad Shanbhag
% sharad.shanbhag@einstein.yu.edu
%------------------------------------------------------------------------
% Created: 30 September, 2009 (SJS)
% 
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%-------------------------------------------------------------------------

hipassOrder = 3;			% Filter Order
lopassOrder = 3;			% lopass filter order

if ~exist('hipassFc', 'var')
	hipassFc = 60;				% hipass corner frequency
end

if ~exist('lopassFc', 'var')
	lopassFc = 15000;			% hipass corner frequency
end

Fnyq = inFs/2;		% Nyquist frequency

% get the coefficients for a butterworth filter, highpass
[hiB, hiA] =  butter(hipassOrder, hipassFc/Fnyq, 'high');

% get the coefficients for a butterworth filter, lopass
[loB, loA] =  butter(lopassOrder, lopassFc/Fnyq, 'low');

in = sin2array(in, 1, inFs);
% hi pass filter the data to remove DC and low frequency crap
out = filtfilt(hiB, hiA, in);
% low pass filter the data to remove RF stuff
out = filtfilt(loB, loA, out);