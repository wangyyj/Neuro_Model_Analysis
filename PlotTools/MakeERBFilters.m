function fcoefs=MakeERBFilters(fs,numChannels,lowFreq)
% 0002 % function [fcoefs]=MakeERBFilters(fs,numChannels,lowFreq)
% 0003 % This function computes the filter coefficients for a bank of
% 0004 % Gammatone filters.  These filters were defined by Patterson and
% 0005 % Holdworth for simulating the cochlea.
% 0006 %
% 0007 % The result is returned as an array of filter coefficients.  Each row
% 0008 % of the filter arrays contains the coefficients for four second order
% 0009 % filters.  The transfer function for these four filters share the same
% 0010 % denominator (poles) but have different numerators (zeros).  All of these
% 0011 % coefficients are assembled into one vector that the ERBFilterBank
% 0012 % can take apart to implement the filter.
% 0013 %
% 0014 % The filter bank contains "numChannels" channels that extend from
% 0015 % half the sampling rate (fs) to "lowFreq".  Alternatively, if the numChannels
% 0016 % input argument is a vector, then the values of this vector are taken to
% 0017 % be the center frequency of each desired filter.  (The lowFreq argument is
% 0018 % ignored in this case.)
% 0019 
% 0020 % Note this implementation fixes a problem in the original code by
% 0021 % computing four separate second order filters.  This avoids a big
% 0022 % problem with round off errors in cases of very small cfs (100Hz) and
% 0023 % large sample rates (44kHz).  The problem is caused by roundoff error
% 0024 % when a number of poles are combined, all very close to the unit
% 0025 % circle.  Small errors in the eigth order coefficient, are multiplied
% 0026 % when the eigth root is taken to give the pole location.  These small
% 0027 % errors lead to poles outside the unit circle and instability.  Thanks
% 0028 % to Julius Smith for leading me to the proper explanation.
% 0029 
% 0030 % Execute the following code to evaluate the frequency
% 0031 % response of a 10 channel filterbank.
% 0032 %    fcoefs = MakeERBFilters(16000,10,100);
% 0033 %    y = ERBFilterBank([1 zeros(1,511)], fcoefs);
% 0034 %    resp = 20*log10(abs(fft(y')));
% 0035 %    freqScale = (0:511)/512*16000;
% 0036 %    semilogx(freqScale(1:255),resp(1:255,:));
% 0037 %    axis([100 16000 -60 0])
% 0038 %    xlabel('Frequency (Hz)'); ylabel('Filter Response (dB)');
% 0039 
% 0040 % Rewritten by Malcolm Slaney@Interval.  June 11, 1998.
% 0041 % (c) 1998 Interval Research Corporation
 
 T = 1/fs;
 if length(numChannels) == 1
     cf = ERBSpace(lowFreq, fs/2, numChannels);
 else
     cf = numChannels(1:end);
     if size(cf,2) > size(cf,1)
         cf = cf';
     end
 end
 
 % Change the followFreqing three parameters if you wish to use a different
 % ERB scale.  Must change in ERBSpace too.
 EarQ = 9.26449;                %  Glasberg and Moore Parameters
 minBW = 24.7;
 order = 1;
 
%  ERB = ((cf/EarQ).^order + minBW^order).^(1/order); %frequency channels that are dependent on the centre of the frequency channel. 
 ERB = abs(unique(diff(cf)))*ones(size(cf));   % This sets the bandwidth of all frequency channels is constant.
 B=1.019*2*pi*ERB;
 
 A0 = T;
 A2 = 0;
 B0 = 1;
 B1 = -2*cos(2*cf*pi*T)./exp(B*T);
 B2 = exp(-2*B*T);
 
 A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
         exp(B*T))/2;
 A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
         exp(B*T))/2;
 A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
         exp(B*T))/2;
 A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
         exp(B*T))/2;
 
 gain = abs((-2*exp(4*i*cf*pi*T)*T + ...
                  2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
                          (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
                           sin(2*cf*pi*T))) .* ...
            (-2*exp(4*i*cf*pi*T)*T + ...
              2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
               (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
                sin(2*cf*pi*T))).* ...
            (-2*exp(4*i*cf*pi*T)*T + ...
              2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
               (cos(2*cf*pi*T) - ...
                sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
            (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
            (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
           (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +  ...
            2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);
     
 allfilts = ones(length(cf),1);
 fcoefs = [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain];
 
 if (0)                        % Test Code
     A0  = fcoefs(:,1);
     A11 = fcoefs(:,2);
     A12 = fcoefs(:,3);
     A13 = fcoefs(:,4);
     A14 = fcoefs(:,5);
     A2  = fcoefs(:,6);
     B0  = fcoefs(:,7);
     B1  = fcoefs(:,8);
     B2  = fcoefs(:,9);
     gain= fcoefs(:,10);    
     chan=1;
     x = [1 zeros(1, 511)];
     y1=filter([A0(chan)/gain(chan) A11(chan)/gain(chan) ...
         A2(chan)/gain(chan)],[B0(chan) B1(chan) B2(chan)], x);
     y2=filter([A0(chan) A12(chan) A2(chan)], ...
             [B0(chan) B1(chan) B2(chan)], y1);
     y3=filter([A0(chan) A13(chan) A2(chan)], ...
             [B0(chan) B1(chan) B2(chan)], y2);
     y4=filter([A0(chan) A14(chan) A2(chan)], ...
             [B0(chan) B1(chan) B2(chan)], y3);
     semilogx((0:(length(x)-1))*(fs/length(x)),20*log10(abs(fft(y4))));
 end