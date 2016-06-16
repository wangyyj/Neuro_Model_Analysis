%------------------------------------------------------------------------
% PowerSpectrum.m
%------------------------------------------------------------------------
% FFT transform
% FFT power spectrum analysis 
%------------------------------------------------------------------------
[data]=uigetfile('*dat', 'Select data for Power Spectrum Analysis','cd');
Fs=input('Input data Fs     ');
L=length(data);
t=(1:L)/Fs;
figure
plot(t,data);
hold on
title('data')
xlabel('Time (s)')
hold off

Y = fft(data,L);
Y_conj = Y.*conj(Y)/L;
f = Fs/L*(0:L/2+2);
figure;
plot(f,Y_conj(1:L/2+3));
%use only half of the points, the rest are symmetric
hold on
title('Power Spectral Density')
xlabel('Frequency (Hz)')
