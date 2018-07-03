clear all;
i=sqrt(-1);

signal_freq=2*pi*[0.1 0.3 0.4+0.005i 0.2 0.1+0.002i]; %% angular frequencies of signal
signal_amp=[1 -0.5 0.2 -0.3 0.01]; %% amplitdues in spectra
signal_length=100; %% length of signal
noise_level=1e-3; %% noise level

signal=signal_amp*exp(i*(signal_freq.')*[0:signal_length-1]); %% signal = sum of periodic functions

%% adding noise
signal=signal+noise_level*rand(1,signal_length);
plot(real(signal))

[freq,amp]=matrixpencil(signal,4); %% matrix pencil method 


%% we print the frequncies and amplitudies in signal
%% both given on beginning and found by matrix pencil 
%% method (sorted by abs(amplitiude))
disp('--freq and amp in signal---')
[absamp,indx]=sort(abs(signal_amp),'DESCEND');
disp(signal_freq(indx(1:length(signal_freq))));
disp(signal_amp(indx(1:length(signal_freq))));

disp('--freq and amp found by MP method---')
[absamp,indx]=sort(abs(amp),'DESCEND');
disp(freq(indx(1:length(freq)))');
disp(amp(indx(1:length(freq)))');

