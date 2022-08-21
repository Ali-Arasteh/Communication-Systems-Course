% Q1)
%%
% a)
%%
% b)
%%
% c)
%%
% 1)
[Audio,Audio_fs] = audioread('Audio.wav');
Audio = Audio.';
time = 0:1/Audio_fs:(length(Audio)-1)/Audio_fs;
figure();
plot(time,Audio);
xlim([0 (length(Audio)-1)/Audio_fs]);
ylim([-1 1]);
xlabel('time');
ylabel('Audio(t)');
title('Audio signal in time domain');
AudioF = 1/Audio_fs*fftshift(fft(Audio));
frequency = -Audio_fs/2:Audio_fs/length(Audio):Audio_fs/2-Audio_fs/length(Audio);
figure();
plot(frequency,abs(AudioF));
xlabel('frequency');
ylabel('Audio(f)');
title('Audio signal in frequency domain');
AudioESD = abs(AudioF).^2;
AudioTotalEnergy = sum(AudioESD);
PartOfAudioEnergy = AudioESD((length(Audio)+1)/2);
for i = 1:length(Audio)/2-1
   if(PartOfAudioEnergy > 0.99*AudioTotalEnergy)
       break;
   end
   PartOfAudioEnergy = PartOfAudioEnergy + AudioESD((length(Audio)+1)/2+i) + AudioESD((length(Audio)+1)/2-i);
end
BWFrequency = i*Audio_fs/length(AudioESD);
%%
% 2)
[Audio,Audio_fs] = audioread('Audio.wav');
Audio = Audio.';
LPF = load('LPF6500Hz44100.mat');
FilteredAudio = filter(LPF.LPF6500Hz,1,Audio);
FilteredAudioF = 1/Audio_fs*fftshift(fft(FilteredAudio));
frequency = -Audio_fs/2:Audio_fs/length(Audio):Audio_fs/2-Audio_fs/length(Audio);
figure();
plot(frequency,abs(FilteredAudioF));
xlabel('frequency');
ylabel('filtered Audio(f)');
title('filtered Audio signal in frequency domain');
sound(FilteredAudio,Audio_fs);
%%
% 3)
L = 10;
[Audio] = audioread('Audio.wav');
Audio = Audio.';
UpsampleAudio = interp(Audio,L);
time = 0:1/(L*Audio_fs):(length(UpsampleAudio)-1)/(L*Audio_fs);
figure();
plot(time,UpsampleAudio);
xlim([0 (length(UpsampleAudio)-1)/(L*Audio_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('upsample Audio(t)');
title('upsample Audio signal in time domain');
%%
% 4,5)
L = 10;
[Audio,Audio_fs] = audioread('Audio.wav');
Audio = Audio.';
AudioF = 1/Audio_fs*fftshift(fft(Audio));
UpsampleAudio = interp(Audio,L);
time = 0:1/(L*Audio_fs):(length(UpsampleAudio)-1)/(L*Audio_fs);
fc = 5e4;
beta = 3;
fm = 6.5e3;
deltaf = beta*fm;
AudioMaximum = max(abs(Audio));
fdelta = deltaf/AudioMaximum;
AudioIntegral = 1/(L*Audio_fs)*cumtrapz(UpsampleAudio);
Xc = cos(2*pi*fc.*time + 2*pi*fdelta.*AudioIntegral);
figure();
plot(time,Xc);
xlim([0 (length(Audio)-1)/Audio_fs]);
ylim([-1 1]);
xlabel('time');
ylabel('Xc(t)');
title('modulated signal in time domain');
frequency = -L*Audio_fs/2:L*Audio_fs/length(Xc):L*Audio_fs/2-L*Audio_fs/length(Xc);
XcF = 1/Audio_fs*fftshift(fft(Xc));
figure();
plot(frequency,abs(XcF));
xlabel('frequency');
ylabel('Xc(f)');
title('modulated signal in frequency domain');
DiffXc = zeros(1,length(Xc));
for i = 2:length(Xc)-2
   DiffXc(i) = L*Audio_fs/2*(Xc(i+1)-Xc(i-1));
end
AbsDiffXc = abs(DiffXc);
LPF = load('LPF6500Hz441000.mat');
Yd = filter(LPF.LPF6500Hz,1,AbsDiffXc);
DemodulatedSignal = Yd - mean(Yd);
DownsampleDemodulatedSignal = downsample(DemodulatedSignal,L)/(4*fdelta);
DownsampleDemodulatedSignalF = 1/Audio_fs*fftshift(fft(DownsampleDemodulatedSignal));
time = 0:1/Audio_fs:(length(DownsampleDemodulatedSignal)-1)/Audio_fs;
frequency = -Audio_fs/2:Audio_fs/length(Audio):Audio_fs/2-Audio_fs/length(Audio);
figure();
plot(time,DownsampleDemodulatedSignal);
xlim([0 (length(Audio)-1)/Audio_fs]);
ylim([-1 1]);
xlabel('time');
ylabel('downsample demodulated signal(t)');
title('downsample demodulated signal in time domain');
figure();
plot(frequency,abs(DownsampleDemodulatedSignalF));
xlabel('frequency');
ylabel('downsample demodulated signal(f)');
title('downsample demodulated signal in frequency domain');
sound(DownsampleDemodulatedSignal,Audio_fs);
MSE = immse(AudioF,DownsampleDemodulatedSignalF);
%%
% Q2)
%%
% a)
fs = 50000;
time = 0:1/fs:0.1-1/fs;
fm = 10;
X = cos(2*pi*fm*time);
fc = 1000;
fdelta = 100;
XIntegral = 1/fs*cumtrapz(X);
V = cos(2*pi*fc.*time + 2*pi*fdelta.*XIntegral);
figure();
subplot(2,1,1);
plot(time,X);
xlim([0 (length(X)-1)/fs]);
ylim([-1 1]);
xlabel('time');
ylabel('x(t)');
title('signal in time domain');
subplot(2,1,2);
plot(time,V);
xlim([0 (length(X)-1)/fs]);
ylim([-1 1]);
xlabel('time');
ylabel('v(t)');
title('modulated signal in time domain');
frequency = -fs/2:fs/length(X):fs/2-fs/length(X);
XF = 1/fs*fftshift(fft(X));
VF = 1/fs*fftshift(fft(V));
figure();
subplot(2,1,1);
plot(frequency,abs(XF));
xlabel('frequency');
ylabel('X(f)');
title('signal in frequency domain');
subplot(2,1,2);
plot(frequency,abs(VF));
xlabel('frequency');
ylabel('V(f)');
title('modulated signal in frequency domain');
VESD = abs(VF).^2;
AudioTotalEnergy = sum(VESD);
PartOfAudioEnergy = VESD(length(VF)/2+length(VF)/50) + VESD(length(VF)/2-length(VF)/50);
for i = 1:length(VF)/50
   PartOfAudioEnergy = PartOfAudioEnergy + VESD(length(VF)/2+length(VF)/50+i) + VESD(length(VF)/2+length(VF)/50-i) + VESD(length(VF)/2+length(VF)/50+i) + VESD(length(VF)/2+length(VF)/50-i);
   if(PartOfAudioEnergy > 0.99*AudioTotalEnergy)
       break;
   end
end
BWFrequency = 2*i*fs/length(VESD);
%%
% b)
fs = 50000;
time = 0:1/fs:0.1-1/fs;
fm = 10;
X = cos(2*pi*fm*time);
fc = 1000;
fdelta = 100;
XIntegral = 1/fs*cumtrapz(X);
V = cos(2*pi*fc.*time + 2*pi*fdelta.*XIntegral);
Y = V.^3;
frequency = -fs/2:fs/length(X):fs/2-fs/length(X);
YF = 1/fs*fftshift(fft(Y));
figure();
plot(frequency,abs(YF));
xlabel('frequency');
ylabel('Y(f)');
title('signal in frequency domain');
BPF = load('BPF1000Hz50000.mat');
Vtilda = filter(BPF.BPF1000Hz50000,1,filter(BPF.BPF1000Hz50000,1,filter(BPF.BPF1000Hz50000,1,Y)));
VtildaF = 1/fs*fftshift(fft(Vtilda));
figure();
plot(frequency,abs(VtildaF));
xlabel('frequency');
ylabel('Vtilda(f)');
title('signal in frequency domain');
%%
% c)
fs = 50000;
time = 0:1/fs:0.1-1/fs;
fm = 10;
X = cos(2*pi*fm*time);
fc = 1000;
fdelta = 100;
XIntegral = 1/fs*cumtrapz(X);
V = cos(2*pi*fc.*time + 2*pi*fdelta.*XIntegral);
DiffV = zeros(1,length(V));
for i = 2:length(V)-2
   DiffV(i) = fs/2*(V(i+1)-V(i-1));
end
EnvDiffLimV = envelope(DiffV);
DemodulatedSignal = EnvDiffLimV - mean(EnvDiffLimV);
DemodulatedSignal = max(DemodulatedSignal).*DemodulatedSignal;
figure();
plot(time,DemodulatedSignal);
xlabel('time');
ylabel('demodulated signal(t)');
title('demodulated signal in time domain');
frequency = -fs/2:fs/length(X):fs/2-fs/length(X);
DemodulatedSignalF = 1/fs*fftshift(fft(DemodulatedSignal));
figure();
plot(frequency,abs(DemodulatedSignalF));
xlabel('frequency');
ylabel('demodulated signal(f)');
title('demodulated signal in frequency domain');
%%
% d)
%%
fs = 50000;
time = 0:1/fs:0.1-1/fs;
fm = 10;
X = cos(2*pi*fm*time);
fc = 1000;
fdelta = 100;
XIntegral = 1/fs*cumtrapz(X);
V = cos(2*pi*fc.*time + 2*pi*fdelta.*XIntegral);
Ai = 0.1;
Wi = 100;
V = V + Ai*cos((2*pi*fc+Wi).*time);
DiffV = zeros(1,length(V));
for i = 2:length(V)-2
   DiffV(i) = fs/2*(V(i+1)-V(i-1));
end
EnvDiffLimV = envelope(DiffV);
DemodulatedSignal = EnvDiffLimV - mean(EnvDiffLimV);
DemodulatedSignal = max(DemodulatedSignal).*DemodulatedSignal;
figure();
plot(time,DemodulatedSignal);
xlabel('time');
ylabel('demodulated signal(t)');
title('demodulated signal in time domain');
frequency = -fs/2:fs/length(X):fs/2-fs/length(X);
DemodulatedSignalF = 1/fs*fftshift(fft(DemodulatedSignal));
figure();
plot(frequency,abs(DemodulatedSignalF));
xlabel('frequency');
ylabel('demodulated signal(f)');
title('demodulated signal in frequency domain');
%%
fs = 50000;
time = 0:1/fs:0.1-1/fs;
fm = 10;
X = cos(2*pi*fm*time);
fc = 1000;
fdelta = 100;
XIntegral = 1/fs*cumtrapz(X);
V = cos(2*pi*fc.*time + 2*pi*fdelta.*XIntegral);
[b,a] = butter(1,10/fs);
V = filter(a,b,V);
Ai = 0.1;
Wi = 100;
V = V + Ai*cos((2*pi*fc+Wi).*time);
V = filter(b,a,V);
DiffV = zeros(1,length(V));
for i = 2:length(V)-2
   DiffV(i) = fs/2*(V(i+1)-V(i-1));
end
EnvDiffLimV = envelope(DiffV);
DemodulatedSignal = EnvDiffLimV - mean(EnvDiffLimV);
DemodulatedSignal = max(DemodulatedSignal).*DemodulatedSignal;
figure();
plot(time,DemodulatedSignal);
xlabel('time');
ylabel('demodulated signal(t)');
title('demodulated signal in time domain');
frequency = -fs/2:fs/length(X):fs/2-fs/length(X);
DemodulatedSignalF = 1/fs*fftshift(fft(DemodulatedSignal));
figure();
plot(frequency,abs(DemodulatedSignalF));
xlabel('frequency');
ylabel('demodulated signal(f)');
title('demodulated signal in frequency domain');