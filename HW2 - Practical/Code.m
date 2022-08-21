% Q1)
%%
% a)
%%
% b)
syms t tau w x(t) c(t) v(t) r(t) HLPF(t) d(t)
fc = 10;
x(t) = exp(-t^2);
figure();
subplot(2,1,1);
fplot(t,x(t));
xlabel('time');
ylabel('x(t)');
title('x signal in time domain');
Xw = fourier(x(t));
subplot(2,1,2);
fplot(w,Xw);
xlabel('omega');
ylabel('X(w)');
title('x signal in frequency domain');
c(t) = cos(2*pi*fc*t);
figure();
subplot(2,1,1);
fplot(t,c(t),[-1 1]);
xlabel('time');
ylabel('c(t)');
title('c signal in time domain');
Cw = fourier(c(t));
subplot(2,1,2);
stem([-20*pi 20*pi],[65 65],'marker','^');
text([-20*pi 20*pi],[75 75],{'\infty','\infty'});
xlim([-40*pi 40*pi]);
ylim([0 80]);
yticks(0:25:50);
xlabel('omega');
ylabel('C(w)');
title('c signal in frequency domain');
v(t) = x(t)*c(t);
figure();
subplot(2,1,1);
fplot(t,v(t));
xlabel('time');
ylabel('v(t)');
title('v signal in time domain');
Vw = fourier(v(t));
subplot(2,1,2);
fplot(w,Vw,[-40*pi 40*pi]);
xlabel('omega');
ylabel('V(w)');
title('v signal in frequency domain');
r(t) = v(t)*c(t);
figure();
subplot(2,1,1);
fplot(t,r(t));
xlabel('time');
ylabel('r(t)');
title('r signal in time domain');
Rw = fourier(r(t));
subplot(2,1,2);
fplot(w,Rw,[-60*pi 60*pi]);
xlim([-60*pi 60*pi]);
ylim([0 1]);
xlabel('omega');
ylabel('R(w)');
title('r signal in frequency domain');
HLPF(t) = sin(20*pi*t)/(pi*t);
HLPFw = rectangularPulse(-20*pi,20*pi,w);
d(t) = int(r(tau)*HLPF(t-tau),tau,-inf,inf);
Dw = Rw*HLPFw;
figure();
subplot(2,1,1);
fplot(t,d(t));
xlim([-5 5]);
ylim([0 1]);
xlabel('time');
ylabel('d(t)');
title('d signal in time domain');
subplot(2,1,2);
fplot(w,Dw);
xlabel('omega');
ylabel('D(w)');
title('d signal in frequency domain');
%%
% Q2)
%%
% a)
[message,message_fs] = audioread('message.wav');
message = message.';
time = 0:1/message_fs:(length(message)-1)/message_fs;
figure();
plot(time,message);
xlim([0 (length(message)-1)/message_fs]);
ylim([-1 1]);
xlabel('time');
ylabel('message');
title('message signal in time domain');
MessageF = 1/message_fs*fftshift(fft(message));
MessageESDF = abs(MessageF).^2;
frequency = -message_fs/2:message_fs/length(message):message_fs/2-message_fs/length(message);
figure();
plot(frequency,MessageESDF);
MessageTotalEnergy = sum(MessageESDF);
PartOfMessageEnergy = MessageESDF(length(message)/2+1);
for i = 1:length(message)/2-1
   if(PartOfMessageEnergy > 0.99*MessageTotalEnergy)
       break;
   end
   PartOfMessageEnergy = PartOfMessageEnergy + MessageESDF(length(message)/2+1+i) + MessageESDF(length(message)/2+1-i);
end
BWFrequency = i*message_fs/length(MessageESDF);
%%
% b)
%%
% c)
L = 14;
[message] = audioread('message.wav');
message = message.';
UpsampleMessage = interp(message,L);
time = 0:1/(L*message_fs):(length(UpsampleMessage)-1)/(L*message_fs);
figure();
plot(time,UpsampleMessage);
xlim([0 (length(UpsampleMessage)-1)/(L*message_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('upsample message');
title('upsample message signal in time domain');
%%
% d)
L = 14;
[message,message_fs] = audioread('message.wav');
message = message.';
UpsampleMessage = interp(message,L);
time = 0:1/(L*message_fs):(length(UpsampleMessage)-1)/(L*message_fs);
UpsampleMessageMultipleInCarrier = UpsampleMessage.*2.*cos(2*pi*1.02e5*time);
BPF = load('BPF102kHz.mat');
FilteredUpsampleMessage = filter(BPF.BPF102kHz.Numerator,1,UpsampleMessageMultipleInCarrier);
figure();
plot(time,FilteredUpsampleMessage);
xlim([0 (length(UpsampleMessage)-1)/(L*message_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('filtered upsample message');
title('filtered upsample message signal in time domain');
FilteredUpsampleMessageF = 1/(L*message_fs)*fftshift(fft(FilteredUpsampleMessage));
frequency = -L*message_fs/2:L*message_fs/length(UpsampleMessage):L*message_fs/2-L*message_fs/length(UpsampleMessage);
figure();
plot(frequency,abs(FilteredUpsampleMessageF));
xlabel('frequency');
ylabel('filtered upsample message');
title('filtered upsample message signal in frequency domain');
%%
% e)
L = 14;
[message,message_fs] = audioread('message.wav');
message = message.';
UpsampleMessage = interp(message,L);
time = 0:1/(L*message_fs):(length(UpsampleMessage)-1)/(L*message_fs);
UpsampleMessageMultipleInCarrier = UpsampleMessage.*2.*cos(2*pi*3e3*time);
BPF = load('BPF3kHz.mat');
FilteredUpsampleMessage = filter(BPF.BPF3kHz.Numerator,1,UpsampleMessageMultipleInCarrier);
figure();
plot(time,FilteredUpsampleMessage);
xlim([0 (length(UpsampleMessage)-1)/(L*message_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('filtered upsample message');
title('filtered upsample message signal in time domain');
FilteredUpsampleMessageF = 1/(L*message_fs)*fftshift(fft(FilteredUpsampleMessage));
frequency = -L*message_fs/2:L*message_fs/length(UpsampleMessage):L*message_fs/2-L*message_fs/length(UpsampleMessage);
figure();
plot(frequency,abs(FilteredUpsampleMessageF));
xlabel('frequency');
ylabel('filtered upsample message');
title('filtered upsample message signal in frequency domain');
UpsampleMessageMultipleInCarrier = FilteredUpsampleMessage.*2.*cos(2*pi*7e3*time);
BPF = load('BPF10kHz.mat');
FilteredUpsampleMessage = filter(BPF.BPF10kHz.Numerator,1,UpsampleMessageMultipleInCarrier);
figure();
plot(time,FilteredUpsampleMessage);
xlim([0 (length(UpsampleMessage)-1)/(L*message_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('filtered upsample message');
title('filtered upsample message signal in time domain');
FilteredUpsampleMessageF = 1/(L*message_fs)*fftshift(fft(FilteredUpsampleMessage));
frequency = -L*message_fs/2:L*message_fs/length(UpsampleMessage):L*message_fs/2-L*message_fs/length(UpsampleMessage);
figure();
plot(frequency,abs(FilteredUpsampleMessageF));
xlabel('frequency');
ylabel('filtered upsample message');
title('filtered upsample message signal in frequency domain');
UpsampleMessageMultipleInCarrier = FilteredUpsampleMessage.*2.*cos(2*pi*92e3*time);
BPF = load('BPF102kHz.mat');
FilteredUpsampleMessage = filter(BPF.BPF102kHz.Numerator,1,UpsampleMessageMultipleInCarrier);
figure();
plot(time,FilteredUpsampleMessage);
xlim([0 (length(UpsampleMessage)-1)/(L*message_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('filtered upsample message');
title('filtered upsample message signal in time domain');
FilteredUpsampleMessageF = 1/(L*message_fs)*fftshift(fft(FilteredUpsampleMessage));
frequency = -L*message_fs/2:L*message_fs/length(UpsampleMessage):L*message_fs/2-L*message_fs/length(UpsampleMessage);
figure();
plot(frequency,abs(FilteredUpsampleMessageF));
xlabel('frequency');
ylabel('filtered upsample message');
title('filtered upsample message signal in frequency domain');
%%
% f)
L = 14;
[message,message_fs] = audioread('message.wav');
message = message.';
UpsampleMessage = interp(message,L);
time = 0:1/(L*message_fs):(length(UpsampleMessage)-1)/(L*message_fs);
UpsampleMessageMultipleInCosCarrier = UpsampleMessage.*cos(2*pi*1.02e5*time);
HilbertTransformOfUpsampleMessage = imag(hilbert(UpsampleMessage));
UpsampleMessageMultipleInSinCarrier = HilbertTransformOfUpsampleMessage.*sin(2*pi*1.02e5*time);
ModulatedMessage = UpsampleMessageMultipleInCosCarrier - UpsampleMessageMultipleInSinCarrier;
figure();
plot(time,ModulatedMessage);
xlim([0 (length(UpsampleMessage)-1)/(L*message_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('modulated message');
title('modulated message with hilbert transform in time domain');
ModulatedMessageF = 1/(L*message_fs)*fftshift(fft(ModulatedMessage));
figure();
plot(frequency,abs(ModulatedMessageF));
xlabel('frequency');
ylabel('modulated message');
title('modulated message with hilbert transform in frequency domain');
%%
% g)
L = 14;
[message,message_fs] = audioread('message.wav');
message = message.';
UpsampleMessage = interp(message,L);
time = 0:1/(L*message_fs):(length(UpsampleMessage)-1)/(L*message_fs);
UpsampleMessageMultipleInCarrier = UpsampleMessage.*2.*cos(2*pi*1.02e5*time);
BPF = load('BPF102kHz.mat');
FilteredUpsampleMessage = filter(BPF.BPF102kHz.Numerator,1,UpsampleMessageMultipleInCarrier);
FilteredUpsampleMessageF = 1/(L*message_fs)*fftshift(fft(FilteredUpsampleMessage));
frequency = -L*message_fs/2:L*message_fs/length(UpsampleMessage):L*message_fs/2-L*message_fs/length(UpsampleMessage);
figure();
plot(frequency,abs(FilteredUpsampleMessageF));
xlabel('frequency');
ylabel('filtered upsample message');
title('filtered upsample message signal in frequency domain');
UpsampleMessageMultipleInCarrier = UpsampleMessage.*2.*cos(2*pi*3e3*time);
BPF = load('BPF3kHz.mat');
FilteredUpsampleMessage = filter(BPF.BPF3kHz.Numerator,1,UpsampleMessageMultipleInCarrier);
UpsampleMessageMultipleInCarrier = FilteredUpsampleMessage.*2.*cos(2*pi*7e3*time);
BPF = load('BPF10kHz.mat');
FilteredUpsampleMessage = filter(BPF.BPF10kHz.Numerator,1,UpsampleMessageMultipleInCarrier);
UpsampleMessageMultipleInCarrier = FilteredUpsampleMessage.*2.*cos(2*pi*92e3*time);
BPF = load('BPF102kHz.mat');
FilteredUpsampleMessage = filter(BPF.BPF102kHz.Numerator,1,UpsampleMessageMultipleInCarrier);
FilteredUpsampleMessageF = 1/(L*message_fs)*fftshift(fft(FilteredUpsampleMessage));
frequency = -L*message_fs/2:L*message_fs/length(UpsampleMessage):L*message_fs/2-L*message_fs/length(UpsampleMessage);
figure();
plot(frequency,abs(FilteredUpsampleMessageF));
xlabel('frequency');
ylabel('filtered upsample message');
title('filtered upsample message signal in frequency domain');
UpsampleMessageMultipleInCosCarrier = UpsampleMessage.*cos(2*pi*1.02e5*time);
HilbertTransformOfUpsampleMessage = imag(hilbert(UpsampleMessage));
UpsampleMessageMultipleInSinCarrier = HilbertTransformOfUpsampleMessage.*sin(2*pi*1.02e5*time);
ModulatedMessage = UpsampleMessageMultipleInCosCarrier - UpsampleMessageMultipleInSinCarrier;
ModulatedMessageF = 1/(L*message_fs)*fftshift(fft(ModulatedMessage));
figure();
plot(frequency,abs(ModulatedMessageF));
xlabel('frequency');
ylabel('modulated message');
title('modulated message with hilbert transform in frequency domain');
%%
% Q3)
%%
% 3.1)
%%
% a)
%%
% b)
%%
% squaring envelope detector
%%
% DSB AM
fs = 1e5;
fc = 5e3;
time = -0.1:1/fs:0.1-1/fs;
InputSignal = sin(2*pi*50*time);
x1 = (1+InputSignal).*cos(2*pi*fc*time);
x2 = x1.^2;
x3 = 2.*x2;
[b,a] = butter(5,5e3/fs);
x4 = filter(b,a,x3);
x5 = real(sqrt(x4));
figure();
plot(time,x5);
xlim([-0.02 0.02]);
ylim([0,2]);
xlabel('time');
ylabel('envelope signal');
title('envelope signal in time domain');
%%
% DSB_SC AM
fs = 1e5;
fc = 5e3;
time = -0.1:1/fs:0.1-1/fs;
InputSignal = sin(2*pi*50*time);
x1 = InputSignal.*cos(2*pi*fc*time);
x2 = x1.^2;
x3 = 2.*x2;
[b,a] = butter(5,5e3/fs);
x4 = filter(b,a,x3);
x5 = real(sqrt(x4));
figure();
plot(time,x5);
xlim([-0.02 0.02]);
ylim([0 1]);
xlabel('time');
ylabel('envelope signal');
title('envelope signal in time domain');
%%
% hilbert transform envelope detector
%%
% DSB AM
fs = 1e5;
fc = 5e3;
time = -0.1:1/fs:0.1-1/fs;
InputSignal = (1+sin(2*pi*50*time)).*cos(2*pi*fc*time);
x1 = imag(hilbert(InputSignal));
x2 = i.*x1;
x3 = x2+InputSignal;
x4 = abs(x3);
figure();
plot(time,x4);
xlim([-0.02 0.02]);
ylim([0,2]);
xlabel('time');
ylabel('envelope signal');
title('envelope signal in time domain');
%%
% DSB_SC AM
fs = 1e5;
fc = 5e3;
time = -0.1:1/fs:0.1-1/fs;
InputSignal = sin(2*pi*50*time).*cos(2*pi*fc*time);
x1 = imag(hilbert(InputSignal));
x2 = i.*x1;
x3 = x2+InputSignal;
x4 = abs(x3);
figure();
plot(time,x4);
xlim([-0.02 0.02]);
ylim([0,1]);
xlabel('time');
ylabel('envelope signal');
title('envelope signal in time domain');
%%
% 3.2)
%%
% a)
fs = 1e5;
fc = 5e3;
time = -0.1:1/fs:0.1-1/fs;
InputSignal = sin(2*pi*50*time).*cos(2*pi*fc*time);
InputSignalModulatedInCarrier = InputSignal.*cos(2*pi*fc*time);
[b,a] = butter(5,5e3/fs);
ModulatedSignal = filter(b,a,InputSignalModulatedInCarrier);
figure();
plot(time,ModulatedSignal);
xlim([-0.02 0.02]);
ylim([-1,1]);
xlabel('time');
ylabel('envelope signal');
title('envelope signal in time domain');
%%
% b)
fs = 1e5;
fc = 5e3;
time = -0.1:1/fs:0.1-1/fs;
InputSignal = sin(2*pi*50*time).*cos(2*pi*fc*time);
InputSignalModulatedInCarrier = InputSignal.*cos(2*pi*fc*time+2*pi*rand());
[b,a] = butter(5,5e3/fs);
ModulatedSignal = filter(b,a,InputSignalModulatedInCarrier);
figure();
plot(time,ModulatedSignal);
xlim([-0.02 0.02]);
ylim([-1,1]);
xlabel('time');
ylabel('envelope signal');
title('envelope signal in time domain');
%%
% Q4)
%%
% a)
[Hello,Hello_fs] = audioread('Hello.mp3');
Hello = (Hello(:,1) + Hello(:,2))/2;
Hello_time = 0:1/Hello_fs:(length(Hello)-1)/Hello_fs;
HelloF = 1/Hello_fs*fftshift(fft(Hello));
Hello_frequency = -Hello_fs/2:Hello_fs/length(Hello):Hello_fs/2-Hello_fs/length(Hello);
figure();
subplot(2,1,1);
plot(Hello_time,Hello);
xlim([0 length(Hello)/Hello_fs]);
ylim([-1 1]);
xlabel('time');
ylabel('Hello');
title('Hello signal in time domain');
subplot(2,1,2);
plot(Hello_frequency,abs(HelloF));
xlabel('frequency');
ylabel('|Hello(f)|');
title('magnitude of Hello signal in frequency domain');
[Goodbye,Goodbye_fs] = audioread('Goodbye.mp3');
Goodbye = (Goodbye(:,1) + Goodbye(:,2))/2;
Goodbye_time = 0:1/Goodbye_fs:(length(Goodbye)-1)/Goodbye_fs;
GoodbyeF = 1/Goodbye_fs*fftshift(fft(Goodbye));
Goodbye_frequency = -Goodbye_fs/2:Goodbye_fs/length(Goodbye):Goodbye_fs/2-Goodbye_fs/length(Goodbye);
figure();
subplot(2,1,1);
plot(Goodbye_time,Goodbye);
xlim([0 length(Goodbye)/Goodbye_fs]);
ylim([-1 1]);
xlabel('time');
ylabel('Goodbye');
title('Goodbye signal in time domain');
subplot(2,1,2);
plot(Goodbye_frequency,abs(GoodbyeF));
xlabel('frequency');
ylabel('|Goodbye(f)|');
title('magnitude of Goodbye signal in frequency domain');
%%
% b)
L = 10;
[Hello,Hello_fs] = audioread('Hello.mp3');
Hello = (Hello(:,1) + Hello(:,2))/2;
UpsampleHello = interp(Hello,L).';
UpsampleHello_time = 0:1/(L*Hello_fs):(length(UpsampleHello)-1)/(L*Hello_fs);
ModulatedUpsampleHello = UpsampleHello.*cos(2*pi*1e5*UpsampleHello_time);
ModulatedUpsampleHelloF = 1/(L*Hello_fs)*fftshift(fft(ModulatedUpsampleHello));
UpsampleHello_frequency = -L*Hello_fs/2:L*Hello_fs/length(UpsampleHello):L*Hello_fs/2-L*Hello_fs/length(UpsampleHello);
figure();
subplot(2,1,1);
plot(UpsampleHello_time,ModulatedUpsampleHello);
xlim([0 length(UpsampleHello)/(L*Hello_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('modulated upsample Hello');
title('modulated upsample Hello with L = 10 signal in time domain');
subplot(2,1,2);
plot(UpsampleHello_frequency,abs(ModulatedUpsampleHelloF));
xlabel('frequency');
ylabel('|modulated upsample Hello (f)|');
title('magnitude of modulated upsample Hello with L = 10 signal in frequency domain');
[Goodbye,Goodbye_fs] = audioread('Goodbye.mp3');
Goodbye = (Goodbye(:,1) + Goodbye(:,2))/2;
UpsampleGoodbye = interp(Goodbye,L).';
UpsampleGoodbye_time = 0:1/(L*Goodbye_fs):(length(UpsampleGoodbye)-1)/(L*Goodbye_fs);
ModulatedUpsampleGoodbye = UpsampleGoodbye.*cos(2*pi*1.3e5*UpsampleGoodbye_time);
ModulatedUpsampleGoodbyeF = 1/(L*Goodbye_fs)*fftshift(fft(ModulatedUpsampleGoodbye));
UpsampleGoodbye_frequency = -L*Goodbye_fs/2:L*Goodbye_fs/length(UpsampleGoodbye):L*Goodbye_fs/2-L*Goodbye_fs/length(UpsampleGoodbye);
figure();
subplot(2,1,1);
plot(UpsampleGoodbye_time,ModulatedUpsampleGoodbye);
xlim([0 length(UpsampleGoodbye)/(L*Goodbye_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('modulated upsample Goodbye');
title('modulated upsample Goodbye with L = 10 signal in time domain');
subplot(2,1,2);
plot(UpsampleGoodbye_frequency,abs(ModulatedUpsampleGoodbyeF));
xlabel('frequency');
ylabel('|modulated upsample Goodbye (f)|');
title('magnitude of modulated upsample Goodbye with L = 10 signal in frequency domain');
%%
% c)
L = 10;
[Hello,Hello_fs] = audioread('Hello.mp3');
Hello = (Hello(:,1) + Hello(:,2))/2;
UpsampleHello = interp(Hello,L).';
UpsampleHello_time = 0:1/(L*Hello_fs):(length(UpsampleHello)-1)/(L*Hello_fs);
ModulatedUpsampleHello = UpsampleHello.*cos(2*pi*1e5*UpsampleHello_time);
[Goodbye,Goodbye_fs] = audioread('Goodbye.mp3');
Goodbye = (Goodbye(:,1) + Goodbye(:,2))/2;
UpsampleGoodbye = interp(Goodbye,L);
UpsampleGoodbye_time = 0:1/(L*Goodbye_fs):(length(UpsampleGoodbye)-1)/(L*Goodbye_fs);
ModulatedUpsampleGoodbye = UpsampleHello.*cos(2*pi*1.3e5*UpsampleGoodbye_time);
RecievedSignal = ModulatedUpsampleHello + ModulatedUpsampleGoodbye;
RecievedSignal_time = 0:1/(L*Hello_fs):(length(UpsampleHello)-1)/(L*Hello_fs);
RecievedSignalF = 1/(L*Goodbye_fs)*fftshift(fft(RecievedSignal));
RecievedSignal_frequency = -L*Hello_fs/2:L*Hello_fs/length(UpsampleHello):L*Hello_fs/2-L*Hello_fs/length(UpsampleHello);
figure();
subplot(2,1,1);
plot(RecievedSignal_time,RecievedSignal);
xlim([0 length(UpsampleHello)/(L*Hello_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('recieved signal');
title('recieved signal in time domain');
subplot(2,1,2);
plot(RecievedSignal_frequency,abs(RecievedSignalF));
xlabel('frequency');
ylabel('|recieved signal (f)|');
title('recieved signal in frequency domain');
RecievedSignalMultipleInCarrier = RecievedSignal.*cos(2*pi*1e5*RecievedSignal_time);
RecievedSignalMultipleInCarrierF = 1/(L*Hello_fs)*fftshift(fft(RecievedSignalMultipleInCarrier));
figure();
subplot(2,1,1);
plot(RecievedSignal_time,RecievedSignalMultipleInCarrier);
xlim([0 length(UpsampleHello)/(L*Hello_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('recieved signal multiple in carrier');
title('recieved signal multiple in carrier in time domain');
subplot(2,1,2);
plot(RecievedSignal_frequency,abs(RecievedSignalMultipleInCarrierF));
xlabel('frequency');
ylabel('|recieved signal multiple in carrier(f)|');
title('recieved signal multiple in carrier in frequency domain');
[b,a] = butter(1,5e3/(L*Hello_fs));
DemodulatedRecievedSignal = filter(b,a,RecievedSignalMultipleInCarrier);
DemodulatedRecievedSignalF = 1/(L*Hello_fs)*fftshift(fft(DemodulatedRecievedSignal));
figure();
subplot(2,1,1);
plot(RecievedSignal_time,DemodulatedRecievedSignal);
xlim([0 length(UpsampleHello)/(L*Hello_fs)]);
ylim([-1 1]);
xlabel('time');
ylabel('demodulated recieved signal');
title('demodulated recieved signal in time domain');
subplot(2,1,2);
plot(RecievedSignal_frequency,abs(DemodulatedRecievedSignalF));
xlabel('frequency');
ylabel('|demodulated recieved signal(f)|');
title('demodulated recieved signal in frequency domain');
%%
% d)
DownsampleDemodulatedRecievedSignal = downsample(DemodulatedRecievedSignal,L);
MSE = immse(DownsampleDemodulatedRecievedSignal,Hello.');
sound(DownsampleDemodulatedRecievedSignal,Hello_fs);
%%
% e)
%%
% n = 1
[b,a] = butter(1,5e3/(L*Hello_fs));
DemodulatedRecievedSignal = filter(b,a,RecievedSignalMultipleInCarrier);
DownsampleDemodulatedRecievedSignal = downsample(DemodulatedRecievedSignal,L);
MSE_1 = immse(DownsampleDemodulatedRecievedSignal,Hello.');
sound(DownsampleDemodulatedRecievedSignal,Hello_fs);
%%
% n = 3
[b,a] = butter(3,5e3/(L*Hello_fs));
DemodulatedRecievedSignal = filter(b,a,RecievedSignalMultipleInCarrier);
DownsampleDemodulatedRecievedSignal = downsample(DemodulatedRecievedSignal,L);
MSE_3 = immse(DownsampleDemodulatedRecievedSignal,Hello.');
sound(DownsampleDemodulatedRecievedSignal,Hello_fs);
[b,a] = butter(5,5e3/(L*Hello_fs));
%%
% n = 5
DemodulatedRecievedSignal = filter(b,a,RecievedSignalMultipleInCarrier);
DownsampleDemodulatedRecievedSignal = downsample(DemodulatedRecievedSignal,L);
MSE_5 = immse(DownsampleDemodulatedRecievedSignal,Hello.');
sound(DownsampleDemodulatedRecievedSignal,Hello_fs);