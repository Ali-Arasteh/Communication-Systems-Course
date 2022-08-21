% Q1)
%%
% a)
Length = 1;
Frequency = 1000;
Sigma = 1;
T = 0.1;
time = 0:1/Frequency:Length;
Y = RandomDigitalWave(Length,Frequency,Sigma,T);
figure();
plot(time,Y);
title('Random Digital Wave');
xlabel('time');
ylabel('x');
%%
% b)
Number = 5000;
Length = 1;
Frequency = 1000;
Sigma = 1;
T = 0.1;
Maxlag = 1000;
shift = -Maxlag/Frequency:1/Frequency:Maxlag/Frequency;
R = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma,T);
figure();
plot(shift,R);
title('Autocorrelation of random digital wave');
xlabel('Shift');
ylabel('Mean Of R');
%%
% c)
Number = 5000;
Length = 1;
Frequency = 1000;
Sigma = 1;
T = 0.1;
Maxlag = 1000;
shift = -Maxlag/Frequency:1/Frequency:Maxlag/Frequency;
% T1 = 0.1;
% T2 = 0.25;
% T3 = 0.5;
% R11 = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma,T1);
% R12 = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma,T2);
% R13 = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma,T3);
% figure();
% plot(shift,R11,'red');
% hold on;
% plot(shift,R12,'green');
% hold on;
% plot(shift,R13,'blue');
% title('Autocorrelation of random digital wave for different T');
% xlabel('Shift');
% ylabel('Mean Of R');
% legend('T = 0.1','T = 0.25','T = 0.5');
Sigma1 = 1;
Sigma2 = 2.5;
Sigma3 = 5;
R21 = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma1,T);
R22 = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma2,T);
R23 = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma3,T);
figure();
plot(shift,R21,'red');
hold on;
plot(shift,R22,'green');
hold on;
plot(shift,R23,'blue');
title('Autocorrelation of random digital wave for different sigma');
xlabel('Shift');
ylabel('Mean Of R');
legend('sigma = 1','sigma = 2.5','sigma = 5');
%%
% d)
Number = 1000;
Length = 1;
Frequency = 1000;
Sigma = 1;
T = 0.1;
Maxlag = 1000;
R = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma,T);
omega = -pi:pi/Maxlag:pi;
Sx = 1/Frequency*fftshift(fft(R));
MagnitudeOfSx = abs(Sx);
figure();
plot(omega,MagnitudeOfSx);
title('Energy Specturm Density');
xlabel('omega');
ylabel('Sx');
%%
% e)
Number = 1000;
Length = 1;
Frequency = 1000;
Sigma = 1;
T = 0.1;
Maxlag = 1000;
R = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma,T);
for i = 1:Number
    X(i,:) = RandomDigitalWave(Length,Frequency,Sigma,T);
end
omega = -pi:pi/Maxlag:pi;
Sx1 = 1/Frequency*fftshift(fft(R));
MagnitudeOfSx1 = abs(Sx1);
for i = 1:Number
    Xf(i,:) = 1/Frequency*fftshift(fft(X(i,:),2*Maxlag+1));
end
for i = 1:2*Maxlag+1
    Sx2(i) = mean(abs(Xf(:,i)).^2);
end
MagnitudeOfSx2 = abs(Sx2);
figure();
subplot(2,1,1);
plot(omega,MagnitudeOfSx1);
title('Energy Specturm Density way 1');
xlabel('omega');
ylabel('Sx way 1');
subplot(2,1,2);
plot(omega,MagnitudeOfSx2);
title('Energy Specturm Density way 2');
xlabel('omega');
ylabel('Sx way 2');
mse = immse(MagnitudeOfSx1,MagnitudeOfSx2);
%%
% f)
Number = 1000;
Length = 1;
Frequency = 1000;
Sigma = 1;
T = 0.1;
Maxlag = 1000;
shift = -Maxlag/Frequency:1/Frequency:Maxlag/Frequency;
time = 0:1/Frequency:Length;
sin_wave = sin(2*pi*200.*time);
for i = 1:Number
    X(i,:) = RandomDigitalWave(Length,Frequency,Sigma,T);
end
for i = 1:Number
    V(i,:) = X(i,:).*sin_wave;
end
Rv = zeros(Number,2*Maxlag+1);
for i = 1:Number
    Rv(i,:) = 1/Frequency*xcorr(V(i,:),Maxlag);
end
for i = 1:2*Maxlag+1
    R(i) = mean(Rv(:,i));
end
figure();
plot(shift,R);
title('Autocorrelation of V');
xlabel('Shift');
ylabel('Mean Of R');
omega = -pi:pi/Maxlag:pi;
Sv = 1/Frequency*fftshift(fft(R));
MagnitudeOfSv = abs(Sv);
figure();
plot(omega,MagnitudeOfSv);
title('Energy Specturm Density');
xlabel('Omega');
ylabel('Sv');
%%
% g)
%%
% Q2
%%
% a)
[impulse_input,fs] = audioread('clip.wav');
time = 1:length(impulse_input);
time = time./fs;
figure();
plot(time,impulse_input.');
title('Input Message Signal');
xlabel('t(s)');
ylabel('x(t)');
Xf = 1/fs*fftshift(fft(impulse_input.'));
omega = -pi:2*pi/length(impulse_input):pi-1/length(impulse_input);
Sx = abs(Xf).^2;
figure();
plot(omega,Sx);
title('Energy Specturm Density Of Input');
xlabel('Omega');
ylabel('Sx');
%%
% b)
Number = 22050*0.1;
fs = 22050;
lpha = 10;
Shift = 22050*0.02;
time = 0:1/fs:(Number-1)/fs;
impulse_input = zeros(1,Number);   
impulse_input(1) = 1;
h = Channel(impulse_input);
figure();
plot(time,h);
title('Impulse response of channel');
xlabel('time');
ylabel('h');
impulse_multiple_input = zeros(1,Number);   
impulse_multiple_input(1) = alpha;       
h_multiple = Channel(impulse_multiple_input);
figure();
plot(time,h_multiple);
title('Response of channel for impulse multiple input');
xlabel('time');
ylabel('h multiple');
shifted_impulse_input = zeros(1,Number);
shifted_impulse_input(Shift+1) = 1;       
shifted_h = Channel(shifted_impulse_input);
figure();
plot(time,shifted_h);
title('Response of channel for shifted impulse input');
xlabel('time');
ylabel('shifted h');
%%
% c)
fs = 22050;
impulse_input = zeros(1,1001);   
impulse_input(1) = fs;
h = Channel(impulse_input); 
Hf = 1/fs*fftshift(fft(h));
figure();
phasedelay(Hf,fs);
title('Phase delay');
figure();
grpdelay(Hf,fs);
title('Group delay');
%%
% d)
[impulse_input,fs] = audioread('clip.wav');
time = 1:length(impulse_input);
time = time./fs;
y = Channel(impulse_input);
figure();
plot(time,y.');
title('Output Message Signal');
xlabel('time');
ylabel('cy');
Yf = 1/fs*fftshift(fft(y.'));
omega = -pi:2*pi/length(y):pi-1/length(y);
Sy = abs(Yf).^2;
figure();
plot(omega,Sy);
title('Energy Specturm Density Of Output');
xlabel('Omega');
ylabel('Sy');
%%
% e)
%%
% f)
[x,fs] = audioread('clip.wav');
Xf = 1/fs*fftshift(fft(x.'));
y = Channel(x);
Yf = 1/fs*fftshift(fft(y.'));
data = iddata(y,x,1/fs);
np = 3;
nz = 2;
tf = tfest(data,np,nz);
pzmap(tf.Denominator,tf.Numerator);
Number = length(x);
fs = 22050;
M = 5000;
impulse_input = zeros(1,Number);   
impulse_input(1) = fs;
h = Channel(impulse_input); 
Hf = 1/fs*fftshift(fft(h))+1e-15;
Hf_Phase = atan(imag(Hf)./real(Hf));
Hf_Magnitude = abs(Hf);
omega = -pi:2*pi/Number:pi-2*pi/Number;
H_1f_Phase = -Hf_Phase;
H_1f_Magnitude = 1./Hf_Magnitude;
H_1f = H_1f_Magnitude.*exp(i.*H_1f_Phase);
Wf = zeros(1,Number);
Wf((Number-M)/2+1:(Number+M)/2) = 1;
H_1f = H_1f.*Wf;
H_1f_Phase = atan(imag(H_1f)./real(H_1f));
H_1f_Magnitude = abs(H_1f);
figure();
plot(omega,H_1f_Phase);
title('Phase of compensator filter');
xlabel('Omega');
ylabel('\angle H^{-1}');
figure();
plot(omega,H_1f_Magnitude);
title('Magnitude of compensator filter');
xlabel('Omega');
ylabel('|H^{-1}|');
CYf = Yf.*H_1f;
cy = fs*ifft(fftshift(CYf));
time = 0:1/fs:(Number-1)/fs;
figure();
plot(time,real(cy));
title('Compensated Received Message Signal');
xlabel('time');
ylabel('cy');
%%
% Q3
%%
% a)
%%
% b)
fs = 22050;
Length = 0.2;
time = -Length/2+1/(2*fs):1/fs:Length/2-1/(2*fs);
h = 1./(pi.*time);
figure();
plot(time,h);
title('Impulse response of channel');
ylim([-100 100]);
xlabel('time');
ylabel('h');
%%
% c)
fs = 22050;
Length = 1;
time = -Length/2+1/(2*fs):1/fs:Length/2-1/(2*fs);
h = 1./(pi.*time);
M = 300;
w = zeros(1,length(time));
w((length(w)-M)/2+1:(length(w)+M)/2) = 1;
hw = h.*w;
HWf = 1/fs*fftshift(fft(hw));
HWf_Phase = atan(imag(HWf)./real(HWf));
HWf_Magnitude = abs(HWf);
omega = -pi:2*pi/(fs*Length):pi-1/(fs*Length);
figure();
plot(omega,HWf_Phase);
title('Phase of Hf for M = 300');
xlabel('omega');
ylabel('Phase of Hf');
figure();
plot(omega,HWf_Magnitude);
xlim([-4 4]);
ylim([0.8 1.2]);
title('Magnitude of Hf for M = 300');
xlabel('omega');
ylabel('Magnitude of Hf');
%%
% d)
fs = 22050;
Length = 1;
time = -Length/2+1/(2*fs):1/fs:Length/2-1/(2*fs);
h = 1./(pi.*time);
M1 = 100;
w1 = zeros(1,length(time));
w1((length(w1)-M1)/2+1:(length(w1)+M1)/2) = 1;
hw1 = h.*w1;
HW1f = 1/fs*fftshift(fft(hw1));
HW1f_Phase = atan(imag(HW1f)./real(HW1f));
HW1f_Magnitude = abs(HW1f);
M2 = 500;
w2 = zeros(1,length(time));
w2((length(w2)-M2)/2+1:(length(w2)+M2)/2) = 1;
hw2 = h.*w2;
HW2f = 1/fs*fftshift(fft(hw2));
HW2f_Phase = atan(imag(HW2f)./real(HW2f));
HW2f_Magnitude = abs(HW2f);
M3 = 1000;
w3 = zeros(1,length(time));
w3((length(w3)-M3)/2+1:(length(w3)+M3)/2) = 1;
hw3 = h.*w3;
HW3f = 1/fs*fftshift(fft(hw3));
HW3f_Phase = atan(imag(HW3f)./real(HW3f));
HW3f_Magnitude = abs(HW3f);
omega = -pi:2*pi/(fs*Length):pi-1/(fs*Length);
figure();
plot(omega,HW1f_Phase,'red');
hold on;
plot(omega,HW2f_Phase,'green');
hold on;
plot(omega,HW3f_Phase,'blue');
title('Phase of Hf for different M');
xlabel('omega');
ylabel('Phase of Hf');
legend('M1 = 100','M2 = 500','M3 = 1000');
figure();
plot(omega,HW1f_Magnitude,'red');
hold on;
plot(omega,HW2f_Magnitude,'green');
hold on;
plot(omega,HW3f_Magnitude,'blue');
xlim([-4 4]);
ylim([0.8 1.2]);
title('Magnitude of Hf for different M');
xlabel('omega');
ylabel('Magnitude of Hf');
legend('M1 = 100','M2 = 500','M3 = 1000');
%%
% e)
%%
% part 1)
fs = 22050;
Length = 1;
time = -Length/2+1/(2*fs):1/fs:Length/2-1/(2*fs);
h = 1./(pi.*time);
M = 300;
w = zeros(1,length(time));
for i = 0:M
   w(i+(length(w)-M)/2+1) = 0.54-0.46*cos(2*pi*i/M);
end
hw = h.*w;
HWf = 1/fs*fftshift(fft(hw));
HWf_Phase = atan(imag(HWf)./real(HWf));
HWf_Magnitude = abs(HWf);
omega = -pi:2*pi/(fs*Length):pi-1/(fs*Length);
figure();
plot(omega,HWf_Phase);
title('Phase of Hf for M = 300');
xlabel('omega');
ylabel('Phase of Hf');
figure();
plot(omega,HWf_Magnitude);
title('Magnitude of Hf for M = 300');
xlabel('omega');
ylabel('Magnitude of Hf');
%%
% part 2)
fs = 22050;
Length = 1;
time = -Length/2+1/(2*fs):1/fs:Length/2-1/(2*fs);
h = 1./(pi.*time);
M1 = 100;
w1 = zeros(1,length(time));
for i = 0:M1
   w1(i+(length(w1)-M1)/2+1) = 0.54-0.46*cos(2*pi*i/M1);
end
hw1 = h.*w1;
HW1f = 1/fs*fftshift(fft(hw1));
HW1f_Phase = atan(imag(HW1f)./real(HW1f));
HW1f_Magnitude = abs(HW1f);
M2 = 500;
w2 = zeros(1,length(time));
for i = 0:M2
   w2(i+(length(w2)-M2)/2+1) = 0.54-0.46*cos(2*pi*i/M2);
end
hw2 = h.*w2;
HW2f = 1/fs*fftshift(fft(hw2));
HW2f_Phase = atan(imag(HW2f)./real(HW2f));
HW2f_Magnitude = abs(HW2f);
M3 = 1000;
w3 = zeros(1,length(time));
for i = 0:M3
   w3(i+(length(w3)-M3)/2+1) = 0.54-0.46*cos(2*pi*i/M3);
end
hw3 = h.*w3;
HW3f = 1/fs*fftshift(fft(hw3));
HW3f_Phase = atan(imag(HW3f)./real(HW3f));
HW3f_Magnitude = abs(HW3f);
omega = -pi:2*pi/(fs*Length):pi-1/(fs*Length);
figure();
plot(omega,HW1f_Phase,'red');
hold on;
plot(omega,HW2f_Phase,'green');
hold on;
plot(omega,HW3f_Phase,'blue');
title('Phase of Hf for different M');
xlabel('omega');
ylabel('Phase of Hf');
legend('M = 100','M = 500','M = 1000');
figure();
plot(omega,HW1f_Magnitude,'red');
hold on;
plot(omega,HW2f_Magnitude,'green');
hold on;
plot(omega,HW3f_Magnitude,'blue');
title('Magnitude of Hf for different M');
xlabel('omega');
ylabel('Magnitude of Hf');
legend('M = 100','M = 500','M = 1000');
%%
% f)
%%
% optional part
%%
% part 1)
fs = 22050;
Length = 1;
time = -Length/2+1/(2*fs):1/fs:Length/2-1/(2*fs);
h = 1./(pi.*time);
M = 300;
w = zeros(1,length(time));
for i = 0:M
   w(i+(length(w)-M)/2+1) = besseli(0,15*pi*sqrt(1-(2*i/M-1)^2))/besseli(0,15*pi);
end
hw = h.*w;
HWf = 1/fs*fftshift(fft(hw));
HWf_Phase = atan(imag(HWf)./real(HWf));
HWf_Magnitude = abs(HWf);
omega = -pi:2*pi/(fs*Length):pi-1/(fs*Length);
figure();
plot(omega,HWf_Phase);
title('Phase of Hf for M = 300');
xlabel('omega');
ylabel('Phase of Hf');
figure();
plot(omega,HWf_Magnitude);
xlim([-4 4]);
ylim([0 1.2]);
title('Magnitude of Hf for M = 300');
xlabel('omega');
ylabel('Magnitude of Hf');
%%
% part 2)
fs = 22050;
Length = 1;
time = -Length/2+1/(2*fs):1/fs:Length/2-1/(2*fs);
h = 1./(pi.*time);
M1 = 100;
w1 = zeros(1,length(time));
for i = 0:M1
   w1(i+(length(w1)-M1)/2+1) = besseli(0,15*pi*sqrt(1-(2*i/M1-1)^2))/besseli(0,15*pi);
end
hw1 = h.*w1;
HW1f = 1/fs*fftshift(fft(hw1));
HW1f_Phase = atan(imag(HW1f)./real(HW1f));
HW1f_Magnitude = abs(HW1f);
M2 = 500;
w2 = zeros(1,length(time));
for i = 0:M2
   w2(i+(length(w2)-M2)/2+1) = besseli(0,15*pi*sqrt(1-(2*i/M2-1)^2))/besseli(0,15*pi);
end
hw2 = h.*w2;
HW2f = 1/fs*fftshift(fft(hw2));
HW2f_Phase = atan(imag(HW2f)./real(HW2f));
HW2f_Magnitude = abs(HW2f);
M3 = 1000;
w3 = zeros(1,length(time));
for i = 0:M3
   w3(i+(length(w3)-M3)/2+1) = besseli(0,15*pi*sqrt(1-(2*i/M3-1)^2))/besseli(0,15*pi);
end
hw3 = h.*w3;
HW3f = 1/fs*fftshift(fft(hw3));
HW3f_Phase = atan(imag(HW3f)./real(HW3f));
HW3f_Magnitude = abs(HW3f);
omega = -pi:2*pi/(fs*Length):pi-1/(fs*Length);
figure();
plot(omega,HW1f_Phase,'red');
hold on;
plot(omega,HW2f_Phase,'green');
hold on;
plot(omega,HW3f_Phase,'blue');
title('Phase of Hf for different M');
xlabel('omega');
ylabel('Phase of Hf');
legend('M = 100','M = 500','M = 1000');
figure();
plot(omega,HW1f_Magnitude,'red');
hold on;
plot(omega,HW2f_Magnitude,'green');
hold on;
plot(omega,HW3f_Magnitude,'blue');
xlim([-4 4]);
ylim([0 1.2]);
title('Magnitude of Hf for different M');
xlabel('omega');
ylabel('Magnitude of Hf');
legend('M = 100','M = 500','M = 1000');
%%
% functions
function Y = RandomDigitalWave(Length,Frequency,Sigma,T)
Delay = T*rand();
Number_Of_Pulse = ceil((Length-Delay))/T;
A = Sigma*randn(1,Number_Of_Pulse);
Y = zeros(1,floor(Length*Frequency)+1);
for i = 1:floor(Delay*Frequency)
    Y(i)=0;
end
for x = 1:Number_Of_Pulse-1
    for i = ceil(Frequency*(Delay+T*(x-1))):floor(Frequency*(Delay+T*(x)))
        Y(i) = A(x);
    end
end
for i = ceil(Frequency*(Delay+T*(Number_Of_Pulse-1))):floor(Frequency*Length)
        Y(i) = A(Number_Of_Pulse);
end
end
function R = AutocorrelationOfRandomDigitalWave(Number,Maxlag,Length,Frequency,Sigma,T)
for i = 1:Number
    X(i,:) = RandomDigitalWave(Length,Frequency,Sigma,T);
end
Rx = zeros(Number,2*Maxlag+1);
for i = 1:Number
    Rx(i,:) = 1/Frequency*xcorr(X(i,:),Maxlag);
end
for i = 1:2*Maxlag+1
    R(i) = mean(Rx(:,i));
end
end
