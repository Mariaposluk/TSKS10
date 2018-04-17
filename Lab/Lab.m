%Laboration i TSKS10 Signaler, kommunikation och information
%Maria Posluk, marpo758, 950310-0829


%%
[y, fs] = audioread('signal-marpo758.wav');
Ts = 1/fs;
Y = fftshift(fft(y));
Y = abs(Y);
length_half = ceil(length(y)/2);
freq = linspace(0, fs/2, length_half);
%%
figure(1)
hold on
plot(freq, Y(1:length_half)); %36 kHz, 93 kHz, 150 kz
title('Amplitudsspektrum'), xlabel('Frekvens [Hz]'), ylabel('|Y|');
hold off

%%
norm_freq = linspace(0, 1, length_half);
figure(2)
hold on
plot(norm_freq, Y(1:length_half)); %
title('Amplitudsspektrum'), xlabel('Normerad frekvens'), ylabel('|Y|');
hold off
%%
[B36,A36] = butter(9, [0.16, 0.20]); %hur v�ljer man ordningen p� filtret?
[B93,A93] = butter(9, [0.44, 0.49]);
[B150,A150] = butter(9, [0.73, 0.77]);
%%
y36 = filter(B36, A36, y);
y93 = filter(B93, A93, y);
y150 = filter(B150, A150, y);
%%
time = length(y)*Ts;
t = (0:Ts:time-Ts);
%%
figure(3)
subplot(3, 1, 1);
plot(t, y36);
subplot(3, 1, 2);
plot(t, y93);
subplot(3, 1, 3);
plot(t, y150);

%xI = y * 2 * cos(2*pi*fc*t);
%xQ = y * (-2) * sin(2*pi*fc*t);




