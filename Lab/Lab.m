%Laboration i TSKS10 Signaler, kommunikation och information
%Maria Posluk, marpo758, 950310-0829


%% Ta ut information om den mottagna signalen y(t)
[y, fs] = audioread('signal-marpo758.wav');
Ts = 1/fs;
Y = fftshift(fft(y));
Y = abs(Y);
length_half = ceil(length(y)/2);
freq = linspace(0, fs/2, length_half);
%% Rita ut amplitudspektrat f�r Y(f)
figure(1)
hold on
plot(freq, Y(1:length_half)); %M�jliga b�rfrekvenser: 36 kHz, 93 kHz, 150 kz
title('Amplitudsspektrum'), xlabel('Frekvens [Hz]'), ylabel('|Y|');
hold off

%% Rita ut amplitudspektrat f�r Y(f) med normerad frekvens f�r att l�sa av bandbredden f�r de tre m�jliga b�rfrekvenserna.
norm_freq = linspace(0, 1, length_half);
figure(2)
hold on
plot(norm_freq, Y(1:length_half)); %
title('Amplitudsspektrum'), xlabel('Normerad frekvens'), ylabel('|Y|');
hold off
%% Filtrera kring de m�jliga b�rfrekvenserna med buttherworth 
[B36,A36] = butter(9, [0.16, 0.20]); %ordning p� filtret??
[B93,A93] = butter(9, [0.44, 0.49]);
[B150,A150] = butter(9, [0.73, 0.77]);

y36 = filter(B36, A36, y);
y93 = filter(B93, A93, y);
y150 = filter(B150, A150, y);
%% R�kna ut en tidsvektor f�r att kunna rita ut y(t)
time = length(y)*Ts;
t = (0:Ts:time-Ts);
%% Rita ut y(t) kring de m�jliga b�rfrekvenserna
figure(3)
subplot(3, 1, 1);
plot(t, y36);
title('Vid frekvens 36 kHz'), xlabel('Tid [s]'), ylabel('y(t)');
subplot(3, 1, 2);
plot(t, y93);
title('Vid frekvens 93 kHz'), xlabel('Tid [s]'), ylabel('y(t)');
subplot(3, 1, 3);
plot(t, y150);
title('Vid frekvens 150 kHz'), xlabel('Tid [s]'), ylabel('y(t)');
%I figur 3 kan man se att det �r 36 kHz som �r b�rfrekvensen.

%% Korrelation

korr = xcorr(y36);

figure(4)
plot(korr); 
title('Autokorrelation vid frekvens 36 kHz'), xlabel('Sampel'), ylabel('Korrelation');
% topp till v�nster: 5.032*10^6 sampel, topp i mitten: 5.2*10^6 sampel, topp till h�ger: 5.368*10^6 sampel.
%Avst�nd till mitten (5.368*10^6-5.032*10^6)/2= 1.68*10^5 sampel
deltasampel = 1.68*10^5;
Tau = deltasampel/fs; %Tau = 0.4200 s = 420 ms

%xI = y * 2 * cos(2*pi*fc*t);
%xQ = y * (-2) * sin(2*pi*fc*t);




