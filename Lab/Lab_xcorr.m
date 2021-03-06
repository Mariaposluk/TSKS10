%Laboration i TSKS10 Signaler, kommunikation och information
%Maria Posluk, marpo758

%%

clear all
close all

%% Ta ut information om den mottagna signalen y(t)
[y, fs] = audioread('signal-marpo758.wav');
Ts = 1/fs;
Y = abs(fft(y));
lengthy = length(y);
length_half = ceil(lengthy/2);
freq = linspace(0, fs/2, length_half);
%% Rita ut amplitudspektrat f�r Y(f)
figure(1)
hold on
plot(freq, Y(1:length_half)); %M�jliga b�rfrekvenser: 36 kHz, 93 kHz, 150 kz.
%Zoomar in d�r det finns en m�rkbar topp men ej n�ra m�jlig b�rfrekvens. =>
%ser att f1 och f2 �r 140500 Hz och 140501 Hz.
title('Amplitudsspektrum'), xlabel('Frekvens [Hz]'), ylabel('|Y|');
hold off

%% Rita ut amplitudspektrat f�r Y(f) med normerad
% frekvens f�r att l�sa av bandbredden f�r de tre
% m�jliga b�rfrekvenserna f�r att anv�nda Matlabs filtrering.
norm_freq = linspace(0, 1, length_half);
figure(2)
hold on
plot(norm_freq, Y(1:length_half));
title('Amplitudsspektrum'), xlabel('Normerad frekvens'), ylabel('|Y|');
hold off

%% Rita ut f1 och f2

freqs = ((freq > 1.40495e5) & (freq < 1.40505e5));
figure(3)
hold on
plot(freq(freqs), Y(freqs));
title('Frekvenserna f1 och f2'), xlabel('Frekvens [Hz]'), ylabel('|Y|');
hold off

%% Filtrera kring de m�jliga b�rfrekvenserna med buttherworth 
[B36,A36] = butter(9, [0.16, 0.20]);
[B93,A93] = butter(9, [0.44, 0.49]);
[B150,A150] = butter(9, [0.73, 0.77]);

y36 = filter(B36, A36, y);
y93 = filter(B93, A93, y);
y150 = filter(B150, A150, y);
%% R�kna ut en tidsvektor f�r att kunna rita ut y(t)
time = length(y)*Ts; %L�ngd p� filen
t = (0:Ts:time-Ts);
%% Rita ut y(t) kring de m�jliga b�rfrekvenserna
figure(4)
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

%% Korrelation och tidsf�rdr�jning

korr = xcorr(y36);

figure(5)
plot(korr); 
title('Autokorrelation vid frekvens 36 kHz'), xlabel('Sampel'), ylabel('Korrelation');
% topp till v�nster: 5.032*10^6 sampel, topp i mitten: 5.2*10^6 sampel, topp till h�ger: 5.368*10^6 sampel.
%Avst�nd till mitten (5.368*10^6-5.032*10^6)/2= 1.68*10^5 sampel
deltasampel = 1.68*10^5;
Tau = deltasampel/fs; %Tau = 0.4200 s = 420 ms

x(1:deltasampel) = y36(1:deltasampel);
for index = deltasampel+1:lengthy
    x(index) = y36(index) - 0.9*x(index - deltasampel);
end



%% I/Q-demodulering
fc = 36000;
phi = 11*pi/20;
xI = zeros(size(y));
xQ = zeros(size(y));
for index = 1:lengthy
    xI(index) = x(index) * 2 * cos(2*pi*fc/fs*index+phi);
    xQ(index) = x(index) * (-2) * sin(2*pi*fc/fs*index+phi);
end

%% Filtrera och lyssna
%L�gpassfiltrera med bredden av fc = 36 kHz som gr�nsfrekvens
[ButterB, ButterA] = butter(9, (0.21-0.15));
filtered_xI = filter(ButterB, ButterA, xI);
filtered_xQ = filter(ButterB, ButterA, xQ);

decimated_xI = decimate(filtered_xI, 10);
decimated_xQ = decimate(filtered_xQ, 10);

%soundsc(decimated_xI, fs/10); %�ven sm� grytor har �ron
soundsc(decimated_xQ, fs/10); %Skrattar b�st som skrattar sist


