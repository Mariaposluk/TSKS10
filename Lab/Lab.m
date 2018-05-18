%Laboration i TSKS10 Signaler, information och kommunikation 
%Maria Posluk, marpo758

%%

%clear all
close all

%% Ta ut information om den mottagna signalen y(t)
[y, fs] = audioread('signal-marpo758.wav');
Ts = 1/fs;
Y = abs(fft(y));
lengthy = length(y);
length_half = ceil(lengthy/2);
freq = linspace(0, fs/2, length_half);
%% Rita ut amplitudspektrat för Y(f)
figure('DefaultAxesFontSize',16);
figure(1)
hold on
plot(freq, Y(1:length_half)); %Möjliga bärfrekvenser: 36 kHz, 93 kHz, 150 kz.
%Zoomar in där det finns en märkbar topp men ej nära möjlig bärfrekvens. =>
%ser att f1 och f2 är 140500 Hz och 140501 Hz.
title('Amplitudsspektrum'), xlabel('Frekvens [Hz]'), ylabel('|Y|');
hold off

%% Rita ut amplitudspektrat för Y(f) med normerad
% frekvens för att läsa av bandbredden för de tre
% möjliga bärfrekvenserna för att använda Matlabs filtrering.
norm_freq = linspace(0, 1, length_half);
figure('DefaultAxesFontSize',16);
figure(2)
hold on
plot(norm_freq, Y(1:length_half));
title('Amplitudsspektrum'), xlabel('Normerad frekvens'), ylabel('|Y|');
hold off

%% Rita ut f1 och f2

freqs = ((freq > 1.40495e5) & (freq < 1.40505e5));
figure('DefaultAxesFontSize',16);
figure(3)
hold on
plot(freq(freqs), Y(freqs));
title('Frekvenserna f1 och f2'), xlabel('Frekvens [Hz]'), ylabel('|Y|');
hold off

%% Filtrera kring de möjliga bärfrekvenserna med buttherworth 
[B36,A36] = butter(9, [0.16, 0.20]);
[B93,A93] = butter(9, [0.44, 0.49]);
[B150,A150] = butter(9, [0.73, 0.77]);

y36 = filter(B36, A36, y);
y93 = filter(B93, A93, y);
y150 = filter(B150, A150, y);

%% Räkna ut en tidsvektor för att kunna rita ut y(t)
time = lengthy*Ts; %Längd på filen
t = (0:Ts:time-Ts);

%% Rita ut y(t) kring de möjliga bärfrekvenserna

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
%I figur 3 kan man se att det är 36 kHz som är bärfrekvensen.

%% Räkna ut deltatau

f1 = 140500;
f2 = 140501;
Yf1 = 4940; %från figur 3
Yf2 = 1254; %från figur 3
w = 0.001*(cos(2*pi*f1*t) + cos(2*pi*f2*t));
W = abs(fft(w));

%Läser av värdena i f1 och f2
figure('DefaultAxesFontSize',16);
figure(5)
hold on
plot(freq, W(1:length_half));
title('W(f)'), xlabel('Frekvens [Hz]'), ylabel('|W|');
hold off

Wf1 = 2600;
Wf2 = 2600;

Hf1 = Yf1/Wf1;
Hf2 = Yf2/Wf2;

smallest_diff1 = inf;
smallest_diff2 = inf;
deltatau = (0.001:0.001:0.5);
diff1 = zeros(500:1);
diff2 = zeros(500:1);
for index = 1:500
    H1 = sqrt(1+0.9^2+1.8*cos(2*pi*f1*deltatau(index)));
    H2 = sqrt(1+0.9^2+1.8*cos(2*pi*f2*deltatau(index)));
    diff1(index) = abs(H1-Hf1);
    diff2(index) = abs(H2-Hf2);
    if (diff1(index) <= smallest_diff1 && diff2(index) <= smallest_diff2)
        tau = deltatau(index); %Resultatet på tau = 0.420 s
        smallest_diff1 = diff1(index);
        smallest_diff2 = diff2(index);
    end
end



%% Rita ut tau

figure('DefaultAxesFontSize',16);
figure(6)
plot(deltatau, diff1, 'r')
hold on
plot(deltatau, diff2, 'b')
legend('|H(f1)| - |Y(f1)|/|W(f1)|', '|H(f2)| - |Y(f2)|/|W(f2)|')
title('Bestämma ??'), xlabel('?? [s]'), ylabel('|H(f)| - |Y(f)|/|W(f)|');
hold off

%Figur 6 ger att deltatau = 0,420 s

%% Ta bort ekot

deltasampel = fs*tau;

x(1:deltasampel) = y36(1:deltasampel);
for index = deltasampel+1:lengthy
    x(index) = y36(index) - 0.9*x(index - deltasampel);
end



%% I/Q-demodulering
fc = 36000;
phi = pi/20;
xI = zeros(size(y));
xQ = zeros(size(y));
for index = 1:lengthy
    xI(index) = x(index) * 2 * cos(2*pi*fc/fs*index+phi);
    xQ(index) = x(index) * (-2) * sin(2*pi*fc/fs*index+phi);  
end


%% Filtrera och lyssna
%Lågpassfiltrera med bredden av fc = 36 kHz som gränsfrekvens
[ButterB, ButterA] = butter(9, (0.20-0.16));
filtered_xI = filter(ButterB, ButterA, xI);
filtered_xQ = filter(ButterB, ButterA, xQ);

decimated_xI = decimate(filtered_xI, 40);
decimated_xQ = decimate(filtered_xQ, 40);

%soundsc(decimated_xI); %Skrattar bäst som skrattar sist
soundsc(decimated_xQ); %Även små grytor har öron


