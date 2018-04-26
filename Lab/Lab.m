%Laboration i TSKS10 Signaler, kommunikation och information
%Maria Posluk, marpo758
%Fråga på labben: 
%fftshift?


%% Ta ut information om den mottagna signalen y(t)
[y, fs] = audioread('signal-marpo758.wav');
Ts = 1/fs;
Y = abs(fft(y));
length_half = ceil(length(y)/2);
freq = linspace(0, fs/2, length_half);
%% Rita ut amplitudspektrat för Y(f)
figure(1)
hold on
plot(freq, Y(1:length_half)); %Möjliga bärfrekvenser: 36 kHz, 93 kHz, 150 kz.
%Zoomar in där det finns en märkbar topp men ej nära möjlig bärfrekvens. =>
%ser att f1 och f2 är 140500 Hz och 140501 Hz.
title('Amplitudsspektrum'), xlabel('Frekvens [Hz]'), ylabel('|Y|');
hold off

%% Rita ut amplitudspektrat för Y(f) med normerad frekvens för att läsa av bandbredden för de tre möjliga bärfrekvenserna.
norm_freq = linspace(0, 1, length_half);
figure(2)
hold on
plot(norm_freq, Y(1:length_half));
title('Amplitudsspektrum'), xlabel('Normerad frekvens'), ylabel('|Y|');
hold off
%% Filtrera kring de möjliga bärfrekvenserna med buttherworth 
[B36,A36] = butter(9, [0.16, 0.20]); %ordning på filtret??
[B93,A93] = butter(9, [0.44, 0.49]);
[B150,A150] = butter(9, [0.73, 0.77]);

y36 = filter(B36, A36, y);
y93 = filter(B93, A93, y);
y150 = filter(B150, A150, y);
%% Räkna ut en tidsvektor för att kunna rita ut y(t)
time = length(y)*Ts; %antalet sampel / sampelfrekvens
t = (0:Ts:time-Ts);
%% Rita ut y(t) kring de möjliga bärfrekvenserna
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
%I figur 3 kan man se att det är 36 kHz som är bärfrekvensen.

%% Korrelation och tidsfördröjning

korr = xcorr(y36);

figure(4)
plot(korr); 
title('Autokorrelation vid frekvens 36 kHz'), xlabel('Sampel'), ylabel('Korrelation');
% topp till vänster: 5.032*10^6 sampel, topp i mitten: 5.2*10^6 sampel, topp till höger: 5.368*10^6 sampel.
%Avstånd till mitten (5.368*10^6-5.032*10^6)/2= 1.68*10^5 sampel
deltasampel = 1.68*10^5;
Tau = deltasampel/fs; %Tau = 0.4200 s = 420 ms

x(1:deltasampel) = y36(1:deltasampel);
for index = deltasampel+1:length(y)   
    % y(t) = x(t-Tau1) + 0.9x(t-Tau2)
    % x(t+(Tau2-Tau1)) = y(t+(Tau2-Tau1)) - 0.9x(t) dubbelkolla sen
    x(index) = y36(index) - 0.9*x(index - deltasampel);
end
    
%% Sigsys metoden (isället för korrelation)

% %H = 1*exp(-1i*2*pi*f*tau1) + 0.9*exp(-1i*2*pi*f*tau2); 
% absH2 = abs(1^2 + 0.9^2 + 2*1*cos(2*pi*f*deltatau)); %(absH)^2
% absY = abs(X)*sqrt(absH2);
% w = 0.001*(cos(2*pi*f1*t)+cos(2*pi*f2*t));
% W = fft(w);


%% I/Q-demodulering
fc = 36000;
phi = 3*pi/5; %testar
xI = zeros(size(y));
xQ = zeros(size(y));
for index = 1:length(y)
    xI(index) = x(index) * 2 * cos(2*pi*fc/fs*index+phi);
    xQ(index) = x(index) * (-2) * sin(2*pi*fc/fs*index+phi);    
end
%xI = y * 2 * cos(2*pi*fc*t);
%xQ = y * (-2) * sin(2*pi*fc*t);

%% Filtrera och lyssna
%Lågpassfiltrera med bredden av fc = 36 kHz som gränsfrekvens
[ButterB, ButterA] = butter(6, (0.2081-0.1524));
filtered_xI = filter(ButterB, ButterA, xI);
filtered_xQ = filter(ButterB, ButterA, xQ);

decimated_xI = decimate(filtered_xI, 10);
decimated_xQ = decimate(filtered_xQ, 10);

soundsc(decimated_xI, fs/10);
soundsc(decimated_xQ, fs/10);





