%% Ta reda p? b?rfrekvensen Fc

clear all
close all

% Fs samplingsfrekvens
[y,Fs] = audioread('signal-marpo758.wav');

% Tidsspecifikationer
dt = 1/Fs;                     % sekunder per sampel
StopTime = 13;                 % fill?ngd i sekunder
t = (0:dt:StopTime-dt)';
N = size(t,1);

% Fouriertransform
Y_unshifted = fft(y)/N;
Y = fftshift(Y_unshifted);
absY = abs(Y)/N;

% Den ska skalas med N pga Matlab
% g?r som den vill om man inte specificerar

% Vi skiftar eftersom matlab l?gger resultatet from fft som 
% 0 1 2 3 ... -3 -2 -1. Efter skiften blir det -3 -2 -1 0 1 2 3

% Frekvensspecifikationer
dF = Fs/N;                      % Frekvensuppl?sning, i hertz
frq = -Fs/2:dF:Fs/2-dF;         % hertz

%% Plotta spektrumet:
figure;
plot(frq,absY);
xlabel('Frekvens (Hz)');
title('Frekvensspektrum f?r den mottagna signalen');
grid on

% B?rfrekvensen ?r n?gon av 36 kHz, 93 kHz eller 150 kHz
%% f1 och f2 ?r den smala spiken.

freqs = ((frq > 1.4e5) & (frq < 1.42e5));
figure
plot(frq(freqs), absY(freqs));

%% Detta ger att :
f1 = 140500;
f2 = 140501;

%% Butterworth-filtrera ut alla tre m?jliga b?rfrekvenser
B = 10000; % Bandbredd
fc1 = 36000;
fc2 = 93000;
fc3 = 150000;

%% fc1
figure
[b,a] = butter(10,[fc1-B fc1+B]/(Fs/2));
yf1 = filter(b, a, y);
subplot(3,1,1)
plot(t,yf1) 
grid on
    
% fc2
[b,a] = butter(10,[fc2-B fc2+B]/(Fs/2));
yf2 = filter(b, a, y);
subplot(3,1,2)
%figure
plot(t,yf2) 
grid on

% fc3
[b,a] = butter(10,[fc3-B fc3+B]/(Fs/2));
yf3 = filter(b, a, y);

subplot(3,1,3)
plot(t,yf3) 
grid on

%% Jag gissar att det ?r fc2 eftersom ploten ser 
% ut som en icke-periodisk signal som inte ?r brus.
fc = fc1;
my_signal = yf1;

%% Skapa w
w = 0.001*(cos(2*pi*f1*t) + cos(2*pi*f2*t));
W_unshifted = fft(w)/N;
W = fftshift(W_unshifted);
absW = abs(W)/N;

figure
plot(frq, absW)

% Eftersom shiften placerar om frekvenserna fr?n -fs/2 till fs/2 s? ?r den 
% frekvensen som jag letar efter placerad p? andra halvan av spektrument. 
% Allts? l?ttare att leta i den oshiftade.
% Man m?ste ocks? dela med dF eftersom varje steg m?ste vara en
% frekvensuppl?sning
Yf1 = abs((Y_unshifted((f1/dF)+1)));
Yf2 = abs((Y_unshifted((f2/dF)+1)));

Wf1 = abs(W_unshifted((f1/dF)+1));
Wf2 = abs(W_unshifted((f2/dF)+1));

Hf1 = Yf1/Wf1;
Hf2 = Yf2/Wf2;


%% Hitta tau2-tau1
threshold = 0.01;

diff1 = zeros(1,500);
diff2 = zeros(1,500);
delta_tau = (0.001:0.001:0.5);

for i = 1:500
    H1 = abs(1 + 0.9*exp(-1i*2*pi*f1*delta_tau(i)));
    H2 = abs(1 + 0.9*exp(-1i*2*pi*f2*delta_tau(i)));
    diff1(i) = abs(H1 - Hf1);
    diff2(i) = abs(H2 - Hf2);
    if (diff1(i) < threshold) && (diff2(i) < threshold)
        tau_diff = delta_tau(i);
    end
end


%% Visa m?jliga v?rden p? tau2-tau1
figure
plot(delta_tau, diff1, 'r')
hold on
plot(delta_tau, diff2, 'b')
legend('abs(H(f1)) - abs(Y(f1))/abs(W(f1))', 'abs(H(f2)) - abs(Y(f2))/abs(W(f2))')


%% Ta bort eko fr?n yf2

% Antal sampel inom tidsf?rdr?jningen
N_tau = tau_diff * Fs;
% G? fr?n d?r ekot b?rjar till slutet.(Plus ett f?r array b?rjar med 1)
for i = N_tau +1:N
    % Ta bort signal en tidsf?rskjutning bak?t
    my_signal(i) = my_signal(i) - 0.9 * my_signal(i-N_tau);
end


%% Demodulera den valda signalen

[b, a] = butter(10, B/(Fs/2), 'low');
d = 4; % Fasvridning pga tidsf?rskjutning, testade olika v?rden
I = 2*cos(2*pi*fc*t + d);
Q = 2*sin(2*pi*fc*t + d);
%%
x_1I = filter(b, a, my_signal.*I); 
x_1Q = -1*filter(b, a, my_signal.*Q);

x_1I = decimate(x_1I, 10);
x_1Q = decimate(x_1Q, 10);


%% Spela upp ljuden
soundsc(x_1I, Fs/10); % Ju fler kockar desto s?mre soppa
soundsc(x_1Q, Fs/10); % Man ska inte g? ?ver ?n efter vatten

