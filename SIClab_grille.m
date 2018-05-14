%% Hitta Fc

[y,Fs] = audioread('signal-chrlu367.wav');
t = linspace(0,length(y)/Fs,length(y));
%plot(t,y);

L = length(y);
T = 1/Fs;

Y = fft(y);
fnorm = linspace(0, 1, L/2);
f = linspace(0, Fs/2, L/2);

plot(f,abs(Y(1:(L/2))));

%kandidater från ploten är att Fc = 131, 36, 93

[b,a] = butter(10, [0.14 0.22]); %nomerade i tidsdomänen

[d,c] = butter(10, [0.44 0.49]);

[f,e] = butter(10, [0.62 0.7]);

y1 = filter(b, a, y);
y2 = filter(d, c, y);
y3 = filter(f, e, y);



%plot(t, y1);
%plot(t, y2);
%plot(t, y3);

% Fc = 36
% f1 och f2 159500 och 159001
%% delta tau

f1 = 159000;
f2 = 159001;
R1 = 1;
R2 = 0.9;

%|Y(f1)| = 4940 se plot i förra section
%|Y(f2)| = 1254 se plot i förra section

[y,Fs] = audioread('signal-chrlu367.wav');

T = 1/Fs;
L = length(y);
t = linspace(0,L/Fs,L);
f = linspace(0, Fs/2, L/2);

w_t = 0.001*(cos(2*pi*f1*t) + cos(2*pi*f2*t));

%HF1ABS = R1^2 + R2^2 + 2*R1*R2*cos(2*pi*f1*tau);
%HF2ABS = R1^2 + R2^2 + 2*R1*R2*cos(2*pi*f2*tau);

W_F = fft(w_t);
fnorm = linspace(0, 1, L/2);

%plot(f ,abs(W_F(1:L/2))); 

% |W(f1)| = 2034 se plot ovan
% |W(f2)| = 1913

% |H(f1)| = 4940/2034  
% |H(f2)| = 1254/1913

Hf1abs = 4940/2034;
Hf2abs = 1254/1913;

Yf1abs = 4940;
Yf2abs = 1252;
diff1 = inf;
diff2 = inf;
close_tau1 = inf;
close_tau2 = inf;

for tau = 0:0.001:0.5
    HF1ABS = sqrt(R1^2 + R2^2 + 2*R1*R2*cos(2*pi*f1*tau));
    diff_cal1 = abs(Hf1abs -  HF1ABS);
    if(diff_cal1 < diff1)
        diff1 = diff_cal1;
        close_tau1 = tau;
    end
    
        HF2ABS = sqrt(R1^2 + R2^2 + 2*R1*R2*cos(2*pi*f2*tau));
    diff_cal2 = abs(Hf2abs -  HF2ABS);
    if(diff_cal2 < diff2)
        diff2 = diff_cal2;
        close_tau2 = tau;
    end
end


%tau1 = 1 ms
%tau2 = 389 ms

%% For real???

% Detta är en kopia från alfred. BYT SKIT!!!!!!
%y(t + tau1) = x(t) + 0.9x(t + tau1 - tau2)
%x(t) = y(t + tau1) - 0.9x(t + tau1 - tau2)

deltatau = 0.420; %Jag har fel delta tau
fc = 36*10^3;

[y,Fs] = audioread('signal-chrlu367.wav');

deltasampel = deltatau*Fs; % delar upp mina sampel i delta sampels för att kansilera ekot. 

% 13sekunder /  deltatau för alla delar. 
x = zeros(5200000, 1);
x(1:deltasampel) = y(1:deltasampel);
for i = 0:28
    raw = x((1 + deltasampel*i) : deltasampel*(i+1));
    echo = y((1 + deltasampel*(i+1)) : deltasampel*(i+2));
    
    x((1 + deltasampel*(i+1)) : (deltasampel*(i+2))) = echo - 0.9*raw;
end

t_scale = (1/Fs)*transpose(0:(5200000-1));


[b,a] = butter(10, [0.14 0.22]); %nomerade i tidsdomänen
x1 = filter(b, a, x); %filter på x 

%I/Q-demodulera
I_filter = 2*cos(2*pi*fc*t_scale + pi/6); %Testa fram fasvriding
Q_filter = -2*sin(2*pi*fc*t_scale + pi/6);
xI = I_filter.*x1;
xQ = Q_filter.*x1;

% Sampla ner till tillåten fs
i_sound=decimate(xI, 50);
q_sound=decimate(xQ, 50);

soundsc(i_sound);

