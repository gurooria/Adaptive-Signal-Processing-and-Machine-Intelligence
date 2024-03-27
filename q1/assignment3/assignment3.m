%% 3.1

clc
clear
close all

x_128 = randn(128, 1);
x_256 = randn(256, 1);
x_512 = randn(512, 1);

pgm_128 = pgm(x_128);
pgm_256 = pgm(x_256);
pgm_512 = pgm(x_512);

% plots
subplot(1,3,1)
x_axis = [1:128];
x_axis = x_axis./128;
plot(x_axis, pgm_128, 'lineWidth', 0.8);
grid minor;
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15)
title('Periodogram of x for N=128', fontSize=15)

subplot(1,3,2)
x_axis = [1:256];
x_axis = x_axis./256;
plot(x_axis, pgm_256, 'lineWidth', 0.8);
grid minor;
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15)
title('Periodogram of x for N=256', fontSize=15)

subplot(1,3,3)
x_axis = [1:512];
x_axis = x_axis./512;
plot(x_axis, pgm_512, 'lineWidth', 0.8);
grid minor;
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15)
title('Periodogram of x for N=512', fontSize=15)

%% 3.1.1

% filter
a = 1;
b = 0.2*[1 1 1 1 1];

% smoothing
smooth_128 = filter(b, a, pgm_128);
smooth_256 = filter(b, a, pgm_256);
smooth_512 = filter(b, a, pgm_512);

% plots
subplot(1,3,1)
x_axis = [1:128];
x_axis = x_axis./128;
plot(x_axis, smooth_128, 'lineWidth', 0.8);
grid minor;
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15)
title('Smoothed periodogram of x for N=128', fontSize=13)

subplot(1,3,2)
x_axis = [1:256];
x_axis = x_axis./256;
plot(x_axis, smooth_256, 'lineWidth', 0.8);
grid minor;
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15)
title('Smoothed periodogram of x for N=256', fontSize=13)

subplot(1,3,3)
x_axis = [1:512];
x_axis = x_axis./512;
plot(x_axis, smooth_512, 'lineWidth', 0.8);
grid minor;
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15)
title('Smoothed periodogram of x for N=512', fontSize=13)

%% 3.1.2

clc
clear
close all

y = randn(1024,1);
PSD = pgm(y);
PSD = reshape(PSD,[128 8]); % 128 rows, 8 columns for 8 sets

x_axis = [1:128];
x_axis = x_axis./128;

% plotting
subplot(2,4,1)
plot(x_axis,PSD(:,1),'lineWidth', 0.8)
title('1st Segment Periodogram', fontSize=12)
xlabel('Normalized Frequency', fontSize=12)
ylabel('Magnitude', fontSize=12)
grid minor

subplot(2,4,2)
plot(x_axis,PSD(:,2),'lineWidth', 0.8)
title('2nd Segment Periodogram', fontSize=12)
xlabel('Normalized Frequency', fontSize=12)
ylabel('Magnitude', fontSize=12)
grid minor

subplot(2,4,3)
plot(x_axis,PSD(:,3),'lineWidth', 0.8)
title('3rd Segment Periodogram', fontSize=12)
xlabel('Normalized Frequency', fontSize=12)
ylabel('Magnitude', fontSize=12)
grid minor

subplot(2,4,4)
plot(x_axis,PSD(:,4),'lineWidth', 0.8)
title('4th Segment Periodogram', fontSize=12)
xlabel('Normalized Frequency', fontSize=12)
ylabel('Magnitude', fontSize=12)
grid minor

subplot(2,4,5)
plot(x_axis,PSD(:,5),'lineWidth', 0.8)
title('5th Segment Periodogram', fontSize=12)
xlabel('Normalized Frequency', fontSize=12)
ylabel('Magnitude', fontSize=12)
grid minor

subplot(2,4,6)
plot(x_axis,PSD(:,6),'lineWidth', 0.8)
title('6th Segment Periodogram', fontSize=12)
xlabel('Normalized Frequency', fontSize=12)
ylabel('Magnitude', fontSize=12)
grid minor

subplot(2,4,7)
plot(x_axis,PSD(:,7),'lineWidth', 0.8)
title('7th Segment Periodogram', fontSize=12)
xlabel('Normalized Frequency', fontSize=12)
ylabel('Magnitude', fontSize=12)
grid minor

subplot(2,4,8)
plot(x_axis,PSD(:,8),'lineWidth', 0.8)
title('8th Segment Periodogram', fontSize=12)
xlabel('Normalized Frequency', fontSize=12)
ylabel('Magnitude', fontSize=12)
grid minor

%% 3.1.3

PSD_mean = mean(PSD,2);
plot(x_axis,PSD_mean,'lineWidth', 0.8)
title('Averaged Periodogram from 8 Segments of 1024 WGN samples', fontSize=12)
xlabel('Normalized Frequency', fontSize=14)
ylabel('Magnitude', fontSize=14)
grid minor

%% 3.2

clc
clear
close all

x = randn(1064,1);
b = 1;
a = [1 0.9];
y = filter(b, a, x);
y = y(41:1064);

% plot
subplot(1,2,1)
plot(x)
axis([0 1064 -4 4])
title('WGN sequence x', fontSize=15)
xlabel('Time', fontSize=15)
ylabel('Magnitude', fontSize=15)
grid minor

subplot(1,2,2)
plot(y)
axis([0 1024 -8 8])
title('Filtered sequence y', fontSize=15)
xlabel('Time', fontSize=15)
ylabel('Magnitude', fontSize=15)
grid minor

%%
%3.2.1
subplot(1,3,1)
[h,w]=freqz([1],[1 0.9],512);
plot(w/(2*pi),abs(h).^2, 'lineWidth',0.8);
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15);
title('PSD of y', fontSize=15)
grid minor

%3.2.2
pgm_y = pgm(y);

x_axis = [1:1024];
x_axis = x_axis./1024;

subplot(1,3,2)
plot(w/(2*pi),abs(h).^2, 'lineWidth',0.8 , 'DisplayName','PSD');
hold on
plot(x_axis, pgm_y, 'lineWidth',0.8, 'DisplayName','Periodogram')
axis([0 0.5 0 300])
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15);
title('PSD and Periodogram of y', fontSize=15)
grid minor
legend show

%3.2.3
subplot(1,3,3)
plot(w/(2*pi),abs(h).^2, 'lineWidth',0.8 , 'DisplayName','PSD');
hold on
plot(x_axis, pgm_y, 'lineWidth',0.8, 'DisplayName','Periodogram')
axis([0.4 0.5 0 300])
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15);
title('PSD and Periodogram of y', fontSize=15)
grid minor
legend show
zoom xon

%% 3.2.4

acf = xcorr(y, 1, 'biased');
a1 = -acf(3)/acf(2);
var = acf(2) + a1 * acf(3);
sd = sqrt(var);
[h1, w1] = freqz(var^(1/2), [1 a1], 512);

% plots
subplot(1,2,1)
plot(x_axis, pgm_y, 'lineWidth',0.8, 'DisplayName','Periodogram');
hold on
plot(w1/(2*pi),abs(h1).^2, 'lineWidth',0.8 , 'DisplayName','Model-based PSD estimate')
hold on
plot(w/(2*pi),abs(h).^2, 'lineWidth',0.8 , 'DisplayName','PSD');
hold on
axis([0 0.5 0 280])
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15);
title('PSD, Periodogram & Model-based PSD estimate of y', fontSize=13)
grid minor
legend show

subplot(1,2,2)
plot(x_axis, pgm_y, 'lineWidth',0.8, 'DisplayName','Periodogram');
hold on
plot(w1/(2*pi),abs(h1).^2, 'lineWidth',1 , 'DisplayName','Model-based PSD estimate')
hold on
plot(w/(2*pi),abs(h).^2, 'lineWidth',1 , 'DisplayName','PSD');
hold on
axis([0.4 0.5 0 280])
xlabel('Normalised Frequency', fontSize=15);
ylabel('Magnitude', fontSize=15);
title('PSD, Periodogram & Model-based PSD estimate of y', fontSize=13)
grid minor
legend show

%%
%3.2.5

clear all
close all
clc
load sunspot.dat

p = [1 2 10];

sun(:, 1) = sunspot(:, 2);
sun(:, 2) = zscore(sun(:, 1)); % zero mean
N = length(sun(:, 1));

h_orders = zeros(N, length(p));

leg_str = cell(1, 1+length(p));
leg_str{1, 1} = char('Periodogram');

legend_str_n = cell(1, 1+length(p));
legend_str_n{1, 1} = char('Periodogram');


for i = 1:length(p)
    [a, var] = aryule(sun(:, 1), p(i)); % Use sun(:,2) for normalised data

    [h(:, i), w] = freqz(var^(1/2), a, length(sun(:, 1))/2); % Use sun(:,2) for normalised data
    
    leg_str{1, i+1} = char(sprintf('Model Order %d', p(i))); 
end

for i = 1:length(p)
    [a_n, var_n] = aryule(sun(:, 2), p(i)); % Use sun(:,2) for normalised data

    [h_n(:, i), w_n] = freqz(var_n^(1/2), a_n, length(sun(:, 2))/2); % Use sun(:,2) for normalised data
    
    legend_str_n{1, i+1} = char(sprintf('Model Order %d', p(i))); 
end

Ps = pgm(sun(:, 1)); % Use sun(:,2) for normalised data
Ps_n = pgm(sun(:, 2));
f = linspace(0, 0.5, N/2+1);

% plots
subplot(1,3,1)
hold on
plot(f, Ps(1:N/2+1), 'lineWidth', 0.8);
for i = 1:length(p)
    plot(w/(2*pi), abs(h(:, i)).^2, 'lineWidth', 0.8);
end
title('PSD Estimate of Sunspot Data', fontSize=13);
xlabel('Normalized Frequencies', fontSize=15);
ylabel('Magnitude', fontSize=15);
legend(leg_str);
grid on
grid minor

subplot(1,3,2)
hold on
plot(f, Ps_n(1:N/2+1), 'lineWidth', 0.8);
for i = 1:length(p)
    plot(w_n/(2*pi), abs(h_n(:, i)).^2, 'lineWidth', 0.8);
end
title('PSD Estimate of zero-mean Data', fontSize=13);
xlabel('Normalized Frequencies', fontSize=15);
ylabel('Magnitude', fontSize=15);
legend(leg_str);
grid on
grid minor

subplot(1,3,3)
hold on
plot(f, Ps_n(1:N/2+1), 'lineWidth', 0.8);
for i = 1:length(p)
    plot(w_n/(2*pi), abs(h_n(:, i)).^2, 'lineWidth', 0.8);
end
axis([0 0.2 0 45])
title('PSD Estimates of zero-mean Data', fontSize=13);
xlabel('Normalized Frequencies', fontSize=15);
ylabel('Magnitude', fontSize=15);
legend(leg_str);
grid on
grid minor

%% 3.3.3 & 3.3.4

close all;
clear;
clc;

load sunspot.dat
x=(sunspot(:,2)); %input data
x=detrend(x);
[rxx,lags]= xcorr(x,'biased');
rxx_pos= rxx(round(length(rxx)/2)+1: length(rxx));

N=length(x);
M=height(rxx);
start=round(M/2);
a_p=zeros(2,1);
E_p=zeros(1,10);
a_ls=zeros(10,10);

for p=1:10
    H=zeros(start-1,p);
    for ind=1:p
        h=rxx(start-(ind-1): M-(ind));
        H(:,ind)=h;
    end
a_p= inv(transpose(H)*H)* transpose(H)*rxx_pos;
a_ls(p,1:p)= a_p';
est= H *a_ls(p,1:p)'; %approximated signal
E_p(p)= (1/(start-1))* (rxx_pos-est)'*(rxx_pos-est);
end

p = (1:10);
MDL = log(E_p) + p.*log(N)./N;
AIC = log(E_p) + 2*p/N;
AICc = AIC + 2.*p.*(p+1)./(N - p - 1);

hold on;
plot(p, MDL, 'LineWidth', 0.8,'DisplayName','MDL');
plot(p, AIC, 'LineWidth', 0.8,'DisplayName','AIC');
plot(p, AICc, 'LineWidth', 0.8,'DisplayName','AICc');
grid minor
xticks(1:1:10);
hold off; legend('show', fontSize=11);
xlabel('Model Order p', fontSize=15); ylabel('Criteria Value', fontSize=15);
title('AR Model Order Selection using LSE', fontSize=15);

%%
%3.3.4
%power spectra
a_p = [ones(10,1) -a_ls];
[a1, e1, r1] = aryule(x, 1);
[a3,e3,r3]= aryule(x,3);
[a10, e10, r10] = aryule(x, 10);

[h1,w1]=freqz(sqrt(e1),a_p(1, 1:2));
[h3,w3]= freqz(sqrt(e3),a_p(3, 1:4));
[h10,w10]=freqz(sqrt(e10),a_p(10, 1:11));
P_xx = pgm(x);

hold on;
x_axis = [1:288];
x_axis = x_axis./288;
plot(w1/(2*pi), abs(h1).^2, 'linewidth', 0.8);
plot(w3/(2*pi), abs(h3).^2, 'linewidth', 0.8);
plot(w10/(2*pi), abs(h10).^2, 'linewidth', 0.8);
plot(x_axis, P_xx,'linewidth', 0.8);
grid minor
hold off;

xlabel('Normalised Frequency',fontSize=15);
ylabel('Magnitude',fontSize=15); title('Power Spectra of AR models for Sunspot data',fontSize=15);
xlim([0 0.5])
legend('AR(1)','AR(3)','AR(10)', 'Periodogram', fontSize=11);

%% 3.3.6
close all;
clear;
clc;

load sunspot.dat
sun = sunspot(:,2);
sun= zscore(sun);
p = 3; % optimal order
MSE=zeros(1,49);
J=zeros(1,49);

for N = 10:5:250
    sun2 = sun(1:N);
    [rxx,lags]= xcorr(sun2, 'biased');
    rxx_pos= rxx(round(length(rxx)/2)+1: length(rxx));
    M=height(rxx);
    start=round(M/2);
    H=zeros(start-1,p);
    for ind=1:p
        h=rxx(start-(ind-1): M-(ind));
        H(:,ind)=h;
    end
    a_p= inv(transpose(H)*H)* transpose(H)*rxx_pos;
    a_ls(p,1:p)= a_p';
    est= H *a_ls(p,1:p)'; %estimate
    J_p=(rxx_pos-est)'*(rxx_pos-est); % cost function
    E_p = (1/(start-1))* (rxx_pos-est)'*(rxx_pos-est); % MSE
    MSE(49-((250-N)/5))=E_p;
    J(49-((250-N)/5))=J_p;
end

x_axis = linspace(10,250,49);

subplot(1,2,1)
plot(x_axis, J, 'lineWidth', 1)
grid minor
title('Cost function against N for AR(3) model', fontSize=15);
xlabel('Data Length N', fontSize=15);
ylabel('Magnitude', fontSize=15);
axis([10 250 0.005 0.04])
xticks([10 50 100 150 200 250])

subplot(1,2,2)
plot(x_axis, MSE, 'lineWidth', 1)
grid minor
title('MSE against N for AR(3) model', fontSize=15);
xlabel('Data Length N', fontSize=15);
ylabel('Magnitude', fontSize=15);
axis([10 250 0 0.003])
xticks([10 50 100 150 200 250])

%% 3.4.1 Dial-tone pad

close all;
clear;
clc;

num = randperm(9,8);
num = [0,2,0,num];

% Dial pad frequencies
tones = [1336, 1209, 1336, 1447, 1209, 1336, 1447, 1209, 1336, 1447;
         941,  697,  697,  697,  770,  770,  770,  852,  852,  852];
x = linspace(0, 0.25, 8192); % gives sampling frequency of f 32768 Hz
dial_seq = []; % output sequence

% generate sequence of 10 digits and 10 gaps
for i=1:10
    dig = sin(2*pi*tones(1, num(i)+1) * x) + sin(2 * pi * tones(2, num(i)+1) * x); 
    gap = zeros(1, 8192);
    dial_seq = [dial_seq,dig,gap];
end

% final digit separate since no gap after
dig = sin(2*pi*tones(1, num(11)+1) * x) + sin(2 * pi * tones(2, num(11)+1) * x); 
dial_seq = [dial_seq,dig];

% plot
subplot(2,1,1)
x_axis = linspace(0,5.25,172032);
plot(x_axis, dial_seq, 'lineWidth', 0.8);
axis([0 5.25 -2 2])
grid minor
xlabel('time(s)', fontSize=15);
ylabel('Amplitude', fontSize=15);
title('Dial Tone signal y[n]', fontSize=15);

subplot(2,1,2)
x_axis = linspace(0,5.25,172032);
plot(x_axis, dial_seq, 'lineWidth', 0.8);
axis([0 0.75 -2 2])
grid minor
xlabel('time(s)', fontSize=15);
ylabel('Amplitude', fontSize=15);
title('Dial Tone signal y[n], difference between 2 tones', fontSize=15);

%% 3.4.2

[s,f,t] = spectrogram(dial_seq, hann(8192), 0, 8192, 32768);

figure(1);hold on;
spectrogram(dial_seq,hann(8192), 0, 8192, 32768, 'yaxis');
ylim([0.2,1.9]);xlim([0,5.25])
title("Spectrogram of dial tone signal y[n]", fontSize=15);
xlabel('time(s)', fontSize=15);
ylabel('Frequency(kHz)', fontSize=15)
c = colorbar;
c.Label.String = 'Power/frequency(dB/Hz)';
c.Label.FontSize = 15;

figure(2)
hold on;
plot(f, mag2db(abs(s(:,1))), 'lineWidth', 0.8); 
plot(f, mag2db(abs(s(:,3))),'lineWidth', 0.8);
title("Magnitude spectrum of FFT segments for 0 and 2", fontSize=15)
xlabel('Frequency(Hz)', fontSize=15);
ylabel('Magnitude(dB)', fontSize=15);
legend('Number 0','Number 2', fontSize=11);
axis([0 2000 -100 80])
grid minor

%% 3.4.3

sd_noise1 = 0.1;
sd_noise2 = 1;
sd_noise3 = 10;
wgn1 = randn(1,length(dial_seq))*sd_noise1;
wgn2 = randn(1,length(dial_seq))*sd_noise2;
wgn3 = randn(1,length(dial_seq))*sd_noise3;
dial_seq1 = dial_seq + wgn1;
dial_seq2 = dial_seq + wgn2;
dial_seq3 = dial_seq + wgn3;

subplot(2,3,1)
plot(5.25/length(dial_seq1):5.25/length(dial_seq1):5.25,dial_seq1-5.25/length(dial_seq1));
xlim([0,5.25]);
grid on;
title("Dial Tone signal y[n] with Noise, var=0.1 ", fontSize=12);
xlabel("time(s)", fontSize=13);
ylabel("Magnitude", fontSize=13);
grid minor;

subplot(2,3,2)
plot(5.25/length(dial_seq2):5.25/length(dial_seq2):5.25,dial_seq2-5.25/length(dial_seq2));
xlim([0,5.25]);
grid on;
title("Dial Tone signal y[n] with Noise, var=1 ", fontSize=12);
xlabel("time(s)", fontSize=13);
ylabel("Magnitude", fontSize=13);
grid minor;

subplot(2,3,3)
plot(5.25/length(dial_seq3):5.25/length(dial_seq3):5.25,dial_seq3-5.25/length(dial_seq3));
xlim([0,5.25]);
grid on;
title("Dial Tone signal y[n] with Noise, var=10 ", fontSize=12);
xlabel("time(s)", fontSize=13);
ylabel("Magnitude", fontSize=13);
grid minor;

subplot(2,3,4);hold on;
spectrogram(dial_seq1,hann(8192), 0, 8192, 32768, 'yaxis');
ylim([0.2,1.9]);xlim([0,5.25])
title("Spectrogram of y[n] with Noise, var=0.1", fontSize=12);
xlabel('time(s)', fontSize=13);
ylabel('Frequency(kHz)', fontSize=13)
c = colorbar;
c.Label.String = 'Power/frequency(dB/Hz)';
c.Label.FontSize = 13;

subplot(2,3,5);hold on;
spectrogram(dial_seq2,hann(8192), 0, 8192, 32768, 'yaxis');
ylim([0.2,1.9]);xlim([0,5.25])
title("Spectrogram of y[n] with Noise, var=1", fontSize=12);
xlabel('time(s)', fontSize=13);
ylabel('Frequency(kHz)', fontSize=13)
c = colorbar;
c.Label.String = 'Power/frequency(dB/Hz)';
c.Label.FontSize = 13;

subplot(2,3,6);hold on;
spectrogram(dial_seq3,hann(8192), 0, 8192, 32768, 'yaxis');
ylim([0.2,1.9]);xlim([0,5.25])
title("Spectrogram of y[n] with Noise, var=10", fontSize=12);
xlabel('time(s)', fontSize=13);
ylabel('Frequency(kHz)', fontSize=13)
c = colorbar;
c.Label.String = 'Power/frequency(dB/Hz)';
c.Label.FontSize = 13;


%% 3.5
fs=4; %sampling frequency of given data
c50=fs*50;
c150=fs*100;

% trial 1 ------------------------------------------------------------
xRRI= RRI_trial1; 
xRRI= detrend(xRRI);
Px=pgm(xRRI); %periodogram
freqs=linspace(0,1,length(Px));

figure(1)
subplot(1,3,1);
plot(freqs(1:length(freqs)/2), Px(1:length(Px)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Periodogram of RRI 1', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);

% window length 50
win_sz50=floor(length(xRRI)/c50); % number of segments
pxx50=zeros(1,c50); % holds the periodogram
for ind=0:win_sz50-1 % going through each segment
    x_50= xRRI(ind*c50+1:(ind+1)*c50); % extracts each segment
    pxx50=pxx50+ pgm(x_50); % accumulate all the segments
end
pxx_avg50= pxx50./win_sz50; % averaged
freqs=linspace(0,1,length(pxx_avg50));
subplot(1,3,2);
plot(freqs(1:length(freqs)/2), pxx_avg50(1:length(pxx_avg50)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Averaged periodogram of RRI 1, 50s Window', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);


win_sz150=floor(length(xRRI)/c150);
pxx150=zeros(1,c150);
for ind=0:win_sz150-1
    x_avg_150= xRRI(ind*c150+1:(ind+1)*c150);
    pxx150=pxx150+ pgm(x_avg_150);
end
pxx_avg150= pxx150./win_sz150;
freqs=linspace(0,1,length(pxx_avg150));
subplot(1,3,3);
plot(freqs(1:length(freqs)/2), pxx_avg150(1:length(pxx_avg150)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Averaged periodogram of RRI 1, 100s Window', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);

% trial 2 ------------------------------------------------------------
xRRI= RRI_trial2; 
xRRI= detrend(xRRI);
Px=pgm(xRRI); %periodogram
freqs=linspace(0,1,length(Px));

figure(2)
subplot(1,3,1);
plot(freqs(1:length(freqs)/2), Px(1:length(Px)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Periodogram of RRI 2', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);

% window length 50
win_sz50=floor(length(xRRI)/c50); % number of segments
pxx50=zeros(1,c50); % holds the periodogram
for ind=0:win_sz50-1 % going through each segment
    x_50= xRRI(ind*c50+1:(ind+1)*c50); % extracts each segment
    pxx50=pxx50+ pgm(x_50); % accumulate all the segments
end
pxx_avg50= pxx50./win_sz50; % averaged
freqs=linspace(0,1,length(pxx_avg50));
subplot(1,3,2);
plot(freqs(1:length(freqs)/2), pxx_avg50(1:length(pxx_avg50)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Averaged periodogram of RRI 2, 50s Window', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);


win_sz150=floor(length(xRRI)/c150);
pxx150=zeros(1,c150);
for ind=0:win_sz150-1
    x_avg_150= xRRI(ind*c150+1:(ind+1)*c150);
    pxx150=pxx150+ pgm(x_avg_150);
end
pxx_avg150= pxx150./win_sz150;
freqs=linspace(0,1,length(pxx_avg150));
subplot(1,3,3);
plot(freqs(1:length(freqs)/2), pxx_avg150(1:length(pxx_avg150)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Averaged periodogram of RRI 2, 100s Window', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);

% trial 3 ------------------------------------------------------------
xRRI= RRI_trial3; 
xRRI= detrend(xRRI);
Px=pgm(xRRI); %periodogram
freqs=linspace(0,1,length(Px));

figure(3)
subplot(1,3,1);
plot(freqs(1:length(freqs)/2), Px(1:length(Px)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Periodogram of RRI 3', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);

% window length 50
win_sz50=floor(length(xRRI)/c50); % number of segments
pxx50=zeros(1,c50); % holds the periodogram
for ind=0:win_sz50-1 % going through each segment
    x_50= xRRI(ind*c50+1:(ind+1)*c50); % extracts each segment
    pxx50=pxx50+ pgm(x_50); % accumulate all the segments
end
pxx_avg50= pxx50./win_sz50; % averaged
freqs=linspace(0,1,length(pxx_avg50));
subplot(1,3,2);
plot(freqs(1:length(freqs)/2), pxx_avg50(1:length(pxx_avg50)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Averaged periodogram of RRI 3, 50s Window', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);


win_sz150=floor(length(xRRI)/c150);
pxx150=zeros(1,c150);
for ind=0:win_sz150-1
    x_avg_150= xRRI(ind*c150+1:(ind+1)*c150);
    pxx150=pxx150+ pgm(x_avg_150);
end
pxx_avg150= pxx150./win_sz150;
freqs=linspace(0,1,length(pxx_avg150));
subplot(1,3,3);
plot(freqs(1:length(freqs)/2), pxx_avg150(1:length(pxx_avg150)/2), 'lineWidth', 0.8);
xlim([0 0.5]);
grid minor
title('Averaged periodogram of RRI 3, 100s Window', fontSize=12);
xlabel('Normalised Frequency', fontSize=15); ylabel('Magnitude', fontSize=15);