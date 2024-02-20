%% Classical and Modern Spectrum Estimation
% Question 1.1

clear all;
close all;
clc;

%% Initialisations

Fs = 100; % sampling frequency
T = 1/Fs; % sampling period
N = 500; % signal length

% plotting axis in time domain
t = (0 : N-1) * T; % signal vector (in seconds, starts from 0)

%% Signal generation

% pulse signal
pulse = zeros(size(t));
pulse(N/2) = 2;
pulseACF = xcorr(pulse, 'biased'); % pulse ACF

% sine signal
sine = sin(2*pi*0.5*t);
sineACF = xcorr(sine, 'biased'); % sine ACF

acfAxis = (-N+1) : (N-1);  

%% PSD estimation

pulsePSD1 = abs(fftshift(fft(pulseACF)));
pulsePSD2 = (abs(fftshift(fft(pulse))).^2)/N;

sinePSD1 = abs(fftshift(fft(sineACF)));
sinePSD2 = (abs(fftshift(fft(sine))).^2)/N;

% plotting axis in frequency domain
f1 = (-Fs/2) : Fs/length(pulseACF) : (Fs/2 - Fs/length(pulseACF));
f2 = (-Fs/2) : Fs/length(pulsePSD2) : (Fs/2 - Fs/length(pulsePSD2));

%% Plot generation

figure(1); 

% pulse wave
subplot(3, 2, 1);
plot(t, pulse, 'LineWidth', 1);
set(gca,'fontsize', 12);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
title('Pulse Signal', 'fontsize', 12);
grid on;
grid minor;

% sine wave
subplot(3, 2, 2);
plot(t, sine, 'LineWidth', 1);
set(gca,'fontsize', 12);
xlabel('Time (s)', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
title('Sine Signal', 'fontsize', 12);
grid on;
grid minor;

% pulse ACF
subplot(3, 2, 3);
plot(acfAxis, pulseACF, 'LineWidth', 1);
set(gca,'fontsize', 12);
xlabel('Time Lag k', 'fontsize', 12);
ylabel('Correlation', 'fontsize', 12);
title('Pulse ACF (Fast-decaying)', 'fontsize', 12);
grid on;
grid minor;

% sine ACF
subplot(3, 2, 4);
plot(acfAxis, sineACF, 'LineWidth', 1);
set(gca,'fontsize', 12);
xlabel('Time Lag k', 'fontsize', 12);
ylabel('Correlation', 'fontsize', 12);
title('Sine ACF (Slow-decaying)', 'fontsize', 12);
grid on;
grid minor;

% pulse PSD
subplot(3, 2, 5);
plot(f1, pulsePSD1, 'r', 'LineWidth', 1);
hold on;
plot(f2, pulsePSD2, '--b', 'LineWidth', 1.6);
hold off;
legend('Def.1', 'Def.2');
set(gca,'fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
title('PSD Estimations for Pulse Signal', 'fontsize', 12);

grid on;
grid minor;

% sine PSD
subplot(3, 2, 6);
plot(f1, sinePSD1, 'r', 'LineWidth', 1);
hold on;
plot(f2, sinePSD2, '--b', 'LineWidth', 1.6);
hold off;
legend('Def.1', 'Def.2');
set(gca,'fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
title('PSD Estimations for Sine Signal', 'fontsize', 12);
grid on;
grid minor;

