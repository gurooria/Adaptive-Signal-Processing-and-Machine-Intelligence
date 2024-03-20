%% 1.3a biased and unbiased ACF estimates

clc;
clear;
close all;

% Initialisations
Fs = 20; % sampling frequency
T = 1/Fs; % sampling period
N = 1024; % signal length
t = (0 : N-1) * T; % signal vector (in seconds, starts from 0)

%% Signal Generation

WGN = randn(1, N);
filtered_WGN = filter([1/4 1/4 1/4 1/4], 1, WGN);
noisy_sine = sin(2*pi*0.5*t) + WGN;

%% Calculate biased + unbiased ACFs

% biased
[WGN_biased, WGN_lag1] = xcorr(WGN, 'biased');
[sine_biased, sine_lag1] = xcorr(noisy_sine, 'biased');
[filtered_biased, filtered_lag1] = xcorr(filtered_WGN, 'biased');

% unbiased
[WGN_unbiased, WGN_lag2] = xcorr(WGN, 'unbiased');
[sine_unbiased, sine_lag2] = xcorr(noisy_sine, 'unbiased');
[filtered_unbiased, filtered_lag2] = xcorr(filtered_WGN, 'unbiased');

%% Calculate correlogram based on ACFs
PSD_wgn_biased = fftshift(real(fft(ifftshift(WGN_biased))));
PSD_sine_biased = fftshift(real(fft(ifftshift(sine_biased))));
PSD_filtered_biased = fftshift(real(fft(ifftshift(filtered_biased))));

PSD_wgn_unbiased = fftshift(real(fft(ifftshift(WGN_unbiased))));
PSD_sine_unbiased = fftshift(real(fft(ifftshift(sine_unbiased))));
PSD_filtered_unbiased = fftshift(real(fft(ifftshift(filtered_unbiased))));

%% Plotting

% WGN ACF
subplot(3, 2, 1)
plot(WGN_lag2, WGN_unbiased, 'LineWidth', 1.2);
hold on
plot(WGN_lag1, WGN_biased, 'LineWidth', 1.2);
hold off
xlabel('Time Lag k', 'fontsize',12);
ylabel('Correlation', 'fontsize',12);
title('ACF of unbiased and biased WGN signals', 'fontsize',12)
legend('Unbiased', 'Biased');
grid on
grid minor
set(gca,'fontsize', 12);
xlim([WGN_lag1(1), WGN_lag1(end)]);

% sine ACF
subplot(3, 2, 3)
plot(sine_lag2, sine_unbiased, 'LineWidth', 1.2);
hold on
plot(sine_lag1, sine_biased, 'LineWidth', 1.2);
hold off
xlabel('Time Lag k', 'fontsize',12);
ylabel('Correlation', 'fontsize',12);
title('ACF of unbiased and biased noisy sine signals', 'fontsize',12)
legend('Unbiased', 'Biased');
grid on
grid minor
set(gca,'fontsize', 12);
xlim([sine_lag1(1), sine_lag1(end)]);

% filtered WGN ACF
subplot(3, 2, 5)
plot(filtered_lag2, filtered_unbiased, 'LineWidth', 1.2);
hold on
plot(filtered_lag1, filtered_biased, 'LineWidth', 1.2);
hold off
xlabel('Time Lag k', 'fontsize',12);
ylabel('Correlation', 'fontsize',12);
title('ACF of unbiased and biased filtered WGN signals', 'fontsize',12)
legend('Unbiased', 'Biased');
grid on
grid minor
set(gca,'fontsize', 12);
xlim([filtered_lag1(1), filtered_lag1(end)]);

% PSD axis
norm_f = linspace(-1,1, length(WGN_lag1));

% WGN PSD
subplot(3, 2, 2)
plot(norm_f, PSD_wgn_unbiased, 'LineWidth', 1.2);
hold on
plot(norm_f, PSD_wgn_biased, 'LineWidth', 1.2);
hold off
xlabel('Normalised Frequency (Cycles/Sample)', 'fontsize',12);
ylabel('Magnitude', 'fontsize',12);
title('Correlograms of unbiased and biased WGN signals', 'fontsize',12)
legend('Unbiased', 'Biased');
grid on
grid minor
set(gca,'fontsize', 12);

% sine PSD
subplot(3, 2, 4)
plot(norm_f, PSD_sine_unbiased, 'LineWidth', 1.2);
hold on
plot(norm_f, PSD_sine_biased, 'LineWidth', 1.2);
hold off
xlabel('Normalised Frequency (Cycles/Sample)', 'fontsize',12);
ylabel('Magnitude', 'fontsize',12);
title('Correlograms of unbiased and biased noisy sine signals', 'fontsize',12)
legend('Unbiased', 'Biased');
grid on
grid minor
set(gca,'fontsize', 12);

% filtered WGN PSD
subplot(3, 2, 6)
plot(norm_f, PSD_filtered_unbiased, 'LineWidth', 1.2);
hold on
plot(norm_f, PSD_filtered_biased, 'LineWidth', 1.2);
hold off
xlabel('Normalised Frequency (Cycles/Sample)', 'fontsize',12);
ylabel('Magnitude', 'fontsize',12);
title('Correlograms of unbiased and biased filtered WGN signals', 'fontsize',12)
legend('Unbiased', 'Biased');
grid on
grid minor
set(gca,'fontsize', 12);



