%% Q1.2a
clear;
close all;
clc;

load sunspot.dat

%% Pre-processing
data = sunspot(:, 2);

data_mean = data - mean(data); % remove mean
data_detrend = detrend(data); % detrend data
data_log = log(data+eps); % log data (with floating pt relative accuracy)
data_log_mean = data_log - mean(data_log); % log data mean removed

%% PSD computation
[psd_data, x1] = periodogram(data, hamming(length(data)), [], 1); % x1 is the axis
[psd_mean, x2] = periodogram(data_mean, hamming(length(data_mean)), [], 1);
[psd_detrend, x3] = periodogram(data_detrend, hamming(length(data_detrend)), [], 1);
[psd_log_mean, x4] = periodogram(data_log_mean, hamming(length(data_log_mean)), [], 1);

%% Graph Plotting
figure(1)

% Time domain plotting
subplot(2, 1, 1);

hold on
grid on
grid minor

plot(1:length(data), data, 'linewidth', 1.2);
plot(1:length(data_mean), data_mean, 'linewidth', 1.2);
plot(1:length(data_detrend), data_detrend, 'linewidth', 1.2);
plot(1:length(data_log_mean), data_log_mean, 'linewidth', 1.2);

hold off
title('Sunspot Time Series with different Pre-Processing', 'fontsize', 12);
set(gca,'fontsize', 12);
legend('Raw', 'Mean Removed', 'Detrend', 'Log + Mean Removed', 'fontsize', 10);
xlabel('Sample Number (Year)', 'fontsize', 12);
ylabel('Sunspot Number', 'fontsize', 12);


% PSD plotting
subplot(2, 1, 2);

hold on
grid on
grid minor

plot(x1, 10*log10(psd_data), 'linewidth', 1.2);
plot(x2, 10*log10(psd_mean), 'linewidth', 1.2);
plot(x3, 10*log10(psd_detrend), 'linewidth', 1.2);
plot(x4, 10*log10(psd_log_mean), 'linewidth', 1.2);

hold off
title('Periodograms of Sunspot Time Series', 'fontsize', 12);
set(gca,'fontsize', 12);
legend('Raw', 'Mean Removed', 'Detrend', 'Log + Mean Removed', 'fontsize', 10);
xlabel('Normalised Frequency (Cycles/Sample)', 'fontsize', 12);
ylabel('Magnitude (dB)', 'fontsize', 12);


%% Q1.2b

clear;
close all;
clc;

load('EEG_Data_Assignment1.mat');

%% Initialisations

sample_length = length(POz); % EEG data
nfft = 5 * fs; % Number of DFT points, recommended to have 5 DFT samples per Hz
POz = POz - mean(POz); % mean removal of data

% fs = sampling frequency -> number of samples per second

% data reshaping for different window lenghts
POz1 = reshape(POz, fs*1, []); % 1s window -> 1200 samples, 80 segments
POz5 = reshape(POz, fs*5, []); % 5s window -> 6000 samples, 16 segments
POz10 = reshape(POz, fs*10, []); % 10s window -> 12000 samples, 8 segments

%% PSD Computation

[psd_standard,x_standard] = periodogram(POz, hamming(sample_length), nfft, fs); 
psd_1 = mean(periodogram(POz1, hamming(length(POz1)), nfft, fs), 2);
psd_5 = mean(periodogram(POz5, hamming(length(POz5)), nfft, fs), 2);
psd_10 = mean(periodogram(POz10, hamming(length(POz10)), nfft, fs), 2);

%% PSD plotting
figure(2)

subplot(2,2,1)
plot(x_standard, 10*log10(psd_standard), 'linewidth', 1.2);
xlim([0, 60]);
set(gca,'fontsize', 12);
title('Standard Periodogram of EEG Recording', 'fontsize', 12);
ylabel('Magnitude (dB)', 'fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12);
grid on
grid minor


subplot(2,2,2)
plot(x_standard, 10*log10(psd_1), 'linewidth', 1.2);
xlim([0, 60]);
set(gca,'fontsize', 12);
title('Averaged Periodogram of EEG with 1s Windowing', 'fontsize', 12);
ylabel('Magnitude (dB)', 'fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12);
grid on
grid minor


subplot(2,2,3)
plot(x_standard, 10*log10(psd_5), 'linewidth', 1.2);
xlim([0, 60]);
set(gca,'fontsize', 12);
title('Averaged Periodogram of EEG with 5s Windowing', 'fontsize', 12);
ylabel('Magnitude (dB)', 'fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12);
grid on
grid minor


subplot(2,2,4)
plot(x_standard, 10*log10(psd_10), 'linewidth', 1.2);
xlim([0, 60]);
set(gca,'fontsize', 12);
title('Averaged Periodogram of EEG with 10s Windowing', 'fontsize', 12);
ylabel('Magnitude (dB)', 'fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12);
grid on
grid minor