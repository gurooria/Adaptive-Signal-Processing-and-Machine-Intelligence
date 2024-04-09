%% Q2.3d Applying ANC to EEG Data
clc
clear
close all
load("EEG_Data_Assignment2.mat");

%%Initialisations
stepSizes = [0.001, 0.01, 0.1];
orders = [5, 10, 20];
var = 0.005;
noise = sqrt(var)*randn(1,length(Cz));
t = (0 : 1/fs : (1/fs)*(length(Cz)-1)); % fs is given
sine = sin(2*pi*50*t);
ref = sine + noise; % corrupted sine signal as reference
windowLength = 4096;
overlap = 2/3;
nfft = 5 * windowLength;

%% Plot spectrogram for the original noisy signal
figure
subplot(1,2,1)
spectrogram(Cz, hamming(windowLength), round(overlap * windowLength), nfft, fs, 'yaxis');
title('Spectrogram of the Noise-Corrupted EEG Data (Cz)', 'fontsize', 12);
xlabel('Time (minute)', 'fontsize', 12)
ylabel('Frequency (Hz)', 'fontsize', 12)
yticks(0:10:60);
ylim([0 60]);
c = colorbar('fontSize', 12);
c.Label.FontSize = 13;
c.Label.String = "Power (dB)";
set(gca, 'fontSize', 12);
set(gcf, 'color', 'w');

subplot(1,2,2)
M = 10; % optimal parameters as observed from above
mu = 0.001;
%Computing ANC algorithm
[xHatANC, ~, ~] = LMS_ANC(Cz, ref, mu, 0, M);
windowLength = 1024;
nfft = 5 * windowLength;
% noisy Cz
[pNoisy, xNoisy] = pwelch(Cz, hamming(windowLength), round(overlap*windowLength), nfft, fs, 'onesided');
pNoisy = 10*log10(pNoisy); % apply log transformation
% denoised Cz
[pDenoised, xDenoised] = pwelch(xHatANC, hamming(windowLength), round(overlap*windowLength), nfft, fs, 'onesided');
pDenoised = 10*log10(pDenoised);

plot(xNoisy, pNoisy, 'b', 'Linewidth', 1.2)
hold on
plot(xDenoised, pDenoised, 'r', 'Linewidth', 1.2)
title('Periodogram of Noisy vs. Denoised EEG data (Cz)', 'Fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12);
xlim([0 70]);
ylabel('Power (dB)', 'fontsize', 12);
legend('Noisy', 'Denoised');
set(gca, 'fontSize', 12);
grid on
grid minor


%% Vary Learning Rates and Model Orders for ANC
count = 1;
figure
hold on

for i = 1 : length(orders)
    for j = 1 : length(stepSizes)
        [xHatANC, ~, ~] = LMS_ANC(Cz, ref, stepSizes(j), 0, orders(i));
        subplot(3,3,count)
        spectrogram(xHatANC, hamming(windowLength), round(overlap*windowLength), nfft, fs, 'yaxis');
        title(sprintf('M = %0.0f, µ = %0.3f', orders(i), stepSizes(j)), 'fontsize', 12)
        ax = gca;
        ax.FontSize = 12;
        xlabel('Time (minute)', 'fontsize', 12);
        ylabel('Frequency (Hz)', 'fontsize', 13);
        yticks(0:10:60);
        ylim([0 60]);
        c = colorbar();
        c.Label.FontSize = 13;
        c.Label.String = "Power (dB)";
        set(gcf,'color','w');
        count = count + 1;
    end
end

%% Periodograms of noisy vs. denoised EEG

M = 10; % optimal parameters as observed from above
mu = 0.001;

%Computing ANC algorithm
[xHatANC, ~, ~] = LMS_ANC(Cz, ref, mu, 0, M);
windowLength = 1024;
nfft = 5 * windowLength;

% noisy Cz
[pNoisy, xNoisy] = pwelch(Cz, hamming(windowLength), round(overlap*windowLength), nfft, fs, 'onesided');
pNoisy = 10*log10(pNoisy); % apply log transformation

% denoised Cz
[pDenoised, xDenoised] = pwelch(xHatANC, hamming(windowLength), round(overlap*windowLength), nfft, fs, 'onesided');
pDenoised = 10*log10(pDenoised);

figure
plot(xNoisy, pNoisy, 'b', 'Linewidth', 1.2)
hold on
plot(xDenoised, pDenoised, 'r', 'Linewidth', 1.2)
title('Periodogram of Noisy vs. Denoised EEG data (Cz)', 'Fontsize', 12);
xlabel('Frequency (Hz)', 'fontsize', 12);
xlim([0 70]);
ylabel('Power (dB)', 'fontsize', 12);
legend('Noisy', 'Denoised');
set(gca, 'fontSize', 12);
grid on
grid minor
