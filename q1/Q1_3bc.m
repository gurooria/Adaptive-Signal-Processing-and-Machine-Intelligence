%% Q1.3b & Q1.3c
clc;
clear;
close all;

%% Initialisations
Fs = 10; % sampling frequency
N = 1024; % signal length
t = 0 : 1/(2*Fs) : Fs; % sample vector
sine = [0.8*sin(2*pi*t*0.5) + 1.5*sin(2*pi*t*1.2) + 1*sin(2*pi*t*2) zeros(1, N-length(t))];
psds = zeros(100, N*2-1);

%% Plotting

% Plot the realisations
figure(1)
hold on
for i = 1:100
    noisySine = sine + wgn(length(sine), 1, 1)';

    % Compute the correlogram using the 'biased' method
    [psd, lags] = xcorr(noisySine, 'biased');
    psds(i,:) = fftshift(real(fft(ifftshift(psd))));
    subplot(2,2,1)
    plot((lags/max(lags))*Fs, real(psds(i,:)), 'color', 'c', 'LineWidth', 1.2);
    hold on
end

% Plot the mean
plot((lags/max(lags))*Fs, mean(real(psds)), 'color', 'b', 'LineWidth',1.2);
xlim([0,3])
ax = gca;
ax.FontSize = 12;
xlabel('Normalised Frequency (Cycles/Sample)')
ylabel('Magnitude')
set(gca,'fontsize', 12)
title('Realisations and Mean of the Noisy Composite Sine Correlograms', 'fontsize', 12)
grid on 
grid minor

% Plot the standard deviation
subplot(2,2,2)
plot((lags/max(lags))*Fs, std(real(psds)),'color', 'b', 'LineWidth', 1.2);
xlim([0,3])
ax = gca;
ax.FontSize = 12;
xlabel('Normalised Frequency (Cycles/Sample)')
ylabel('Magnitude')
set(gca, 'fontsize', 12)
title('Standard Deviation of the Noisy Composite Sine Correlograms','fontsize',12)
grid on
grid minor
set(gcf, 'color', 'w')

% 1.3c

% Plot the realisations in dB
hold on
for i = 1:100
    noisySine = sine + wgn(length(sine), 1, 1)';

    % Compute the correlogram using the 'biased' method
    [psd, lags] = xcorr(noisySine, 'biased');
    psds(i,:) = fftshift(real(fft(ifftshift(psd))));
    subplot(2,2,3)
    plot((lags/max(lags))*Fs,  10*log10(real(psds(i,:))), 'color', 'c', 'LineWidth', 1.2);
    hold on
end

% Plot the mean
plot((lags/max(lags))*Fs, 10*log10(mean(real(psds))), 'color','b', 'LineWidth',1.2);
xlim([0,3])
ax = gca;
ax.FontSize = 12;
xlabel('Normalised Frequency (Cycles/Sample)')
ylabel('Magnitude (dB)')
set(gca,'fontsize', 12)
title('Realisations and Mean of the Noisy Composite Sine Correlograms','fontsize', 12)
grid on 
grid minor

% Plot the standard deviation
subplot(2,2,4)
plot((lags/max(lags))*Fs, 10*log10(std(real(psds))), 'color', 'b', 'LineWidth', 1.2);
xlim([0,3])
ax = gca;
ax.FontSize = 12;
xlabel('Normalised Frequency (Cycles/Sample)')
ylabel('Magnitude (dB)')
set(gca, 'fontsize', 12)
title('Standard Deviation of the Noisy Composite Sine Correlograms','fontsize', 12)
grid on
grid minor
set(gcf, 'color', 'w')