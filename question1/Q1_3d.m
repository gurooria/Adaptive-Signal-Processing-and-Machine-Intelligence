%% 1.3d: Generating Complex Exponential Signals

clc;
clear;
close all;

% Initialisations
N = 1024; % signal length
n1 = [20 35 40]; % signal length vector
n2 = [45 50 100];
colors = {'b','r','c'};

%%
figure(1)
hold on

subplot(1,2,1)
for i = 1: length(n1) % go through each signal length
    % given by coursework for complex exponentials
    n = 0 : n1(i);
    noise = 0.2/sqrt(2) * (randn(size(n)) + 1j * randn(size(n)));
    noisySignal = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    zeroPadding = zeros(1, N-length(noisySignal));
    psd = abs((fft([noisySignal zeroPadding]))./length(n));
    x = (0:N-1)/N;
    plot(x, 10*log10(psd),'color', colors{i}, 'LineWidth', 1.2)
    hold on
end
grid on
grid minor
ax = gca;
ax.FontSize = 12;
legend('N = 20','N = 35','N = 40')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Periodograms of Noisy Complex Exponential Signals', 'fontsize', 12)
xlim([0.1 0.5])

subplot(1,2,2)
for i = 1: length(n2) % go through each signal length
    % given by coursework for complex exponentials
    n = 0 : n2(i);
    noise = 0.2/sqrt(2) * (randn(size(n)) + 1j * randn(size(n)));
    noisySignal = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    zeroPadding = zeros(1, N-length(noisySignal));
    psd = abs((fft([noisySignal zeroPadding]))./length(n));
    x = (0:N-1)/N;
    plot(x, 10*log10(psd),'color', colors{i}, 'LineWidth', 1.2)
    hold on
end
grid on
grid minor
ax = gca;
ax.FontSize = 12;
legend('N = 45','N = 50', 'N = 100')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Periodograms of Noisy Complex Exponential Signals', 'fontsize', 12)
xlim([0.1 0.5])
