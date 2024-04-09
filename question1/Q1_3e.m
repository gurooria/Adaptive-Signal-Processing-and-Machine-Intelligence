%% Q1.3e MUSIC
clc;
clear;
close all;

% Initialisations
N = 256;
n = 0 : (30-1);
psds = zeros(100, 256);

%%
figure(1)
hold on

subplot(3,2,1)
p = 1; % order
for i = 1:100
    noise = 0.2/sqrt(2) * (randn(size(n)) + 1j*randn(size(n)));
    noisySignal = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    [~, R] = corrmtx(noisySignal, 14, 'mod');
    [psds(i, :), x] = pmusic(R, p, [ ], 1);
    plot(x, psds(i, :), 'color', 'c', 'LineWidth', 1.2)
    hold on;
end

% mean
plot(x, mean(psds), 'color', 'b', 'LineWidth', 1.2)
ax = gca;
ax.FontSize = 12;
ylabel('Magnitude')
xlabel('Normalised Frequency (Cycles/Sample)')
xlim([0.2 0.4])
title('Realisations and Mean of MUSIC Pseudospectrum, p = 1', 'fontsize', 12)
grid on
grid minor

% standard deviation
subplot(3,2,2)
plot(x, std(psds), 'color','b','LineWidth', 1.2)
ax = gca;
ax.FontSize = 12;
ylabel('Magnitude')
xlabel('Normalised Frequency (Cycles/Sample)')
xlim([0.2 0.4])
title('Standard Deviation of MUSIC Pseudospectrum, p = 1', 'fontsize', 12)
grid on
grid minor
set(gcf,'color','w')


subplot(3,2,3)
p = 2; % order
for i = 1:100
    noise = 0.2/sqrt(2) * (randn(size(n)) + 1j*randn(size(n)));
    noisySignal = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    [~, R] = corrmtx(noisySignal, 14, 'mod');
    [psds(i, :), x] = pmusic(R, p, [ ], 1);
    plot(x, psds(i, :), 'color', 'c', 'LineWidth', 1.2)
    hold on;
end

% mean
plot(x, mean(psds), 'color', 'b', 'LineWidth', 1.2)
ax = gca;
ax.FontSize = 12;
ylabel('Magnitude')
xlabel('Normalised Frequency (Cycles/Sample)')
xlim([0.2 0.4])
title('Realisations and Mean of MUSIC Pseudospectrum, p = 2', 'fontsize', 12)
grid on
grid minor

% standard deviation
subplot(3,2,4)
plot(x, std(psds), 'color','b','LineWidth', 1.2)
ax = gca;
ax.FontSize = 12;
ylabel('Magnitude')
xlabel('Normalised Frequency (Cycles/Sample)')
xlim([0.2 0.4])
title('Standard Deviation of MUSIC Pseudospectrum, p = 2', 'fontsize', 12)
grid on
grid minor
set(gcf,'color','w')


subplot(3,2,5)
p = 8; % order
for i = 1:100
    noise = 0.2/sqrt(2) * (randn(size(n)) + 1j*randn(size(n)));
    noisySignal = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    [~, R] = corrmtx(noisySignal, 14, 'mod');
    [psds(i, :), x] = pmusic(R, p, [ ], 1);
    plot(x, psds(i, :), 'color', 'c', 'LineWidth', 1.2)
    hold on;
end

% mean
plot(x, mean(psds), 'color', 'b', 'LineWidth', 1.2)
ax = gca;
ax.FontSize = 12;
ylabel('Magnitude')
xlabel('Normalised Frequency (Cycles/Sample)')
xlim([0.2 0.4])
title('Realisations and Mean of MUSIC Pseudospectrum, p = 8', 'fontsize', 12)
grid on
grid minor

% standard deviation
subplot(3,2,6)
plot(x, std(psds), 'color','b','LineWidth', 1.2)
ax = gca;
ax.FontSize = 12;
ylabel('Magnitude')
xlabel('Normalised Frequency (Cycles/Sample)')
xlim([0.2 0.4])
title('Standard Deviation of MUSIC Pseudospectrum, p = 8', 'fontsize', 12)
grid on
grid minor
set(gcf,'color','w')
