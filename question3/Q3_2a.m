%% 3.2a AR Power Spectra of Frequency Modulated Signal
clc
clear
close all

% Initialisations
N = 1500; % specified by the question
n = 1 : N;
fs = 2000;
var = 0.05;
noise = sqrt(var).*randn(1,N) + 1j*sqrt(var).*randn(1,N); % circular complex-valued white noise

% generate f(n) and phi(n)
f = [100*ones(1,500), 100 + ((501:1000)-500)/2, 100 + (((1001:1500)-1000)/25).^2];
phi = cumtrapz(f);

% plot fequency
figure
subplot(1,2,1)
plot(f, 'b', 'LineWidth', 1.2)
ylim([0 500])
xlabel('Time Step n', 'fontsize', 12)
ylabel('Frequency (Hz)', 'fontsize', 12)
title('Frequency f(n) of Frequency Modulated Signal','fontsize', 12)
set(gca, 'fontSize', 12)
grid on
grid minor

% generate y(n)
y = exp(1j*((2*pi)/fs)*phi) + noise;

% using aryule to find the AR coefficient
orders = [1, 5, 10];
colors = {[0 0 1],[1 0 0],[1 0 1]};
count = 1;

subplot(1,2,2)
hold on
for i = orders
    a = aryule(y, i);
    [h, w] = freqz(1, a, N, fs);
    psd = 10*log10(abs(h).^2);
    plot(w, psd, 'color', colors{count}, 'LineWidth', 1.2)
    xlabel('Frequency (Hz)', 'fontsize', 12)
    ylabel('Power (dB)', 'fontsize', 12)
    title('Power Spectra of Frequency Modulated Signal', 'fontsize', 12)
    ax = gca;
    ax.FontSize = 12;
    grid on
    grid minor
    count = count+1;
end
legend('AR(1)', 'AR(5)', 'AR(10)')
set(gcf,'color','w')


%% Splitting into 3 segments

figure
hold on
for i = 1:3
    a = aryule(y(500*(i-1)+1:500*i), 1);
    [h, w] = freqz(1, a, N/3, fs);
    psd = 10*log10(abs(h).^2);
    subplot(1,3,i)
    plot(w, psd, 'b', 'LineWidth', 1.2)
    xlabel('Frequency (Hz)', 'FontSize', 12)
    ylabel('Power (dB)', 'fontsize', 12)
    title(sprintf('AR(1) Power Spectrum of Segment %0.0f', i), 'fontsize', 12)
    ax = gca;
    ax.FontSize = 12;
    grid on
    grid minor
end
set(gcf,'color','w')