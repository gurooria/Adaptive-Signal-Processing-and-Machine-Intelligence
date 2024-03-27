%% Q1.4b & c PSD of ARMA Process
clc;
clear;
close all;

% 1.4b Initialisations
n = 1000;
removed = 500; % discard the first 500 samples to remove the transient filter effects
a = [1, -2.76, 3.81, -2.65, 0.92]; % filter coefficients
x =  filter(1, a, randn(n, 1));
x = x(removed+1 : end);

figure(1)
hold on

% true power spectrum
[h, ~] = freqz(1, a, length(x));
psd = abs(h).^2;
orders = [2, 4, 11];
MSEsB = [];

for i = 1: length(orders)
    % The Yule-Walker method provides an autoregressive power spectral density estimate
    [aEstimate, e] = aryule(x, orders(i));
    [h, w] = freqz(e^(1/2), aEstimate, length(x));
    psdEstimate = abs(h).^2;

    subplot(2,3,i)
    plot(w/pi, 10*log10(psd), 'LineWidth', 1.2);
    hold on
    plot(w/pi, 10*log10(psdEstimate), 'r', 'LineWidth', 1.2)
    MSEsB = [MSEsB, mean((10*log10(psd) - 10*log10(psdEstimate)).^2)];
    ax = gca;
    ax.FontSize = 12;
    axis tight
    xlabel('Normalised Frequency (Cycles/Sample)')
    ylabel('Magnitude (dB)')
    title('Autoregressive PSD Estimation, n = 1000', 'fontsize', 12)
    switch i 
        case 1
            legend('True', 'AR(2)')
        case 2
            legend('True', 'AR(4)')
        case 3
            legend('True', 'AR(11)')
    end
    grid on 
    grid minor
    set(gcf, 'color','w');
    hold on
end


% 1.4c
n = 10000;
a = [1, -2.76, 3.81, -2.65, 0.92];
x =  filter(1, a, randn(n,1));
x = x(removed+1 : end);

% true power spectrum
[h, ~] = freqz(1,a,length(x));
psd = abs(h).^2;
orders = [2, 4, 11];
MSEsC = [];
for i = 1: length(orders)
    [aEstimate, e] = aryule(x, orders(i));
    [h, w] = freqz(1, aEstimate, length(x));
    psdEstimate = abs(h).^2;
    subplot(2,3,i+3)

    plot(w/pi, 10*log10(psd), 'LineWidth', 1.2);
    hold on
    plot(w/pi, 10*log10(psdEstimate), 'r', 'LineWidth', 1.2)
    MSEsC = [MSEsC,mean((10*log10(psd) - 10*log10(psdEstimate)).^2)];
    ax = gca;
    ax.FontSize = 12;
    axis tight
    xlabel('Normalised Frequency (Cycles/Sample)')
    ylabel('Magnitude (dB)')
    title('Autoregressive PSD Estimation, n = 10000', 'fontsize', 12)
    switch i 
        case 1
            legend('True', 'AR(2)')
        case 2
            legend('True', 'AR(4)')
        case 3
            legend('True', 'AR(11)')
    end
    grid on 
    grid minor
    set(gcf, 'color','w');
    hold on
end
