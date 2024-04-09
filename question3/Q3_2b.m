%% Q3.2b CLMS Power Spectra of Frequency Modulated Signal
clc
clear
close all

% Intialisations
N = 1500;
L = 1024;
n = 1 : N;
fs = 2000;
var = 0.05;
noise = sqrt(var).*randn(1,N) + 1j*sqrt(var).*randn(1,N); % circular complex-valued white noise

% generate f(n)
f = [100*ones(1,500), 100 + ((501:1000)-500)/2, 100 + (((1001:1500)-1000)/25).^2];
phi = cumtrapz(f);

% generate y(n)
y = exp(1j*((2*pi)/fs)*phi) + noise;
x = [0, y(1:end-1)];
stepSizes = [0.001, 0.05, 0.1, 0.5];
count = 1;

figure
hold on
for stepSize = stepSizes
    coeff = complex(zeros(1, N));
    [coeff, ~] = CLMS_fm(x, y, stepSize, 1);
    H = zeros(L, N);

    % code given by coursework
    for n = 1:N
        [h, w] = freqz(1, [1; -conj(coeff(n))], L);
        H(:, n) = abs(h).^2;
    end

    % Remove outliers
    medianH = 50 * median(median(H));
    H(H > medianH) = medianH;

    subplot(2,2,count)
    surf(1:N, ((w*fs)/(2*pi)), H, 'LineStyle','none');
    view(2)
    xlabel('Time Step n', 'fontsize', 12)
    ylabel('Frequency (Hz)', 'fontsize', 12)
    c = colorbar();
    c.Label.FontSize = 13;
    c.Label.String = "Power (dB)";
    title(sprintf('CLMS Time-Frequency Spectrum Estimation, \\mu = %0.3f', stepSize), 'fontsize', 12)
    ax = gca;
    ax.FontSize = 12;
    grid on
    grid minor
    count = count+1;
    hold on
end
set(gcf,'color','w')