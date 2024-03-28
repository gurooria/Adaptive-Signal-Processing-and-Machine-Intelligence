%% Q3.3c Implementation of the DFT-CLMS algorithm - FM signal
clc
clear
close all

% Intialisations
N = 1500;
fs = 2000;
var = 0.05;
noise = sqrt(var).*randn(1,N) + 1j*sqrt(var).*randn(1,N); % circular complex-valued white noise
gammas = [0, 0.01, 0.1, 0.4];
stepSize = 1;
count = 1;

% generate f(n)
f = [100*ones(1,500), 100 + ((501:1000)-500)/2, 100 + (((1001:1500)-1000)/25).^2];
phi = cumtrapz(f);

% generate y(n)
y = exp(1j*((2*pi)/fs)*phi) + noise;
L = 1024;
w = (0:(L-1)) .* (fs / L);
x = (1/L) * exp(1j * 2 * (0:N-1)' * pi * (0:(L-1))/L).';

figure
hold on
for gamma = gammas
    coeff = complex(zeros(1, N)); 
    [coeff, ~] = CLMS_dft(x, y, stepSize, gamma, L);
    H = abs(coeff);
    medianH = 50 * median(median(H));
    H(H > medianH) = medianH;
    
    % Plotting
    subplot(2,2,count)
    surf(1:N, (0:(L-1)).*(fs/L), H, 'LineStyle','none');
    view(2)
    c = colorbar();
    c.Label.FontSize = 13;
    c.Label.String = "Power (dB)";
    xlabel('Time Step n', 'fontsize', 12);
    ylabel('Frequency (Hz)', 'fontsize', 12);
    ylim([0 700])
    title(sprintf('CLMS-DFT Time-Frequency Spectrum Estimation, \\gamma = %0.3f', gamma), 'fontsize', 12)
    ax = gca;
    ax.FontSize = 12;
    grid on
    grid minor
    hold on
    count = count+1;
end

set(gcf,'color','w')