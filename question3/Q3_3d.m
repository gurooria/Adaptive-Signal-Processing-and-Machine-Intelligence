%% Q3.3d Implementation of the DFT-CLMS algorithm - EEG signal
clc
clear
close all

% Intialisations
load('EEG_Data_Assignment1.mat')
N = 1200;
fs = 1200;
timeStart = 1000;
y = POz(timeStart : timeStart+N-1);
y = detrend(y);
L = 1024;
x = (1/L)*exp(1j*2*(0:N-1)'*pi*(0:(L-1))/L).';
stepSize = 1;
count = 1;
gamma = 0;

figure
coeff = complex(zeros(1, N));
[coeff, ~] = CLMS_dft(x, y, stepSize, gamma, L);
H = abs(coeff);
medianH = 50*median(median(H));
H(H > medianH) = medianH;
surf(1:N, (0:(L-1)).*(fs/L), H, 'LineStyle','none');
view(2)
c = colorbar();
c.Label.FontSize = 13;
c.Label.String = "Power (dB)";
xlabel('Time Step n', 'fontsize', 12)
ylabel('Frequency (Hz)', 'fontsize', 12)
title('CLMS-DFT Time-Frequency Spectrum Estimation', 'fontsize', 12)
ax = gca;
ax.FontSize = 12;
grid on
grid minor
ylim([0 100])
set(gcf,'color','w')