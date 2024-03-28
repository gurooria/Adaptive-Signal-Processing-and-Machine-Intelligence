%% 3.1e Frequency Estimation using CLMS and ACLMS
clc
clear
close all

% Intialisations
N = 500;
n = 1 : N; % time vector
fo = 100; % system frequency
fs = 10000; % sampling frequency
clarkeMatrix = sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
stepSize = 0.02;

%% Balanced system
vMagnitudes = ones(1, 3); % peak voltage magnitudes
deltas = zeros(2, 1); % phase distortions

% compute voltages
vAll = zeros(3, N);
vAll(1, :) = vMagnitudes(:, 1) * cos(2*pi*(fo/fs)*n);
vAll(2, :) = vMagnitudes(:, 2) * cos(2*pi*(fo/fs)*n + deltas(1) - 2*pi/3);
vAll(3, :) = vMagnitudes(:, 3) * cos(2*pi*(fo/fs)*n + deltas(2) + 2*pi/3);
vProjected = clarkeMatrix * vAll; % V0 = 0, Valpha, Vbeta
v = complex(vProjected(2, :), vProjected(3, :)); % change to complex clarke voltage using Valpha and Vbeta

% CLMS & ACLMS
[h, ~] = CLMS_voltage(v, stepSize, 1);
foCLMS = (fs/(2*pi)) * atan(imag(h) ./ real(h));
[g, h, ~] = ACLMS_voltage(v, stepSize, 1);
foACLMS = ((fs/(2*pi)) * atan(sqrt((imag(h)).^2 - abs(g).^2)./real(h)));

% Plotting
figure
hold on
subplot(1,2,1)
plot(abs(foCLMS), 'LineWidth', 1.2)
hold on
plot(abs(foACLMS), 'r', 'LineWidth', 1.2)
hold on
yline(fo, 'k--', 'LineWidth', 1.2)
title("Frequency Estimation for Balanced System", 'fontsize', 12);
xlabel("Sample Number n", 'FontSize', 12);
ylabel('Frequency fo (Hz)', 'FontSize', 12);
ylim([0 150])
legend('CLMS','ACLMS','True Value','Interpreter','latex')
ax = gca;
ax.FontSize = 12;
grid on
grid minor
set(gcf,'color','w')

%% Unbalanced System
vMagnitudes = [1, 1, 1]; % peak voltage magnitudes
deltas = [0, 1]; % phase distortions

% compute voltages
vAll = zeros(3, N);
vAll(1, :) = vMagnitudes(:, 1) * cos(2*pi*(fo/fs)*n);
vAll(2, :) = vMagnitudes(:, 2) * cos(2*pi*(fo/fs)*n + deltas(1) - 2*pi/3);
vAll(3, :) = vMagnitudes(:, 3) * cos(2*pi*(fo/fs)*n + deltas(2) + 2*pi/3);
vProjected = clarkeMatrix * vAll; % V0 = 0, Valpha, Vbeta
v = complex(vProjected(2, :), vProjected(3, :)); % change to complex clarke voltage using Valpha and Vbeta

% CLMS & ACLMS
[h, ~] = CLMS_voltage(v, stepSize, 1);
foCLMS = (fs/(2*pi)) * atan(imag(h) ./ real(h));
[g, h, ~] = ACLMS_voltage(v, stepSize, 1);
foACLMS = ((fs/(2*pi)) * atan(sqrt((imag(h)).^2 - abs(g).^2)./real(h)));

% Plotting
hold on
subplot(1,2,2)
plot(abs(foCLMS), 'LineWidth', 1.2)
hold on
plot(abs(foACLMS), 'r', 'LineWidth', 1.2)
hold on
yline(fo, 'k--', 'LineWidth', 1.2)
title("Frequency Estimation for Unbalanced System", 'fontsize', 12);
xlabel("Sample Number n", 'FontSize', 12);
ylabel('Frequency fo (Hz)', 'FontSize', 12);
ylim([0 150])
legend('CLMS','ACLMS','True Value','Interpreter','latex')
ax = gca;
ax.FontSize = 12;
grid on
grid minor
set(gcf,'color','w')