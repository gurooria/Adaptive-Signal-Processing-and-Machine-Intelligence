%% 3.1c Balanced & Unbalanced systems
clc
clear
close all

% Intialisations
N = 1000;
n = 1 : N; % time vector
fo = 100; % system frequency
fs = 10000; % sampling frequency
clarkeMatrix = sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];

%% Balanced - V = [1, 1, 1], delta = [0, 0]
vMagnitudes = ones(1, 3); % peak voltage magnitudes
deltas = zeros(2, 1); % phase distortions

% compute voltages
vAll = zeros(3, 1000);
vAll(1, :) = vMagnitudes(:, 1) * cos(2*pi*(fo/fs)*n);
vAll(2, :) = vMagnitudes(:, 2) * cos(2*pi*(fo/fs)*n + deltas(1) - 2*pi/3);
vAll(3, :) = vMagnitudes(:, 3) * cos(2*pi*(fo/fs)*n + deltas(2) + 2*pi/3);

% clarke transform
vProjected = clarkeMatrix * vAll; % V0 = 0, Valpha, Vbeta
v = complex(vProjected(2, :), vProjected(3, :)); % change to complex clarke voltage using Valpha and Vbeta
circularity = abs(mean((v).^2)/mean(abs(v).^2));

% Plotting
figure
subplot(1,3,1)
hold on
scatter(real(v), imag(v), 12, 'filled');
xlim([-2 2])
ylim([-2 2])
xlabel('Real', 'fontSize', 12)
ylabel('Imaginiary', 'fontSize', 12)
set(gca, 'Fontsize', 12)
title(sprintf('Balanced system, ρ = 0'), 'FontSize', 12)
grid on
grid minor

%% Unbalanced - varying voltages
vMagnitudes = [0,1,1; 1,0,1; 1,1,0];
v = zeros(3, N); % to hold 3 conditions
deltas = zeros(2, 1); % phase distortion still normal

subplot(1,3,2)
hold on
circularity1 = [];

% Compute the voltages
for i = 1 : size(v, 1)
    % set the voltages
    vAll = zeros(3, 1000);% Va, Vb, Vc
    vAll(1, :) = vMagnitudes(i, 1) * cos(2*pi*(fo/fs)*n);
    vAll(2, :) = vMagnitudes(i, 2) * cos(2*pi*(fo/fs)*n + deltas(1) - 2*pi/3);
    vAll(3, :) = vMagnitudes(i, 3) * cos(2*pi*(fo/fs)*n + deltas(2) + 2*pi/3);
    vProjected = clarkeMatrix * vAll; % V0 = 0, Valpha, Vbeta
    v(i, :) = complex(vProjected(2, :), vProjected(3, :));
    scatter(real(v(i,:)), imag(v(i,:)), 12, 'filled')
    circularity1 = [circularity1, abs(mean((v(i, :)).^2)/mean(abs(v(i, :)).^2));];
end

xlim([-2 2])
ylim([-2 2])
xlabel('Real', 'fontSize', 12)
ylabel('Imaginiary', 'fontSize', 12)
set(gca, 'Fontsize', 12)
legend('$V_{a, b, c}$ = [0, 1, 1]','$V_{a, b, c}$ = [1, 0, 1]','$V_{a, b, c}$ = [1, 1, 0]')
title('Unbalanced system: Voltage, ρ = 0.8', 'FontSize', 12)
grid on
grid minor

%% Unbalanced - varying phases
vMagnitudes = ones(1, 3);
deltas = [pi/3,0; -pi/3,0; 0,pi/3; 0,-pi/3];
v = zeros(size(deltas,1), N);

subplot(1,3,3)
hold on
circularity2 = [];

% Compute the voltages
for i = 1 : size(v, 1)
    vAll = zeros(3, 1000);% Va, Vb, Vc
    vAll(1, :) = vMagnitudes(:, 1) * cos(2*pi*(fo/fs)*n);
    vAll(2, :) = vMagnitudes(:, 2) * cos(2*pi*(fo/fs)*n + deltas(i, 1) - 2*pi/3);
    vAll(3, :) = vMagnitudes(:, 3) * cos(2*pi*(fo/fs)*n + deltas(i, 2) + 2*pi/3);
    vProjected = clarkeMatrix * vAll; % V0 = 0, Valpha, Vbeta
    v(i, :) = complex(vProjected(2, :), vProjected(3, :));
    scatter(real(v(i,:)), imag(v(i,:)), 12, 'filled')
    circularity2 = [circularity2, abs(mean((v(i, :)).^2)/mean(abs(v(i, :)).^2));];
end

xlim([-2 2])
ylim([-2 2])
xlabel('Real', 'fontSize', 12)
ylabel('Imaginiary', 'fontSize', 12)
set(gca, 'Fontsize', 12)
legend('$\Delta_{b, c}$ = [$\frac{\pi}{3}$, 0]','$\Delta_{b, c}$ = [$\frac{-\pi}{3}$, 0]', '$\Delta_{b, c}$ = [0, $\frac{\pi}{3}$]', '$\Delta_{b, c}$ = [0, $\frac{-\pi}{3}$]')
title('Unbalanced system: Phase, ρ = 0.635', 'FontSize', 12)
grid on
grid minor