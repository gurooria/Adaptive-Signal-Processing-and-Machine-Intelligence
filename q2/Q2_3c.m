%% Q2.3c
clc
clear
close all

% Initialisations
N = 1000;
deltas = 3;
order = 5; 
stepSize = 0.01; 
iterations = 100;
a = 1;
b = [1 0 0.5];
x = sin(0.01*pi.*(0:N-1)); 
corruptedSignals = zeros(iterations, N);
xHatsALE = zeros(iterations, N);
mspeALE = zeros(iterations, 1);
xHatsANC = zeros(iterations, N);
mspeANC = zeros(iterations, 1);

% Generate results
for iteration = 1: iterations
    wgn = randn(1, N); % has unit variance
    coloredNoise = filter(b, a, wgn);
    s = x + coloredNoise;
    u = 1.2 * coloredNoise + 0.1;
    corruptedSignals(iteration, :) = s;

    % ALE
    [xHatALE, ~, ~] = LMS_ALE(s, stepSize, 0, order, deltas);
    xHatsALE(iteration, :) = xHatALE;
    mspeALE(iteration) = mean((x(100:end) - xHatALE(100:end)').^2);
    
    % ANC
    [xhatANC, ~, ~] = LMS_ANC(s, u, stepSize, 0, order);
    xHatsANC(iteration,:) = xhatANC;
    mspeANC(iteration) = mean((x(100:end)-xhatANC(100:end)').^2);

end

%% Plotting
figure
hold on

% plot ALE
subplot(1,3,1)
sPlot = plot(corruptedSignals', 'b', 'LineWidth', 1.2, 'DisplayName','$s(n)$');
hold on
xHatALEPlot = plot(xHatsALE', 'r', 'LineWidth', 1.2, 'DisplayName', '$\hat{x}(n)$');
hold on
xPlot = plot(x, 'y', 'LineWidth', 1.5, 'DisplayName', '$x(n)$');
meanMSPE = mean(mspeALE);
title(sprintf('Prediction using ALE, MSPE = %0.3f', meanMSPE), 'fontsize', 12);
xlabel('Sample Number n', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
set(groot,'defaultLegendInterpreter','latex');
set(gca, 'fontSize', 12);
allPlots = [sPlot(1), xHatALEPlot(1), xPlot];
legend(allPlots)
grid on
grid minor

% plot ANC
subplot(1,3,2)
sPlot = plot(corruptedSignals', 'b', 'LineWidth', 1.2, 'DisplayName', '$s(n)$');
hold on
xHatANCPlot = plot(xHatsANC', 'r', 'LineWidth', 1.2, 'DisplayName', '$\hat{x}(n)$');
hold on
xPlot = plot(x, 'y', 'LineWidth', 1.5, 'DisplayName', '$x(n)$');
meanMSPE = mean(mspeANC);
title(sprintf('Prediction using ANC, MSPE = %0.3f', meanMSPE), 'fontsize', 12);
xlabel('Sample Number n', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
set(gca, 'fontSize', 12);
set(groot, 'defaultLegendInterpreter', 'latex');
allPlots = [sPlot(1), xHatANCPlot(1), xPlot];
legend(allPlots)
grid on
grid minor 
set(gcf,'color','w')

% plot realisation ensemble
aleEnsemble = mean(xHatsALE);
ancEnsemble = mean(xHatsANC);
subplot(1,3,3)
plot(aleEnsemble, 'b', 'LineWidth', 1.2);
hold on
plot(ancEnsemble, 'r', 'LineWidth', 1.2);
hold on
plot(x, 'y', 'LineWidth', 1.2);
title('Comparison between the Ensemble Means of ALE vs ANC', 'fontsize', 12);
xlabel('Sample Number n', 'fontsize', 12);
ylabel('Magnitude', 'fontsize', 12);
set(gca, 'fontSize', 12);
set(groot, 'defaultLegendInterpreter', 'latex');
legend('ALE', 'ALC', '$x(n)$')
grid on
grid minor 
set(gcf,'color','w')
