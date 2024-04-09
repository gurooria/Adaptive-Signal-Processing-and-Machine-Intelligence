%% Q4.1a Predicting Time Series with LMS
clc
clear
close all

% Intialisations
load('time-series.mat')
stepSize = 0.00001;
gamma = 0;
y = y - mean(y); % zero-mean
order = 4; % AR(4)

%% LMS & Plotting
[yHat, ~, error] = LMS(y, stepSize, gamma, order);

figure
subplot(1,2,1)
plot(y, 'LineWidth', 1.2)
hold on
plot(yHat, 'r', 'LineWidth', 1.2)
xlabel('Time Step n', 'fontsize', 12)
ylabel('Magnitude', 'fontsize', 12)
title('One-Step Ahead LMS Prediction of Time Series Data', 'fontsize', 12)
legend('Zero-Mean Time Series','LMS Prediction')
ax = gca;
ax.FontSize = 12;
grid on
grid minor
set(gcf,'color','w')

subplot(1,2,2)
plot(750:1000, y(750:end), 'LineWidth', 1.2)
hold on
plot(750:1000, yHat(750:end), 'r', 'LineWidth', 1.2)
xlim([750,1000])
xlabel('Time Step n','fontsize', 12)
ylabel('Magnitude', 'fontsize', 12)
title('AR(4) Prediction of Time Series Data', 'fontsize', 12)
legend('Original Time Series','AR(4) Prediction')
ax = gca;
ax.FontSize = 12; 
grid on
grid minor
set(gcf,'color','w')

% Calculate MSE
MSE = 10*log10(mean(error.^2));
predictionGain = 10*log10(var(yHat)/var(error));

MSEend = 10*log10(mean(error(750:1000).^2));
predictionGainEnd = 10*log10(var(yHat(750:1000))/var(error(750:1000)));