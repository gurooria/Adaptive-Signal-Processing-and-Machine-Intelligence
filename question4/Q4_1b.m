%% Q4.1b Dynamical Perceptron and Activation Function
clc
close
clear all

% Initialisations
load('time-series.mat')
stepSize = 0.00001;
gamma = 0;
y = y - mean(y);
order = 4;
alpha = 1;

%% LMS & Plotting
[yHat, ~, error] = LMS_dp(y, stepSize, order, alpha, 0, [0;0;0;0], 0);

figure
subplot(1,2,1)
plot(y, 'LineWidth', 1.2)
hold on
plot(yHat, 'r', 'LineWidth', 1.2)
xlabel('Time Step n', 'fontsize', 12)
ylabel('Magnitude', 'fontsize', 12)
title('One-Step Ahead LMS Prediction of Time Series Data', 'fontsize', 12)
legend('Zero-Mean Time Series','LMS Prediction (+ Dyn. Perceptron)')
ax = gca;
ax.FontSize = 12; 
grid on
grid minor

subplot(1,2,2)
plot(750:1000, y(750:end), 'LineWidth', 1.2)
hold on
plot(750:1000, yHat(750:end), 'r', 'LineWidth', 1.2)
xlim([750,1000])
xlabel('Time Step n', 'fontsize', 12)
ylabel('Magnitude', 'fontsize', 12)
title('One-Step Ahead LMS Prediction of Time Series Data', 'fontsize', 12)
legend('Zero-Mean Time Series','LMS Prediction (+ Dyn. Perceptron)')
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
