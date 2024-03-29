%% Q4.1c Scaled Activation Function
clc
clear
close all

% Intialisations
load('time-series.mat')
y = y - mean(y);
order = 4;
gamma = 0;
alphaStart = 40;
alphaEnd = 100;
alphas = alphaStart:0.1:alphaEnd;
stepSize = 0.00001;
MSEs= [];
predictiveGains = [];

% LMS
for alpha = alphas
    [yHat, ~, error] = LMS_dp(y, stepSize, order, alpha, 0, [0;0;0;0], 0);
    MSEs = [MSEs, 10*log10(mean(abs(error).^2))]; % mean squared & in dB
    predictiveGains = [predictiveGains, 10*log10(var(yHat)/var(error))];
end

% Plotting
figure
subplot(2,2,1)
plot(alphas, MSEs, 'LineWidth', 1.2)
hold on
[minError, index] = min(MSEs); % find the minimum MSE value
plot(alphaStart + index*0.1, minError,'r*', 'MarkerSize', 12, 'LineWidth', 0.9)
xlabel('\alpha', 'fontsize', 12)
ylabel('MSE (dB)', 'fontsize', 12)
title('Finding Optimal \alpha with MSE, \mu = 0.00001', 'fontSize', 12)
ax = gca;
ax.FontSize = 12; 
grid on 
grid minor

subplot(2,2,2)
plot(alphas, predictiveGains, 'LineWidth', 1.2)
hold on
[maxGain, index] = max(predictiveGains);
plot(alphaStart + index*0.1, maxGain, 'r*', 'MarkerSize', 12, 'LineWidth', 0.9)
xlabel('\alpha', 'fontsize', 12)
ylabel('Prediction Gain (dB)', 'fontsize', 12)
title('Finding Optimal \alpha with Predictive Gain, \mu = 0.00001', 'fontSize', 12)
ax = gca;
ax.FontSize = 12; 
grid on 
grid minor
set(gcf,'color','w')

% repeat for a smaller step size
stepSize = 0.0000001;
MSEs= [];
predictiveGains = [];

% LMS
for alpha = alphas
    [yHat, ~, error] = LMS_dp(y, stepSize, order, alpha, 0, [0;0;0;0], 0);
    MSEs = [MSEs, 10*log10(mean(abs(error).^2))]; % mean squared & in dB
    predictiveGains = [predictiveGains, 10*log10(var(yHat)/var(error))];
end

subplot(2,2,3)
plot(alphas, MSEs, 'LineWidth', 1.2)
hold on
[minError, index] = min(MSEs); % find the minimum MSE value
plot(alphaStart + index*0.1, minError,'r*', 'MarkerSize', 12, 'LineWidth', 0.9)
xlabel('\alpha', 'fontsize', 12)
ylabel('MSE (dB)', 'fontsize', 12)
title('Finding Optimal \alpha with MSE, \mu = 0.0000001', 'fontSize', 12)
ax = gca;
ax.FontSize = 12; 
grid on 
grid minor

subplot(2,2,4)
plot(alphas, predictiveGains, 'LineWidth', 1.2)
hold on
[maxGain, index] = max(predictiveGains);
plot(alphaStart + index*0.1, maxGain, 'r*', 'MarkerSize', 12, 'LineWidth', 0.9)
xlabel('\alpha', 'fontsize', 12)
ylabel('Prediction Gain (dB)', 'fontsize', 12)
title('Finding Optimal \alpha with Predictive Gain, \mu = 0.0000001', 'fontSize', 12)
ax = gca;
ax.FontSize = 12; 
grid on 
grid minor
set(gcf,'color','w')

% Plot the actual prediction
% calculating the MSE and Gains at the end of the prediction
mseEnd = 10*log10(mean(error(750:1000).^2));
predictiveGainEnd = 10*log10(var(yHat(750:1000))/var(error(750:1000)));

figure
subplot(1,2,1)
plot(y, 'LineWidth', 1.2)
hold on
plot(yHat, 'r', 'LineWidth', 1.2)
xlabel('Time Step n', 'fontsize', 12)
ylabel('Magnitude', 'fontsize', 12)
title('One-Step Ahead LMS Prediction of Time Series Data', 'fontsize', 12)
legend('Zero-Mean Time Series','LMS Prediction (+ Scaled Dyn. Perceptron)')
ax = gca;
ax.FontSize = 12; 
grid on
grid minor
set(gcf,'color','w')

subplot(1,2,2)
plot(750:1000, y(750:1000), 'LineWidth', 1.2)
hold on
plot(750:1000, yHat(750:1000), 'r', 'LineWidth', 1.2)
xlim([750,1000])
xlabel('Time Step n', 'fontsize', 12)
ylabel('Magnitude', 'fontsize', 12)
title('One-Step Ahead LMS Prediction of Time Series Data (Zoomed In)', 'fontsize', 12)
legend('Zero-Mean Time Series','LMS Prediction (+ Scaled Dyn. Perceptron)')
ax = gca;
ax.FontSize = 12; 
grid on
grid minor
set(gcf,'color','w')
