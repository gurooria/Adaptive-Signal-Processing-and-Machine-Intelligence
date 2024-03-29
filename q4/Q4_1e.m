%% Q4.1e Pre-Trained Weights
clc
close all
clear

% Initialisations
load('time-series.mat')
stepSize = 0.0000001;
gamma = 0;
order = 4;
bias = 1;
alphaStart =  40;
alphaEnd = 100;
alphas = alphaStart:0.1:alphaEnd;

% Pre-training
a = 67.5; % temporary alpha
weightPretrained = zeros(order+bias, 1);
segment = 20; % train on first 20 samples
ySegment = y(1:segment); % for pretraining
epochs = 100;
for e = 1:epochs
    [~, weights, ~] = LMS_dp(y, stepSize, order, a, 0, weightPretrained, bias);
    weightPretrained = weights(:, end);
end

% Find optimal alpha
MSEs = [];
predictiveGains = [];
for alpha = alphas
    [yHat, ~, error] = LMS_dp(y, stepSize, order, alpha, 0, weightPretrained, bias);
    MSEs = [MSEs, 10*log10(mean(abs(error).^2))]; % mean squared & in dB
    predictiveGains = [predictiveGains, 10*log10(var(yHat)/var(error))];
end

%Plotting
figure
subplot(2,2,1)
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

subplot(2,2,2)
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

%% LMS with the optimal alpha
alpha = 68.3;
[yHat, ~, error] = LMS_dp(y, stepSize, order, alpha, 0, weightPretrained, bias);
MSE = 10*log10(mean(abs(error).^2)); % mean squared & in dB
predictiveGain = 10*log10(var(yHat)/var(error));

% Plotting
subplot(2,2,3)
plot(y, 'LineWidth', 1.2)
hold on
plot(yHat, 'r', 'LineWidth', 1.2)
xlabel('Time Step n', 'fontsize', 12)
ylabel('Magnitude', 'fontsize', 12)
title('One-Step Ahead LMS Prediction of Time Series Data', 'fontsize', 12)
legend('Time Series','LMS Prediction (+ Pre-trained)')
ax = gca;
ax.FontSize = 12; 
grid on
grid minor
set(gcf,'color','w')

subplot(2,2,4)
plot(750:1000,y(750:end), 'LineWidth', 1.2)
hold on
plot(750:1000, yHat(750:end), 'r', 'LineWidth', 1.2)
xlim([750,1000])
xlabel('Time Step n', 'fontsize', 12)
ylabel('Magnitude', 'fontsize', 12)
title('One-Step Ahead LMS Prediction of Time Series Data (Zoomed In)', 'fontsize', 12)
legend('Time Series','LMS Prediction (+ Pre-trained)')
ax = gca;
ax.FontSize = 12; 
grid on
grid minor
set(gcf,'color','w')