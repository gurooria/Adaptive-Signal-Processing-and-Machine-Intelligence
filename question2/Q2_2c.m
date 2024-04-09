%% Q2.2c GNGD Algorithm
clc
clear
close all

% Intialisations
N =  1000;
realisations = 1000;
rho = 0.001;
alpha = 0.8;

b = [1, 0.9]; % MA model
order = length(b) - 1;
a = 1;

var = 0.5;
wgn = zeros(realisations, N);
for i = 1 : realisations
    wgn(i, :) = sqrt(var) * randn(1, N);
end

N = N/2;
weightsBen = zeros(realisations, N);
weightsGNGD = weightsBen;
mu = 0.1;
epsilonStart = 1/mu;

for i = 1: realisations
   x = filter(b, a, wgn(i, :));
   x =  x(501 : end);
   [~, weightB, ~, ~] = LMS_GASS(wgn(i, 501:end), x, 0.1, 0.002, alpha, order, 3);
   [~, weightG, ~] = GNGD(wgn(i, 501:end), x, 1, epsilonStart, order, 0.005);
   weightsBen(i, :) = weightB; 
   weightsGNGD(i, :) = weightG;
end

%% Plotting
figure

subplot(1,2,1)
hold on
errorsBen = b(2) * ones(realisations, N) - weightsBen;
errorsGNGD = b(2) * ones(realisations, N) - weightsGNGD;
plot(mean(errorsBen), 'b', 'LineWidth', 1.2);
hold on
plot(mean(errorsGNGD), 'r', 'LineWidth', 1.2);
ax = gca;
ax.FontSize = 12;
legend('Benveniste', 'GNGD', 'fontsize', 12);
xlabel('Time Step', 'fontsize', 12)
ylabel('Weight Error', 'fontsize', 12)
title('Weight Error Curves of Benveniste vs. GNGD', 'FontSize', 12)
grid on
grid minor

subplot(1,2,2)
plot(mean(10*log10(errorsBen.^2)), 'b', 'LineWidth', 1.2);
hold on
plot(mean(10*log10(errorsGNGD.^2)), 'r', 'LineWidth', 1.2);
ax = gca;
ax.FontSize = 12;
legend('Benveniste','GNGD', 'fontsize', 12);
xlabel('Time Step','fontsize', 12)
ylabel('Squared Weight Error (dB)', 'fontsize', 12)
xlim([0 100])
title('Squared Weight Error Curves of Benveniste vs. GNGD', 'FontSize', 12)
grid on
grid minor