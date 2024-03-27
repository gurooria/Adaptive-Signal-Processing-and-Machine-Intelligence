%% Q2.2a GASS algorithms
clc
clear
close all

% Initialisations
N =  1000; % sample length
realisations = 1000;
rho = 0.001;
alpha = 0.8;

b = [1, 0.9]; % the given MA model
order = length(b) - 1;
a = 1;

% Generate 1000 realisations of wgn
var = 0.5;
wgn = zeros(realisations, N);
for i = 1 : realisations
    wgn(i, :) = sqrt(var) * randn(1, N);
end

%% Plotting

figure
hold on
ax = gca;
ax.FontSize = 12;
colors = {'b', 'r', 'm', 'c', [0.5, 0, 0.5]};

stepSizes = [0.01, 0.1, 0.2, 0.2, 0.2];
subplot(1,3,1)
for i = 1:5
    weights = zeros(realisations,N);
    xHats = weights;
    MSPE = zeros(realisations,1);
    
    for realisation = 1: realisations
       x = filter(b, a, wgn(realisation, :)); % shape noise with MA filter
       [xHat, weight, error, endStepSize] = LMS_GASS(wgn(realisation, :), x, stepSizes(i), rho, alpha, order, i);
       weights(realisation, :) = weight; 
    end
      errors = b(2) * ones(realisations,N) - weights;
      plot(mean(errors), 'color', colors{i}, 'LineWidth', 1.2);
      hold on
end

ax = gca;
ax.FontSize = 12;
set(gca,'fontsize', 12);
legend('\mu = 0.01','\mu = 0.1', 'Benveniste', 'Ang and Farhang', 'Matthews and Xie', 'fontsize', 12);
xlabel('Time Step', 'fontsize', 12)
ylabel('Weight Error', 'fontsize', 12)
title('Weight Error Curves, GASS intial \mu = 0.2')
grid on
grid minor

stepSizes = [0.01, 0.1, 0.1, 0.1, 0.1];
subplot(1,3,2)
for i = 1:5
    weights = zeros(realisations, N);
    xHats = weights;
    MSPE = zeros(realisations, 1);
    
    for realisation = 1 : realisations
       x = filter(b, a, wgn(realisation, :));
       [xHat, weight, error, endStepSize] = LMS_GASS(wgn(realisation, :), x, stepSizes(i), rho, alpha, order, i);
       weights(realisation, :) = weight; 
    end
      errors = b(2) * ones(realisations, N) - weights;
      plot(mean(errors), 'color', colors{i}, 'LineWidth', 1.2);
      hold on
end
ax = gca;
ax.FontSize = 12;
legend('\mu = 0.01','\mu = 0.1', 'Benveniste', 'Ang and Farhang', 'Matthews and Xie', 'fontsize', 12);
xlabel('Time Step', 'fontsize', 12)
ylabel('Weight Error', 'fontsize', 12)
title('Weight Error Curves, GASS intial \mu = 0.1')
grid on
grid minor
set(gcf,'color','w')
set(gca,'fontsize', 12);

stepSizes = [0.01, 0.1, 0.1, 0.1, 0.1];
subplot(1,3,3)
for i = 1 : 5
    weights = zeros(realisations, N);
    xHats = weights;
    MSPE = zeros(realisations, 1);
    
    for realisation = 1 : realisations
       x = filter(b, a, wgn(realisation, :));
       [xHat, weight, error, endStepSize] = LMS_GASS(wgn(realisation, :), x, stepSizes(i), rho, alpha, order, i);
       weights(realisation, :) = weight; 
    end
      errors = b(2) * ones(realisations, N) - weights;
      plot(10*log10(mean(errors.^2)), 'color', colors{i}, 'LineWidth', 1.2);
      hold on
end
ax = gca;
ax.FontSize = 12;
legend('\mu = 0.01','\mu = 0.1', 'Benveniste', 'Ang and Farhang', 'Matthews and Xie', 'fontsize', 12);
xlabel('Time Step', 'fontsize', 12)
ylabel('Squared Weight Error (dB)', 'fontsize', 12)
title('Squared Weight Error Curves, GASS intial \mu = 0.1')
grid on
grid minor
set(gcf,'color','w')
set(gca,'fontsize', 12);