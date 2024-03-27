%Q2.1b Implementation of an LMS Adaptive Predictor
clc
clear
close all

% Initialisations
N = 1000; % specified 1000 samples
stepSizes = [0.01, 0.05]; % specified step sizes
realisations = 100; % 100 realisations
a0 = 1;
a1 = -0.1;
a2 = -0.8;
a = [a0, a1, a2]; % initial parameters
var = 0.25; % variance of noise
wgn = sqrt(var) .* randn(N,1); % white gaussian noise with specified std
x =  filter(1 , a , wgn); % shape the noise using the filter
leak = 0;

%% Plotting

% squared prediction error for 1 realisation
figure (1)
hold on

colors = {'b','r'};
for stepSize = 1: length(stepSizes) % loop through each step size
    [xPredicted, ~, error] =  LMS(x, stepSizes(stepSize), leak, length(a)-1); % exclude a0
    subplot(1,2,1)
    plot(10*log10(error.^2), 'color', colors{stepSize}, 'LineWidth', 1.2)
    hold on
end
title("Squared Prediction Error e^{2}(n)  with 1 Realisation", 'fontsize', 12);
xlabel("Time Step", 'fontsize', 12);
ylabel("Squared Prediction Error (dB)", 'fontsize', 12);
legend('\mu = 0.01','\mu = 0.05')
grid on; grid minor;
set(gca,'fontsize', 12);

% squared prediction error for 100 realisations
empiricalMisadjust = []; % *for Q2_1c
errors = zeros(N, realisations); % plot error over 100 realisations
for stepSize = 1: length(stepSizes)  
    for i = 1: realisations
        wgn = sqrt(var) .* randn(N,1); % regenerate wgn for every realisation
        x =  filter(1, a, wgn);
        [xPredicted, ~, error] =  LMS(x, stepSizes(stepSize), leak, length(a)-1); % excluding a0
        errors (:, i)= error;
    end
    subplot(1,2,2)
    plot(10*log10(mean(errors.^2,2)), 'color', colors{stepSize}, 'LineWidth', 1.2)
    hold on
    % Q2.1c
    % assume steady state after t = 600
    empiricalMisadjust = [empiricalMisadjust, (mean(mean(errors(600:end, :) .^ 2, 1), 2) / 0.25 -1)];
end

title("Squared Prediction Error e^{2}(n)  with 100 Realisations", 'fontsize', 12);
xlabel("Time Step", 'fontsize', 12);
ylabel("Squared Prediction Error (dB)", 'fontsize', 12);
legend('\mu = 0.01','\mu = 0.05')
grid on; grid minor;
set(gcf,'color','w')
set(gca,'fontsize', 12);
