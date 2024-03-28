%Q2.1d Estimation of the Steady State Values of the Adaptive Filter Coefficients
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
leaks = [0.1, 0.5, 0.8];

%% Plotting

hold on
colors = {[0 0 1],[1 0 0],[1 0 1]};

% step size = 0.01
coefficients = zeros(length(a)-1, N, realisations);
for leak = 1: length(leaks)
    for realisation = 1: realisations
        wgn = sqrt(var) .* randn(N, 1);
        x =  filter(1, a, wgn);
        [~, weights, ~] =  LMS(x, stepSizes(1), leaks(leak), length(a)-1); % excluding a0
        coefficients(:, :, realisation) = weights;
    end
    
    subplot(2,2,1)
    plot(mean(coefficients(1, :, :), 3), 'color', colors{leak}, 'LineWidth', 1.2);
    hold on
    if leak == length(leaks)
        yline(0.1, 'LineWidth', 1.2)
    end

    title("Time Evolution of a_1 over 100 Realisations, \mu = 0.01", 'fontsize', 12);
    xlabel("Time Step" , 'fontsize', 12);
    ylabel("Weight Value", 'fontsize', 12);
    legend('$\gamma$ = 0.1', '$\gamma$ = 0.5', '$\gamma$ = 0.8')
    grid on; grid minor;
    set(gca,'fontsize', 12);
    
    subplot(2,2,2)
    plot(mean(coefficients(2, :, :), 3), 'color', colors{leak}, 'LineWidth', 1.2);
    hold on
    if leak == length(leaks)
        yline(0.8, 'LineWidth', 1.2)
    end
    title("Time Evolution of a_2 over 100 Realisations, \mu = 0.01", 'fontsize', 12);
    xlabel("Time Step", 'fontsize', 12);
    ylabel("Weight Value", 'fontsize', 12);
    legend('$\gamma$ = 0.1', '$\gamma$ = 0.5', '$\gamma$ = 0.8')
    grid on; grid minor;
    set(gcf,'color','w');
    set(gca,'fontsize', 12);
end

% step size = 0.05
coefficients = zeros(length(a)-1, N, realisations);
for leak = 1: length(leaks)
    for realisation = 1: realisations
        wgn = sqrt(var) .* randn(N, 1);
        x =  filter(1, a, wgn);
        [~, weights, ~] =  LMS(x, stepSizes(2), leaks(leak), length(a)-1); % excluding a0
        coefficients(:, :, realisation) = weights;
    end
    
    subplot(2,2,3)
    plot(mean(coefficients(1, :, :), 3), 'color', colors{leak}, 'LineWidth', 1.2);
    hold on
    if leak == length(leaks)
        yline(0.1, 'LineWidth', 1.2)
    end

    title("Time Evolution of a_1 over 100 Realisations, \mu = 0.05", 'fontsize', 12);
    xlabel("Time Step" , 'fontsize', 12);
    ylabel("Weight Value", 'fontsize', 12);
    legend('$\gamma$ = 0.1', '$\gamma$ = 0.5', '$\gamma$ = 0.8')
    grid on; grid minor;
    set(gca,'fontsize', 12);
    
    subplot(2,2,4)
    plot(mean(coefficients(2, :, :), 3), 'color', colors{leak}, 'LineWidth', 1.2);
    hold on
    if leak == length(leaks)
        yline(0.8, 'LineWidth', 1.2)
    end
    title("Time Evolution of a_2 over 100 Realisations, \mu = 0.05", 'fontsize', 12);
    xlabel("Time Step", 'fontsize', 12);
    ylabel("Weight Value", 'fontsize', 12);
    legend('$\gamma$ = 0.1', '$\gamma$ = 0.5', '$\gamma$ = 0.8')
    grid on; grid minor;
    set(gcf,'color','w');
    set(gca,'fontsize', 12);
end