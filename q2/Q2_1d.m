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
leaks = [0];

%%
coefficients = zeros(length(a)-1, N, realisations);

% loop through 100 realisations
figure(1)
hold on
for stepSize = 1: length(stepSizes)
    for leak = 1: length(leaks)

        % record the weights over 100 realisations
        for realisation = 1: realisations
            wgn = sqrt(var) .* randn(N,1);
            x =  filter(1, a, wgn);
            [~, weight, ~] =  LMS(x, stepSizes(stepSize), leaks(leak), length(a)-1); % exclude a0
            coefficients(:, :, realisation) = weight;
        end
        
        subplot(1,2,stepSize)
        h(1) = plot(mean(coefficients(1, :, :), 3), 'b', 'LineWidth', 1.2);
        hold on 
        if leak == length(leaks)
            yline(0.1, 'LineWidth', 1.2)
        end
        h(2) = plot(mean(coefficients(2,:,:), 3), 'r', 'LineWidth', 1.2);
        hold on
        if leak == length(leaks)
            yline(0.8, 'LineWidth', 1.2)
        end

        title(sprintf('Time Evolution of the coefficients over 100 Realisations, µ = %0.2f', stepSizes(stepSize)), 'fontsize', 12);
        xlabel("Time Step" , 'fontsize', 12);
        ylabel("Weight Value", 'fontsize', 12);
        legend([h(1), h(2)], '$\hat{a}_{1}$', '$\hat{a}_{2}$', 'interpreter','latex', 'fontsize', 12)
        grid on;
        grid minor;
        set(gcf, 'color', 'w')
        set(gca,'fontsize', 12);
    end
end