%% Q2.3a Minimum Value of Delay
clc
clear
close all

% Intialisations
N = 1000;
deltas = [1, 2, 3, 4];
stepSize = 0.01; 
iterations = 100;
M = 5; % filter order 
a = 1;
b = [1 0 0.5];
x = sin(0.01*pi.*(0:N-1)); % REAL x

% Plotting
figure
hold on
for i = 1: length(deltas)
    xHats = zeros(iterations,N);
    corruptedSignals =  zeros(iterations,N);
    MSPE = zeros(iterations,1);
    
    for j = 1: iterations
        wgn = randn(1, N); 
        coloredNoise = filter(b, a, wgn);
        s = x + coloredNoise;
        corruptedSignals(j, :) = s;
        [xHat, ~, ~] = LMS_ALE(s, stepSize, 0, M, deltas(i));
        xHats(j, :) = xHat;
        MSPE(j) = mean((x(100:end) - xHat(100:end)').^2);
    end

    subplot(1, 4, i) % one plot for each delta
    sPlot = plot(corruptedSignals', 'b', 'LineWidth', 1.2, 'DisplayName', 's(n)');
    hold on
    xHatPlot = plot(xHats', 'r', 'LineWidth', 1.2, 'DisplayName', 'x_{hat}(n)');
    hold on
    xPlot = plot(x, 'y', 'LineWidth', 1.5, 'DisplayName', 'x(n)');
    meanMSPE = mean(MSPE);
    title(sprintf('MSPE = %0.3f, Delay = %0.0f', meanMSPE, deltas(i)), 'fontsize', 12)
    ax = gca;
    ax.FontSize = 12;
    xlabel('Sample Number (n)','fontsize', 12)
    ylabel('Magnitude', 'fontsize', 12)
    allPlots = [sPlot(1), xHatPlot(1), xPlot];
    legend(allPlots)
    grid on
    grid minor
end

set(gcf,'color','w')