%% Q2.3b
clc
clear
close all

%% Initialisations for MSPE vs. delays
N = 1000;
deltas = 1:25;
orders = [5, 10, 15, 20]; % filter order
stepSize = 0.01; 
iterations = 100;
colors = {'b', 'r', 'm', 'c'};
a = 1;
b = [1 0 0.5];
x = sin(0.01*pi.*(0:N-1)); % REAL x
mspeDelays = zeros(length(orders), length(deltas));

% Plotting
figure
subplot(1,2,1)
hold on
for order = 1 : length(orders)
    for i = 1 : length(deltas)
        xHats = zeros(iterations, N);
        corruptedSignals =  zeros(iterations, N);
        MSPE = zeros(iterations, 1);
        
        for j = 1: iterations
            wgn = randn(1, N); 
            coloredNoise = filter(b, a, wgn);
            s = x + coloredNoise;
            corruptedSignals(j, :) = s;
            [xHat, ~, ~] = LMS_ALE(s, stepSize, 0, orders(order), deltas(i));
            xHats(j, :) = xHat;
            MSPE(j) = mean((x(100:end) - xHat(100:end)').^2);
        end
        mspeDelays(order, i) = mean(MSPE);
    end
    plot(mspeDelays(order, :), 'Color', colors{order}, 'LineWidth', 1.2)
end

xlabel('Delay', 'FontSize', 12);
ylabel('MSPE', 'FontSize', 12);
legend('M = 5','M = 10','M = 15','M = 20')
title('Dependence between Delay and MSPE', 'FontSize', 12);
grid on
grid minor
set(gca, 'fontSize', 12);
set(gcf,'color','w')

%% Initialisations for Order vs. MSPE
N = 1000;
deltas = 3;
stepSize = 0.01; 
iterations = 100;
orders = 1:20; 
a = 1;
b = [1 0 0.5]; 
x = sin(0.01*pi.*(0:N-1)); 
mspeOrders = zeros(1, length(deltas));

subplot(1,2,2)
hold on
for i = 1 : length(orders)
    xHats = zeros(iterations, N);
    corruptedSignals = zeros(iterations, N);
    MSPE = zeros(iterations, 1);
    
    for j = 1 : iterations
        wgn = randn(1,N); % has unit variance
        coloredNoise = filter(b, a, wgn);
        s = x + coloredNoise;
        corruptedSignals(j, :) = s;
        [xHat, ~, ~] = LMS_ALE(s, stepSize, 0, orders(i), deltas);
        xHats(j, :) = xHat;
        MSPE(j) = mean((x(100:end) - xHat(100:end)').^2);
    end
    mspeOrders(i) =  mean(MSPE); 
end

plot(mspeOrders,'Color', 'b', 'LineWidth', 1.2)
set(groot,'defaultLegendInterpreter', 'latex');
xlabel('Model Order M', 'fontsize', 12);
ylabel('MSPE', 'fontsize', 12);
title('Effect of Filter Order M on MSPE, Delay = 3', 'fontsize', 12)
grid on
grid minor
set(gca, 'fontSize', 12)
set(gcf, 'color', 'w')
