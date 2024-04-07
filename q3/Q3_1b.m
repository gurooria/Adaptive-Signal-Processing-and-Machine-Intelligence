%% Q3.1b Bivariate Wind Data
clc
clear
close all
lowWindDat = load('low-wind.mat');
medWindDat = load('medium-wind.mat');
highWindDat = load('high-wind.mat');

%% Initialisations
v_low = complex(lowWindDat.v_east,lowWindDat.v_north);
v_med = complex(medWindDat.v_east,medWindDat.v_north);
v_high = complex(highWindDat.v_east,highWindDat.v_north);

% circularity
cirLow = abs(mean((v_low).^2)/mean(abs(v_low).^2));
cirMed = abs(mean((v_med).^2)/mean(abs(v_med).^2));
cirHigh =  abs(mean((v_high).^2)/mean(abs(v_high).^2));

% centre of mass
centerLow = [mean(real(v_low)),mean(imag(v_low))];
centreMed = [mean(real(v_med)),mean(imag(v_med))];
centreHigh = [mean(real(v_high)),mean(imag(v_high))];

%% Plotting 
figure

% low speed
subplot(2,3,1)
scatter(real(v_low), imag(v_low), 12, 'filled')
hold on
scatter(centerLow(1), centerLow(2), 20, 'filled', 'r')
xlabel("Real (East)", 'fontsize', 12)
ylabel("Imaginary (North)",'fontsize',12)
title(sprintf('Low-speed Data with ρ = %0.2f',cirLow), 'fontsize', 12)
ax = gca;
ax.FontSize = 12;
grid on
grid minor

% med speed
subplot(2,3,2)
scatter(real(v_med), imag(v_med) ,12, 'filled')
hold on
scatter(centreMed(1), centreMed(2), 20, 'filled', 'r')
xlabel("Real (East)", 'fontsize', 12)
ylabel("Imaginary (North)",'fontsize',12)
title(sprintf('Medium-speed Data with ρ = %0.2f',cirMed), 'fontsize', 12)
ax = gca;
ax.FontSize = 12;
grid on
grid minor

% high speed
subplot(2,3,3)
scatter(real(v_high), imag(v_high), 10, 'filled')
hold on
scatter(centreHigh(1), centreHigh(2), 20, 'filled', 'r')
xlabel("Real (East)", 'fontsize', 12)
ylabel("Imaginary (North)",'fontsize',12)
title(sprintf('High-speed Data with ρ = %0.2f', cirHigh), 'fontsize', 12)
ax = gca;
ax.FontSize = 12;
grid on
grid minor
set(gcf,'color','w')
%% CLMS and ACLMS
mu = [0.1; 0.01; 0.001];
N = length(v_high);
filterLength = 25;
allSpeeds = [v_low, v_med, v_high];

hold on
for j = 1:3
    errorCLMS = complex(zeros(length(filterLength), N));
    errorACLMS = errorCLMS;
    speed = allSpeeds(:,j);
    
    for i = 1:filterLength
        x = [0; speed(1:end-1)]';
        [~, errorCLMS(i, :)] = CLMS_wind(x, speed', mu(j), i);
        [~, ~, errorACLMS(i, :)] = ACLMS_wind(x, speed', mu(j), i);
    end
    
    eCLMS = squeeze(mean(abs(errorCLMS).^2, 2));
    eACLMS = squeeze(mean(abs(errorACLMS).^2, 2));
    
    subplot(2,3,j+3)
    plot(10*log10(eCLMS), 'LineWidth', 1.2)
    hold on
    plot(10*log10(eACLMS), 'r', 'LineWidth', 1.2)
    if j == 1
    	title('Learning Curves for Low Speed', 'fontsize', 12)
    elseif j == 2
        title('Learning Curves for Medium Speed', 'fontsize', 12)
    else
        title('Learning Curves for High Speed', 'fontsize', 12)
    end
    xlabel('Filter Order')
    ylabel('MSPE (dB)')
    legend('CLMS','ACLMS')
    ax = gca;
    ax.FontSize = 12;
    grid on
    grid minor
    axis tight
end

set(gcf,'color','w')