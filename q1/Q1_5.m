%% Q1.5 Real world signals: Respiratory sinus arrhythmia from RRI
clear
clc
close all

%% Data Processing
RAW = readmatrix('ecg_s2.csv');
data = RAW(:,3);
plot(data)

int1_finish = 123399;
int2_start = 147394;
int2_finish = 249358;
int3_start = 270200;
int3_finish = 379678;
fs = 500; % sampling frequency

Trial_1 = data(1:int1_finish);
Trial_2 = data(int2_start+1:int2_finish);
Trial_3 = data(int3_start+1:int3_finish);

%%
[RRI_trial1, RRI_fs1] = ECG_to_RRI(Trial_1, fs);
%%
[RRI_trial2, RRI_fs2] = ECG_to_RRI(Trial_2, fs);
%%
[RRI_trial3, RRI_fs3] = ECG_to_RRI(Trial_3, fs);

%% Q1.5a & b Initialisations
RRI1 = detrend(RRI_trial1);
RRI2 = detrend(RRI_trial2);
RRI3 = detrend(RRI_trial3);
RRIs = {RRI1; RRI2; RRI3};
fs = 4; % 4 Hz

figure(1)
hold on

for i = 1 : 3 % loop through each trial
    RRI = RRIs{i};
    N = length(RRI);

    [standardPSD, standardX] = pwelch(RRI, hamming(N), 0, 1024*2, 4, 'onesided');
    [pSD150, f150] = pwelch(RRI, hamming(150*fs), 0, 1024, 4, 'onesided');
    [pSD100, f100] = pwelch(RRI, hamming(100*fs), 0, 1024, 4, 'onesided');
    [pSD50, f50] = pwelch(RRI, hamming(50*fs), 0, 1024, 4, 'onesided');

    subplot(3,2,2*i-1)
    plot(standardX, 10*log10(standardPSD), 'LineWidth', 1.2) 
    ax = gca;
    ax.FontSize = 12;
    xlabel('Normalised Frequency (Cycles/Sample)')
    ylabel('Magnitude (dB)')
    title(sprintf('Standard Periodogram of RRI%d', i), 'fontsize', 12)
    
    grid on 
    grid minor
    subplot(3,2,2*i)
    plot(f150, 10*log10(pSD150), 'LineWidth', 1.2)
    hold on
    plot(f100, 10*log10(pSD100), 'LineWidth', 1.2)
    hold on
    plot(f50, 10*log10(pSD50), 'LineWidth', 1.2)
    ax = gca;
    ax.FontSize = 12;
    legend('150 s window','100 s window','50 s window')
    set(gcf,'color','w')
    
    xlabel('Normalised Frequency (Cycles/Sample)')
    ylabel('Magnitude (dB)')
    title(sprintf('Averaged Periodogram of RRI%d', i), 'fontsize', 12)
    grid on
    grid minor
end 


%% Q1.5c
figure(2)
hold on
orders = [5, 10, 7]; % order for each trial

for i = 1:3
    RRI = RRIs{i};
    N = length(RRI);
    [standardPSD, standardX] = pwelch(RRI, hamming(N), 0, 1024*2, 4, 'onesided');

    subplot(1,3,i)
    plot(standardX, 10*log10(standardPSD), 'LineWidth', 1.2)
    hold on

    % Test the order of AR models
    for p = orders(i)
        [h, w] = pyulear(RRI, p, 1024, 4);
        plot(w, 10*log10(h), 'r', 'LineWidth', 1.2)
        hold on
    end
    ax = gca;
    ax.FontSize = 12;
    xlabel('Normalised Frequency (Cycles/Sample)')
    ylabel('Magnitude (dB)')
    xlim([0.01 1.99])
    title(sprintf('Autoregressive PSD Estimate of RRI%d',i),'fontsize', 12)
    switch i
        case 1
            legend('Standard', 'AR(5)')
        case 2
            legend('Standard', 'AR(10)')
        case 3
            legend('Standard', 'AR(8)')
    end
    grid on
    grid minor
end 

set(gcf,'color','w')