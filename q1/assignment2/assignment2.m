%%2.1

x = randn(1, 1000);
xcor = xcorr(x, 'unbiased');
x_axis = [-999:999];

figure(1)
stem(x_axis, xcor, 'lineWidth', 0.8);
axis([-50 50 -0.8 1.3])
grid minor;
zoom xon;
xlabel('Time Lag \tau', fontSize=15); ylabel('Autocorrelation', fontSize=15)
title('Stem plot of the unbiased estimate of ACF for x', fontSize=15);

%%
%2.1.4

x = randn(1,1000);
y = filter(ones(4,1),[1],x);
ycor = xcorr(y, 'unbiased');
x_axis = [-999:999];

figure(2)
stem(x_axis, ycor,'lineWidth', 0.8);
axis([-20 20 -2 10]);
grid minor;
zoom xon;
xlabel('Time Lag \tau', fontSize=15); ylabel('Autocorrelation', fontSize=15)
title('Unbiased estimate of ACF for filtered samples z, using MA filter of order 4', fontsize=12);

%%
%2.2.1
y = filter(ones(9,1),[1],x);
xycorr = xcorr(x,y,'unbiased');

figure(3)
stem(x_axis, xycorr,'lineWidth', 0.8);
axis([-20 20 -0.5 1.5]);
grid minor;
zoom xon;
xlabel('Time Lag \tau', fontSize=15); ylabel('Cross-correlation', fontSize=15)
title('Unbiased estimate of CCF between x and y', fontsize=15);

%%
%2.2.2
y = filter(ones(4,1),[1],x);
xycorr = xcorr(x,y,'unbiased');

figure(3)
stem(x_axis, xycorr,'lineWidth', 0.8);
axis([-20 20 -0.5 1.5]);
grid minor;
zoom xon;
xlabel('Time Lag \tau', fontSize=15); ylabel('Cross-correlation', fontSize=15)
title('Unbiased estimate of CCF between x and z', fontsize=15);

%%
%2.3.1

a1 = -2.5 + 5.*rand(100,1);
a2 = -1.5 + 3.*rand(100,1);

for i= 1:100
     if (a2(i) + a1(i) < 1) && (a2(i) - a1(i) < 1) && abs(a2(i)) < 1
     hold on
     plot(a1(i),a2(i),'*b', 'display','off');
     end
end

title('Stable (a1, a2) coefficient pairs', fontSize=15)
xlabel('Coefficient a1', fontSize=15);
ylabel('Cofficient a2', fontSize=15);
grid minor;
axis([-2.5 2.5 -1.5 1.5])

%%
%2.3.2
load sunspot.dat
N=5;
[a,b] = xcorr(sunspot(1:N,2),'unbiased');
zero_mean = (sunspot(1:N,2)) - mean(sunspot(1:N,2));
[c,d] = xcorr(zero_mean,'unbiased');
figure(1)
stem(b,a, 'lineWidth', 0.75);
hold on
stem(d,c, 'lineWidth', 0.75);
grid minor;
xlabel('Year Lag', fontSize=15);
ylabel('Autocorrelation', fontSize=15);
title(['ACF estimate of Sunspot Data for N=' num2str(N)], fontSize=15);
axis([-N N -400 600]);
legend('Non-centered','Zero-mean');

%%
%2.3.3
load sunspot.dat
p=10; % Order 
sun = sunspot(:,2); % Given that sunspot data is in the second column.
sun_norm = zscore(sun); % Standardizes values to 0 mean and unit variance.

parc_len = 10;
parc = zeros(parc_len,1);
parc_norm = zeros(parc_len,1);

for i = 1:parc_len
    coefs = aryule(sun, i);
    parc(i) = -coefs(i+1);
end

for i = 1:parc_len
    coefs_norm = aryule(sun_norm, i);
    parc_norm(i) = -coefs_norm(i+1);
end

% generate x axis
x = linspace(1, parc_len, parc_len);
stem(x,parc,'DisplayName','PCF', 'lineWidth',0.9);
hold on;
stem(x,parc_norm,'DisplayName','Zero-Mean PCF', 'lineWidth',0.9);
grid on;
grid minor;
xlabel('Time Lag \tau', fontSize=15);
ylabel('Partial Correlation', fontSize=15);
title('PCF of sunspot time series', fontSize=15);

line([0,10], [0.095, 0.095],...
     'Color', 'r',...
     'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Upper 95% bound');
 
line([0,10], [-0.095, -0.095],...
 'Color', 'r',...
 'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Lower 95% bound');

legend('show');

%%
%2.3.4
clc
close all
clear
load sunspot.dat

p=10;
sun = sunspot(:,2);
sun_n = zscore(sun); % standardise
N = length(sun_n);
for order=1:p
    [ar_coeffs,loss_fuction] = aryule(sun_n,order);
    MDL(order) = log(loss_fuction) + (order/N)*log(N);
    AIC(order) = log(loss_fuction) + (2*order)/N;
    AIC_c(order) = AIC(order) + 2*order*(order+1)/(N-order-1);
end

figure (4)
title('AR Model Order Selection for sunspot series', fontSize=15)
xlabel('Model Order p',fontSize=15)
ylabel('Criteria Value', fontSize=15)
grid minor;
axis([1 10 -1.9 -0.9])
hold on
plot(MDL, 'lineWidth', 0.8, 'DisplayName','MDL')
plot(AIC, 'lineWidth', 0.8, 'DisplayName','AIC')
plot(AIC_c, 'lineWidth', 0.8, 'DisplayName','AIC_c')
legend('show')

%%
%2.3.5
close all;
clear;
clc;

load sunspot.dat
sun=sunspot(:, 2); % selects sunspot data
sun = zscore(sun); % standardise
N = 150; % sample N

% order 1
pred_11 = zeros(N,1);
a = ar(sun(1:N), 1, 'yw'); % Polynomial model with identifiable parameters
pred_11 = predict(a, sun(1:N), 1); % m = 1

pred_12 = zeros(N,1);
b = ar(sun(1:N), 1, 'yw');
pred_12 = predict(b, sun(1:N), 2); % m = 2

pred_15 = zeros(N,1);
c = ar(sun(1:N), 1, 'yw');
pred_15 = predict(c, sun(1:N), 5); % m = 5

pred_110 = zeros(N,1);
d = ar(sun(1:N), 1, 'yw');
pred_110 = predict(d, sun(1:N), 10); % m =10

% order 2
pred_21 = zeros(N,1);
e = ar(sun(1:N), 2, 'yw');
pred_21 = predict(e, sun(1:N), 1);

pred_22 = zeros(N,1);
f = ar(sun(1:N), 2, 'yw');
pred_22 = predict(f, sun(1:N), 2);

pred_25 = zeros(N,1);
g = ar(sun(1:N), 2, 'yw');
pred_25 = predict(g, sun(1:N), 5);

pred_210 = zeros(N,1);
h = ar(sun(1:N), 2, 'yw');
pred_210 = predict(h, sun(1:N), 10);

% order 10

pred_101 = zeros(N,1);
i = ar(sun(1:N), 10, 'yw');
pred_101 = predict(i, sun(1:N), 1);

pred_102 = zeros(N,1);
j = ar(sun(1:N), 10, 'yw');
pred_102 = predict(j, sun(1:N), 2);

pred_105 = zeros(N,1);
k = ar(sun(1:N), 10, 'yw');
pred_105 = predict(k, sun(1:N), 5);

pred_1010 = zeros(N,1);
l = ar(sun(1:N), 10, 'yw');
pred_1010 = predict(l, sun(1:N), 10);

% plotting
x = linspace(1, N, N);

subplot(3,4,1)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_11, 'lineWidth', 0.75);
title('AR(1), m=1');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,2)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_12, 'lineWidth', 0.75);
title('AR(1), m=2');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,3)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_15, 'lineWidth', 0.75);
title('AR(1), m=5');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,4)
plot(x, sun(1:N),'DisplayName','Real Data', 'lineWidth', 0.75);
hold on
plot(x, pred_110,'DisplayName','Predicted', 'lineWidth', 0.75);
title('AR(1), m=10');
xlabel('Sample number N');
ylabel('Signal Values');
legend('show')
grid minor

subplot(3,4,5)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_21, 'lineWidth', 0.75);
title('AR(2), m=1');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,6)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_22, 'lineWidth', 0.75);
title('AR(2), m=2');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,7)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_25, 'lineWidth', 0.75);
title('AR(2), m=5');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,8)
plot(x, sun(1:N),'DisplayName','Real Data', 'lineWidth', 0.75);
hold on
plot(x, pred_210,'DisplayName','Predicted', 'lineWidth', 0.75);
title('AR(2), m=10');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,9)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_101, 'lineWidth', 0.75);
title('AR(10), m=1');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,10)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_102, 'lineWidth', 0.75);
title('AR(10), m=2');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,11)
plot(x, sun(1:N), 'lineWidth', 0.75);
hold on
plot(x, pred_105, 'lineWidth', 0.75);
title('AR(10), m=5');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

subplot(3,4,12)
plot(x, sun(1:N),'DisplayName','Real Data', 'lineWidth', 0.75);
hold on
plot(x, pred_1010,'DisplayName','Predicted', 'lineWidth', 0.75);
title('AR(10), m=10');
xlabel('Sample number N');
ylabel('Signal Values');
grid minor

legend('FontSize',11)

%% 2.4.1
clear all
close all
clc

load NASDAQ.mat
prices = NASDAQ.Close;
prices_norm = zscore(prices);
time = NASDAQ.Date;

p=10; % Order
N = length(prices_norm);
x = zeros(N,10);

for order=1:p
    [ar_coeffs,loss_fuction] = aryule(prices_norm,order);
    MDL(order) = log(loss_fuction) + (order/N)*log(N);
    AIC(order) = log(loss_fuction) + (2*order)/N;
    AIC_c(order) = AIC(order) + 2*order*(order+1)/(N-order-1);
end

figure(1)
hold on;
plot(AIC,'DisplayName','AIC', 'lineWidth', 0.8);
plot(MDL,'DisplayName','MDL', 'lineWidth', 0.8);
grid minor;
xlabel('Model Order p', fontSize=15);
ylabel('Criteria Value', fontSize=15);
title(['AIC and MDL for AR model of order p = 0 to ' num2str(p)], fontSize=15);
legend('show');

p=10; % Order 

parc_len = 10;
parc = zeros(parc_len,1);
parc_norm = zeros(parc_len,1);

for i = 1:parc_len
    coefs = aryule(prices, i);
    parc(i) = -coefs(i+1);
end

for i = 1:parc_len
    coefs_norm = aryule(prices_norm, i);
    parc_norm(i) = -coefs_norm(i+1);
end

% plot & x axis
x = linspace(1, parc_len, parc_len);

figure(2)
stem(x,parc,'DisplayName','Empirical PCF','lineWidth',0.9);
hold on;
stem(x,parc_norm,'DisplayName','Zero Mean PCF','lineWidth',0.9);
line([0,10], [0.07, 0.07],...
     'Color', 'r',...
     'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Upper 95% bound');
 
line([0,10], [-0.07, -0.07],...
 'Color', 'r',...
 'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Lower 95% bound');
grid on;
grid minor;
xlabel('Time Lag \tau', fontSize=15);
ylabel('Partial Correlation', fontSize=15);
title('PCF of NASDAQ Data', fontSize=15);
legend('show');

%% 2.4.1 c)

clear all
close all
clc

n = 1:50:1001;
var = 1:50:1001;
load NASDAQ.mat
prices = NASDAQ.Close;
r_xx = xcorr(zscore(prices),'unbiased');

[N,VAR] = meshgrid(n,var);
crlb_var = 2*(VAR.^2)./N;
crlb_a = (VAR)./(N*r_xx(924));

figure(1)
h = heatmap(n, var, crlb_var);
h.Colormap = parula;
h.ColorScaling = 'log';
title('CRLB Heatmap for estimated variance of driving noise')
xlabel('N')
ylabel('True variance of driving noise')
set(struct(h).Axes.Title,'FontSize',13)
set(struct(h).Axes.XLabel,'FontSize',13)
set(struct(h).Axes.YLabel,'FontSize',13)

figure(2)
h = heatmap(n, var, crlb_a);
h.Colormap = parula;
h.ColorScaling = 'log';
title('CRLB Heatmap for a1');
xlabel('N')
ylabel('True variance of driving noise')
set(struct(h).Axes.Title,'FontSize',13)
set(struct(h).Axes.XLabel,'FontSize',13)
set(struct(h).Axes.YLabel,'FontSize',13)

%% 2.5

clc
clear
close all

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
[RRI_trial1,RRI_fs1] = ECG_to_RRI(Trial_1,fs);
%%
[RRI_trial2,RRI_fs2] = ECG_to_RRI(Trial_2,fs);
%%
[RRI_trial3,RRI_fs3] = ECG_to_RRI(Trial_3,fs);

%% 2.5 averaged h[n]
alpha1 = 1;
alpha2 = 0.6;
hd = RRI_trial1(1:980);
hr = 60./hd;
nbins = 15;

h_avg = reshape(hr,[10 98]);

h_avg_a = zeros(1,size(h_avg,2));
for i = 1:size(h_avg_a,2)
    h_avg_a(i) = (1/10)*alpha1*sum(h_avg(:,i));
end

h_avg_alpha2 = zeros(1,size(h_avg,2));
for i = 1:size(h_avg_alpha2,2)
    h_avg_alpha2(i) = (1/10)*alpha2*sum(h_avg(:,i));
end

subplot(1,3,1)
histogram(hr, nbins, 'Normalization', 'pdf',...
          'FaceColor', [0.4660 0.6740 0.1880]);
xlabel('Heart Rate', fontSize=13)
ylabel('Probability', fontSize=13)
title('PDE of original heart rates', fontSize=13)
grid on

subplot(1,3,2)
histogram(h_avg_a, nbins, 'Normalization', 'pdf',...
          'FaceColor', [0.4660 0.6740 0.1880]);
title(['PDE of averaged heart rates for \alpha = ' num2str(alpha1)], fontSize=13);
ylabel('Probability', fontSize=13)
xlabel('Heart Rate', fontSize=13)
grid on

subplot(1,3,3)
histogram(h_avg_alpha2, nbins, 'Normalization', 'pdf',...
          'FaceColor', [0.4660 0.6740 0.1880]);
title(['PDE of averaged heart rates for \alpha = ' num2str(alpha2)], fontSize=13)
ylabel('Probability', fontSize=13)
xlabel('Heart Rate', fontSize=13)
grid on

%% 2.5 acf

h_data_1(:, 1) = RRI_trial1;
h_data_1(:, 1) = detrend(h_data_1(:, 1));
corr_1(:, 1) = xcorr(h_data_1(:, 1), 'biased');

h_data_2(:, 1) = RRI_trial2;
h_data_2(:, 1) = detrend(h_data_2(:, 1));
corr_2(:, 1) = xcorr(h_data_2(:, 1), 'biased');

h_data_3(:, 1) = RRI_trial3;
h_data_3(:, 1) = detrend(h_data_3(:, 1));
corr_3(:, 1) = xcorr(h_data_3(:, 1), 'biased');

subplot(1,3,1)
N = length(h_data_1(:, 1));
x = linspace(-N+1, N-1, 2*N-1);
plot(x, corr_1(:, 1), 'lineWidth', 0.8);
axis([-986 986 -0.0025 0.0075]);
grid minor
title('ACF of Trial 1 RRI', fontSize=15);
xlabel('Time Lag \tau', fontSize=15);
ylabel('Autocorrelation', fontSize=15);

subplot(1,3,2)
N = length(h_data_2(:, 1));
x = linspace(-N+1, N-1, 2*N-1);
plot(x, corr_2(:, 1), 'lineWidth', 0.8);
axis([-877 877 -0.001 0.0025]);
grid minor
title('ACF of Trial 2 RRI', fontSize=15);
xlabel('Time Lag \tau', fontSize=15);
ylabel('Autocorrelation', fontSize=15);

subplot(1,3,3)
N = length(h_data_3(:, 1));
x = linspace(-N+1, N-1, 2*N-1);
plot(x, corr_3(:, 1), 'lineWidth', 0.8);
axis([-874 874 -0.009 0.012]);
grid minor
title('ACF of Trial 3 RRI', fontSize=15);
xlabel('Time Lag \tau', fontSize=15);
ylabel('Autocorrelation', fontSize=15);

%% 2.5 d

h_data_n1 = zscore(h_data_1(:,1));
h_data_n2 = zscore(h_data_2(:,1));
h_data_n3 = zscore(h_data_3(:,1));

p=10; % Order 
parc_len = 10;
parc1 = zeros(parc_len,1);
parc_norm1 = zeros(parc_len,1);
parc2 = zeros(parc_len,1);
parc_norm2 = zeros(parc_len,1);
parc3 = zeros(parc_len,1);
parc_norm3 = zeros(parc_len,1);

for i = 1:parc_len
    coefs1 = aryule(h_data_1, i);
    parc1(i) = -coefs1(i+1);
end
for i = 1:parc_len
    coefs_norm1 = aryule(h_data_n1, i);
    parc_norm1(i) = -coefs_norm1(i+1);
end

for i = 1:parc_len
    coefs2 = aryule(h_data_2, i);
    parc2(i) = -coefs2(i+1);
end
for i = 1:parc_len
    coefs_norm2 = aryule(h_data_n2, i);
    parc_norm2(i) = -coefs_norm2(i+1);
end

for i = 1:parc_len
    coefs3 = aryule(h_data_3, i);
    parc3(i) = -coefs3(i+1);
end
for i = 1:parc_len
    coefs_norm3 = aryule(h_data_n3, i);
    parc_norm3(i) = -coefs_norm3(i+1);
end

% trial 1
subplot(2,3,1)
x = linspace(1, parc_len, parc_len);
hold on;
stem(x,parc_norm1,'DisplayName','Zero Mean PCF', 'lineWidth',0.9);
line([0,10], [0.07, 0.07],...
     'Color', 'r',...
     'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Upper 95% bound');
 line([0,10], [-0.07, -0.07],...
     'Color', 'r',...
     'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Lower 95% bound');
grid minor;
xlabel('Time Lag \tau', fontSize=15);
ylabel('Partial Correlation', fontSize=15);
title('PCF of RRI 1', fontSize=15);
legend('show');

% trial 2
subplot(2,3,2)
x = linspace(1, parc_len, parc_len);
hold on;
stem(x,parc_norm2,'DisplayName','Zero Mean PCF', 'lineWidth',0.9);
line([0,10], [0.07, 0.07],...
     'Color', 'r',...
     'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Upper 95% bound');
 line([0,10], [-0.07, -0.07],...
     'Color', 'r',...
     'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Lower 95% bound');
grid minor;
xlabel('Time Lag \tau', fontSize=15);
ylabel('Partial Correlation', fontSize=15);
title('PCF of RRI 2', fontSize=15);
legend('show');

% trial 3
subplot(2,3,3)
x = linspace(1, parc_len, parc_len);
hold on;
stem(x,parc_norm3,'DisplayName','Zero Mean PCF', 'lineWidth',0.9);
line([0,10], [0.07, 0.07],...
     'Color', 'r',...
     'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Upper 95% bound');
 line([0,10], [-0.07, -0.07],...
     'Color', 'r',...
     'LineWidth', 1, 'LineStyle', '--', 'DisplayName','Lower 95% bound');
grid minor;
xlabel('Time Lag \tau', fontSize=15);
ylabel('Partial Correlation', fontSize=15);
title('PCF of RRI 3', fontSize=15);
legend('show');

% 2.5 MDL & AIC & AICc
p=10; % Order
N1 = length(h_data_n1);
N2 = length(h_data_n2);
N3 = length(h_data_n3);

for order=1:p
    [ar_coeffs1,loss_fuction1] = aryule(h_data_n1,order);
    MDL1(order) = log(loss_fuction1) + (order/N1)*log(N1);
    AIC1(order) = log(loss_fuction1) + (2*order)/N1;
    AIC_c1(order) = AIC1(order) + 2*order*(order+1)/(N1-order-1);
end

for order=1:p
    [ar_coeffs2,loss_fuction2] = aryule(h_data_n2,order);
    MDL2(order) = log(loss_fuction2) + (order/N2)*log(N2);
    AIC2(order) = log(loss_fuction2) + (2*order)/N2;
    AIC_c2(order) = AIC2(order) + 2*order*(order+1)/(N2-order-1);
end

for order=1:p
    [ar_coeffs3,loss_fuction3] = aryule(h_data_n3,order);
    MDL3(order) = log(loss_fuction3) + (order/N3)*log(N3);
    AIC3(order) = log(loss_fuction3) + (2*order)/N3;
    AIC_c3(order) = AIC3(order) + 2*order*(order+1)/(N3-order-1);
end

subplot(2,3,4)
hold on;
plot(AIC1,'DisplayName','AIC', 'lineWidth', 0.8);
plot(MDL1,'DisplayName','MDL', 'lineWidth', 0.8);
plot(AIC_c1,'DisplayName','AICc', 'lineWidth', 0.8);
grid minor;
xlabel('Model Order p', fontSize=15);
ylabel('Criteria Value', fontSize=15);
title('AIC, MDL & AICc for RRI 1', fontSize=15);
legend('show');

subplot(2,3,5)
hold on;
plot(AIC2,'DisplayName','AIC', 'lineWidth', 0.8);
plot(MDL2,'DisplayName','MDL', 'lineWidth', 0.8);
plot(AIC_c2,'DisplayName','AICc', 'lineWidth', 0.8);
grid minor;
xlabel('Model Order p', fontSize=15);
ylabel('Criteria Value', fontSize=15);
title('AIC, MDL & AICc for RRI 2', fontSize=15);
legend('show');

subplot(2,3,6)
hold on;
plot(AIC3,'DisplayName','AIC', 'lineWidth', 0.8);
plot(MDL3,'DisplayName','MDL', 'lineWidth', 0.8);
plot(AIC_c3,'DisplayName','AICc', 'lineWidth', 0.8);
grid minor;
xlabel('Model Order p', fontSize=15);
ylabel('Criteria Value', fontSize=15);
title('AIC, MDL & AICc for RRI 3', fontSize=15);
legend('show');

