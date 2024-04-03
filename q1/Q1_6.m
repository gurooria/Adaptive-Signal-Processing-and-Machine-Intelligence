%% Q1.6 Robust Regression
clc
clear
close all
load("PCAPCR.mat")

%% 1.6a
% Obtain the singular value decomposition of X and Xnoise

% SVD of X
[~, sX, ~] = svd(X);
[~, sN, ~] = svd(Xnoise);

figure(1)

subplot(1,3,1)
stem(sum(sX), 'LineWidth', 1.2)
axis tight
ax = gca;
ax.FontSize = 12;
xlabel("Singular Value Index", 'fontsize', 12)
ylabel("Magnitude", 'fontsize', 12)
title("SVD Singular Values of X", 'fontsize', 12)
grid on
grid minor

subplot(1,3,2)
stem(sum(sN), 'LineWidth', 1.2)
axis tight
ax = gca;
ax.FontSize = 12;
xlabel("Singular Value Index", 'fontsize', 12)
ylabel("Magnitude", 'fontsize', 12)
title("SVD Singular Values of X_{noise}", 'fontsize', 12)
grid on 
grid minor

squaredError = sum((sX-sN).^2);
subplot(1,3,3)
stem(squaredError, 'LineWidth', 1.2)
axis tight
ax = gca;
ax.FontSize = 12;
xlabel("Singular Value Index", 'fontsize', 12)
ylabel("Squared Error", 'fontsize', 12)
title('Squared Error between Singular Values of X and X_{noise}', 'fontsize', 12)
grid on
grid minor
set(gcf, 'color', 'w')

%% 1.6b
[uNew, sNew, vNew] = svds(Xnoise, 3);
noiseNew = uNew*sNew*vNew';
error1 = sum((X-Xnoise).^2);
error2 = sum((X-noiseNew).^2);

stem(error1, 'LineWidth', 1.2)
hold on
stem(error2, 'r', 'LineWidth', 1.2)
axis tight
ax = gca;
ax.FontSize = 12;
xlabel('Singular Value Index', 'fontsize', 12)
ylabel('Squared Error', 'fontsize', 12)
title('Squared Error between X and X_{noise} vs. X and low rank approximation of X_{noise}', 'fontsize', 12)
legend('X and $X_{noise}$','X and $\tilde{X}_{noise}$')
grid on 
grid minor
set(gcf,'color','w')

%% 1.6c
% training set
bOls = inv(Xnoise' * Xnoise) * Xnoise' * Y;
bPcr = vNew * inv(sNew) * uNew' * Y;
yOls = Xnoise * bOls;
yPcr = Xnoise * bPcr;
eOls =  sum((Y-yOls).^2);
ePcr =  sum((Y-yPcr).^2);

%test set
yOlsTest = Xtest * bOls;
yPcrTest = Xtest * bPcr;
eOlsTest =  sum((Ytest-yOlsTest).^2);
ePcrTest =  sum((Ytest-yPcrTest).^2);

figure(1)

subplot(1,3,1)
stem(ePcr, 'LineWidth', 1.2)
hold on
stem(eOls, 'color', 'red', 'LineWidth', 1.2)
axis tight
ax = gca;
ax.FontSize = 12;
xlabel("Singular Value Index ", 'fontsize', 12)
ylabel("Squared Error", 'fontsize', 12)
title("In-sample Squared Error", 'fontsize', 12)
legend('PCR', 'OLS')
grid on
grid minor

subplot(1,3,2)
stem(ePcrTest, 'LineWidth', 1.2)
hold on
stem(eOlsTest, 'color', 'red',  'LineWidth', 1.2)
axis tight
ax = gca;
ax.FontSize = 12;
xlabel("Singular Value Index", 'fontsize', 12)
ylabel("Squared Error", 'fontsize', 12)
title("Out-sample Squared Error", 'fontsize', 12)
legend('PCR', 'OLS')
grid on
grid minor

% 1.6d

errorsOls = [];
errorsPcr = [];

for i = 1:10 % ensemble of 10
    
    [yOls, yOlsh] = regval(bOls);
    [yPcr, yPcrh] = regval(bPcr);
    
    errorsOls= [errorsOls; sum((yOls-yOlsh).^2)];
    errorsPcr = [errorsPcr; sum((yPcr-yPcrh).^2)];
    
end

subplot(1,3,3)
stem(mean(errorsPcr), 'LineWidth', 1.2)
hold on
stem(mean(errorsOls), 'color', 'red', 'LineWidth', 1.2)

axis tight
ax = gca;
ax.FontSize = 12;
xlabel('Singular Value Index', 'fontsize', 12)
ylabel('Mean Squared Error', 'fontsize', 12)
title('Mean Squared Error over an Ensemble of Test Data ', 'fontsize', 12)
legend('PCR', 'OLS')
grid on 
grid minor
set(gcf,'color','w')
