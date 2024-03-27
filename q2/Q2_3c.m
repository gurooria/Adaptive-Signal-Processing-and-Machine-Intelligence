%% Q2.3c
clc
clear
close all

sampNo = 1000;
deltas = [3];% only one delay value to be tested
mu = 0.01; 
its = 100;
M = 5; % filter order (minimum one that is tested in part b)
a = 1;
b = [1 0 0.5]; % since we have eta(n) = v(n) + 0.5v(n-2) (+0v(n-1))
% generating the clean component of s(n)
x = sin(0.01*pi.*(0:sampNo-1)); 

figure
hold on

s_all =  zeros(its,sampNo);
xhat_allALE = zeros(its,sampNo);
MSPE_ALE = zeros(its,1);

xhat_allANC = zeros(its,sampNo);
MSPE_ANC = zeros(its,1);

for j = 1: its
    j
    v = randn(1,sampNo); % has unit variance
    % generating the coloured noise component of s(n)
    eta = filter(b,a,v);
    % constructing the input signal to ALE LMS
    s = x + eta;
    u = 1.2*eta+0.1;
    s_all(j,:) = s;

    % ALE
    [xhatALE,wALE,errALE] = ALE_LMS(s,mu,0,M,deltas);
    xhat_allALE(j,:) = xhatALE;
    MSPE_ALE(j) = mean((x(100:end)-xhatALE(100:end)').^2);
    % ANC
    [noiseEst,wANC,xhatANC] = ANC_LMS(s,u,mu,0,M);
    xhat_allANC(j,:) = xhatANC;
    MSPE_ANC(j) = mean((x(100:end)-xhatANC(100:end)').^2);

end
% plot ALE
subplot(1,3,1)
p1 = plot(s_all','b','LineWidth',1.5,'DisplayName','$s(n)$')
hold on
p2 = plot(xhat_allALE','r','LineWidth',1.5,'DisplayName','$\hat{x}(n)$')
hold on
p3 = plot(x,'y','LineWidth',2,'DisplayName','$x(n)$')
MSPE4title = mean(MSPE_ALE);
title(sprintf('MSPE = %0.3f, $\\Delta$ = 3, M =5',MSPE4title),'Interpreter','latex','fontsize',15)
ax = gca;
ax.FontSize = 15;
xlabel('Sample Index (n)','fontsize',15)
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
allps = [p1(1),p2(1),p3];
legend(allps)
grid on
grid minor 
% sgtitle('Empirical Justification of Minimum Delay','fontsize',18)
% plot ALE
subplot(1,3,2)
p1 = plot(s_all','b','LineWidth',1.5,'DisplayName','$s(n)$')
hold on
p2 = plot(xhat_allANC','r','LineWidth',1.5,'DisplayName','$\hat{x}(n)$')
hold on
p3 = plot(x,'y','LineWidth',2,'DisplayName','$x(n)$')
MSPE4title = mean(MSPE_ANC);
title(sprintf('MSPE = %0.3f, $\\Delta$ = 3, M =5',MSPE4title),'Interpreter','latex','fontsize',15)
ax = gca;
ax.FontSize = 15;
xlabel('Sample Index (n)','fontsize',15)
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
allps = [p1(1),p2(1),p3];
legend(allps)
grid on
grid minor 
set(gcf,'color','w')

% ensemble of realisations
ANC_ensemble = mean(xhat_allANC);
ALE_ensemble = mean(xhat_allALE);
subplot(1,3,3)
plot(ALE_ensemble,'b','LineWidth',2);
hold on
plot(ANC_ensemble,'r','LineWidth',2);
hold on
plot(x,'y','LineWidth',2);
title('Ensemble Means: ALE vs ANC','Interpreter','latex','fontsize',15)
ax = gca;
ax.FontSize = 15;
xlabel('Sample Index (n)','fontsize',15)
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
legend('ALE','ALC','$x(n)$')
grid on
grid minor 
set(gcf,'color','w')
