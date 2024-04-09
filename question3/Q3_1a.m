%% Q3.1a Complex LMS and Widely Linear Modelling
clear
close all
clc

% Initialisations
N = 1000;
iterations = 100;
mu = 0.1;
var = 1;
b = [1.5 + 1i, 2.5 - 0.5i];
a = 1;
order = length(b);
y = complex(zeros(1, N));

hCLMS = complex(zeros(iterations, order, N));
eCLMS = complex(zeros(iterations, N));
hACLMS = complex(zeros(iterations, order, N));
gACLMS = complex(zeros(iterations, order, N));
eACLMS = complex(zeros(iterations, N));

for i = 1 : iterations
    x = wgn(1, N, pow2db(var), 'complex');
    for j = 2:  N
        y(j) = a*x(j) + b(1)*x(j-1) + b(2)*conj(x(j-1));
    end
    [hCLMS(i,:,:), eCLMS(i,:)] = CLMS(y, x, mu, order);
    [hACLMS(i,:,:), gACLMS(i,:,:), eACLMS(i,:)] = ACLMS(y, x, mu, order);
end

% Plotting
figure
subplot(1,2,1)
hold on
scatter(real(y), imag(y), 12, 'filled')
scatter(real(x),imag(x), 12, 'filled')
xlabel('Real')
ylabel('Imaginary')
grid on
grid minor
set(gca, 'Fontsize', 12)
legend('y - WLMA(1)','x - WGN')
title('Visualisations of the Distribution of x and y', 'Fontsize', 12)
xlim([-10 10])

eCLMS = squeeze(mean(abs(eCLMS).^2, 1));
eACLMS = squeeze(mean(abs(eACLMS).^2, 1));

subplot(1,2,2)
hold on
plot(pow2db(eCLMS), 'linewidth', 1.2)
plot(pow2db(eACLMS), 'linewidth', 1.2)
xlabel('Sample Number n')
ylabel('Square Error (dB)')
set(gca, 'Fontsize', 12)
legend('CLMS', 'ACLMS')
title('Averaged Learning Curves of CLMS and ACLMS', 'Fontsize', 12)
grid on
grid minor