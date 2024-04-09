function [g, h, error] = ACLMS_voltage(xInput, mu, order)
    N = length(xInput);
    g = complex(zeros(order,N));
    h = complex(ones(order,N));
    error = complex(zeros(1,N));
    lag = zeros(order, length(xInput)); % x(n-k)

    for i = 1 : order
        lag(i, :) = [zeros(1, i), xInput(1:length(xInput)-i)]; 
    end
    
    for k = 1 : length(xInput)
        xOutput(k) = h(:, k)' * lag(:, k) + g(:, k)' * conj(lag(:, k));
        error(k) = xInput(k) - xOutput(k);
        h(:, k+1) = h(:, k) + mu * conj(error(k)) * lag(:, k);
        g(:, k+1) = g(:, k) + mu * conj(error(k)) * conj(lag(:, k));
    end
    g = g(:, 2:end);
    h =  h(:, 2:end);
end