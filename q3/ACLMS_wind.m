function [g, h, error] = ACLMS_wind(xInput, y, mu, order)
    N = length(xInput);
    g = complex(zeros(order,N));
    h = g;
    error = complex(zeros(1,N));
    lag = zeros(N+order-1,1);
    lag(order:N+order-1) = xInput;
    
    for k = 1: length(xInput)
        xOutput(k) = h(:, k)' * lag(k:k+order-1) + g(:,k)' * conj(lag(k:k+order-1));
        error(k) = y(k) - xOutput(k);
        h(:, k+1) = h(:, k) + mu*conj(error(k))*lag(k:k+order-1);
        g(:, k+1) = g(:, k) + mu*conj(error(k))*conj(lag(k:k+order-1));
    end
    g = g(:, 2:end);
    h =  h(:, 2:end);
end