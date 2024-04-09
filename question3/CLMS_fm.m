function [h, error] = CLMS_fm(xInput, y, mu, order)
    N = length(xInput);
    h = complex(zeros(order,N));
    error = complex(zeros(1,N));
    lag = zeros(N+order-1,1);
    lag(order:N+order-1) = xInput;

    for k = 1: length(xInput)
        xOutput(k) = h(:, k)' * lag(k:k+order-1);
        error(k) = y(k) - xOutput(k);
        h(:, k+1) = h(:, k) + mu * conj(error(k)) * lag(k:k+order-1);
    end
    h =  h(:, 2:end);
end