function [h, error] = CLMS_voltage(xInput, mu, order)
    N = length(xInput);
    h = complex(ones(order,N));
    error = complex(zeros(1,N));
    lag = zeros(order,length(xInput)); % x(n-k)
    
    for i = 1 : order
        lag(i, :) = [zeros(1,i), xInput(1: length(xInput)-i)]; 
    end
    
    for k = 1: length(xInput)
        xOutput(k) = h(:, k)' * lag(:, k);
        error(k) = xInput(k) - xOutput(k);
        h(:, k+1) = h(:, k) + mu * conj(error(k)) * lag(:, k);
    end
    h =  h(:, 2:end);
end