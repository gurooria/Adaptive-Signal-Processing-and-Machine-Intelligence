function [xPredicted, weight, error] = LMS_ALE(xInput, step, leak, order, delay)
    xPredicted = zeros(length(xInput), 1);
    error = xPredicted;
    weight = zeros(order, length(xInput)+1);
    lag = zeros(order, length(xInput));

    for i = 1 : order
        lag(i, :) = [zeros(1, i+delay-1), xInput(1:length(xInput) - (i+delay-1))]; 
    end
    
    for j = 1 : length(xInput)
        xPredicted(j) = weight(:, j)' * lag(:, j);
        error(j) = xInput(j) - xPredicted(j);
        weight(:, j+1) = (1 - step*leak) .* weight(:, j) + step * error(j) * lag(:, j);
    end
    weight =  weight(:, 2:end);
end