function [xPredicted, noiseEstimation, weight] = LMS_ANC(xInput, ref, step, leak,order)
    noiseEstimation = zeros(length(xInput), 1);
    xPredicted = noiseEstimation;
    weight = zeros(order, length(xInput)+1);
    lag = zeros(order, length(xInput));
    
    for i = 1 : order
        lag(i, :) = [zeros(1, i-1), ref(1 : length(xInput)-(i-1))]; 
    end
    
    for k = 1 : length(xInput)
        noiseEstimation(k) = weight(:, k)' * lag(:, k);
        xPredicted(k) = xInput(k) - noiseEstimation(k);
        weight(:, k+1) = (1-step*leak) .* weight(:, k) + step * xPredicted(k) * lag(:, k);
    end
    weight = weight(:, 2:end);
end