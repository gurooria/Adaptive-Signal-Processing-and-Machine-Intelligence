function [xPredicted, weight, error] = GNGD(noise, toPredict, step, epsilonStart, order, rho)
    xPredicted = zeros(length(noise), 1);
    error = xPredicted;
    weight = zeros(order, length(noise)+1);
    lag = zeros(order, length(noise)); 
    epsilon = (1/step) * ones(length(noise)+1, 1);
    epsilon(1) = epsilonStart;
    beta=1;

    for i = 1 : order
        lag(i, :) = [zeros(1,i), noise(1 : length(noise)-i)]; 
    end

    for j = 1 : length(noise)
        xPredicted(j) = weight(:, j)' * lag(:, j);
        error(j) = toPredict(j) - xPredicted(j);
        weight(:, j+1) = weight(:, j) + ((beta * error(j)) / (epsilon(j) + (lag(:, j))' * lag(:, j))) * lag(:, j);
        if j > 1
            epsilon(j+1) = epsilon(j) - (rho*step) * ((error(j) * error(j-1) * (lag(:, j)') ...
                * lag(:, j-1)) / ((epsilon(j-1) + (lag(:, j-1)') * lag(:, j-1)))^2);
        end
    end
    weight = weight(:, 2:end);
    epsilon = epsilon(:, 2:end);
end
