function [xOutput, weight, error, alphas] = LMS_dp(xInput, stepSize, order, alpha, alphaUpdate, weight0, bias)
    % Includes a tanh activation function
    xOutput = zeros(length(xInput), 1);
    error = xOutput;
    weight = zeros(order+bias, length(xInput)+1);
    weight(:, 1) = weight0;
    lag = zeros(order+bias, length(xInput));
    alphas = alpha;

    for i = 1: order
        lag(i, :) = [zeros(1,i), xInput(1:length(xInput)-i)']; 
    end
    
    for k = 1: length(xInput)
        xOutput(k) = alpha * tanh(weight(:,k)' * lag(:,k));
        error(k) = xInput(k) - xOutput(k);
        activationFunction = alpha * (1-(xOutput(k)/alpha)^2);
        weight(:, k+1) = weight(:,k) + (stepSize*activationFunction*error(k)) .* lag(:,k);
        
        if alphaUpdate
            alpha = alpha + 0.3*error(k)*(xOutput(k)/alpha);
            alphas = [alphas, alpha];
        else
            % nothing, alphas remains as alpha
        end
    end
    weight =  weight(:, 2:end);
end