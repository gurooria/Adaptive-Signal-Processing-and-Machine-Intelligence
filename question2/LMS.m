function [xPredicted, w, error] = LMS(xInput, step, leak, order)
    
    % Initialisations
    xPredicted = zeros(length(xInput), 1);
    error = xPredicted;
    w = zeros(order, length(xInput)+1); % since order can be > 1 (here 2)
    xLagged = zeros(order, length(xInput)); % x(n-k)
    
    % creating two vectors with a lag i, i.e. the order length
    for i = 1: order
        xLagged(i,:) = [zeros(1,i), xInput(1: length(xInput)-i)']; 
    end
    
    for j = 1: length(xInput)
        % calculate the prediction
        xPredicted(j) = w(:, j)' * xLagged(:, j);
        error(j) = xInput(j) - xPredicted(j);
        % update weights
        w(:, j+1) = (1 - step*leak) .* w(:, j) + step * error(j) * xLagged(:, j);
    end

    w =  w(:, 2:end);
end