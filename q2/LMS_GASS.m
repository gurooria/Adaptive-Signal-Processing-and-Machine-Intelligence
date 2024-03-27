function [xPredicted, weight, error, stepSize] = LMS_GASS(noise, toPredict, step, rho, alpha, order, mode)
    xPredicted = zeros(length(noise), 1);
    error = xPredicted;
    weight = zeros(order, length(noise)+1);
    lag = zeros(order, length(noise));
    psi = weight;
    leak = 0;
    
    if mode < 3 % normal LMS
        stepSize = step * ones(1, length(noise)+1);
    else
        stepSize = zeros(1, length(noise)+1);
        stepSize(1) = step;
    end

    for i = 1: order
        lag(i,:) = [zeros(1, i), noise(1:length(noise)-i)]; 
    end
    
    for j = 1: length(noise)
        xPredicted(j) = weight(:, j)' * lag(:, j);
        error(j) = toPredict(j) - xPredicted(j);
        weight(:,j+1) = (1 - step*leak).*weight(:,j) + stepSize(j)*error(j)*lag(:,j);        
        
        switch mode % GASS
            case 1
                % nothing
            case 2
                % nothing
            case 3 % benveniste
                stepSize(j+1) = stepSize(j) + rho * error(j) * lag(:, j)' * psi(:, j);
                psi(:, j+1) = (eye(order, order) - (stepSize(j) * lag(:, j)' * (lag(:,j)))) * psi(:, j) + error(j) * lag(:, j);
            case 4 % ang and farhang
                stepSize(j+1) = stepSize(j) + rho * error(j) * lag(:, j)' * psi(:, j);
                psi(:, j+1) = alpha * psi(:, j) + error(j) * lag(:, j);
            otherwise % matthews and xie
                stepSize(j+1)= stepSize(j) + rho * error(j) * lag(:, j)' * psi(:, j);
                psi(:, j+1) = error(j) * lag(:, j);
        end
    end

    weight =  weight(:, 2:end);
    stepSize = stepSize(:, 2:end);
end