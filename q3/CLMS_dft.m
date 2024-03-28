function [h, error] = CLMS_dft(xInput, y, mu, gamma, L)
    h = complex(zeros(L, length(y)));
    error = complex(zeros(1, length(y)));

    for k = 1: length(xInput)
        xOutput = h(:, k)' * xInput(:, k);
        error(k) = y(k) - xOutput;
        h(:, k+1) = (1 - gamma*mu) * h(:, k) + mu * conj(error(k)) * xInput(:, k);
    end
    h =  h(:, 2:end);
end