function PGM = pgm(x)
    N = length(x);
    xdft = fft(x);
    PGM = (1/N).*(abs(xdft).^2);
end

