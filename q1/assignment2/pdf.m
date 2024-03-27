function [y] = pdf(v,bins)
    [freq,x] = hist(v,bins); % frequency & center of each bin
    freq = freq/trapz(x,freq); % normalises the pdf
    bar(x,freq);
end
