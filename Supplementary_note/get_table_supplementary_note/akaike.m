function aicc = akaike(chi, n, k)
%% This function returns the corrected Akaike index if chi is the minimum chi squared, n in the number of data and k the number of parameters

aicc = 2 * k + n * log(chi / n) + (2 * k^2 + 2 * k) / (n - k - 1);

end