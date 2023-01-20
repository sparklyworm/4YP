function f = normalpdf(x, mu, sigma)
    % x is column vector
    f = 1/ (sigma*sqrt(2*pi)) * exp(-0.5 * ((x - mu)/sigma).^2);

end