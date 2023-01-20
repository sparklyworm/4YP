function [mu, cov] = posteriorGP(x_test, x_train, y_train, mu_test, mu_train, sigma_noise, sigma_signal, lengthScale)
% mu and cov are parameters of this posterior distribution P(y_test| x_test, x_train, y_train)
% x_test is test data
% x_train and y_train are training data
    
    if nargin == 6
        noisy = true;
    else
        noisy = false;
    end
   
    K_test = covRBF(x_test, [], sigma_signal, lengthScale);
    K_test_train = covRBF(x_test, x_train, sigma_signal, lengthScale);
    K_train_test = K_test_train';
    K_train = covRBF(x_train, [], sigma_signal, lengthScale);
    K_train = addJitter(K_train); % to ensure positive semidefiniteness (for good inverse)

    if noisy
        K_train = K_train + sigma_noise^2 * eye(length(x_train));
    end

    mu = mu_test + K_test_train * (K_train \ (y_train - mu_train));
    cov = K_test - K_test_train * (K_train \ K_train_test);



end