clc, clear
% some data to setup 
N = 100;
x_domain = 10;
x_data = linspace(0, x_domain, N)';
y_data = sin(2* pi* 1/x_domain * 2 * x_data) + cos (4/x_domain * x_data) ;
y_data = flip(y_data);

% hyperparameters for GP model and prior mean 0 
sigma_noise = 0.5; % noise variance 
sigma_signal = 1; % signal variance (how fast changing/squiggly do we want our curve to be)
lengthScale = 1; % bandwidth (how much covariance)
mu_test = zeros(N, 1);
%% plot
figure
plot(x_data, y_data, 'DisplayName', 'True Data'); 

%% GP Model
index_obs = [50];
x_obs = x_data(index_obs);
y_obs = y_data(index_obs);
mu_train = zeros(length(y_obs), 1);

[mu, cov] = posteriorGP(x_data, x_obs, y_obs, mu_test, mu_train, sigma_noise, sigma_signal, lengthScale);
%% plot
hold on 
plot(x_obs, y_obs, 'rx')
hold on 
plot(x_data, mu, 'DisplayName', 'GP Model')
legend 
%% acquisition 
% calculate value of the acquisition function at each x 
best_so_far = max(mu);
epsilon = 0.03; % how much exploration
variance = diag(cov);
A_PI =  acquisitionPI(best_so_far + epsilon, mu, x_data, variance);
A_EI = acquisitionEI(best_so_far + epsilon, mu, x_data, variance);
[next_PI, index_PI] = max(A_PI);
[next_EI, index_EI] = max(A_EI);
%% plot
figure
hold on
plot(x_data, A_PI, 'g')
hold on
plot(x_data(index_PI), next_PI, 'go', 'DisplayName','Probability of Improvement')
hold on 
plot(x_data, A_EI, 'm')
hold on
plot(x_data(index_EI), next_EI, 'mo', 'DisplayName', 'Expected Improvement')
legend 

%% step through each loop to visualize 
index_obs = [];
next_index = 2;
max_iter = 20;
for i = 1:max_iter

    % update GP model
    index_obs = [index_obs next_index];
    x_obs = x_data(index_obs);
    y_obs = y_data(index_obs);
    mu_train = zeros(length(y_obs), 1);
    [mu, cov] = posteriorGP(x_data, x_obs, y_obs, mu_test, mu_train, sigma_noise, sigma_signal, lengthScale);
    
    % acquisition
    best_so_far = max(mu);
    epsilon = 0; % how much exploration
    variance = diag(cov);
    A_EI = acquisitionEI(best_so_far + epsilon, mu, x_data, variance);
    [next_EI, next_index] = max(A_EI);

    % plot GP and data
    figure
    subplot(2, 1, 1)
    plot(x_data, y_data, 'DisplayName', 'True data'); % true curve
    hold on 
    plot(x_obs, y_obs, 'rx', 'DisplayName', 'Observed points') % observed points
    hold on 
    plot(x_data, mu, 'DisplayName','GP model') % posterior GP
    legend 

    subplot(2, 1, 2)
    plot(x_data, A_EI, 'm')
    hold on
    plot(x_data(next_index), next_EI, 'mo', 'DisplayName', 'Expected improvement')
    
end
