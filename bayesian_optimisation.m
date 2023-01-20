%% bayesian optimisation approach (to find point of max probability)
clc, clear
load("interesting-landscape\interesting-landscape-5.mat")

%% add path
addpath(genpath('machineLearning'))

%% our observed data
selective_pressure = 1;
true_mu = 5 * 10^-8;
numGen = 600;
[model_fitness, ~] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen, selective_pressure);
true_fitness = zeros(numGen, 1);
for i = 1:numGen
    true_fitness(i) = model_fitness(i) + 0.05*(-1 + 2*rand(1));
end
plot(1:numGen, true_fitness);

%% specify our domain 
mu_candidates = linspace(1*10^-8, 1*10^-7, 50);
num_candidates = length(mu_candidates);
mean_prior_gp = zeros(num_candidates, 1);

%% hyperparameters for GP model
sigma_noise = 1; % noise variance 
sigma_signal = 1; % signal variance (the amplitude of cov)
lengthScale = 3 * 10^-8; % how fast changing/squiggly do we want our curve to be
next_index = 1;
index_obs = [];
y_obs = [];
max_iter = 10;
for j = 1:max_iter

    index_obs = [index_obs next_index];
    num_obs =length(index_obs);
    mu_obs_new = mu_candidates(next_index);
    mu_obs = mu_candidates(index_obs);
    
    sse = evaluateSSE(true_fitness, gene_to_fitness, N, mu_obs_new, numGen, selective_pressure);
    y_obs_new = 1/sse ;
    y_obs = [y_obs; y_obs_new];
    
    mean_obs = zeros(num_obs, 1);
    [mean_posterior, cov_posterior] = posteriorGP(mu_candidates, mu_obs, y_obs, mean_prior_gp, mean_obs, sigma_noise, sigma_signal, lengthScale);
   
   
    % acquisition 
    best_so_far = max(mean_posterior);
    epsilon = 0.01; % how much exploration
    variance = diag(cov_posterior);
    exploration = 10;
    A_EI = acquisitionEI_explore(best_so_far + epsilon, mean_posterior, mu_candidates, variance, exploration);
    [next_EI, next_index] = max(A_EI);
    
    % stopping criteria: when next_index is alrd visited 
    if ~isempty(find(index_obs == next_index, 1))
        break
    end

%     figure
%     plot(mu_candidates, mean_posterior)
%     hold on 
%     plot(mu_obs, y_obs, 'rx')
%     hold on 
%     plot(mu_candidates(next_index), 0, 'go')
%     
%     figure
%     plot(1:1:length(mu_candidates), A_EI)

    
end
%%
figure
plot(mu_candidates, mean_posterior)
hold on 
plot(mu_obs, y_obs, 'rx')
hold on 
plot(mu_candidates(next_index), 0, 'go')

figure
plot(1:1:length(mu_candidates), A_EI)

%% normalise with area under curve so that it is a probability distribution 
% first find where the last non negative index is 
last_non_neg_ind = find(mean_posterior >=0, 1, 'last');
first_non_neg_ind = find(mean_posterior >= 0, 1);

y_values = mean_posterior(first_non_neg_ind:last_non_neg_ind);
x_domain = mu_candidates(first_non_neg_ind:last_non_neg_ind);
area_under_curve = trapz(x_domain, y_values);
normalised_mean_posterior = max(mean_posterior/area_under_curve,0);

y_pdf_prior = normpdf(mu_candidates, 6*10^-8, 2*10^-8)'; % our Gauss prior belief of distibution of mu values!
mean_post_post = normalised_mean_posterior .* y_pdf_prior; % the posterior
auc_post = trapz(mu_candidates, mean_post_post);
mean_post_post = mean_post_post /auc_post; % normalised posterior 

% some variables to help plot the dotted mean line of posterior
[m, i] = max(mean_post_post);
mu_i = mu_candidates(i);
x_m = [mu_i mu_i];
y_m = [0 m];

figure
plot(mu_candidates, normalised_mean_posterior, 'r', DisplayName='Likelihood (non Gauss)', LineWidth=2)
hold on
plot(mu_candidates, y_pdf_prior, 'b', DisplayName='Prior (Gauss)', LineWidth=2)
hold on
plot(mu_candidates, mean_post_post, 'g', DisplayName='Posterior (non Gauss)', LineWidth=2)
hold on
plot(x_m, y_m, 'g--', 'HandleVisibility','off')
xlabel("Mutation rate")
ylabel("Probability density")
legend boxoff

label_mean = sprintf("%.5g", mu_i);
word = '\mu = ' + label_mean;
text(mu_i, 1*10^6, word, 'FontSize', 12)



%% option 1: fit a gaussian pdf to our mean_posterior so that our curve can be decribed by parameters, a mean and a variance
[~, mean_mle_ind] = max(normalised_mean_posterior);
mean_mle_1 = mu_candidates(mean_mle_ind);
min_F = @(var)sum((normalised_mean_posterior - normpdf(mu_candidates, mean_mle_1, var)').^2);
options = optimset('TolX', 10^-11); % so that optimisation algo doesn't stop prematurely as we are working with v small orders oif magnitude here
var_mle_1 = fminbnd(min_F, 0, 3*10^-8, options);

% use the fitted gaussian to compute posterior, which will then be a gaussian!
% (gauss x gauss = gauss)
fitted_gaus_prior = normpdf(mu_candidates, mean_mle_1, var_mle_1)';
fitted_gaus_post = fitted_gaus_prior .* y_pdf_prior;
auc_fitted_post = trapz(mu_candidates, fitted_gaus_post);
fitted_gaus_post = fitted_gaus_post / auc_fitted_post;

% some variables to help plot the dotted mean line of posterior
[m, i] = max(fitted_gaus_post);
mu_i = mu_candidates(i);
x_m = [mu_i mu_i];
y_m = [0 m];

% some variables to help plot the dotted mean line of gauss fitted
% likelihood
[m_mle, ~] = max(fitted_gaus_prior);
x_mle = [mean_mle_1 mean_mle_1];
y_mle = [0 m_mle];

figure
plot(mu_candidates, normalised_mean_posterior, 'k--', DisplayName='Original likelihood (non Gauss)')
hold on
plot(mu_candidates, fitted_gaus_prior, 'r', DisplayName='Likelihood (Gauss fitted)', LineWidth=2)
hold on
plot(mu_candidates, y_pdf_prior, 'b', DisplayName='Prior (Gauss)', LineWidth=2)
hold on
plot(mu_candidates, fitted_gaus_post, 'g', DisplayName='Posterior (Gauss)', LineWidth=2)
hold on
plot(x_m, y_m, 'g--','HandleVisibility','off')
hold on
plot(x_mle, y_mle, 'r--', 'HandleVisibility','off')
xlabel("Mutation rate")
ylabel("Probability density")
legend boxoff

label_mean = sprintf("%.5g", mu_i);
word = '\mu = ' + label_mean;
text(mu_i, 1*10^6, word, 'FontSize', 12)

label_mle = sprintf("%.5g", mean_mle_1);
word_mle = '\mu = ' + label_mle;
text(mean_mle_1, 1*10^6, word_mle, 'FontSize', 12, 'HorizontalAlignment','right')

%% option 2: fit a gaussian curve straight to mean_post_post
[~, mean_mle_ind] = max(mean_post_post);
mean_mle_2 = mu_candidates(mean_mle_ind);
min_F = @(var)sum((mean_post_post - normpdf(mu_candidates, mean_mle_2, var)').^2);
options = optimset('TolX', 10^-11); % so that optimisation algo doesn't stop prematurely as we are working with v small orders oif magnitude here
var_mle_2 = fminbnd(min_F, 0, 3*10^-8, options);
fit_gaus_post_post = normpdf(mu_candidates, mean_mle_2, var_mle_2)';

% some variables to help plot the dotted mean line of posterior
[m, i] = max(fit_gaus_post_post);
mu_i = mu_candidates(i);
x_m = [mu_i mu_i];
y_m = [0 m];

figure
plot(mu_candidates, normalised_mean_posterior, 'r', DisplayName='Original likelihood (non Gaus)', LineWidth=2)
hold on
plot(mu_candidates, y_pdf_prior, 'b', DisplayName='Prior (Gauss)', LineWidth=2)
hold on
plot(mu_candidates, fit_gaus_post_post, 'g', DisplayName='Posterior (Gauss fitted)', LineWidth=2)
hold on
plot(mu_candidates, mean_post_post, 'k--', DisplayName='Original posterior (non Gauss)')
hold on
plot(x_m, y_m, 'g--', 'HandleVisibility','off')
xlabel("Mutation rate")
ylabel("Probability density")
legend boxoff

label_mean = sprintf("%.5g", mu_i);
word = '\mu = ' + label_mean;
text(mu_i, 1*10^6, word, 'FontSize', 12)



%% plot the noisy true data and our predicted curve (specified by mu_mean)
[predicted_fitness, ~] = adaptiveWalk(gene_to_fitness, N, mu_i, numGen, selective_pressure);
figure
plot(1:numGen, true_fitness);
hold on
plot(1:numGen, predicted_fitness, LineWidth=2, Color='#FFC300');

%%
figure
plot(1:numGen, true_fitness, HandleVisibility="off");
hold on
plot(1:numGen, predicted_fitness_1, LineWidth=2, Color='r', DisplayName="MLE");
hold on
plot(1:numGen, predicted_fitness, LineWidth=2, Color='#FFC300', DisplayName="MAP");
legend boxoff

%%
% some variables to help plot the dotted mean line of gauss fitted
% likelihood
[m, ~] = max(fitted_gaus_prior);
x_m = [mean_mle_1 mean_mle_1];
y_m = [0 m];
hold on
plot(x_m, y_m, 'r--', 'HandleVisibility','off')

label_mean = sprintf("%.5g", mean_mle_1);
word = '\mu = ' + label_mean;
text(mean_mle_1, 1*10^6, word, 'FontSize', 12, 'HorizontalAlignment','right')