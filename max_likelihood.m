clc, clear
load("interesting-landscape\interesting-landscape-5.mat")


%% our observed data
selective_pressure = 1;
true_mu = 5 * 10^-8;
numGen = 600;
[model_fitness, ~] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen, selective_pressure);
true_fitness = zeros(numGen, 1);
for i = 1:numGen
    true_fitness(i) = model_fitness(i) + 0.05*(-1 + 2*rand(1));
end
%%
figure
plot([1:numGen], true_fitness);

%% max likelihood

mu_candidates = [linspace(5*10^-9, 10^-8, 5), linspace(10^-8, 10^-7, 5), linspace(10^-7, 10^-6, 5), linspace(10^-6, 10^-5, 5)];
%mu_candidates = linspace(1*10^-8, 1*10^-7, 30);
num_candidates = length(mu_candidates);
num_samples = 10;
sum_squared_errors = zeros(num_candidates, 1);

% for each candidate mu 
for c = 1:num_candidates

    mu_c = mu_candidates(c);
    samples = zeros(numGen, num_samples);

    % generate 10 samples
    for s = 1:num_samples
        [fitness, ~] = adaptiveWalk(gene_to_fitness, N, mu_c, numGen, selective_pressure);
        samples(:, s) = fitness;
    end
    avg_fitness = sum(samples, 2)/num_samples;
    sum_squared_errors(c) = sum((avg_fitness - true_fitness).^2);
    
    % some prints to track progress
    if c == 1
        fprintf("\n================================================ \n")
        fprintf("Total number of candidates: %i \n", num_candidates)
        fprintf("Start... \n")
    end
    fprintf("%i ", c)
    if c == num_candidates
        fprintf("\nOperation complete! \n")
    end
end

%% plot
figure
plot(mu_candidates, sum_squared_errors)
figure
probability = 1./sum_squared_errors;
plot(mu_candidates, probability)


