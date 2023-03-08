% week 5
% to find range for threshold values
% 3 landscapes, 5 experiments each
clc,clear

% true data
true_a = 0.2;
true_mu = 1 * 10 ^-6;
true_n = 8; 
% population parameters
N = 10^8;
numGen = 600;
selective_pressure = 1;

mu_logthreshold = log10(3);
n_threshold = 1;
numExp = 5;
numLandscape =3;
numSamples = 30;

std_n = 0.8;
std_time = 115;

fdiff_average = zeros(numLandscape, 1); % average difference in evolution curve
fdiff_stds = zeros(numLandscape, 1);
tdiff_average = zeros(numLandscape, 1); % average difference in t_01_99
tdiff_stds = zeros(numLandscape, 1);
n_average = zeros(numLandscape, 1); % average difference in n
n_stds = zeros(numLandscape, 1);
obj_average = zeros(numLandscape, 1); % average size of objective function
obj_stds = zeros(numLandscape, 1);

for k = 1:numLandscape

    % generate a few replicate experiments
    gene_to_fitness =  assignFitness(true_n, true_a);  
    experiments = cell(numExp, 1);
    figure
    for i = 1:numExp
        [true_fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen, selective_pressure);
        experiments{i} = true_fitness;
        hold on
        plot(1:numGen, true_fitness)
    end
    title_text = sprintf("mutation rate = %.2g, n = %i", true_mu, true_n);
    title(title_text)
    
    % get distribution for each experiment
    mu_distributions = cell(numExp, 1);
    n_distributions = cell(numExp, 1);
    f_diff_distributions = cell(numExp, 1);
    t_diff_distributions = cell(numExp, 1);
    obj_distributions = cell(numExp, 1);

    
    % initialise prior
    mu_prior = makedist('Loguniform','Lower',10^-9,'Upper',10^-4);
    n_domain =1:1:10;

    for j = 1:numExp
        fprintf("====== Experiment %i ======\n", j)
        true_fitness = experiments{j};
        true_t_01 = computeTimeToPercentageMax(true_fitness, 0.01);
        true_t_99 = computeTimeToPercentageMax(true_fitness, 0.99);
        true_t_01_to_99 = true_t_99 - true_t_01;
       
        mu_samples = zeros(numSamples, 1);
        n_samples = zeros(numSamples, 1);
        f_samples = zeros(numGen, numSamples);
        f_distances = zeros(numSamples, 1);
        t_diffs = zeros(numSamples, 1);
        obj_funcs = zeros(numSamples, 1);
        for i = 1:numSamples
            
            mu_logdiff  = mu_logthreshold + 1;
            n_diff = n_threshold + 1;
    
            while ~(mu_logdiff <= mu_logthreshold && n_diff <= n_threshold)
                mu_s = random(mu_prior,1);
                n_s = randsample(n_domain,1);

                mu_logdiff = abs(log10(true_mu) - log10(mu_s));
                n_diff = abs(true_n - n_s);
            end
            mu_samples(i) = mu_s;
            n_samples(i) = n_s;
            
            [sim_fitness, ~, ~] = simulator(mu_s, n_s, true_a, false);
            inflpts =  findInflectionPoints2ndOrder(sim_fitness);
            num_levels = round((length(inflpts) - 1)/2);
            f_samples(:, i) = sim_fitness;
    
            t_99 = computeTimeToPercentageMax(sim_fitness, 0.99);
            t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
            t_01_to_99 = t_99 - t_01;
            t_diff = abs(t_01_to_99 - true_t_01_to_99);
           
            f_distances(i) = sqrt(sum((sim_fitness(t_01:t_99) - true_fitness(true_t_01:true_t_01 + t_01_to_99)).^2)/t_01_to_99);
            t_diffs(i) = t_diff;
            obj_funcs(i) = 1/std_n * abs(num_levels - true_n) + 1/std_time * t_diff;

            % disp progress
            if mod(i, numSamples/10) ==0
                fprintf("sample %i \n", i)
            end
        end
        
        % record data
        mu_distributions{j} = mu_samples;
        n_distributions{j} = n_samples;
        f_diff_distributions{j} = f_distances;
        t_diff_distributions{j} = t_diffs;
        obj_distributions{j} = obj_funcs;
    end
    
    
    % combine all experiments
    numSamples = length(mu_distributions{1});
    n_distribution = zeros(numSamples*numExp, 1);
    f_diff_distribution =  zeros(numSamples*numExp, 1);
    t_diff_distribution =  zeros(numSamples*numExp, 1);
    obj_distribution = zeros(numSamples*numExp, 1);
    for j = 1:numExp
        start = (j -1)*numSamples + 1;
        finish = start + numSamples -1;
        n_distribution(start:finish) = n_distributions{j};
        f_diff_distribution(start:finish) = f_diff_distributions{j};
        t_diff_distribution(start:finish) = t_diff_distributions{j};
        obj_distribution(start:finish) = obj_distributions{j};
    end
    
    fdiff_average(k) = mean(f_diff_distribution)
    fdiff_stds(k) = std(f_diff_distribution)
    tdiff_average(k) = mean(t_diff_distribution)
    tdiff_stds(k) = std(t_diff_distribution)
    n_average(k) = mean(n_distribution)
    n_stds(k) = std(n_distribution)
    obj_average(k) = mean(obj_distribution)
    obj_stds(k) = std(obj_distribution)
end