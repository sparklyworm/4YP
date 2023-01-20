clc, clear

n = 4; % number of genes
a = 0.5; % random tune param
N = 10^8; % population size
mu = 10^-8; % mutation rate

gene_to_fitness =  assignFitness(n, a);    % maps genotype to fitness value
genotypes = gene_to_fitness.keys();        % string array of genotypes
plotFitnessLandscape(gene_to_fitness, n, true);

%%
genotype_count = dictionary(genotypes, zeros(2^n, 1)); % count of individuals of each genotype
genotype_p = dictionary(genotypes, zeros(2^n, 1));     % frequency """
gene2genes = hammingMapFitness(gene_to_fitness);       % maps genotype to genotypes 1 mutation step away & are of higher fitness 

wildtype = genotypes(1);
genotype_count(wildtype) = N;
genotype_p(wildtype) = 1;



%%
tic
numGen = 600;
fitness = zeros(numGen, 1);
genotype_count_gen = cell(numGen, 1);
for gen = 1:numGen

    disp(gen)

    r_mean = calculateMeanFitness(genotype_count, gene_to_fitness);

    fitness(gen) = r_mean;                    % record for plotting
    genotype_count_gen{gen} = genotype_count; % record for plotting 
    
    p = updateFrequency(genotype_p, gene_to_fitness, r_mean); % apply selective pressure
    count = sampleNewGeneration(N, p); % multivariate normal distribution 
    genotype_count = dictionary(genotypes, count);
    genotype_count_new = genotype_count;

    for i = 1: length(genotypes)
        
        gene = genotypes(i);
        if genotype_count(gene) == 0
            continue
        end

        g = gene2genes(gene);
        % mutate_g is the array of genotypes that are 1 mutation step away from gene
        if length(g{1}) == 1
            mutate_g = g{1};
        else
            mutate_g = string(g{1}); % string array ["0110", "1011"]
        end

        lambda = mu*genotype_count(gene);
        for k = 1:length(mutate_g)
            new_gene = mutate_g(k);
            new_ind = poissrnd(lambda);
            
            genotype_count_new(new_gene) = genotype_count_new(new_gene) + new_ind;
            genotype_count_new(gene) = genotype_count_new(gene) - new_ind;
        end

    end
    genotype_count = genotype_count_new;
    p = genotype_count_new.values()/sum(genotype_count_new.values());
    genotype_p = dictionary(genotypes, p);


end

time_taken = toc;

figure
plot([1:numGen], fitness);
save("draft_adaptive_walk_.mat")

%%
threshold = 10^7; logScale = false;
plotAdaptiveWalk(genotype_count_gen, mu, threshold, logScale);

