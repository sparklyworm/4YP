clc, clear
%%
n = 4; % number of genes !!!!!!!!!!!!!!!!!!!!!! 
a = 0.3; % random tune param
N = 10^8; % population size
mu = 1.1 *10^-8; % mutation rate

gene_to_fitness =  assignFitness(n, a);    % maps genotype to fitness value
genotypes = gene_to_fitness.keys();        % string array of genotypes

plotFitnessLandscape(gene_to_fitness, n, true);
% save("interesting-landscape.mat")
%%
avg_fitness = zeros(n, 1);
mutated_gene_count = zeros(n, 1);
for i = 2:2^n
    gene = genotypes(i);
    num_mutated = length(find(convertStringsToChars(gene) == '1'));
    avg_fitness(num_mutated) = avg_fitness(num_mutated) + gene_to_fitness(gene);
    mutated_gene_count(num_mutated) =  mutated_gene_count(num_mutated) + 1;
end
avg_fitness = avg_fitness ./ mutated_gene_count;
%%
selective_pressure = 1;
MU = [1] * 10^-8; % [ 1 5 10] * 10^-8
numSimulation = 1;
numGen = 600; % 600
numPlot = 1;
figure
for mu = MU
    disp(numPlot)
    for i = 1:numSimulation
        [fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);
        subplot(1, 3, numPlot)
        plot([1:numGen], fitness);
        title(sprintf("mutation rate %g", mu))
        hold on
    end
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'off';
    numPlot = numPlot + 1;
end

%%
figure
plot([1:numGen], fitness);

for i = 1:n
    hold on
    plot(1:numGen, avg_fitness(i)*ones(numGen), 'r')
end