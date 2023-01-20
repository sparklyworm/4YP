% -1: no mutation (wildtype)
%  1: mutation

k = 1;              % each locus interacts with 1 other locus ( 1 + 1 = 2)
N = 4;              % total number of loci
num_interact = k + 1;

% I = interaction matrix: defines which loci interacts with each other
% I = [ 1 1 0; 1 0 1; 0 1 1] for k = 1 (pairwise interaction) and N = 3 
zloc = nchoosek(1:N, num_interact);
m = size(zloc,1);
I = zeros(m,N);
I(sub2ind(size(I),repmat((1:m)',[1,num_interact]),zloc)) = 1;

% all possible combinations for the k+1 interacting loci
% 0 wildtype, 1 mutate
% if k+1 = 2, X = [0 0; 0 1; 1 0; 1 1]
% then replacing all zeros with -1
% to differentiate from interaction matrix I
X = dec2bin(0:2^num_interact-1)-'0';
X(X==0) = -1 ;

% build epistatic matrix Z
% all combinations for each of the sets of k+1 loci
Z = zeros(2^num_interact * nchoosek(N,num_interact), N);
j = 1; % index counter
for i = 1:m
    z = zeros(2^num_interact, N);
    ind = find(I(i, :) == 1);
    z(:, ind) = X;
    Z(j:j+2^num_interact-1, :) = z;
    j = j + 2^num_interact;
end

% fitness contribution from epistatically interacting loci
fitness_epi = rand(1,2^num_interact * nchoosek(N,num_interact));

% different genotypes (wildtype -1, mutation 1)
G0 = dec2bin(0:2^N-1)-'0';
G = dec2bin(0:2^N-1)-'0';
G(G==0) = -1 ;

% fitness value of each genotype
F = zeros(1,2^N);
for i = 1: 2^N
    genotype = G(i, :);
    epi_mat = genotype .* I;
    [~, r] = ismember(epi_mat, Z, 'rows');
    fitness_value = sum(fitness_epi(r))/nchoosek(N,num_interact);
    F(i) = fitness_value;
end

visualise = [G0, F']



