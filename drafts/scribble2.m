clc, clear
n = 3;
a = 0.01;

gene_to_fitness =  assignFitness(n, a);
gene_arr = gene_to_fitness.keys();
gene2genes = dictionary(gene_arr, cell(2^n, 1));

%%
g0 = "000";
c = gene2genes(g0)
if isempty(c{1})
    disp("empty")
end
length(c{1})
%%
gene1 = "101";
gene2genes(g0) = {gene1};
%%
gene2 = "110";
c = gene2genes(g0)
newc = [c, {gene2}]
gene2genes(g0) = {newc}


%%
gene3 = "111";
c = gene2genes(g0)
newc = [c{1} {gene3}] % unpack current cell, concacenate with new cell
gene2genes(g0) = {newc} % pack it up as one again 

%%
gene2genes = hammingMapFitness(gene_to_fitness)

