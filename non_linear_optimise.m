clc, clear
load("interesting-landscape\interesting-landscape-5.mat")
%%
numGen = 600;
selective_pressure = 1;

mu = 1 * 10^-7;
[fitness, ~] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);
plot([1:numGen], fitness);
%%
weights = [repmat(1,160,1); repmat(1,440,1)];
func_min = @(x)sum(((adaptiveWalkParamOpt(x) - fitness).*weights).^2);
options = optimoptions('fmincon', 'StepTolerance', 10^-10);

fmincon(func_min, 1*10^-8, [], [], [], [], 1*10^-9, 1*10^-4)

%%
