%clc, clear
% to "guess" what was the underlying means+std that generated this data 

% r = the ranks of the ordered data (the first=earliest=smallest)
% n = total number of data
% observed = the observed values of our data as they come in 

load('simulation_mu_1.1e-08_selective_0.9.mat', 'plateauSimMat', 'plateauMeans')
%sort(plateauSimMat(:, 2))

%%
n = 100;
%observed = [14; 17; 18; 19; 20];index = 2;
observed = [43; 43; 44; 44; 44; 44; 45; 45; 45]; index = 1;
r = [100:-1:(100-length(observed)+1)];
sigma_best = 1;
mu_best = observed(end);

num_iter = 20;
mu_options = optimoptions(@fmincon,'StepTolerance', 0.01);
sigma_options = optimoptions(@fmincon,'StepTolerance', 0.01);
mu_arr = zeros(1,num_iter);
sigma_arr = zeros(1, num_iter);
%%
for i = 1: num_iter
    
    mu_old = mu_best;
    sigma_old = sigma_best;
    
    % differentiable: differentiate find min point of sum squares 
    syms mu
    f = sum((expectedNormalOrderArray(r, n, 0, sigma_best) + mu - observed).^2);
    mu_best = double(solve(diff(f)==0));
    
    % round up to speed up 
    if round(mu_best) > mu_best
        mu_best = round(mu_best);
    end
   
    tic
    % not differentiable: use fmincon to optimise
    func_min_sigma = @(sigma)sum((expectedNormalOrderArray(r, n, 0, sigma) + mu_best - observed).^2);
    sigma_best = fmincon(func_min_sigma, sigma_best,  [], [],[],[],[],[],[], sigma_options);
    toc

    
    disp(i)
    disp(sigma_best)
    disp(mu_best)
    mu_arr(i) = mu_best;
    sigma_arr(i) = sigma_best;

    if abs(mu_best - mu_old) < 0.5 && abs(sigma_best - sigma_old) < 0.05
        break
    end

end


%func_min_mu = @(mu)sum((expectedNormalOrderArray(r, n, 0, sigma_best) + mu - observed).^2);
%mu_best = fmincon(func_min_mu, mu_best, [], [],[],[],[],[],[], mu_options);
   
%%
figure
plot(1:i, mu_arr(1:i), 'DisplayName','estimate \mu')
hold on
plot(1:i, repmat(plateauMeans(index), 1, i), 'r', 'DisplayName','true \mu')
xlabel("Iteration")
ylabel("Estimate of Mean \mu")
legend box off

figure
plot(1:i, sigma_arr(1:i), 'DisplayName','estimate \sigma')
std_true = std(plateauSimMat(:, index));
hold on
plot(1:i, repmat(std_true, 1, i), 'r', 'DisplayName','true \sigma')
xlabel("Iteration")
ylabel("Estimate of Standard Deviation \sigma")
legend boxoff
%%
x = 43:0.1:53;
y = normpdf(x, mu_best, sigma_best);
hold on
plot(x, 96*y, LineWidth=2)