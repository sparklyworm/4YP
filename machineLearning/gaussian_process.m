%% gaussian process
clc, clear
x_start = 0;
x_end = 10;
num_test = 50;
x_test = linspace(x_start, x_end, num_test)';
mu_test = zeros(num_test,1);

%% prior
% sigma and lengthScale are specifications of our prior!!!!
sigma = 1;
lengthScale = 2;
K = covRBF(x_test, [],sigma, lengthScale); % covariance function (RBF cov)
figure
imagesc(K)
%%
% prior of our curves (our initial believe in how the curves are : smooth
% and curvy)

num_samples = 5;
rng("default")
R = mvnrnd(mu_test,K,num_samples); % R is a 5 x 50 matrix, where each row holds the y-values of curves

figure
for i = 1: num_samples
    plot(x_test, R(i,:))
    hold on
end
r = mvnrnd(mu_test,K);
plot(x_test, r, 'k.')
plot(x_test, mu_test, 'Color', 'black', 'LineWidth', 1.5 )
hold on 
top_sd = mu_test + sigma;
bottom_sd = mu_test - sigma;
patch('Faces', linspace(1,2*num_test,2*num_test),'Vertices',[[x_test; flip(x_test)], [top_sd; flip(bottom_sd)]], 'FaceColor', 'magenta','EdgeColor', 'none','FaceAlpha', '0.2')
ylim([-6,6])

%% posterior

%training data 
x_train = [1; 3; 5];
y_train = [1; 9; 6];
n_train = length(x_train);

% hyperparameter
sigma_noise = 0.5; % noise variance 
sigma_signal = 1; % signal variance (how fast changing/squiggly do we want our curve to be)
lengthScale = 1; % bandwidth (how much covariance - how much 2 rv's that are far apart can influence each other)

for i = 1:n_train
    
    % calculate posterior
    mu_train = zeros(i, 1);
    [mu, cov] = posteriorGP(x_test, x_train(1:i), y_train(1:i), mu_test, mu_train, sigma_noise, sigma_signal, lengthScale);
    figure
    imagesc(cov);

    num_samples = 5;
    rng("default")
    posterior_samples = mvnrnd(mu,cov,num_samples);
    sd = sqrt(diag(cov));
    
    figure
    for s = 1: num_samples
        f = posterior_samples(s,:)';
        %top_sd = f + sd;
        %bottom_sd = f - sd;
        %patch('Faces', linspace(1,2*num_test,2*num_test),'Vertices',[[x_test; flip(x_test)], [top_sd; flip(bottom_sd)]], 'FaceColor', 'magenta','EdgeColor', 'none','FaceAlpha', '0.2')
        plot(x_test, posterior_samples(s,:))
        hold on
    end

    % plot mean
    hold on 
    plot(x_test, mu, 'Color','#144272', 'LineWidth', 1.5);

    % plot cov 
    top_sd = mu + diag(cov);
    bottom_sd = mu - diag(cov);
    patch('Faces', linspace(1,2*num_test,2*num_test),'Vertices',[[x_test; flip(x_test)], [top_sd; flip(bottom_sd)]], 'FaceColor', 'magenta','EdgeColor', 'none','FaceAlpha', '0.2')

    % plot training points
    hold on 
    plot(x_train(1:i), y_train(1:i), 'bo', 'LineWidth', 2, 'MarkerSize', 7);
    
    ylim([min(y_train)-5, max(y_train)+5])
    

end


%% 3d visualisation
i = 3;
figure
for s = 1: num_samples
    f = posterior_samples(s,:)';
    top_sd = f + sd;
    bottom_sd = f - sd;
    patch('Faces', linspace(1,2*num_test,2*num_test),'Vertices',[[x_test; flip(x_test)], [top_sd; flip(bottom_sd)]], 'FaceColor', 'magenta','EdgeColor', 'none','FaceAlpha', '0.2')
    plot3(x_test, posterior_samples(s,:), zeros(num_test,1))
    hold on
end

% plot mean
hold on 
plot3(x_test, mu, zeros(num_test,1),'Color','#144272', 'LineWidth', 1.5);

% plot training points
hold on 
plot3(x_train(1:i), y_train(1:i), zeros(i,1), 'bo', 'LineWidth', 2, 'MarkerSize', 7);

ylim([min(y_train)-5, max(y_train)+5])

for i = 1:num_test
    mu_i = mu(i);
    sd_i = sd(i);
    y = linspace(mu_i - 4*sd_i, mu_i + 4*sd_i, 100)';
    z = normpdf(y, mu_i, sd_i);
    x = x_test(i) * ones(100, 1);
    plot3(x, y, z)
    hold on
end
