clc, clear
mu_x1 = 0;
mu_x2 = 0;
sigma_x1 = 0.25;
sigma_x2 = 1;
sigma_x1x2 = 0.3;
x1 = -3:0.1:3;
x2 = -3:0.1:3;

Mu = [mu_x1 mu_x2];
Sigma = [sigma_x1 sigma_x1x2; sigma_x1x2 sigma_x2];
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];
y_ = mvnpdf(X,Mu,Sigma);
y = reshape(y_,length(x2),length(x1));
surf(x1,x2,y, 'EdgeColor','none')
% shading interp % (same effect as EdgeColor none)
clim([min(y(:))-0.5*range(y(:)),max(y(:))])
imagesc(x1, x2, y)
axis([-3 3 -3 3 0 0.4])
xlabel('x1')
ylabel('x2')
zlabel('Probability Density')

%%
x = 1:0.1:10;
y = normpdf(x, 5, 1);
x2 = x * 2;
y2 = normpdf(x2, 10, 2);
plot(x, y, 'blue')
hold on
plot(x2, y2, 'red')

%%
[y,z] = meshgrid(linspace(0,10,40));
for off=50:50:200    
    x = off + zeros(size(z));
    % My standin for your vorticity data
    c = cos((x+y)/5) .* cos((x+z)/5);
    surf(x,y,z,c)
    shading interp
    hold on
end
hold off
xlim([0 200])
%%
means = linspace(1, 10, 10);
x = linspace(1, 10, 100);
for i = 1:10
    mu = means(i);
    y = i * ones(100, 1);
    z = normpdf(x, mu, 1);
    plot3(x, y, z)
    hold on
end

hold on
z = zeros(10, 1);
plot3(means, means, z)




