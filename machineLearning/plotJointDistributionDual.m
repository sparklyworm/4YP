function plotJointDistributionDual(X, Mu, Cov)
    % X is n x 2 matrix
    % first column is x1
    % second column is x2
    x1 = X(:, 1);
    x2 = X(:, 2);
    [X1,X2] = meshgrid(x1,x2);
    X = [X1(:) X2(:)];
    y = mvnpdf(X,Mu',Cov);
    y = reshape(y,length(x2),length(x1));
    
    figure
    imagesc(x1, x2, y)
    xlabel('w0')
    ylabel('w1')
    set(gca,'YDir','normal')
end