%%
a = 5;                 % true mean (it's winter...)
s = 4;                 % true noise std 
n = 4;                 % # data points
y = a + s*randn(n,1);  % Simulates fake data

va  = linspace(-15,15,1000);
N   = length(va);
like = exp(-.5*sum((repmat(y,1,N)-repmat(va,n,1)).^2,1)/s^2);
lik = like/sum(like); % normalise

figure,hold on,grid on
plot(y,0,'.','color','g','markersize',20);
hli=plot(va,lik,'r','linewidth',2);
%%


a   = 2; % true a
sig = 5; % true sigma
t   = linspace(0,10,20)'; % regressor (known)

% generate data
y = a*t + sig*randn(size(t));

% plot data
figure
plot(t,y,'.','linewidth',2)

% grid
va = linspace(0,5,100);     % values of a
vs = linspace(.01,80,100);  % values of s^2

posterior = zeros(length(vs),length(va));

for i=1:length(vs)
for j=1:length(va)
S  = va(j)*t;                                        % prediction
li = vs(i).^(-n/2)*exp(-sum((y-S).^2)/2/vs(i));      % likelihood
pr = 1/vs(i) * normpdf(va(i),0,1000);                % prior
posterior(i,j) = li*pr;                              % posterior
end
end
posterior = posterior / sum(posterior(:));

figure
imagesc(va,vs,posterior);
xlabel('a');
ylabel('s^2')

% gibbs sampling
niter = 10000;      % number of iterations
a_s=zeros(niter,1); % container for samples for a
v_s=zeros(niter,1); % container for samples for s^2

% initial samples (you can change this to the maximum likelihood):
a_s(1)=2;
v_s(1)=10;

for i=2:niter
    % samples a (Gaussian)
    a_s(i) = normrnd(sum(t.*y)./sum(t.*t),sqrt(v_s(i-1)/sum(t.*t)));
    
    % sample sig^2 (inverse-Gamma)
    v_s(i) = 1/gamrnd( n/2, 2/sum( (y-a_s(i-1)*t).^2 ) );
end


figure
imagesc(va,vs,posterior);
xlabel('a');
ylabel('s^2');
hold on
plot(a_s(1:20:end),v_s(1:20:end),'.k')  % plot every 20th sample

% plot marginals
figure
subplot(1,2,1),hold on
[n,x]=hist(a_s,va);n=n/sum(n);            % histogram from samples
bar(x,n,'r');
plot(va,sum(posterior,1)','linewidth',2); % sum one dimension of grid

subplot(1,2,2),hold on
[n,x]=hist(v_s,vs);n=n/sum(n);
bar(x,n,'r');
plot(vs,sum(posterior,2),'linewidth',2);