function newGenCount = sampleNewGeneration(N, p)
    
    % calculate covariance matrix
    cov = N* diag(p);
    
    % check symmetry 
    if ~issymmetric(cov)
        disp(cov)
    end

    % sample from multivariate normal distribution; floor to get integers 
    count = floor(mvnrnd(N*p, cov));
    count = abs(count');
    
    
    newGenCount = count;
end

%{
    p
    % find entries which are non zero 
    nonZeroIndex = find(p)
    % get the p which are non zero 
    nonZeroP = p(nonZeroIndex);

    if length(nonZeroP) == 1
        newGenCount = zeros(length(p), 1);
        newGenCount(nonZeroP) = N;
        return
    else
       [min_p, min_i] = min(nonZeroP);
       nonZeroP_new = nonZeroP(find(nonZeroP ~= min_p));
    end
  
    % calculate covariance matrix
    %cov = calculateMultinomialCovariance(nonZeroP_new, N)
    cov = N* diag(nonZeroP_new)
    % sample from multivariate normal distribution; floor to get integers 
    count = floor(mvnrnd(N*nonZeroP_new, cov));
    count = count'
    
    % reinsert the element with min_p
    if min_i == 1
        count = [0; count];
    elseif min_i == length(nonZeroP)
        count = [count; 0];
    elseif min_i == length(nonZeroP) -1
        count = [count(1: min_i -1); 0; count(end)];
    else
        count = [count(1: min_i -1); 0; count(min_i: end)];
    end

    
    % make sure everything sums up to N
    if sum(count) ~= N
        diff = N - sum(count);
        [~, max_i] = max(count);
        if diff > 0
            count(min_i) = min(floor(normrnd(N*min_p, sqrt(N*min_p*(1-min_p)))), diff);
            remain = diff - floor(normrnd(N*min_p, sqrt(N*min_p*(1-min_p))));
            count(max_i) = count(max_i) + remain;
        else
             count(max_i) = count(max_i) - diff;
        end
    end

    newGenCount = zeros(length(p), 1)
    newGenCount(nonZeroIndex) = count





%}