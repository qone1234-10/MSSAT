function [sigma2hat, spike_count] = BEMA_combined(data, alpha)

    if nargin < 2
        alpha = 0.2;
    end

    [n, p] = size(data);
    S = (data' * data) / n;
    lambda = eigs(S,p);          

    gamma = min(p, n) / max(p, n);

    k_min = floor(min(p, n) * alpha);
    k_max = floor(min(p, n) * (1 - alpha));
    k = k_min:k_max;

   
    predictor = zeros(length(k), 1);
    for i = 1:length(k)
        predictor(i) = qmp(k(i) / min(p, n), max(n, p), min(n, p)) * max(p, n) / n;
    end

   
    l_reversed = flip(lambda(k)) ;

    
    sigma2hat =  predictor\l_reversed;

    
    cutoff = sigma2hat * ((1 + sqrt(gamma))^2 ...
              + qtw(0.9) * max(p, n)^(-2/3) * gamma^(-1/6) * (1 + sqrt(gamma))^(4/3)) ...
              * max(p, n) / n;

   
    spike_count = sum(lambda > cutoff);
end
