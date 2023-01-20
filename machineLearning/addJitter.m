function jitter_K = addJitter(K)
    jitter = 1*10^-5 * eye(length(K));
    jitter_K = K + jitter;

end