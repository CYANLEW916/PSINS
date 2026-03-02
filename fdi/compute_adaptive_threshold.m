function [T_adaptive, T_adaptive_proj, diagnostics] = compute_adaptive_threshold( ...
    V_w, H_w, sigma_vec, K, P_D_star, j)
% COMPUTE_ADAPTIVE_THRESHOLD  Non-central chi-square adaptive thresholds.
%   Implements eqs.(72)-(80), (89)-(92) for SWGLT energy/projection tests.
%
%   Inputs:
%     V_w       - (n-m) x n whitened parity matrix
%     H_w       - n x m whitened configuration matrix (W*H)
%     sigma_vec - n x 1 or 1 x n noise standard deviations
%     K         - tolerance ratio
%     P_D_star  - target detection probability
%     j         - sliding window length
%
%   Outputs:
%     T_adaptive      - threshold for energy test FD_j ~ chi2(n-m)
%     T_adaptive_proj - threshold for max-projection test FI_j ~ chi2(1)
%     diagnostics     - struct with per-sensor intermediate variables

    [n, m] = size(H_w);
    n_parity = n - m;
    sigma_vec = sigma_vec(:);

    HTH = H_w' * H_w;
    inv_HTH = inv(HTH);
    gamma_0 = trace(inv_HTH);
    inv_HTH_sq = inv_HTH * inv_HTH;

    f_t = zeros(n, 1);
    lambda_j = zeros(n, 1);
    beta_i = zeros(n, 1);
    c_norm_sq = zeros(n, 1);

    for i = 1:n
        c_i = V_w(:, i);
        h_i = H_w(i, :)';

        beta_i(i) = h_i' * inv_HTH_sq * h_i;
        c_norm_sq(i) = c_i' * c_i;

        f_t(i) = sigma_vec(i) * sqrt(K * gamma_0 / beta_i(i) + 1 / c_norm_sq(i));
        lambda_j(i) = j * (K * gamma_0 / beta_i(i) * c_norm_sq(i) + 1);
    end

    [lambda_min, weakest_idx] = min(lambda_j);
    T_adaptive = ncx2inv(1 - P_D_star, n_parity, lambda_min);

    lambda_min_proj = lambda_min;
    T_adaptive_proj = ncx2inv(1 - P_D_star, 1, lambda_min_proj);

    T_fixed = chi2inv(1 - 0.01, n_parity);
    if T_adaptive >= T_fixed
        fprintf('Notice: T_adaptive >= T_fixed under current tuning.\n');
    end

    P_FA_adapt = 1 - chi2cdf(T_adaptive, n_parity);
    P_FA_proj = n * (1 - chi2cdf(T_adaptive_proj, 1));

    fprintf('Adaptive threshold (energy): T=%.4f (P_FA=%.6f)\n', T_adaptive, P_FA_adapt);
    fprintf('Adaptive threshold (proj):   T=%.4f (P_FA_bonf=%.6f)\n', ...
        T_adaptive_proj, P_FA_proj);
    fprintf('Fixed threshold (P_FA=0.01): T=%.4f\n', T_fixed);
    fprintf('Weakest sensor: %d (lambda_min=%.4f)\n', weakest_idx, lambda_min);

    diagnostics.f_t = f_t;
    diagnostics.lambda_j = lambda_j;
    diagnostics.beta_i = beta_i;
    diagnostics.c_norm_sq = c_norm_sq;
    diagnostics.gamma_0 = gamma_0;
    diagnostics.lambda_min = lambda_min;
    diagnostics.weakest_sensor = weakest_idx;
    diagnostics.T_fixed = T_fixed;
    diagnostics.P_FA_adapt = P_FA_adapt;
    diagnostics.P_FA_proj = P_FA_proj;
end
