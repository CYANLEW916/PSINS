function results = fdi_swglt(Z, H, sigma_vec, cfg)
% FDI_SWGLT  Sliding Window Whitened GLT with robust matched-filter isolation.
%   results = FDI_SWGLT(Z, H, sigma_vec, cfg) computes detection/isolation
%   statistics using sliding-window parity-mean tests and robust isolation
%   logic (margin ratio, Bayesian posterior, and dwell confirmation).
%
%   Inputs:
%     Z         - N x n measurement matrix
%     H         - n x m configuration matrix
%     sigma_vec - n x 1 or 1 x n sensor noise standard deviations
%     cfg       - struct with j/window_length, K_tolerance, P_D_star,
%                 N_dwell, rho_threshold, P_isol
%
%   Output:
%     results   - struct containing FD/FI curves, decisions, thresholds,
%                 posteriors, and adaptive-threshold diagnostics
%
% See also  compute_parity_matrix, compute_adaptive_threshold, fdi_wglt.

    [N, n] = size(Z);
    m = size(H, 2);
    n_parity = n - m;

    if isfield(cfg, 'j')
        j = cfg.j;
    else
        j = cfg.window_length;
    end

    [~, V_w, W, H_w] = compute_parity_matrix(H, sigma_vec);
    VW = V_w * W;

    [T_adapt, T_adapt_proj, diag_info] = compute_adaptive_threshold( ...
        V_w, H_w, sigma_vec, cfg.K_tolerance, cfg.P_D_star, j);

    c_tilde = V_w;
    c_norm_sq = sum(c_tilde.^2, 1)';

    P_tilde_all = (VW * Z')';

    FD_energy = zeros(N, 1);
    FD_max = zeros(N, 1);
    FI = zeros(N, n);
    isolated = zeros(N, 1);
    f_hat = zeros(N, 1);
    iso_status = zeros(N, 1);
    posterior = zeros(N, n);

    buf = zeros(j, n_parity);
    buf_sum = zeros(1, n_parity);
    buf_idx = 0;
    filled = 0;

    dwell_candidate = 0;
    dwell_count = 0;

    for k = 1:N
        old_idx = mod(buf_idx, j) + 1;
        if filled >= j
            buf_sum = buf_sum - buf(old_idx, :);
        end
        buf(old_idx, :) = P_tilde_all(k, :);
        buf_sum = buf_sum + P_tilde_all(k, :);
        buf_idx = buf_idx + 1;
        filled = min(filled + 1, j);

        if filled < j
            continue;
        end

        P_bar = (buf_sum' / j);
        FD_energy(k) = j * (P_bar' * P_bar);

        for i = 1:n
            ci = c_tilde(:, i);
            FI(k, i) = j * (ci' * P_bar)^2 / c_norm_sq(i);
        end
        FD_max(k) = max(FI(k, :));

        if FD_energy(k) <= T_adapt
            isolated(k) = 0;
            iso_status(k) = 0;
            dwell_candidate = 0;
            dwell_count = 0;
            continue;
        end

        [FI_sorted, sort_idx] = sort(FI(k, :), 'descend');
        k_star = sort_idx(1);

        if FI_sorted(2) > 0
            rho = FI_sorted(1) / FI_sorted(2);
        else
            rho = Inf;
        end

        log_lik = FI(k, :) / 2;
        log_lik_shifted = log_lik - max(log_lik);
        exp_lik = exp(log_lik_shifted);
        posterior(k, :) = exp_lik / sum(exp_lik);

        if k_star == dwell_candidate
            dwell_count = dwell_count + 1;
        else
            dwell_candidate = k_star;
            dwell_count = 1;
        end

        margin_ok = (rho >= cfg.rho_threshold);
        posterior_ok = (posterior(k, k_star) >= cfg.P_isol);
        dwell_ok = (dwell_count >= cfg.N_dwell);

        if margin_ok && posterior_ok && dwell_ok
            isolated(k) = k_star;
            iso_status(k) = 2;
            ci_star = c_tilde(:, k_star);
            f_hat(k) = sigma_vec(k_star) * (ci_star' * P_bar) / c_norm_sq(k_star);
        elseif margin_ok && posterior_ok
            isolated(k) = k_star;
            iso_status(k) = 1;
            ci_star = c_tilde(:, k_star);
            f_hat(k) = sigma_vec(k_star) * (ci_star' * P_bar) / c_norm_sq(k_star);
        else
            isolated(k) = -1;
            iso_status(k) = 1;
        end
    end

    results.FD_energy = FD_energy;
    results.FD_max = FD_max;
    results.FI = FI;
    results.isolated = isolated;
    results.f_hat = f_hat;
    results.iso_status = iso_status;
    results.posterior = posterior;
    results.T_adaptive = T_adapt;
    results.T_adaptive_proj = T_adapt_proj;
    results.diagnostics = diag_info;
end
