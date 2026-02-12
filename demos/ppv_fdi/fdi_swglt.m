function [FD, FI, T_adaptive] = fdi_swglt(Z, sigma_vec, cfg)
% FDI_SWGLT  PPV-aided Sliding Window Weighted GLT for fault detection and isolation.
%   [FD, FI, T_adaptive] = FDI_SWGLT(Z, sigma_vec, cfg) applies the complete
%   PPV-aided SWGLT algorithm with PCA denoising and adaptive thresholds.
%
%   Inputs:
%     Z         - N x 9 measurement matrix
%     sigma_vec - 1 x 9 noise standard deviation vector
%     cfg       - configuration structure from config.m
%
%   Outputs:
%     FD         - N x 1 detection function values
%     FI         - N x 9 isolation function values
%     T_adaptive - scalar adaptive detection threshold
%
% See also  compute_parity_matrix, compute_ppv, compute_adaptive_threshold, config.

    [N, n] = size(Z);
    H = cfg.H;
    j = cfg.window_length;           % PCA window length (25 samples)
    pca_thresh = cfg.pca_threshold;  % cumulative variance ratio threshold

    % Compute weighted parity matrix
    [~, V_w, W_inv_half] = compute_parity_matrix(H, sigma_vec);
    n_parity = size(V_w, 1);        % parity space dimension (6)

    % Compute adaptive threshold
    T_adaptive = compute_adaptive_threshold(H, V_w, sigma_vec, cfg.K_tolerance);

    % Pre-allocate outputs
    FD = zeros(N, 1);
    FI = zeros(N, n);

    % Compute all weighted parity vectors first
    P_w_all = zeros(n_parity, N);
    for k = 1:N
        z = Z(k, :)';
        P_w_all(:, k) = V_w * W_inv_half * z;
    end

    % Sliding window PCA processing
    for k = 1:N
        if k < j
            % Not enough samples for full window: use available samples
            P_buffer = P_w_all(:, 1:k);
            if k < 3
                % Too few samples for meaningful PCA, use raw parity vector
                P_r = P_w_all(:, k);
            else
                P_r = compute_ppv(P_buffer, pca_thresh);
            end
        else
            % Full window available
            P_buffer = P_w_all(:, (k-j+1):k);
            P_r = compute_ppv(P_buffer, pca_thresh);
        end

        % Detection function: FD = P_r' * P_r
        FD(k) = P_r' * P_r;

        % Isolation function (when detection occurs)
        Vi_w = V_w * W_inv_half;
        for i = 1:n
            Vi = Vi_w(:, i);
            FI(k, i) = abs(P_r' * Vi) / (Vi' * Vi);
        end
    end

end
