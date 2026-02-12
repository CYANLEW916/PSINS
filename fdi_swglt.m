function [FD, FI, T_adaptive] = fdi_swglt(Z, sigma_vec, cfg)
% FDI_SWGLT  PPV-aided Sliding Window Weighted GLT for fault detection/isolation.
%   [FD, FI, T_adaptive] = FDI_SWGLT(Z, sigma_vec, cfg) applies the full
%   PPV-aided SWGLT algorithm combining weighted parity analysis, PCA-based
%   denoising via sliding windows, and adaptive fault-tolerant thresholds.
%
%   Algorithm flow:
%     1. Compute weighted parity vectors P_w for all time steps
%     2. Sliding window of j=25 samples: apply PCA to get principal parity vector P_r
%     3. Detection: FD = P_r' * P_r, compare against adaptive threshold T_df
%     4. Isolation: FI(i) = P_r' * V_w(:,i) / (V_w(:,i)' * V_w(:,i))
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
    j = cfg.window_length;           % PCA window length (25)
    pca_thresh = cfg.pca_threshold;  % cumulative variance ratio (0.85)

    % Compute weighted parity matrix
    [~, V_w, W_inv_half] = compute_parity_matrix(H, sigma_vec);
    n_parity = size(V_w, 1);        % 6

    % Compute adaptive threshold
    T_adaptive = compute_adaptive_threshold(H, V_w, sigma_vec, cfg.K_tolerance);

    % Pre-allocate
    FD = zeros(N, 1);
    FI = zeros(N, n);

    % Step 1: compute all weighted parity vectors
    P_w_all = zeros(n_parity, N);
    for k = 1:N
        z = Z(k, :)';
        P_w_all(:, k) = V_w * W_inv_half * z;
    end

    % Step 2-6: sliding window PCA + detection + isolation
    for k = 1:N
        if k < j
            % Not enough samples for full window
            if k < 3
                % Too few for PCA, use raw parity vector
                P_r = P_w_all(:, k);
            else
                P_buffer = P_w_all(:, 1:k);
                P_r = compute_ppv(P_buffer, pca_thresh);
            end
        else
            % Full window: [k-j+1, k]
            P_buffer = P_w_all(:, (k-j+1):k);
            P_r = compute_ppv(P_buffer, pca_thresh);
        end

        % Step 3: detection function
        FD(k) = P_r' * P_r;

        % Step 6: isolation function
        for i = 1:n
            Vi = V_w(:, i);   % i-th column of V_w
            FI(k, i) = abs(P_r' * Vi) / (Vi' * Vi);
        end
    end

end
