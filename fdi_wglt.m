function [FD, FI, V_w] = fdi_wglt(Z, sigma_vec, cfg)
% FDI_WGLT  Weighted Generalized Likelihood Test for fault detection/isolation.
%   [FD, FI, V_w] = FDI_WGLT(Z, sigma_vec, cfg) applies the weighted GLT
%   that handles heterogeneous sensor noise by pre-whitening measurements.
%
%   Detection:  P_w = V_w * W * Z,  FD = P_w' * P_w
%   Isolation:  FI(i) = (P_w'*V_w(:,i))^2 / (V_w(:,i)'*V_w(:,i))
%   Threshold:  T_D = chi2inv(1-PFA, n-m) ≈ 16.81  (same as GLT)
%
%   Inputs:
%     Z         - N x 9 measurement matrix
%     sigma_vec - 1 x 9 noise standard deviation vector
%     cfg       - configuration structure from config.m
%
%   Outputs:
%     FD  - N x 1 detection function values
%     FI  - N x 9 isolation function values
%     V_w - (n-m) x n weighted parity matrix
%
% See also  compute_parity_matrix, fdi_glt, fdi_swglt, config.

    [N, n] = size(Z);
    H = cfg.H;

    % Compute weighted parity matrix
    [~, V_w, W, ~] = compute_parity_matrix(H, sigma_vec);

    % Pre-allocate
    FD = zeros(N, 1);
    FI = zeros(N, n);

    for k = 1:N
        z = Z(k, :)';                 % 9x1

        % Weighted parity vector: P_w = V_w * W * z
        P_w = V_w * W * z;   % 6x1

        % Detection function (already normalized by weighting)
        FD(k) = P_w' * P_w;

        % Isolation function using V_w columns
        for i = 1:n
            Vi = V_w(:, i);           % i-th column of V_w
            FI(k, i) = (P_w' * Vi)^2 / (Vi' * Vi);
        end
    end

end
