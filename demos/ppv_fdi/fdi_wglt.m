function [FD, FI, V_w] = fdi_wglt(Z, sigma_vec, cfg)
% FDI_WGLT  Weighted Generalized Likelihood Test for fault detection and isolation.
%   [FD, FI, V_w] = FDI_WGLT(Z, sigma_vec, cfg) applies the weighted GLT
%   algorithm that accounts for heterogeneous sensor noise characteristics.
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
    [~, V_w, W_inv_half] = compute_parity_matrix(H, sigma_vec);

    % Pre-allocate outputs
    FD = zeros(N, 1);
    FI = zeros(N, n);

    % Process each time step
    for k = 1:N
        z = Z(k, :)';                 % 9 x 1 measurement vector

        % Weighted parity vector
        P_w = V_w * W_inv_half * z;   % (n-m) x 1

        % Detection function: FD = P_w' * P_w (already normalized by weighting)
        FD(k) = P_w' * P_w;

        % Isolation function for each sensor
        Vi_w = V_w * W_inv_half;      % weighted parity matrix applied to whitened space
        for i = 1:n
            Vi = Vi_w(:, i);
            FI(k, i) = (P_w' * Vi)^2 / (Vi' * Vi);
        end
    end

end
