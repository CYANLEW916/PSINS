function [FD, FI, V] = fdi_glt(Z, cfg)
% FDI_GLT  Traditional Generalized Likelihood Test for fault detection and isolation.
%   [FD, FI, V] = FDI_GLT(Z, cfg) applies the standard GLT algorithm assuming
%   uniform noise variance across all sensors.
%
%   Inputs:
%     Z   - N x 9 measurement matrix
%     cfg - configuration structure from config.m
%
%   Outputs:
%     FD  - N x 1 detection function values
%     FI  - N x 9 isolation function values
%     V   - (n-m) x n parity matrix
%
% See also  compute_parity_matrix, fdi_wglt, fdi_swglt, config.

    [N, n] = size(Z);
    H = cfg.H;

    % Compute unweighted parity matrix
    [V, ~, ~] = compute_parity_matrix(H);

    % Traditional GLT uses uniform noise variance (assume INS level)
    % This causes false alarms when ISIS noise is much larger
    sigma_uniform = cfg.sigma_ins_gyro;  % default to gyro; works for both types
    % since the key issue is treating all sensors equally
    sigma2 = sigma_uniform^2;

    % Pre-allocate outputs
    FD = zeros(N, 1);
    FI = zeros(N, n);

    % Process each time step
    for k = 1:N
        z = Z(k, :)';             % 9 x 1 measurement vector

        % Parity vector
        P = V * z;                % (n-m) x 1

        % Detection function: FD = P'*P / sigma^2
        FD(k) = (P' * P) / sigma2;

        % Isolation function for each sensor
        for i = 1:n
            Vi = V(:, i);         % i-th column of V
            FI(k, i) = (P' * Vi)^2 / (sigma2 * (Vi' * Vi));
        end
    end

end
