function [FD, FI, V] = fdi_glt(Z, sigma_uniform, cfg)
% FDI_GLT  Traditional Generalized Likelihood Test for fault detection/isolation.
%   [FD, FI, V] = FDI_GLT(Z, sigma_uniform, cfg) applies the standard GLT
%   algorithm assuming uniform noise variance sigma_uniform^2 for all sensors.
%
%   Using a uniform sigma leads to false alarms when heterogeneous sensors
%   (INS vs ISIS) have very different noise levels.
%
%   Detection:  FD = P'*P / sigma^2,  where P = V * Z
%   Isolation:  FI(i) = (P'*V(:,i))^2 / (sigma^2 * V(:,i)'*V(:,i))
%   Threshold:  T_D = chi2inv(1-PFA, n-m) ≈ 16.81  (chi^2(6) at PFA=0.01)
%
%   Inputs:
%     Z             - N x 9 measurement matrix
%     sigma_uniform - scalar uniform noise std (INS-level, applied to all sensors)
%     cfg           - configuration structure from config.m
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

    sigma2 = sigma_uniform^2;

    % Pre-allocate
    FD = zeros(N, 1);
    FI = zeros(N, n);

    for k = 1:N
        z = Z(k, :)';             % 9x1

        % Parity vector
        P = V * z;                % 6x1

        % Detection function
        FD(k) = (P' * P) / sigma2;

        % Isolation function for each sensor
        for i = 1:n
            Vi = V(:, i);         % i-th column of V
            FI(k, i) = (P' * Vi)^2 / (sigma2 * (Vi' * Vi));
        end
    end

end
