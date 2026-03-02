function [FD_raw, FI_raw, T_raw] = fdi_raw_isis(Z, sigma_vec)
% FDI_RAW_ISIS  Naive ISIS cross-check baseline for fault detection.
%   [FD_raw, FI_raw, T_raw] = FDI_RAW_ISIS(Z, sigma_vec) detects faults by
%   directly comparing each INS unit against the ISIS reference per-axis,
%   without using parity space or whitening. This represents the simplest
%   possible fault detection approach using heterogeneous redundant sensors.
%
%   Residual model:  r_a = Z_INS_a - Z_ISIS_a  ~ N(0, sigma_INS^2 + sigma_ISIS^2)
%   Detection:       FD = sum_a (r_a^2 / sigma_r_a^2)  ~ chi2(3) under H0
%
%   Inputs:
%     Z         - N x 9 measurement matrix (INS1 xyz, INS2 xyz, ISIS xyz)
%     sigma_vec - 1 x 9 noise standard deviation vector
%
%   Outputs:
%     FD_raw - N x 2 detection statistic [INS1_vs_ISIS, INS2_vs_ISIS]
%     FI_raw - N x 6 per-axis normalized residuals squared
%              (columns 1-3: INS1 X/Y/Z, columns 4-6: INS2 X/Y/Z)
%     T_raw  - scalar detection threshold chi2inv(0.99, 3)
%
% See also  fdi_glt, fdi_wglt, fdi_swglt.

    N = size(Z, 1);
    sigma_vec = sigma_vec(:)';

    % Residual variances per axis: var(Z_INS_a - Z_ISIS_a) = sigma_INS_a^2 + sigma_ISIS_a^2
    sigma_r_sq = zeros(1, 6);
    for a = 1:3
        sigma_r_sq(a)   = sigma_vec(a)^2   + sigma_vec(6+a)^2;  % INS1 axis a vs ISIS axis a
        sigma_r_sq(3+a) = sigma_vec(3+a)^2 + sigma_vec(6+a)^2;  % INS2 axis a vs ISIS axis a
    end

    % Compute residuals: INS - ISIS per axis
    r_ins1 = Z(:, 1:3) - Z(:, 7:9);   % N x 3: INS1 - ISIS
    r_ins2 = Z(:, 4:6) - Z(:, 7:9);   % N x 3: INS2 - ISIS

    % Normalized squared residuals (per-axis isolation statistics)
    FI_raw = zeros(N, 6);
    for a = 1:3
        FI_raw(:, a)   = r_ins1(:, a).^2 / sigma_r_sq(a);
        FI_raw(:, 3+a) = r_ins2(:, a).^2 / sigma_r_sq(3+a);
    end

    % Detection statistics: sum over 3 axes per INS unit ~ chi2(3) under H0
    FD_raw = zeros(N, 2);
    FD_raw(:, 1) = sum(FI_raw(:, 1:3), 2);   % INS1 vs ISIS
    FD_raw(:, 2) = sum(FI_raw(:, 4:6), 2);   % INS2 vs ISIS

    % Threshold for chi2(3) at PFA = 0.01
    T_raw = chi2inv(0.99, 3);

end
