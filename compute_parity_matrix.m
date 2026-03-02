function [V, V_w, W, H_w] = compute_parity_matrix(H, sigma_vec)
% COMPUTE_PARITY_MATRIX  Compute unweighted/whitened parity-space matrices.
%   [V, V_w, W, H_w] = COMPUTE_PARITY_MATRIX(H, sigma_vec) computes:
%     V    : unweighted parity matrix, V*H = 0 and V*V' = I
%     W    : whitening matrix R^{-1/2} (eq.13), NOT W^{-1/2}
%     H_w  : whitened configuration matrix H_w = W*H (eq.16)
%     V_w  : whitened parity matrix, V_w*H_w = 0 and V_w*V_w' = I
%
%   Inputs:
%     H         - n x m configuration matrix
%     sigma_vec - 1 x n or n x 1 noise standard deviations
%
%   Outputs:
%     V   - (n-m) x n unweighted parity matrix
%     V_w - (n-m) x n whitened parity matrix
%     W   - n x n whitening matrix R^{-1/2}
%     H_w - n x m whitened configuration matrix
%
% See also  fdi_glt, fdi_wglt, fdi_swglt.

    [n, m] = size(H);

    % Unweighted parity matrix via SVD
    r = rank(H);
    if r < n
        [U_full, ~, ~] = svd(H);
        V = U_full(:, r+1:end)';
    else
        V = zeros(0, n);
    end

    assert(norm(V * H, 'fro') < 1e-10, 'FAIL: V*H != 0');
    assert(norm(V * V' - eye(n - m), 'fro') < 1e-10, 'FAIL: V*V'' != I');

    if nargin < 2 || isempty(sigma_vec)
        V_w = [];
        W = [];
        H_w = [];
        return;
    end

    sigma_vec = sigma_vec(:);
    W = diag(1 ./ sigma_vec);  % eq.(13)
    H_w = W * H;               % eq.(16)

    % Whitened parity matrix using orthonormal null-space basis
    V_w = null(H_w')';         % eq.(17)

    % --- Offline Verification (Theory §7.1) ---
    assert(norm(V_w * H_w, 'fro') < 1e-10, 'FAIL: V_w * H_w != 0');
    assert(norm(V_w * V_w' - eye(n - m), 'fro') < 1e-10, 'FAIL: V_w * V_w'' != I');

    HTH_inv = inv(H_w' * H_w);
    proj_sum = H_w * HTH_inv * H_w' + V_w' * V_w;
    assert(norm(proj_sum - eye(n), 'fro') < 1e-10, 'FAIL: projection complementarity');

    for i = 1:n
        ei = zeros(n, 1);
        ei(i) = 1;
        c_i = V_w * ei;
        h_i = H_w' * ei;
        lhs = h_i' * HTH_inv * h_i + norm(c_i)^2;
        assert(abs(lhs - 1) < 1e-10, ...
            sprintf('FAIL: identity eq.26 for sensor %d', i));
    end
    fprintf('All offline verifications passed.\n');
end
