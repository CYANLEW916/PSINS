function [V, V_w, W_inv_half] = compute_parity_matrix(H, sigma_vec)
% COMPUTE_PARITY_MATRIX  Compute parity matrices for GLT and WGLT.
%   [V, V_w, W_inv_half] = COMPUTE_PARITY_MATRIX(H, sigma_vec) computes
%   the unweighted parity matrix V and the weighted parity matrix V_w.
%
%   Unweighted V satisfies: V*H = 0, V*V' = I
%   Weighted V_w satisfies: V_w*(W^{-1/2}*H) = 0, V_w*V_w' = I
%
%   Inputs:
%     H          - n x m configuration matrix (9 x 3)
%     sigma_vec  - 1 x n noise standard deviation vector for each sensor
%
%   Outputs:
%     V          - (n-m) x n unweighted parity matrix
%     V_w        - (n-m) x n weighted parity matrix
%     W_inv_half - n x n diagonal matrix W^{-1/2}
%
% See also  fdi_glt, fdi_wglt, fdi_swglt.

    [n, ~] = size(H);

    %% Unweighted parity matrix via SVD
    [U, ~, ~] = svd(H);
    r = rank(H);  % r = 3
    V = U(:, r+1:end)';  % (n-r) x n

    % Verify orthogonality
    assert(norm(V * H) < 1e-10, 'V*H ~= 0');
    assert(norm(V * V' - eye(n - r)) < 1e-10, 'V*V'' ~= I');

    %% Weighted parity matrix
    if nargin >= 2 && ~isempty(sigma_vec)
        W_inv_half = diag(1 ./ sigma_vec);         % W^{-1/2}
        H_w = W_inv_half * H;                      % weighted configuration matrix

        [U_w, ~, ~] = svd(H_w);
        r_w = rank(H_w);
        V_w = U_w(:, r_w+1:end)';                 % (n-r) x n weighted parity matrix

        % Verify: V_w * H_w = V_w * W^{-1/2} * H ≈ 0
        assert(norm(V_w * H_w) < 1e-10, 'V_w * W^{-1/2} * H ~= 0');
    else
        V_w = [];
        W_inv_half = [];
    end

end
