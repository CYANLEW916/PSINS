function P_r = compute_ppv(P_buffer, pca_threshold)
% COMPUTE_PPV  Compute Principal Parity Vector using PCA.
%   P_r = COMPUTE_PPV(P_buffer, pca_threshold) applies PCA to a window of
%   parity vectors and returns the principal parity vector (PPV).
%
%   Inputs:
%     P_buffer      - d x j matrix of j consecutive parity vectors (each column is d x 1)
%     pca_threshold - cumulative variance ratio threshold (e.g., 0.85)
%
%   Output:
%     P_r           - d x 1 principal parity vector
%
% See also  fdi_swglt, config.

    [d, j] = size(P_buffer);

    % Step A2: zero-mean, unit-variance normalization per row
    mu = mean(P_buffer, 2);
    sigma = std(P_buffer, 0, 2);

    % Handle zero-variance rows (constant values)
    sigma(sigma < 1e-15) = 1;

    X_norm = (P_buffer - mu) ./ sigma;

    % Step A3: covariance matrix and eigenvalue decomposition
    C = (X_norm * X_norm') / (j - 1);  % d x d covariance matrix

    [Q, Lambda] = eig(C);
    eigvals = diag(Lambda);

    % Sort eigenvalues in descending order
    [eigvals_sorted, sort_idx] = sort(eigvals, 'descend');
    Q = Q(:, sort_idx);

    % Step A4: select principal components by cumulative variance ratio
    total_var = sum(eigvals_sorted);
    if total_var < 1e-15
        % Degenerate case: no variance
        P_r = mean(P_buffer, 2);
        return;
    end

    cum_ratio = cumsum(eigvals_sorted) / total_var;
    s = find(cum_ratio >= pca_threshold, 1, 'first');
    if isempty(s)
        s = d;  % keep all components
    end

    % Step A5: project onto principal component subspace and reconstruct
    Q_s = Q(:, 1:s);                      % top-s eigenvectors
    X_proj = Q_s * (Q_s' * X_norm);       % reconstruct in PC subspace

    % Denormalize and compute mean
    X_recon = X_proj .* sigma + mu;
    P_r = mean(X_recon, 2);               % principal parity vector (d x 1)

end
