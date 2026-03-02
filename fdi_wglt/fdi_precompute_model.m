function model = fdi_precompute_model(R, H, K, PdStar, j)
%FDI_PRECOMPUTE_MODEL Precompute whitened parity-space model and adaptive threshold.
% Inputs:
%   R      - Measurement noise covariance (n x n).
%   H      - Redundant configuration matrix (n x m).
%   K      - MSE degradation tolerance factor.
%   PdStar - Target detection probability at tolerable fault level.
%   j      - Sliding-window length in samples.
% Output:
%   model  - Struct containing precomputed matrices, vectors, and thresholds.

    n = size(H, 1);
    m = size(H, 2);
    if n <= m
        error('H must be redundant with n > m.');
    end

    if ~isequal(size(R), [n, n])
        error('R must be n-by-n and consistent with H.');
    end

    L = chol(R, 'lower');
    W = L \ eye(n);
    Ht = W * H;

    [U, ~, ~] = svd(Ht, 'full');
    Vt = U(:, m + 1:end)';

    sigma = sqrt(diag(R));
    G = inv(Ht' * Ht);
    G2 = G * G;
    gamma0 = trace(G);

    c = zeros(n - m, n);
    h = zeros(m, n);
    beta = zeros(n, 1);
    c2 = zeros(n, 1);
    ft = zeros(n, 1);
    lambda = zeros(n, 1);

    for i = 1:n
        ei = zeros(n, 1);
        ei(i) = 1;

        c(:, i) = Vt * ei;
        h(:, i) = Ht' * ei;
        beta(i) = h(:, i)' * G2 * h(:, i);
        c2(i) = c(:, i)' * c(:, i);

        ft(i) = sigma(i) * sqrt(K * gamma0 / beta(i) + 1 / c2(i));
        lambda(i) = j * (K * gamma0 * c2(i) / beta(i) + 1);
    end

    lambdaMin = min(lambda);
    Td = chi2inv(0.99, n - m);
    Tadapt = ncx2inv(1 - PdStar, n - m, lambdaMin);
    PfaAdapt = 1 - chi2cdf(Tadapt, n - m);

    model.n = n;
    model.m = m;
    model.R = R;
    model.H = H;
    model.W = W;
    model.Ht = Ht;
    model.Vt = Vt;
    model.sigma = sigma;
    model.G = G;
    model.c = c;
    model.h = h;
    model.beta = beta;
    model.gamma0 = gamma0;
    model.c2 = c2;
    model.ft = ft;
    model.lambda = lambda;
    model.lambdaMin = lambdaMin;
    model.Td = Td;
    model.Tadapt = Tadapt;
    model.PfaAdapt = PfaAdapt;
    model.j = j;
end
