function out = fdi_run_sliding_wglt(Z, model)
%FDI_RUN_SLIDING_WGLT Run sliding-window WGLT detection/isolation on measurement data.
% Inputs:
%   Z      - Measurements (N x n).
%   model  - Precomputed model from fdi_precompute_model.
% Output:
%   out    - Struct with FD/FI statistics, alarms, isolation index, and fault estimate.

    N = size(Z, 1);
    n = model.n;
    d = n - model.m;
    j = model.j;

    if size(Z, 2) ~= n
        error('Z must have n columns matching model.n.');
    end

    P = zeros(N, d);
    FDj = nan(N, 1);
    alarm = false(N, 1);
    kstar = zeros(N, 1);
    fhat = nan(N, 1);
    FIj = nan(N, n);

    buf = zeros(d, j);
    cnt = 0;

    for k = 1:N
        p = model.Vt * model.W * Z(k, :)';
        P(k, :) = p';

        cnt = cnt + 1;
        buf(:, mod(cnt - 1, j) + 1) = p;
        if cnt < j
            continue;
        end

        pbar = mean(buf, 2);
        FDj(k) = j * (pbar' * pbar);
        if FDj(k) <= model.Tadapt
            continue;
        end

        alarm(k) = true;
        for i = 1:n
            num = model.c(:, i)' * pbar;
            FIj(k, i) = j * num * num / model.c2(i);
        end

        [~, idx] = max(FIj(k, :));
        kstar(k) = idx;
        fhat(k) = model.sigma(idx) * (model.c(:, idx)' * pbar) / model.c2(idx);
    end

    out.P = P;
    out.FDj = FDj;
    out.FIj = FIj;
    out.alarm = alarm;
    out.kstar = kstar;
    out.fhat = fhat;
end
