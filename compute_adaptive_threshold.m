function T_df = compute_adaptive_threshold(H, V_w, sigma_vec, K)
% COMPUTE_ADAPTIVE_THRESHOLD  Compute adaptive fault-tolerant threshold.
%   T_df = COMPUTE_ADAPTIVE_THRESHOLD(H, V_w, sigma_vec, K) computes the
%   minimum adaptive detection threshold across all sensors.
%
%   For sensor i, tolerable fault amplitude:
%     f_t(i) = sigma_i * sqrt( K * h_i * (H'H) * h_i' + 1/||V_w(:,i)||^2 )
%   Adaptive threshold for sensor i:
%     T_df(i) = ||V_w(:,i)||^2 * f_t(i)^2
%   Final threshold = min over all sensors.
%
%   Inputs:
%     H         - n x m configuration matrix (9 x 3)
%     V_w       - (n-m) x n weighted parity matrix
%     sigma_vec - 1 x n noise standard deviation vector
%     K         - tolerance factor (e.g. 0.01)
%
%   Output:
%     T_df      - scalar adaptive detection threshold
%
% See also  fdi_swglt, config.

    n = size(H, 1);
    HTH = H' * H;   % m x m

    T_df_all = zeros(n, 1);

    for i = 1:n
        sigma_i = sigma_vec(i);
        h_i = H(i, :);                          % 1 x m
        Vi = V_w(:, i);                          % (n-m) x 1
        Vi_norm_sq = Vi' * Vi;                   % ||V_w(:,i)||^2

        % Tolerable fault amplitude
        f_t_i = sigma_i * sqrt(K * (h_i * HTH * h_i') + 1 / Vi_norm_sq);

        % Adaptive threshold for this sensor
        T_df_all(i) = Vi_norm_sq * f_t_i^2;
    end

    % Most sensitive sensor determines detection threshold
    T_df = min(T_df_all);

end
