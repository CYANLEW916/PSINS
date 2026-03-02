function cfg = config()
% CONFIG  Parameter configuration for SWGLT fault detection.
%   cfg = CONFIG() returns a structure containing all simulation parameters
%   including sensor noise, fault injection, algorithm settings, and plotting.
%
% See also  main_simulation, fdi_glt, fdi_wglt, fdi_swglt.

    %% Simulation parameters
    cfg.ts = 0.02;           % sampling period (s), 50 Hz
    cfg.rng_seed = 42;       % random seed for reproducibility

    %% Sensor noise characteristics (from paper Table 1)
    dph2rps = pi / 180 / 3600;  % 1 deg/h in rad/s
    g0 = 9.8;                   % gravitational acceleration (m/s^2)

    % Gyroscope bias instability (deg/h -> rad/s)
    cfg.sigma_ins_gyro = 0.01 * dph2rps;   % INS1 & INS2 gyro noise std (rad/s)
    cfg.sigma_isis_gyro = 0.1 * dph2rps;   % ISIS gyro noise std (rad/s)

    % Accelerometer bias instability (g -> m/s^2)
    cfg.sigma_ins_acc = 1e-4 * g0;         % INS1 & INS2 accel noise std (m/s^2)
    cfg.sigma_isis_acc = 5e-3 * g0;        % ISIS accel noise std (m/s^2)

    % Noise standard deviation vectors (9 sensors: INS1 xyz, INS2 xyz, ISIS xyz)
    cfg.sigma_gyro = [cfg.sigma_ins_gyro * ones(1,3), ...
                      cfg.sigma_ins_gyro * ones(1,3), ...
                      cfg.sigma_isis_gyro * ones(1,3)];

    cfg.sigma_acc = [cfg.sigma_ins_acc * ones(1,3), ...
                     cfg.sigma_ins_acc * ones(1,3), ...
                     cfg.sigma_isis_acc * ones(1,3)];

    %% Configuration matrix H (9x3)
    cfg.H = repmat(eye(3), 3, 1);
    cfg.n_sensors = 9;
    cfg.n_dof = 3;
    cfg.n_parity = cfg.n_sensors - cfg.n_dof;

    %% Fault injection parameters
    cfg.fault1.sensor_idx = 2;
    cfg.fault1.type = 'gyro';
    cfg.fault1.intervals = [120 135; 440 455; 595 610];
    cfg.fault1.magnitudes = [0.5; 1.0; 2.0] * dph2rps;

    cfg.fault2.sensor_idx = 4;
    cfg.fault2.type = 'acc';
    cfg.fault2.intervals = [100 115; 420 435; 575 590];
    cfg.fault2.magnitudes = [0.005; 0.01; 0.02] * g0;

    cfg.fault3.sensor_idx = 3;
    cfg.fault3.type = 'gyro';
    cfg.fault3.interval = [260 310];
    cfg.fault3.rate = 0.02 * dph2rps;
    cfg.fault3.t_start_ref = cfg.fault3.interval(1);

    %% GLT parameters
    cfg.PFA = 0.01;
    cfg.T_D = chi2inv(1 - cfg.PFA, cfg.n_parity);

    %% SWGLT parameters
    cfg.window_length = 25;              % Sliding window length for SWGLT (samples)
    cfg.K_tolerance = 0.01;              % tolerance factor for adaptive threshold
    cfg.P_D_star = 0.93;                 % target detection probability (eq.80)
    cfg.N_dwell = 2;                     % consecutive-window dwell for isolation
    cfg.rho_threshold = 1.3;             % FI(k1)/FI(k2) margin threshold
    cfg.P_isol = 0.65;                   % Bayesian posterior threshold for isolation
    cfg.min_isolation_votes = 2;         % min passed checks in {margin, posterior, dwell}

end
